#!/bin/bash
# ============================================================
# pipeline.sh
# 三层并行 WGS 联合变异检测流程
#   Layer 1 (样本级): fastp -> bwa-mem2 -> sort -> markdup -> HaplotypeCaller
#   Layer 2 (染色体级): GenomicsDBImport
#   Layer 3 (染色体级): GenotypeGVCFs
#   收尾: MergeVcfs
# 用法：
#   bash pipeline.sh all            # 跑全流程
#   bash pipeline.sh samples        # 只跑样本级
#   bash pipeline.sh genomicsdb     # 只跑GenomicsDB导入
#   bash pipeline.sh genotype       # 只跑分型
#   bash pipeline.sh merge          # 只跑最终合并
# ============================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

# ------------------------------------------------------------
# 断点续跑标记机制
# ------------------------------------------------------------
is_step_done() {
    local marker="$1"
    [[ -f "${marker}.done" ]]
}

mark_step_done() {
    local marker="$1"
    touch "${marker}.done"
}

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
}

# ------------------------------------------------------------
# 0. 参考基因组索引自动检测/构建（缺失才建，避免返工）
# ------------------------------------------------------------
prepare_reference() {
    log "检查参考基因组索引 ..."

    if [[ ! -f "${REF}.fai" ]]; then
        log "构建 .fai ..."
        samtools faidx "${REF}"
    fi

    local dict="${REF%.fa}.dict"
    if [[ ! -f "${dict}" ]]; then
        log "构建 sequence dictionary ..."
        gatk CreateSequenceDictionary -R "${REF}" -O "${dict}"
    fi

    if [[ ! -f "${REF}.bwt.2bit.64" ]]; then
        log "构建 bwa-mem2 索引（较慢，仅首次执行）..."
        bwa-mem2 index "${REF}"
    fi

    # 生成染色体列表（默认排除 decoy/alt/random/unplaced，如需全部纳入请手动编辑 chr_list.txt）
    if [[ ! -f "${CHR_LIST_FILE}" ]]; then
        log "生成染色体列表 ${CHR_LIST_FILE} ..."
        awk '{print $1}' "${REF}.fai" \
            | grep -viE 'decoy|alt|random|Un|HLA|chrEBV' \
            > "${CHR_LIST_FILE}"
        log "染色体列表内容（请核对是否符合预期，尤其是否需要保留 decoy/alt 用于分型）："
        cat "${CHR_LIST_FILE}"
    fi
}

# ------------------------------------------------------------
# Layer 1a: fastp -> bwa-mem2 -> sort -> markdup（这几步扩展性好，纵向堆线程划算）
# 供 GNU parallel 调用，一行一个样本
# ------------------------------------------------------------
run_sample_align() {
    local sample_id="$1"
    local r1="$2"
    local r2="$3"

    local clean_r1="${WORKDIR}/01_clean/${sample_id}_R1.clean.fastq.gz"
    local clean_r2="${WORKDIR}/01_clean/${sample_id}_R2.clean.fastq.gz"
    local fastp_json="${WORKDIR}/01_clean/${sample_id}.fastp.json"
    local fastp_html="${WORKDIR}/01_clean/${sample_id}.fastp.html"
    local sorted_bam="${WORKDIR}/02_bam/${sample_id}.sorted.bam"
    local dedup_bam="${WORKDIR}/02_bam/${sample_id}.dedup.bam"
    local dedup_metrics="${WORKDIR}/02_bam/${sample_id}.dedup_metrics.txt"

    local marker="${WORKDIR}/logs/${sample_id}"

    # ---- Step 1: fastp 质控清洗 ----
    if ! is_step_done "${marker}.fastp"; then
        log "[${sample_id}] fastp 清洗 ..."
        fastp \
            -i "${r1}" -I "${r2}" \
            -o "${clean_r1}" -O "${clean_r2}" \
            --thread "${THREADS_PER_SAMPLE}" \
            --json "${fastp_json}" --html "${fastp_html}" \
            --detect_adapter_for_pe \
            2> "${WORKDIR}/logs/${sample_id}.fastp.log"
        mark_step_done "${marker}.fastp"
    rm -f "${r1}" "${r2}"
	log "[${sample_id}] 原始 FASTQ 已删除"
    else
        log "[${sample_id}] fastp 已完成，跳过"
    fi

    # ---- Step 2: bwa-mem2 比对 + samtools sort，管道直出，不落地 SAM/未排序BAM ----
    if ! is_step_done "${marker}.align"; then
        log "[${sample_id}] bwa-mem2 比对 + 排序 ..."
        bwa-mem2 mem -M \
            -R "@RG\tID:${sample_id}\tSM:${sample_id}\tPL:MGI\tLB:${sample_id}" \
            -t "${THREADS_PER_SAMPLE}" \
            "${REF}" "${clean_r1}" "${clean_r2}" \
        | samtools sort -@ "${THREADS_PER_SAMPLE}" -m 6G -o "${sorted_bam}" - \
            2> "${WORKDIR}/logs/${sample_id}.align.log"
        samtools index "${sorted_bam}"
        mark_step_done "${marker}.align"
	rm -f "${clean_r1}" "${clean_r2}"
    else
        log "[${sample_id}] 比对+排序 已完成，跳过"
    fi

    # ---- Step 3: MarkDuplicates（自动建索引，避免返工） ----
    if ! is_step_done "${marker}.markdup"; then
        log "[${sample_id}] MarkDuplicates ..."
        gatk MarkDuplicates \
            -I "${sorted_bam}" \
            -O "${dedup_bam}" \
            -M "${dedup_metrics}" \
            --CREATE_INDEX true \
            2> "${WORKDIR}/logs/${sample_id}.markdup.log"
        mark_step_done "${marker}.markdup"
        # 清理排序中间BAM，节省磁盘（可选，按需保留）
         rm -f "${sorted_bam}" "${sorted_bam}.bai"
    else
        log "[${sample_id}] MarkDuplicates 已完成，跳过"
    fi

    log "[${sample_id}] align阶段完成"
}
export -f run_sample_align is_step_done mark_step_done log

# ------------------------------------------------------------
# Layer 1b: HaplotypeCaller 单独拉出（限速步骤，PairHMM线程收益递减，
# 策略是少线程、多并发，横向铺开而非纵向堆线程）
# ------------------------------------------------------------
run_sample_hc() {
    local sample_id="$1"
    local dedup_bam="${WORKDIR}/02_bam/${sample_id}.dedup.bam"
    local gvcf="${WORKDIR}/03_gvcf/${sample_id}.g.vcf.gz"
    local marker="${WORKDIR}/logs/${sample_id}"

    if ! is_step_done "${marker}.hc"; then
        log "[${sample_id}] HaplotypeCaller (GVCF模式, ${HC_THREADS}线程, FASTEST_AVAILABLE PairHMM) ..."
        gatk --java-options "-Xmx${HC_JAVA_MEM}" HaplotypeCaller \
            -R "${REF}" \
            -I "${dedup_bam}" \
            -O "${gvcf}" \
            --emit-ref-confidence GVCF \
            --native-pair-hmm-threads "${HC_THREADS}" \
            -pairHMM FASTEST_AVAILABLE \
            --do-not-run-physical-phasing true \
            --create-output-variant-index false \
            -stand-call-conf 0 \
            2> "${WORKDIR}/logs/${sample_id}.hc.log"

        # --create-output-variant-index false 省了HC阶段的建索引时间，
        # 但 GenomicsDBImport 要求输入gVCF必须有索引，这里单独补建（很快，不是瓶颈）
        log "[${sample_id}] 补建gVCF索引 ..."
        gatk IndexFeatureFile -I "${gvcf}" \
            2>> "${WORKDIR}/logs/${sample_id}.hc.log"

        mark_step_done "${marker}.hc"
    else
        log "[${sample_id}] HaplotypeCaller 已完成，跳过"
    fi
}
export -f run_sample_hc

run_all_samples() {
    log "=== Layer 1a: 样本级并行 - fastp+bwa-mem2+markdup (${PARALLEL_JOBS} jobs, 每样本 ${THREADS_PER_SAMPLE} 线程) ==="

    if [[ ! -f "${SAMPLE_SHEET}" ]]; then
        log "错误：找不到样本清单 ${SAMPLE_SHEET}"
        log "格式要求（无表头，tab分隔）：sample_id<TAB>R1.fastq.gz<TAB>R2.fastq.gz"
        exit 1
    fi

    parallel -j "${PARALLEL_JOBS}" --colsep '\t' \
        run_sample_align {1} {2} {3} \
        :::: "${SAMPLE_SHEET}"

    log "=== Layer 1a 完成 ==="

    log "=== Layer 1b: 样本级并行 - HaplotypeCaller (${HC_PARALLEL_JOBS} jobs, 每样本 ${HC_THREADS} 线程, 堆内存${HC_JAVA_MEM}) ==="
    log "    （HC的PairHMM线程收益在4线程后递减，此处用更高并发换总吞吐，而非单任务堆线程）"

    parallel -j "${HC_PARALLEL_JOBS}" --colsep '\t' \
        run_sample_hc {1} \
        :::: "${SAMPLE_SHEET}"

    log "=== Layer 1b 完成 ==="

    # 生成 GenomicsDBImport 需要的 sample_map.tsv
    log "生成 sample_map.tsv ..."
    : > "${SAMPLE_MAP_FILE}"
    while IFS=$'\t' read -r sample_id r1 r2; do
        echo -e "${sample_id}\t${WORKDIR}/03_gvcf/${sample_id}.g.vcf.gz" >> "${SAMPLE_MAP_FILE}"
    done < "${SAMPLE_SHEET}"

    log "=== Layer 1 全部完成 ==="
}

# ------------------------------------------------------------
# Layer 2: 按染色体独立构建 GenomicsDB（替代 CombineGVCFs）
# ------------------------------------------------------------
run_genomicsdb_for_chr() {
    local chr="$1"
    local db_path="${WORKDIR}/04_genomicsdb/gdb_${chr}"
    local marker="${WORKDIR}/logs/genomicsdb_${chr}"

    if is_step_done "${marker}"; then
        log "[GenomicsDB ${chr}] 已完成，跳过"
        return 0
    fi

    log "[GenomicsDB ${chr}] 导入中 ..."
    # GenomicsDBImport 要求 workspace 目录不能预先存在
    rm -rf "${db_path}"

    gatk --java-options "-Xmx12g" GenomicsDBImport \
        --genomicsdb-workspace-path "${db_path}" \
        --batch-size 50 \
        -L "${chr}" \
        --sample-name-map "${SAMPLE_MAP_FILE}" \
        --reader-threads 2 \
        2> "${WORKDIR}/logs/genomicsdb_${chr}.log"

    mark_step_done "${marker}"
    log "[GenomicsDB ${chr}] 完成"
}
export -f run_genomicsdb_for_chr

run_all_genomicsdb() {
    log "=== Layer 2: 染色体级并行 GenomicsDBImport (${GENOMICSDB_PARALLEL_JOBS} jobs) ==="

    if [[ ! -f "${SAMPLE_MAP_FILE}" ]]; then
        log "错误：找不到 ${SAMPLE_MAP_FILE}，请先跑完 Layer 1 或手动准备该文件"
        exit 1
    fi

    parallel -j "${GENOMICSDB_PARALLEL_JOBS}" \
        run_genomicsdb_for_chr {} \
        :::: "${CHR_LIST_FILE}"

    log "=== Layer 2 完成 ==="
}

# ------------------------------------------------------------
# Layer 3: 按染色体独立 GenotypeGVCFs（内存最吃紧，并发数要低）
# ------------------------------------------------------------
run_genotypegvcf_for_chr() {
    local chr="$1"
    local db_path="${WORKDIR}/04_genomicsdb/gdb_${chr}"
    local out_vcf="${WORKDIR}/05_genotype_vcf/joint_${chr}.vcf.gz"
    local marker="${WORKDIR}/logs/genotype_${chr}"

    if is_step_done "${marker}"; then
        log "[Genotype ${chr}] 已完成，跳过"
        return 0
    fi

    log "[Genotype ${chr}] 分型中 ..."
    gatk --java-options "-Xmx40g" GenotypeGVCFs \
        -R "${REF}" \
        -V "gendb://${db_path}" \
        -O "${out_vcf}" \
        2> "${WORKDIR}/logs/genotype_${chr}.log"

    mark_step_done "${marker}"
    log "[Genotype ${chr}] 完成"
}
export -f run_genotypegvcf_for_chr

run_all_genotype() {
    log "=== Layer 3: 染色体级并行 GenotypeGVCFs (${GENOTYPE_PARALLEL_JOBS} jobs，内存敏感，谨慎调大) ==="

    parallel -j "${GENOTYPE_PARALLEL_JOBS}" \
        run_genotypegvcf_for_chr {} \
        :::: "${CHR_LIST_FILE}"

    log "=== Layer 3 完成 ==="
}

# ------------------------------------------------------------
# 收尾: 按染色体顺序合并所有 VCF
# ------------------------------------------------------------
run_merge() {
    log "=== 最终合并 MergeVcfs ==="

    local final_vcf="${WORKDIR}/06_final/joint_genotyped.vcf.gz"
    local input_args=()

    while read -r chr; do
        local f="${WORKDIR}/05_genotype_vcf/joint_${chr}.vcf.gz"
        if [[ ! -f "${f}" ]]; then
            log "错误：缺少 ${f}，请先完成 Layer 3"
            exit 1
        fi
        input_args+=("-I" "${f}")
    done < "${CHR_LIST_FILE}"

    gatk MergeVcfs \
        "${input_args[@]}" \
        -O "${final_vcf}" \
        2> "${WORKDIR}/logs/merge.log"

    tabix -p vcf "${final_vcf}"
    log "最终联合分型 VCF: ${final_vcf}"
}

# ------------------------------------------------------------
# 主入口
# ------------------------------------------------------------
main() {
    local mode="${1:-all}"

    prepare_reference

    case "${mode}" in
        all)
            run_all_samples
            run_all_genomicsdb
            run_all_genotype
            run_merge
            ;;
        samples)      run_all_samples ;;
        genomicsdb)   run_all_genomicsdb ;;
        genotype)     run_all_genotype ;;
        merge)        run_merge ;;
        *)
            echo "用法: bash pipeline.sh [all|samples|genomicsdb|genotype|merge]"
            exit 1
            ;;
    esac

    log "流程结束（模式: ${mode}）"
}

main "$@"
