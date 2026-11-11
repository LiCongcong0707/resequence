#!/bin/bash
# ============================================================
# 07_population_filter.sh
# vcftools 群体级过滤：MAF、biallelic、平均深度、缺失率
# vcftools本身单线程，按染色体拆分后用parallel并行跑
#
# 用法：bash 07_population_filter.sh
# ============================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }
is_step_done() { [[ -f "$1.done" ]]; }
mark_step_done() { touch "$1.done"; }

INPUT_VCF="${WORKDIR}/07_filter/final_pass_variants.vcf.gz"
POP_FILTER_DIR="${WORKDIR}/08_population_filter"
SLICE_DIR="${POP_FILTER_DIR}/chr_slices"
mkdir -p "${POP_FILTER_DIR}" "${SLICE_DIR}"

VCFTOOLS_PARALLEL_JOBS=25   # 按染色体数并行，跟前面Layer2/3思路一致

if [[ ! -f "${INPUT_VCF}" ]]; then
    log "错误：找不到 ${INPUT_VCF}，请先完成SNP/Indel硬过滤+合并"
    exit 1
fi

if [[ ! -f "${INPUT_VCF}.tbi" ]]; then
    log "缺少索引，补建 ..."
    tabix -f -p vcf "${INPUT_VCF}"
fi

# ------------------------------------------------------------
# 单染色体处理函数：切片 -> vcftools过滤 -> 输出
# ------------------------------------------------------------
run_vcftools_for_chr() {
    local chr="$1"
    local slice_vcf="${SLICE_DIR}/${chr}.vcf.gz"
    local out_prefix="${POP_FILTER_DIR}/${chr}_popfilt"
    local marker="${WORKDIR}/logs/popfilter_${chr}"

    if is_step_done "${marker}"; then
        log "[${chr}] 已完成，跳过"
        return 0
    fi

    log "[${chr}] 切片 ..."
    bcftools view -r "${chr}" "${INPUT_VCF}" -Oz -o "${slice_vcf}"
    tabix -f -p vcf "${slice_vcf}"

    log "[${chr}] vcftools 群体过滤中 ..."
    vcftools --gzvcf "${slice_vcf}" \
        --maf 0.05 \
        --min-alleles 2 --max-alleles 2 \
        --min-meanDP 5 \
        --max-missing 0.75 \
        --recode --recode-INFO-all \
        --out "${out_prefix}" \
        2> "${WORKDIR}/logs/popfilter_${chr}.log"

    # vcftools输出的是纯文本.recode.vcf，压缩+建索引方便后续拼接
    bgzip -f "${out_prefix}.recode.vcf"
    tabix -f -p vcf "${out_prefix}.recode.vcf.gz"

    # 切片用完即删，不占地方
    rm -f "${slice_vcf}" "${slice_vcf}.tbi"

    mark_step_done "${marker}"
    log "[${chr}] 完成"
}
export -f run_vcftools_for_chr is_step_done mark_step_done log
export INPUT_VCF SLICE_DIR POP_FILTER_DIR WORKDIR

log "=== 按染色体并行跑vcftools群体过滤 (${VCFTOOLS_PARALLEL_JOBS} jobs) ==="
parallel -j "${VCFTOOLS_PARALLEL_JOBS}" \
    run_vcftools_for_chr {} \
    :::: "${CHR_LIST_FILE}"

# ------------------------------------------------------------
# 拼接所有染色体结果
# ------------------------------------------------------------
log "=== 拼接染色体结果 ==="
FINAL_POP_VCF="${POP_FILTER_DIR}/final_population_filtered.vcf.gz"

CONCAT_LIST="${POP_FILTER_DIR}/concat_list.txt"
: > "${CONCAT_LIST}"
while read -r chr; do
    f="${POP_FILTER_DIR}/${chr}_popfilt.recode.vcf.gz"
    if [[ ! -f "${f}" ]]; then
        log "警告：${chr} 没有找到过滤结果 ${f}，跳过"
        continue
    fi
    echo "${f}" >> "${CONCAT_LIST}"
done < "${CHR_LIST_FILE}"

bcftools concat -f "${CONCAT_LIST}" -Oz -o "${FINAL_POP_VCF}"
tabix -f -p vcf "${FINAL_POP_VCF}"

log "=== 完成 ==="
log "最终群体过滤VCF: ${FINAL_POP_VCF}"

# ------------------------------------------------------------
# 简单统计对比
# ------------------------------------------------------------
log "=== 过滤前后对比 ==="
echo "过滤前变异数: $(zcat "${INPUT_VCF}" | grep -vc '^#')"
echo "过滤后变异数: $(zcat "${FINAL_POP_VCF}" | grep -vc '^#')"
bcftools stats "${FINAL_POP_VCF}" > "${POP_FILTER_DIR}/final_population_filtered.stats.txt"
log "详细统计: ${POP_FILTER_DIR}/final_population_filtered.stats.txt"
