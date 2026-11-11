#!/bin/bash
# ============================================================
# 06_extract_filter_snp_indel.sh
# 从联合分型VCF中拆分SNP/Indel，分别套用GATK硬过滤标准阈值
# （因为HC阶段用了 -stand-call-conf 0，真正的质量把关就靠这一步）
#
# 用法：bash 06_extract_filter_snp_indel.sh
# ============================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

JOINT_VCF="${WORKDIR}/06_final/joint_genotyped.vcf.gz"
FILTER_DIR="${WORKDIR}/07_filter"
mkdir -p "${FILTER_DIR}"

# 是否先把多等位基因位点拆成biallelic（多数下游工具，包括PCA/关联分析工具，
# 都要求biallelic输入；硬过滤阈值也是针对单个等位基因设计的，多等位基因混杂位点
# 会让过滤判断不准确）。默认开启，如果你有特殊理由要保留多等位基因原样，改成false
SPLIT_MULTIALLELIC=false

if [[ ! -f "${JOINT_VCF}" ]]; then
    log "错误：找不到 ${JOINT_VCF}，请先跑完 pipeline.sh merge"
    exit 1
fi

# ------------------------------------------------------------
# Step 0（可选）: 拆分多等位基因位点为biallelic
# ------------------------------------------------------------
INPUT_VCF="${JOINT_VCF}"
if [[ "${SPLIT_MULTIALLELIC}" == "true" ]]; then
    log "=== Step 0: 拆分多等位基因位点(bcftools norm -m -any) ==="
    NORM_VCF="${FILTER_DIR}/joint_genotyped.norm.vcf.gz"
    if [[ ! -f "${NORM_VCF}" ]]; then
        bcftools norm -f "${REF}" -m -any "${JOINT_VCF}" -Oz -o "${NORM_VCF}"
        tabix -p vcf "${NORM_VCF}"
    else
        log "已存在，跳过"
    fi
    INPUT_VCF="${NORM_VCF}"
fi

# ------------------------------------------------------------
# Step 1: 拆分SNP和Indel
# ------------------------------------------------------------
log "=== Step 1: 拆分SNP/Indel ==="

RAW_SNP="${FILTER_DIR}/raw_snps.vcf.gz"
RAW_INDEL="${FILTER_DIR}/raw_indels.vcf.gz"

if [[ ! -f "${RAW_SNP}" ]]; then
    log "提取SNP ..."
    gatk SelectVariants \
        -R "${REF}" \
        -V "${INPUT_VCF}" \
        --select-type-to-include SNP \
        -O "${RAW_SNP}" \
        2> "${WORKDIR}/logs/select_snp.log"
fi

if [[ ! -f "${RAW_INDEL}" ]]; then
    log "提取Indel ..."
    gatk SelectVariants \
        -R "${REF}" \
        -V "${INPUT_VCF}" \
        --select-type-to-include INDEL \
        -O "${RAW_INDEL}" \
        2> "${WORKDIR}/logs/select_indel.log"
fi

# ------------------------------------------------------------
# Step 2: SNP 硬过滤（GATK官方推荐阈值）
# ------------------------------------------------------------
log "=== Step 2: SNP硬过滤 ==="

FILTERED_SNP="${FILTER_DIR}/filtered_snps.vcf.gz"

gatk VariantFiltration \
    -R "${REF}" \
    -V "${RAW_SNP}" \
    --filter-expression "QD < 2.0"              --filter-name "QD2" \
    --filter-expression "FS > 60.0"              --filter-name "FS60" \
    --filter-expression "MQ < 40.0"              --filter-name "MQ40" \
    --filter-expression "MQRankSum < -12.5"      --filter-name "MQRankSum-12.5" \
    --filter-expression "ReadPosRankSum < -8.0"  --filter-name "ReadPosRankSum-8" \
    --filter-expression "SOR > 3.0" --filter-name "SOR3" \
    -O "${FILTERED_SNP}" \
    2> "${WORKDIR}/logs/filter_snp.log"
    # 注：未加入SOR>3.0（GATK较新版本官方推荐会加此项，用于弥补FS在高深度下效力下降的问题）
    # 此处按公司报告给的参数对齐，如需加回SOR，追加一行：
    #   --filter-expression "SOR > 3.0" --filter-name "SOR3" \

# ------------------------------------------------------------
# Step 3: Indel 硬过滤（GATK官方推荐阈值，注意跟SNP阈值不一样）
# ------------------------------------------------------------
log "=== Step 3: Indel硬过滤 ==="

FILTERED_INDEL="${FILTER_DIR}/filtered_indels.vcf.gz"

gatk VariantFiltration \
    -R "${REF}" \
    -V "${RAW_INDEL}" \
    --filter-expression "QD < 2.0"                --filter-name "QD2" \
    --filter-expression "FS > 200.0"              --filter-name "FS200" \
    --filter-expression "ReadPosRankSum < -20.0"  --filter-name "ReadPosRankSum-20" \
    --filter-expression "SOR > 10.0" --filter-name "SOR10" \
    -O "${FILTERED_INDEL}" \
    2> "${WORKDIR}/logs/filter_indel.log"
    # 注：未加入SOR>10.0，同上，按公司报告参数对齐。如需加回：
    #   --filter-expression "SOR > 10.0" --filter-name "SOR10" \

# ------------------------------------------------------------
# Step 4: 合并SNP+Indel过滤结果，产出最终PASS-only高置信度变异集
# ------------------------------------------------------------
log "=== Step 4: 合并并提取PASS位点 ==="

MERGED_FILTERED="${FILTER_DIR}/merged_filtered.vcf.gz"
gatk MergeVcfs \
    -I "${FILTERED_SNP}" \
    -I "${FILTERED_INDEL}" \
    -O "${MERGED_FILTERED}" \
    2> "${WORKDIR}/logs/merge_filtered.log"
tabix -p vcf "${MERGED_FILTERED}"

FINAL_PASS="${FILTER_DIR}/final_pass_variants.vcf.gz"
gatk SelectVariants \
    -R "${REF}" \
    -V "${MERGED_FILTERED}" \
    --exclude-filtered \
    -O "${FINAL_PASS}" \
    2> "${WORKDIR}/logs/select_pass.log"

log "=== 完成 ==="
log "SNP+Indel过滤后合并结果（含FILTER标记，未剔除）: ${MERGED_FILTERED}"
log "最终PASS-only高置信度变异集: ${FINAL_PASS}"

# ------------------------------------------------------------
# Step 5: 简单QC统计，看一眼过滤前后变异数量变化
# ------------------------------------------------------------
log "=== Step 5: 过滤前后统计对比 ==="
echo "原始SNP数量:   $(zcat "${RAW_SNP}" | grep -vc '^#')"
echo "原始Indel数量: $(zcat "${RAW_INDEL}" | grep -vc '^#')"
echo "PASS后总变异数: $(zcat "${FINAL_PASS}" | grep -vc '^#')"
bcftools stats "${FINAL_PASS}" > "${FILTER_DIR}/final_pass_variants.stats.txt"
log "详细统计已写入: ${FILTER_DIR}/final_pass_variants.stats.txt"
