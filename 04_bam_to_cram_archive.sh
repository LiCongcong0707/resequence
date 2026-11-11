#!/bin/bash
# ============================================================
# 04_bam_to_cram_archive.sh
# 批量BAM -> CRAM，验证完整性后删除原BAM，通常省30~60%空间
# 用法：bash 04_bam_to_cram_archive.sh bam_list.txt
#   bam_list.txt: 每行一个BAM绝对路径
# ============================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

BAM_LIST="${1:?用法: bash 04_bam_to_cram_archive.sh bam_list.txt}"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

convert_one() {
    local bam="$1"
    local cram="${bam%.bam}.cram"

    if [[ -f "${cram}" ]]; then
        log "[$(basename "${bam}")] CRAM已存在，跳过转换"
    else
        log "[$(basename "${bam}")] 转换为CRAM ..."
        samtools view -@ 4 -C -T "${REF}" -o "${cram}" "${bam}"
        samtools index "${cram}"
    fi

    # ---- 完整性校验：比对record数是否一致，不一致绝不删原文件 ----
    local bam_count cram_count
    bam_count=$(samtools view -@ 4 -c "${bam}")
    cram_count=$(samtools view -@ 4 -c -T "${REF}" "${cram}")

    if [[ "${bam_count}" == "${cram_count}" ]]; then
        log "[$(basename "${bam}")] 校验通过 (${bam_count} reads)，删除原BAM"
        rm -f "${bam}" "${bam}.bai"
    else
        log "[$(basename "${bam}")] 警告：record数不一致(BAM=${bam_count}, CRAM=${cram_count})，保留原BAM，不删除！"
    fi
}
export -f convert_one log
export REF

# 并行转换：CPU充裕，这一步IO密集，并发数可以适当高于HC那种计算密集步骤
parallel -j 20 convert_one {} :::: "${BAM_LIST}"

log "全部转换完成，请检查上方日志是否有'警告'字样"
