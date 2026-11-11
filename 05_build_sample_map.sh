#!/bin/bash
# ============================================================
# 05_build_sample_map.sh
# 自动发现两批命名格式不同的gVCF，统一生成 sample_map.tsv
#   公司老样本: {sample}.gvcf.gz         (在 03_gvcf/gvcf/ 子目录下)
#   新样本:     {sample}.g.vcf.gz        (直接在 03_gvcf/ 下)
# 顺便补建缺失的tabix索引（GenomicsDBImport要求索引必须存在且是最新的）
#
# 用法：bash 05_build_sample_map.sh
# ============================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

COMPANY_GVCF_DIR="${WORKDIR}/03_gvcf/gvcf"      # 公司老样本目录
NEW_GVCF_DIR="${WORKDIR}/03_gvcf"                # 新样本目录（注意不能递归进company子目录）

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

# ------------------------------------------------------------
# 从文件名推断 sample_id，兼容两种命名格式
# ------------------------------------------------------------
extract_sample_id() {
    local fname
    fname="$(basename "$1")"
    if [[ "${fname}" == *.g.vcf.gz ]]; then
        echo "${fname%.g.vcf.gz}"
    elif [[ "${fname}" == *.gvcf.gz ]]; then
        echo "${fname%.gvcf.gz}"
    else
        echo ""   # 不认识的命名格式，调用处会跳过
    fi
}

# ------------------------------------------------------------
# 确保索引存在且比gvcf本身新，缺失/过期就重建
# ------------------------------------------------------------
ensure_index() {
    local gvcf="$1"
    local tbi="${gvcf}.tbi"

    # 只判断存在性，不比较mtime——文件传输(scp/rsync/浏览器下载等)经常打乱修改时间，
    # 会导致明明有效的索引被误判为"过期"而重新生成。
    # 如果担心索引确实损坏/不匹配，可以用 tabix -l 做一次轻量完整性校验（见下方注释）。
    if [[ -f "${tbi}" ]]; then
        return 0
    fi

    log "  建索引: $(basename "${gvcf}")"
    tabix -f -p vcf "${gvcf}"
}

# 如果想在"只看存在性"之外再加一层保险，可以用下面这段替换上面的函数体，
# 用tabix尝试读取索引里的染色体列表，读取失败说明索引确实损坏/不匹配，才重建：
#
 ensure_index() {
     local gvcf="$1"
     local tbi="${gvcf}.tbi"
     if [[ -f "${tbi}" ]] && tabix -l "${gvcf}" &>/dev/null; then
         return 0
     fi
     log "  建索引: $(basename "${gvcf}")"
     tabix -f -p vcf "${gvcf}"
 }

# ------------------------------------------------------------
# 主流程
# ------------------------------------------------------------
: > "${SAMPLE_MAP_FILE}"
declare -A seen_ids   # 用于检测两批之间sample_id是否重名

log "=== 扫描新样本目录: ${NEW_GVCF_DIR} (仅顶层 *.g.vcf.gz，不递归) ==="
while IFS= read -r -d '' f; do
    sample_id="$(extract_sample_id "${f}")"
    if [[ -z "${sample_id}" ]]; then
        log "  警告：无法识别命名格式，跳过 ${f}"
        continue
    fi
    if [[ -n "${seen_ids[${sample_id}]:-}" ]]; then
        log "  错误：sample_id '${sample_id}' 重复出现（${seen_ids[${sample_id}]} 与 ${f}），请检查命名冲突！"
        exit 1
    fi
    seen_ids[${sample_id}]="${f}"
    ensure_index "${f}"
    echo -e "${sample_id}\t${f}" >> "${SAMPLE_MAP_FILE}"
    log "  [新样本] ${sample_id} -> ${f}"
done < <(find "${NEW_GVCF_DIR}" -maxdepth 1 -name "*.g.vcf.gz" -print0)

log "=== 扫描公司老样本目录: ${COMPANY_GVCF_DIR} (*.gvcf.gz) ==="
if [[ ! -d "${COMPANY_GVCF_DIR}" ]]; then
    log "错误：找不到 ${COMPANY_GVCF_DIR}"
    exit 1
fi
while IFS= read -r -d '' f; do
    sample_id="$(extract_sample_id "${f}")"
    if [[ -z "${sample_id}" ]]; then
        log "  警告：无法识别命名格式，跳过 ${f}"
        continue
    fi
    if [[ -n "${seen_ids[${sample_id}]:-}" ]]; then
        log "  错误：sample_id '${sample_id}' 重复出现（${seen_ids[${sample_id}]} 与 ${f}），请检查命名冲突！"
        exit 1
    fi
    seen_ids[${sample_id}]="${f}"
    ensure_index "${f}"
    echo -e "${sample_id}\t${f}" >> "${SAMPLE_MAP_FILE}"
    log "  [公司样本] ${sample_id} -> ${f}"
done < <(find "${COMPANY_GVCF_DIR}" -maxdepth 1 -name "*.gvcf.gz" -print0)

total=$(wc -l < "${SAMPLE_MAP_FILE}")
log "=== 完成，共 ${total} 个样本，清单: ${SAMPLE_MAP_FILE} ==="
log "内容预览："
cat "${SAMPLE_MAP_FILE}"