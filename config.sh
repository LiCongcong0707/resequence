#!/bin/bash
# ============================================================
# config.sh
# 全局配置：按需修改，不要改动 pipeline.sh 主逻辑
# ============================================================

# ---------- 基础路径 ----------
export REF="/home/h3043/reseq_work/ref/reference.fasta"                     # 参考基因组（务必与公司一致）
#export DBSNP="/data/ref/dbsnp.vcf.gz"                      # BQSR/已知位点（如需要）
#export MILLS_INDELS="/data/ref/Mills_and_1000G_gold_standard.indels.vcf.gz"

export RAWDATA_DIR="/home/h3043/reseq_work/rawdata"                         # 原始 fastq 目录
export WORKDIR="/home/h3043/reseq_work/workdir/CV_batch01"                  # 工作目录（所有中间结果都在这）
export SAMPLE_SHEET="${WORKDIR}/samples.tsv"                # 样本清单：sample_id\tR1\tR2

mkdir -p "${WORKDIR}"/{01_clean,02_bam,03_gvcf,04_genomicsdb,05_genotype_vcf,06_final,logs}

# ---------- 并行度设置（按你给的架构） ----------
export THREADS_PER_SAMPLE=8          # 单样本内部（fastp/bwa/HC）使用的线程数
export PARALLEL_JOBS=10              # 样本级并行数：fastp/bwa/HC 同时跑多少个样本
export GENOMICSDB_PARALLEL_JOBS=29   # 染色体级并行：GenomicsDBImport 同时跑多少条染色体
export GENOTYPE_PARALLEL_JOBS=29      # 染色体级并行：GenotypeGVCFs 同时跑多少条染色体（内存最吃紧，并发要低）

# ---------- HC 单独拉出来的并行策略 ----------
# HaplotypeCaller 的 --native-pair-hmm-threads 超过4线程后收益递减明显，
# 与其纵向堆线程，不如横向铺并发数，把104个逻辑核用在"更多样本同时算"上
export HC_THREADS=4                  # 单任务线程数，4是PairHMM并行效率的甜点区，超过收益很小
export HC_PARALLEL_JOBS=26           # 4×26=104，逻辑核打满；内存充裕不构成限制
export HC_JAVA_MEM="8g"              # 单样本WGS的GVCF模式一般8G堆内存够用，不需要之前的16G

# 经验公式提示（按你的机器实际内存/核数调整，不要照抄）：
#   PARALLEL_JOBS * THREADS_PER_SAMPLE            <= 总CPU核数
#   GENOTYPE_PARALLEL_JOBS * 单染色体内存峰值(~8-15G) <= 总内存

# ---------- 染色体列表 ----------
# 自动从参考基因组 .fai 提取主要染色体（排除 decoy/alt/random，除非你要求全部纳入分型）
export CHR_LIST_FILE="${WORKDIR}/chr_list.txt"

# ---------- gVCF 样本清单（用于 GenomicsDBImport） ----------
export SAMPLE_MAP_FILE="${WORKDIR}/sample_map.tsv"    # 格式：sample_id\tgvcf_path
