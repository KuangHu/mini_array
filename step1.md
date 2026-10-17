# Step 1｜CRISPRidentify 完整 array 回收与 repeat 种子库建立

**用途：给本地 coding / bioinformatics agent 的实施文档。**  
**版本：v1.0，2026-10-13。**  
**范围：只实施 Step 1；本文件不表示软件已经安装或真实数据已经运行。**

## 0. 给执行 agent 的任务

在当前项目中实现一条可恢复、可审计的流程：

```text
已有 genome / MAG / phage / plasmid FASTA
    → CRISPRidentify：独立发现 arrays
    → 坐标、repeat/spacer 单元与来源校验
    → 高质量、局部边界完整的 array 种子集合
    → 提取 repeat consensus
    → 双链等价精确去重 + sequence-family 分组
    → 独立的 Cas operon 注释与 subtype 证据关联
    → repeat 库、array 坐标表、family × subtype 支持表
```

这套库将用于后续全基因组搜索分散的 `[repeat][spacer][partial repeat]`。本阶段**不实施**两-repeat array 救回、短 repeat 搜索、15–60 nt 间隔筛选、PAM 推断、未知蛋白关联或 DNA rearrangement 分析。

### 不得改变的设计原则

1. **array 的保留不依赖附近是否有 Cas。** 无邻近 Cas、subtype 未知或分类冲突，都不是丢弃高质量 array 的理由。
2. **sequence family 与 Cas subtype 是两个独立维度。** 同一 repeat / family 可以关联多个 subtype；不得强制一对一。
3. **后续搜索库保留所有高质量 NR100 consensus，而不只保留每个 family 的代表。**
4. **高质量种子与原始候选分开保存。** Possible、边界不确定及低分候选不进入主种子库，但不得从审计记录中消失。
5. 此处“完整”指**局部 array 不明显受 contig 边界或缺口截断**，不要求宿主是 complete genome，也不代表已经用 reads 排除了组装错误。
6. CRISPRidentify 是主 array caller。不要在它之前强制使用 MinCED 命中作为通行条件，也不要以 CRISPRCasTyper 的 array 检出结果取交集。

**执行顺序：环境与接口核对 → 小型测试 → 最多 100 个输入的 pilot → 汇报。默认不自动启动全库。**

---

## 1. 输入与项目布局

优先复用当前项目已有 FASTA 和 manifest。先检查当前项目已知目录，不要擅自下载新数据库、猜测集群路径或扫描整个文件系统。没有可用输入时，报告缺失的 manifest / FASTA 路径。

```text
project/
  config/step1.yaml
  input/genomes.tsv
  input/staged/
  software/CRISPRidentify/
  software/CRISPRCasTyper/
  envs/identify/
  envs/driver/
  scripts/step1.py
  tests/
  work/
  raw/identify/
  raw/cctyper/
  results/step1/
  provenance/
  logs/
```

### `input/genomes.tsv`

UTF-8 TSV；一个 FASTA 输入实体一行。文件路径相对 manifest 所在目录解析，运行时转成绝对路径。

| 字段 | 要求 |
|---|---|
| `genome_id` | 必填、唯一；代表一个 assembly / MAG / 独立 phage 或 plasmid 输入 |
| `fasta_path` | 必填；允许 `.fna/.fa/.fasta` 及对应 gzip 文件 |
| `source_type` | 必填：`isolate`, `MAG`, `phage`, `plasmid`, `unknown` |
| `taxon` | 可选；原有物种或更高分类标签，不在本阶段重新做分类 |
| `sample_id` | 可选；保留已有采样来源 |
| `independent_group_id` | 可选；已有去冗余簇、样本独立组或其他明确定义的独立单元 |
| `completeness`, `contamination` | 可选；只复用已有数值，缺失写 `NA` |

单个细菌 assembly 内的 chromosome 和 plasmid 可保留在同一个 FASTA；不要把每条 contig 误当成独立基因组。独立提交的 phage/plasmid 也不要仅凭相似物种标签推断宿主。

### FASTA 预处理

- 不修改原文件；解压、清理 header 后写入 `input/staged/`。
- 为输入分配稳定的安全 `genome_uid`，保留它到原始 `genome_id` 的映射。后文命令中的 `G` 使用该 UID。
- 给每条记录分配全局唯一、仅含字母数字的 `contig_uid`，例如基于 `genome_id + 原始记录索引 + sequence SHA256` 的 hash。
- 输出 `contigs.tsv`：`genome_id, genome_uid, contig_uid, original_id, original_header, record_index, length, sequence_sha256, topology`。
- 不直接使用可能含空格、点号、竖线的原始 header 作为工作文件名。CRISPRidentify 的输入拆分代码会改写部分 header，必须提前消除碰撞风险。[R5]
- 序列统一大写；保留并记录 IUPAC 模糊字符，不擅自替换为某个碱基；非法字符导致输入校验失败。
- 不按 N50、最短 contig 或是否含 Cas 预先删除整个输入。局部 array 质量在后面判断。
- `topology` 未提供则为 `unknown`；不得从文件名推断 circular。
- 记录输入原文件和 staged FASTA 的 hash。完全重复的输入可复用计算，但保留所有来源映射；不得算成新的独立支持。

`independent_group_id` 缺失时，只能报告去重输入数或 assembly 数，不能称为“独立谱系数”。

---

## 2. 软件安装与版本锁定

### 工具分工

| 工具 | 本阶段用途 |
|---|---|
| CRISPRidentify | array 检出、分数、repeat/spacer 单元、consensus |
| CRISPRCasTyper | Cas genes / operons 的独立注释和 subtype 分类 |
| Python + edlib | 解析、序列校验、短 repeat 的全长编辑距离比较 |
| pytest | 单元测试与回归测试 |

CRISPRidentify 当前上游环境固定了 Python 3.7.6 及旧版机器学习依赖；当前 CRISPRCasTyper 文档要求 Python ≥3.10。**必须分环境，不能把两者硬装到一个环境中。** [R3, R6]

以下是安装骨架；`PROJECT_ROOT` 由 agent 根据实际项目设置。初次下载后立即记录 commit，测试通过后锁定。后续复跑不得再次追踪浮动的 `master`。

```bash
set -euo pipefail
P="${PROJECT_ROOT:?Set PROJECT_ROOT to the existing project directory}"
mkdir -p "$P/software" "$P/envs" "$P/provenance"

git clone --recurse-submodules \
  https://github.com/BackofenLab/CRISPRidentify.git \
  "$P/software/CRISPRidentify"
git clone \
  https://github.com/Russel88/CRISPRCasTyper.git \
  "$P/software/CRISPRCasTyper"

git -C "$P/software/CRISPRidentify" rev-parse HEAD \
  > "$P/provenance/CRISPRidentify.commit"
git -C "$P/software/CRISPRCasTyper" rev-parse HEAD \
  > "$P/provenance/CRISPRCasTyper.commit"

mamba env create -p "$P/envs/identify" \
  -f "$P/software/CRISPRidentify/environment.yml"

mamba create -y -p "$P/envs/driver" python=3.11 pip
conda run -p "$P/envs/driver" python -m pip install \
  "$P/software/CRISPRCasTyper" edlib pyyaml pytest biopython

conda list -p "$P/envs/identify" --explicit \
  > "$P/provenance/identify.conda.explicit.txt"
conda list -p "$P/envs/driver" --explicit \
  > "$P/provenance/driver.conda.explicit.txt"
conda run -p "$P/envs/driver" python -m pip freeze \
  > "$P/provenance/driver.pip.freeze.txt"
```

已安装时先核对版本和路径，不重复 clone，不覆盖现有环境。若旧环境无法求解，记录 blocker；不得默默升级 scikit-learn / TensorFlow 后假定旧模型兼容。

### 运行前必须完成

保存两个工具的 `--help`、commit、数据库/模型路径与 SHA256、环境锁文件。检查 CRISPRidentify 的主分类模型、orientation 所需资源及底层可执行程序可用。

本阶段设 `--cas False`，不依赖 CRISPRidentify 内嵌 CRISPRcasIdentifier 的额外模型下载。Cas 注释走独立的 CRISPRCasTyper。

**以锁定 checkout 的 CLI 和实际输出为准。** 上游 README、不同版本及源码可能不完全一致；发现差异时在 `provenance/interface_notes.md` 中记录，不猜测参数或列名。

---

## 3. 运行 CRISPRidentify

### 默认策略

对 staged FASTA **按输入实体独立运行**；每个任务使用独立 cwd、结果目录和 JSON 目录。先用少量 CPU 测试，不申请 GPU，不让每个任务默认占用全部 CPU。

下列原生参数已在上游 CLI 核对；实际实现必须由配置生成命令，而不是在多处维护不同参数。[R2]

```bash
set -euo pipefail
P="${PROJECT_ROOT:?Set PROJECT_ROOT}"
G="${GENOME_UID:?Set GENOME_UID to a validated safe identifier}"
THREADS="${THREADS:-4}"
RUN="$P/work/identify/$G"

mkdir -p "$RUN/cwd" "$RUN/json" "$P/raw/identify" "$P/logs"
cd "$RUN/cwd"

export CUDA_VISIBLE_DEVICES=""
export OMP_NUM_THREADS="$THREADS"
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

conda run --no-capture-output -p "$P/envs/identify" \
  python "$P/software/CRISPRidentify/CRISPRidentify.py" \
  --file "$P/input/staged/$G.fna" \
  --result_folder "$P/raw/identify/$G" \
  --json_report "$RUN/json" \
  --fasta_report True \
  --model ALL \
  --cas False \
  --is_element False \
  --strand True \
  --cpu "$THREADS" \
  --parallel True \
  --fast_run False \
  --min_repeats 3 \
  --min_len_rep 21 \
  --max_len_rep 55 \
  --min_len_spacer 18 \
  --max_len_spacer 78 \
  > "$P/logs/identify.$G.log" 2>&1
```

注意：

- `--json_report` 的值是目录，不是 `True`。
- 本阶段使用工具默认的长度搜索范围；**不要把后续 mini-array 的 15–60 nt spacer 窗口提前施加到这里**。这些上游长度参数包含平均长度约束，不能理解成每个单元都必须落在范围内。
- 保留默认的 degenerate-repeat enhancement。当前源码部分参数使用 Python `type=bool`；不要假设给 `--degenerated` 或 `--enhancement_*` 传字符串 `False` 一定等同于关闭。[R2]
- 所有类别的原始结果保留。`--fasta_report` 的导出集合不一定等同于全部候选，不能只读导出的 FASTA 就声称已解析所有类别。[R4]
- 不以“进程退出码为 0”作为唯一成功条件；还要验证输入记录覆盖、输出格式和子程序错误。分别记录 `success`, `success_zero_arrays`, `failed`, `parse_failed`。

---

## 4. 解析与坐标归一化：先保证没有错位

实现版本特异的 adapter，不使用没有校验的正则表达式拼接整张 summary。

### 权威数据与检查

读取每条 contig 的 JSON（有 repeat 序列和起点），同时读取 CSV / 文本用于分数和类别交叉验证。保留 raw category，另映射成统一标签：

```text
bona_fide
alternative
possible
possible_discarded
low_score
```

当前源码中的 JSON 字段包括 `list_repeats`, `list_repeat_starts`, `list_spacers`, `consensus_repeat`, `strand`, `certainty_score`。CSV 的字段及大小写用实际文件头读取。[R4]

### 统一坐标

所有项目表和 BED 使用 **0-based、左闭右开 `[start0, end0)`**；原始字段单独保存。

**不要将 JSON 顶层 `start/end/length` 不加检查地当成 BED。** 上游 JSON/CSV 输出路径的坐标处理并非完全相同，而且负链表示会反转 repeat 列表。[R4]

按锁定版本的逐 repeat 起点构建区间，并用原始 FASTA 逐碱基验算：

1. 确认 `list_repeat_starts` 的坐标语义；本次查看源码中它保留内部 0-based 起点。
2. 每个 repeat 与对应起点成对处理。`Reversed` 输出的 repeat 需要反向互补后才对应 genomic-forward 序列。
3. 将成对的单元按 genomic-forward 起点排序。
4. `repeat_end0 = repeat_start0 + len(repeat_sequence_forward)`。
5. 对每个 repeat 断言 `reference[start0:end0] == sequence_forward`。
6. spacer 用相邻 repeat 之间的原始 FASTA 区间提取，再与工具 spacer 列表进行方向一致的核对。
7. array 区间由所有 repeat 的最小起点、最大终点确定。
8. 断言 `n_spacers == n_repeats - 1`，并且 `R1+S1+...+Rn` 与 reference array 区间完全一致。

`Forward (Orientation was not computed)` 表示输出采用正向表示，**不等于转录方向已经确定**。分别保存 `sequence_representation` 和 `orientation_prediction`；后者缺失写 `unknown`。

坐标/序列不能自洽时隔离该调用并报错，不能偷偷移动一个碱基“修好”。对已知的上游坐标展示差异使用明确 adapter，并记录差异。

### 零结果、重复表示与 consensus

- `Summary.csv` 缺失或汇总文件出现 `No arrays found` 可以是正常零结果；必须结合该任务完成状态和可解析的空候选结果判断。[R4]
- 一个 locus 的 Alternative representation 不算新的独立 array。保留 `locus_id` 与 `call_id`，主种子只选 bona-fide representation。`locus_id` 对应同一 contig 的上游 region/array index；`call_id` 再区分类别与具体表示。发布的 `array_id` 用来源、校验后的区间和 consensus 生成稳定 ID，不使用跨文件不唯一的数字编号。
- 保留 `consensus_reported` 和方向信息。主库使用 CRISPRidentify 的 consensus；不在多个 array 或一个 family 内重新投票生成“新 consensus”。
- 如果 consensus 含非 A/C/G/T 符号，先进入 reserve；不得随意补碱基或删除字符。
- 未能通过坐标验证的结果不得进入任何高质量导出。

---

## 5. 高质量种子定义与保留层级

以下是**本项目的保守默认值**，不是已验证的通用最优阈值。CRISPRidentify 的 certainty score 是分类器分数；`0.90` 不能解释为真实世界中“90% 为真”。上游 bona-fide 判定对应的分数边界约为 `0.75`。[R1]

### `seed_hq`：同时满足

| 条件 | 默认规则 |
|---|---|
| 调用类别与分数 | primary `bona_fide`，`score >= 0.90` |
| 单元结构 | 所有坐标/序列验证通过；至少 3 个完整、相似的 repeat，至少 2 个 spacer |
| 完整 repeat 支持 | 与本 array consensus 的全长归一化编辑相似度 ≥0.90，长度比 ≥0.90；至少 3 个 repeat 达标 |
| Spacer 非简单重复 | 至少 2 个不同的 spacer 序列；这里只用精确去重，不用预设 target |
| 局部边界 | array 两端各有至少 100 bp 可检查的 contig 序列 |
| 缺口和模糊碱基 | array 本体及两端各 100 bp 内无非 A/C/G/T 字符 |
| Repeat query | consensus 非空，且只含 A/C/G/T |

完整 repeat 的相似度定义与第 7 节相同；copy 与 consensus 已统一方向时，不需要改变其生物学方向来“挽救”匹配。

### 保留但不进入主库

- `seed_reserve`：例如 bona-fide 分数 0.75–0.90、Possible、完整 repeat 支持不足或局部边界不确定；保留具体原因，允许后续单独分析。
- `context_only`：Alternative / Possible Discarded / Low score 等原始调用，保留位置和审计信息，不作为独立种子。
- 解析失败或输入异常不是生物学低置信候选；放入错误记录，不混入 reserve 统计。

**正常末端退化不是整条 array 的自动否决项。** 如果仍有至少 3 个完整支持 repeat 且局部边界可靠，可以进入 HQ，并标注 `terminal_repeat_degenerate=true`。不要从退化末端另外生成短 consensus。

已知 circular contig 在线性化起点附近的 array，本版标记 `circular_origin_unresolved` 并进入 reserve；不假装已经完成跨起点恢复。不要给混合 FASTA 中的每条记录一律加 circular 属性。

“没有明显边界/缺口问题”不等于“没有 assembly artifact”。没有 reads 或 assembly graph 支持时，`assembly_validation = not_assessed`。

### 必须保留的对照统计

在同一批 pilot 输出中，比较 `score >=0.75` 与 `>=0.90` 在其余 QC 不变时的 array 数、NR100 repeat 数和新增 family 数。不要只报告严格过滤后的成功数字。

---

## 6. 提取 repeat consensus 与 NR100 库

对于每个 HQ array：

```text
r = uppercase(consensus)
canonical(r) = lexicographically_min(r, reverse_complement(r))
repeat_id = "rep_" + SHA256(canonical(r))[:20]
```

检测 hash 碰撞；发生时使用完整 SHA256。反向互补等价只用于去重，不能用 `canonical` 的方向冒充真实转录方向。

保存完整来源关系：

```text
array_id → repeat_id
repeat_id → canonical_sequence
repeat_id → genome_id / contig_uid / array_id / source_strand
```

同一个 repeat 出现在 100 个近克隆 assembly 中，不应报告为 100 次独立起源。分别计算 `n_arrays`, `n_assemblies`, `n_samples`, `n_independent_groups`；未知维度写 `NA`。

**主搜索库：`repeats.hq.nr100.fna`。**  
这是全部 HQ consensus 的 exact/RC 去重库，而不是 RF90 family 的代表库。

---

## 7. 按 sequence family 分组：不丢失成员

本阶段的 family 是**操作性相似序列分组**，不宣称等同于进化家族或功能兼容组。

默认采用确定性的、全长比对的 greedy-centroid clustering。对短 repeat 使用编辑距离，避免把很短的局部 motif 命中当成全长同源。

```text
d(a,b) = min(
    global_Levenshtein_distance(a, b),
    global_Levenshtein_distance(a, reverse_complement(b))
)

similarity(a,b) = 1 - d(a,b) / max(len(a), len(b))
length_ratio(a,b) = min(len(a), len(b)) / max(len(a), len(b))

pass ⇔ similarity >= 0.90 AND length_ratio >= 0.90
```

这是“归一化编辑相似度”，**不要在报告中改称 BLAST percent identity**。

### 聚类约束

- 输入为 NR100 序列；按序列长度降序、canonical sequence 字典序排序，保证重跑确定性。
- 新序列与所有长度范围允许的已建 centroid 比较；多个 centroid 通过时选距离最小者，平局按 centroid ID。
- 无 centroid 通过则建立新 family；成员始终必须满足其 centroid 的阈值。
- 不使用单连接传递合并：`A~B`、`B~C` 并不保证 `A~C`。
- 即使已 canonicalize，近似匹配时仍比较另一个序列的两个方向；lexicographic canonicalization 不能保证近似序列选择同一方向。
- 不做 circular permutation，不把不同边界的任意短重叠强行并组。
- 保存 `family_id, representative_repeat_id, member_repeat_id, edit_distance, similarity, aligned_orientation`。
- `family_id` 根据 centroid 序列及算法版本生成；同一输入重跑稳定。新增样本导致重新聚类时升库版本，不承诺 family ID 永远不变。

edlib 可用于精确 global edit-distance 计算；先按长度分桶。pilot 实测 unique-repeat 规模和聚类耗时，不生成不必要的全距离矩阵。需要加速时必须保留上述最终判定规则，并通过小规模 exhaustive 比较验证，不能无声改成短 k-mer 命中聚类。

导出 `repeat_families.tsv` 和 `repeat_family_members.tsv`。**family representative 仅用于摘要；不能替代 NR100 搜索库。**

---

## 8. Cas subtype 注释：独立证据，而不是保留门槛

CRISPRCasTyper 在同一份 staged FASTA 上运行，负责识别 Cas operons。它自己的 array calling 结果仅作审计参考；不得覆盖 CRISPRidentify 的坐标或决定 HQ array 是否保留。

对 isolate 默认 `--prodigal single`；对 MAG / phage / plasmid 默认 `meta`。`unknown` 默认 `meta` 并记录。此为项目运行选择，可在 pilot 后按输入情况调整。[R6]

```bash
set -euo pipefail
P="${PROJECT_ROOT:?Set PROJECT_ROOT}"
G="${GENOME_UID:?Set GENOME_UID}"
MODE="${PRODIGAL_MODE:?Set PRODIGAL_MODE to single or meta}"
THREADS="${THREADS:-4}"

mkdir -p "$P/raw/cctyper" "$P/logs"
conda run --no-capture-output -p "$P/envs/driver" \
  cctyper "$P/input/staged/$G.fna" "$P/raw/cctyper/$G" \
  --threads "$THREADS" \
  --prodigal "$MODE" \
  --ccd 10000 \
  --skip_blast \
  --no_plot \
  > "$P/logs/cctyper.$G.log" 2>&1
```

`--skip_blast` 关闭它的已知 repeat 补充搜索；本阶段不需要这条补充检出路径。上述 flags 已在原生 CLI 核对。[R7]

### 读取什么、如何关联

优先使用 `cas_operons.tab` 中的 Cas-based `Prediction`、operon 区间和支持基因；原始 putative / ambiguous / hybrid 结果也保留。**不要把 `CRISPR_Cas.tab` 的综合 Prediction 当作完全独立于 repeat 序列的证据。** [R6]

CRISPRCasTyper 的部分结果表仅在有相应结果时生成；确认任务成功且确实无 Cas 结果时，缺少 `cas_operons.tab` 可以是正常空结果。仅凭文件缺失不能判断无 Cas，也不能一律判为任务失败。[R6]

将 Cas 坐标按锁定版本归一化，并测试 off-by-one。对同一 contig 的 array `[a0,a1)` 与 operon `[c0,c1)`：

```text
gap_bp = max(0, c0-a1, a0-c1)
```

默认仅关联 `gap_bp <= 10000` 的 operon；10000 bp 是操作窗口，不是功能归属证明。

| 情形 | `cas_link_status` 与处理 |
|---|---|
| 窗口内恰有一个可信 operon | `proximal_operon`；记录 subtype、距离、operon ID |
| 多个可信 operon，但已判定 subtype 相同 | `multiple_operons_same_subtype`；保留全部 operon ID |
| 多个不同 subtype，或含无法解析的 hybrid / ambiguous | `ambiguous_local`；保留候选集合，不强行最近者定型 |
| 只有低置信 operon | `putative_only`；不当作确定 subtype |
| 无窗口内 operon | `no_proximal_operon`；保留 HQ array；不是“全基因组没有 Cas”的结论 |
| Cas 任务失败 | `annotation_failed`；不能伪装为 orphan / unknown |

不得将另一条 contig 上的 Cas，仅凭属于同一个 MAG，就赋给该 array。不得根据同 family 的其他成员替未分类 array “补齐”局部 subtype。

### 可选 RepeatTyper

不是完成 Step 1 的必要条件。启用时，单独输出 `repeat_predicted_subtype` 和 `repeat_prediction_probability`；不能覆盖 `cas_subtype_local`，也不能据此过滤未知 repeat。它的输入是每行一条 repeat 的纯文本，不是 FASTA。[R6]

出现预测与局部 Cas 不一致时，保留 `subtype_conflict`。即使同一个 exact repeat 对应多个 Cas subtype，也保留多对多关系，不用多数投票抹掉少数关联。

---

## 9. 交付数据与字段

全部 TSV 带 header；缺失值统一为 `NA`；有列表的字段使用 JSON 数组，不使用 Python `eval` 解析。稳定排序后导出。

| 文件 | 内容 |
|---|---|
| `contigs.tsv` | 原始 ID、staged ID、长度、hash、来源和 topology |
| `arrays.all.tsv` | 所有成功解析的候选、类别、分数、QC 标签、locus/call ID |
| `arrays.hq.tsv` | HQ primary arrays，包含 repeat_id、family_id、局部 Cas 证据 |
| `arrays.reserve.tsv` | 保留候选与未进 HQ 的原因 |
| `array_units.tsv.gz` | 每个成功解析调用的 repeat/spacer 坐标、顺序和 genomic-forward 序列 |
| `repeats.hq.nr100.fna` | 后续全基因组搜索使用的主 query 库 |
| `repeat_sources.tsv` | repeat_id → array_id → genome/contig/方向的全部来源关系 |
| `repeat_family_members.tsv` | 每个 NR100 repeat 的 centroid、距离和方向 |
| `repeat_families.tsv` | 每个 family 的代表、成员数、各层级来源计数 |
| `array_cas_links.tsv` | 每个 array 与局部 Cas operon 的证据，不限一个 operon |
| `family_subtype_summary.tsv` | family × subtype × evidence_status 的支持数 |
| `known_array_regions.bed` | 坐标可靠的 bona-fide arrays，包括未达 HQ seed 阈值者 |
| `array_context_regions.tsv` | Alternative 边界、Possible 和端部状态等软上下文，供后续排除假 mini-array |
| `run_status.tsv` | 每个输入、每个阶段的成功、零结果、失败、耗时与资源记录 |
| `report.md` | pilot 结果、测试、过滤损失、资源需求和剩余问题 |
| `provenance.json` | 输入/参数/软件/模型/输出 hash，运行时间与版本 |

`known_array_regions.bed` 使用全局 `contig_uid`。保留原 array ID，不随意向外扩展几 kb；Alternative 和 Possible 是独立的软上下文，不一律并入硬 mask。

### `arrays.all.tsv` / `arrays.hq.tsv` 最小字段

```text
array_id, call_id, locus_id, genome_id, contig_uid,
start0, end0, category_raw, category, score,
n_repeats, n_full_repeats, n_spacers, n_unique_spacers,
consensus_reported, consensus_canonical, repeat_id, family_id,
orientation_prediction, coordinate_check_pass,
left_flank_bp, right_flank_bp, ambiguous_bases_in_context,
terminal_repeat_degenerate, qc_tier, qc_reasons,
cas_subtype_local, cas_link_status, cas_operon_ids,
assembly_validation
```

`family_subtype_summary.tsv` 至少包含：

```text
family_id, subtype, evidence_status,
n_arrays, n_repeat_sequences, n_assemblies,
n_samples, n_independent_groups
```

这里的数组数使用独立 primary arrays，不能把 Alternative 重复计数。

---

## 10. 默认配置

由 agent 保存为 `config/step1.yaml`。除工具默认搜索参数外，QC/聚类阈值均是本项目设计值，必须随结果记录。

```yaml
schema_version: 1
manifest: input/genomes.tsv
output_dir: results/step1
seed: 42

execution:
  max_genomes: 100
  allow_full_run: false
  workers: 2
  threads_per_job: 4
  resume: true

identify:
  model: ALL
  min_repeats: 3
  min_len_rep: 21
  max_len_rep: 55
  min_len_spacer: 18
  max_len_spacer: 78
  fast_run: false
  strand: true
  cas: false
  is_element: false

seed_qc:
  min_score: 0.90
  min_full_repeats: 3
  min_unique_spacers: 2
  context_bp: 100
  max_ambiguous_bases_in_context: 0
  min_copy_normalized_edit_similarity: 0.90
  min_copy_length_ratio: 0.90
  require_cas: false

families:
  method: deterministic_greedy_centroid
  min_normalized_edit_similarity: 0.90
  min_length_ratio: 0.90
  compare_both_strands: true
  export_all_nr100_queries: true

cas_annotation:
  enabled: true
  max_local_gap_bp: 10000
  cross_contig_assignment: false
  repeat_typer_optional: false
```

配置中的相对路径按项目根目录解析；manifest 内的 FASTA 相对路径按 manifest 所在目录解析。

同一个输出目录不得混合不同参数/模型版本。缓存键至少含 input hash、工具 commit、模型/数据库 hash、参数 hash 和 adapter 版本。

---

## 11. 实现接口与失败处理

由 agent 新建或接入项目现有 runner；不要求为此重写整个项目，也不要求新增复杂工作流框架。至少包含这些可独立测试的组件：

```text
prepare_inputs
run_identify
parse_identify
apply_seed_qc
build_repeat_catalog
cluster_repeat_families
run_cctyper
link_array_cas
export_and_report
```

以下是**待实现的项目 CLI 契约，不是 CRISPRidentify 自带命令**：

```bash
python scripts/step1.py doctor --config config/step1.yaml
python scripts/step1.py run --config config/step1.yaml --limit 100 --resume
python scripts/step1.py validate --config config/step1.yaml
```

每个任务使用临时输出，校验通过后原子提交。中断重跑不得复用半截输出；输入、参数或版本变化必须触发相应重算。并发资源满足 `workers × threads_per_job` 不超过实际分配 CPU。

单任务失败可以继续处理其他输入，但最终报告必须写明分母、失败清单及 `complete` / `partial` 状态。失败不是“未检测到 array”。

---

## 12. 必须包含的测试

### A. 单元 / 回归测试

| 测试 | 验收条件 |
|---|---|
| Header 与来源碰撞 | 两个 genome 均含 `contig_1`、带版本点号的 accession 不互相覆盖 |
| 正反链解析 | synthetic/固定 raw fixture 的 repeat、spacer 均能逐碱基映射回 reference |
| 坐标边界 | 专门覆盖 JSON/CSV 的端点差异；人为 ±1 错位被发现，不静默修复 |
| 真零结果 vs 失败 | 合法空输出给 `success_zero_arrays`；崩溃或缺损输出不能假装空 |
| HQ vs reserve | 同一合法结构只改变 score / 端部缺口时，分层按配置改变 |
| 末端退化 | 3 个完整 repeat + 合理退化末端不会仅因末端退化被拒；2 个 repeat 不进入本阶段 HQ |
| Spacer 规则 | 完全相同的 spacer 串列不能过 diversity gate；长度 >60 nt 不因 Step 2 窗口被预先删除 |
| RC 去重 | `r` 与 `RC(r)` 只有一个 repeat_id，保留两条来源与方向 |
| 非传递聚类 | 构造 A 先成为 centroid 的 A~B、B~C、A!~C fixture；C 不能仅因匹配 B 而并入 A |
| Orphan 保留 | 移除 Cas 链接不改变已合格 array 的 HQ 状态 |
| Subtype 多对多 | 同 repeat 关联两个 subtype 时两种证据均保留 |
| 跨 contig 隔离 | 同 MAG 不同 contig 的 Cas 不形成局部链接 |
| 可恢复与确定性 | 原输入重跑不改变结果；修改输入 hash 后重算；重复运行无重复行 |

### B. 生物学 smoke test

使用至少两个有独立注释来源的完整 CRISPR array 阳性实例，来自不同 repeat 序列；保存 accession、来源和预期区间。可复用上游测试文件，但**上游输出本身不是独立 ground truth**。

额外构造同一阳性 locus 的反向互补输入，以及保留 array/足够两侧序列但不含 Cas operon 的片段。后者只是人工 orphan 测试，不能写成发现天然 orphan system。

检查：
- 预期 array 被回收，repeat consensus 与参考在合理边界/方向下对应；
- 边界一致的正反链调用得到同一 NR100 repeat；若 caller 改变了 repeat 边界，单独记录并解释，不用强行 circular permutation 消除差异；
- Cas 缺失不导致 seed 被过滤；
- 小型 tandem-repeat / composition-matched 随机对照的输出与 QC 去向被完整报告。

若阳性失败，区分 detector 漏检、边界差异、过严 QC 和 parser 错误。不能只靠重新定义“正确结果”让测试通过。小型 smoke test 不产生“全数据库 sensitivity / specificity”的声明。

### C. Pilot

最多 100 个真实输入，按已有 taxon / source_type 做确定性抽样；尽量涵盖完整和碎片化组装。数据不足时使用现有数量并报告，不为凑样本下载新大库。

记录每个输入的 wall time、峰值内存（可获取时）、总碱基数和各级输出数。以输入大小和 source_type 分层估计全库 CPU-hours，而不是直接将最好的一条记录耗时乘总数。

---

## 13. 验收与最终汇报

pilot 完成后交付代码、配置、测试、结果表及一份 `report.md`，至少回答：

1. 成功 / 真零结果 / 失败的输入分别多少？实际分析了多少条 contig、多少 bp？
2. raw、bona-fide、HQ、reserve 分别多少？每条额外过滤规则损失多少 array 和多少 unique repeat？
3. `score>=0.75` 与 `>=0.90` 得到的种子多样性差多少？
4. HQ NR100 repeat 和操作性 family 数分别多少？来源最多的是哪些 taxon / source_type？
5. 局部 subtype 可关联、无邻近 Cas、低置信、多 subtype 冲突和注释失败分别多少？
6. 每个 repeat 是否能追溯到原始 FASTA、准确区间与 raw 调用？所有 HQ 坐标验证是否通过？
7. 哪些测试已通过、失败或尚未执行？当前结果是否完整？
8. 全库预计计算成本、磁盘需求和当前 blocker 是什么？

**可交付条件：** 所有发布的 HQ arrays 通过序列/坐标校验；来源和分组关系无丢失；失败任务没有混为阴性；所有规则和版本可复现。没有原始 reads 时明确保留 assembly validation 的局限。

达到条件后停止并汇报。`allow_full_run: false` 时不得自行开启全集群扫描，也不得顺手开始 Step 2。

---

## 14. 已核对的上游依据

访问核对日期：2026-10-13。下面是公开论文、官方 README 或源码；执行时还需记录**实际安装 commit**。项目自定的 QC、聚类与邻接规则不得归称为这些工具的原始默认功能。

- **[R1]** Mitrofanov et al., *CRISPRidentify: identification of CRISPR arrays using machine learning approach*, NAR，DOI `10.1093/nar/gkaa1158`。用于类别/分数含义，不用于证明本项目自定的 HQ 阈值已达到某个精度。
- **[R2]** CRISPRidentify CLI：输入参数、JSON 路径、CPU 与布尔参数定义。  
  `https://github.com/BackofenLab/CRISPRidentify/blob/master/CRISPRidentify.py`
- **[R3]** CRISPRidentify 环境：  
  `https://github.com/BackofenLab/CRISPRidentify/blob/master/environment.yml`
- **[R4]** CRISPRidentify 输出及内部坐标实现：  
  `https://github.com/BackofenLab/CRISPRidentify/blob/master/components/components_output_maker.py`  
  `https://github.com/BackofenLab/CRISPRidentify/blob/master/components/components_detection_refinement.py`
- **[R5]** CRISPRidentify FASTA/header 处理：  
  `https://github.com/BackofenLab/CRISPRidentify/blob/master/components/components_helpers.py`
- **[R6]** CRISPRCasTyper 官方安装、输出字段及 RepeatTyper 说明：  
  `https://github.com/Russel88/CRISPRCasTyper`
- **[R7]** CRISPRCasTyper CLI：  
  `https://github.com/Russel88/CRISPRCasTyper/blob/master/bin/cctyper`

**本阶段最终产物不是“新系统候选列表”，而是一个可信、保留未知类别、可追溯到原始 array 的 repeat 种子库。**
