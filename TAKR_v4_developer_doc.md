# TAKR v4 开发者文档

**文档版本**: v1.1
**日期**: 2026-04-26
**状态**: 进行中
**角色**: 程序员实施日志

---

## 1. 项目概述

**项目名称**: TAKR v4 事件驱动祖先核型重建
**代码位置**: `/media/40T/wlx/zrg/users/zhangrenang/OrthoIndex/soi/`
**设计文档**: `/media/40T/wlx/zrg/users/zhangrenang/OrthoIndex/TAKR_v4_event_driven_design.md`
**Python版本**: 3.11
**依赖**: OR-Tools >= 9.4, networkx, ete3, pulp, numpy

### 1.1 代码结构

```
soi/
├── AK.py                       # 主入口 (3765行), 保留v3路径
│                                # use_v4=True 时调用 takr_event_driven.py
├── takr_event_driven.py        # [新] 事件驱动重建主流程 (调用 ColoredGraph)
├── takr_colored_graph.py       # [新] ColoredGraph 类 — 颜色标签+环检测+路径覆盖
├── takr_events.py              # 统一事件类型 + 树解析 + TAKREvent (已存在)
├── evolution_simulator_ak.py   # 模拟器, 生成真值数据
├── eval_ak.py                 # [新] 评估模块: 加载→匹配→算指标
├── chromosome_path_cover.py    # CP-SAT路径覆盖 (v3使用)
├── telomere_constraint.py      # 端粒约束管理器
└── validator.py                # 染色体-事件验证器
```

---

## 2. 大计划

| 序号 | 大功能 | 状态 | 完成日期 |
|------|--------|------|----------|
| 1 | 统一事件类型（模拟器、AKR、eval） | 已完成 | 2026-04-25 |
| 2 | 统一树解析（复用AK.py的Newick解析） | 已完成 | 2026-04-25 |
| 3 | 检查模拟器事件输出格式 | 已完成 | 2026-04-25 |
| 4 | 取消单独translocation，添加URT | 已完成 | 2026-04-25 |
| 5 | 更新TAKR_v4设计文档 | 已完成 | 2026-04-26 |
|| 6 | takr_event_driven.py v1 实现（共识祖先+分支级检测） | 已完成 | 2026-04-26 |
|| 7 | 模拟器事件输出改造 (TAKREvent格式) | 已完成 | 2026-04-26 |
|| 8 | AKR事件输出改造 (branch级) | 已完成 | 2026-04-26 |
|| 9 | eval_ak.py 评估模块 | 已完成 | 2026-04-26 |
|| 10 | **ColoredGraph 类实现 (v2 核心算法)** | 规划中 | - |
|| 11 | **ColoredGraph indel/loss/dup 集成** | 规划中 | - |
|| 12 | **ColoredGraph 结构重排检测 + 路径覆盖** | 规划中 | - |
|| 13 | **端到端测试 + v1→v2 F1对比** | 规划中 | - |
|| 14 | **真实数据测试** | 规划中 | - |

---

## 3. 小计划（当前迭代）

### 迭代 10.0: ColoredGraph 类实现

**新文件**: `soi/takr_colored_graph.py` (~250行)

**目标**: 实现 `ColoredGraph` 类，封装无向图 + 颜色标签 + 基础操作。

#### 接口定义

```python
class ColoredGraph:
    def __init__(self, hog_level: str):
        self.graph = nx.Graph()   # edge attr: 'colors' = set of (child_id, chrom_idx)
        self.hog_level = hog_level
        self.events: List[TAKREvent] = []

    # 构建
    def add_edge(self, h1, h2, child_id, chrom_idx):
        # 如果边已存在，只添加颜色；否则创建新边
    def add_child(self, source_id, child_graph):
        # 遍历 child_graph 所有染色体，对连续 HOG 调用 add_edge

    # 查询
    def get_colors(self, h1, h2) -> Set[Tuple[str, int]]:
    def shared_edges(self) -> List:
    def unique_edges(self) -> List:

    # 边管理
    def remove_edge_color(self, h1, h2, color):
        # 移除一个颜色，如果无剩余颜色则移除边
    def remove_edge_colors(self, h1, h2, colors_to_remove):
```

#### 实现细节

存储格式:
```python
# 每条边存储颜色集
G.add_edge(h1, h2, colors={(child_A, 0), (child_B, 1)})

# 查询颜色
G[h1][h2]['colors']  # Set[Tuple[str, int]]

# 移除一个颜色
G[h1][h2]['colors'].discard((child_A, 0))
if not G[h1][h2]['colors']:
    G.remove_edge(h1, h2)
```

#### 验收标准

- [ ] `add_edge()` / `add_child()` 正确添加颜色标签
- [ ] `shared_edges()` 返回多颜色边，`unique_edges()` 返回单颜色边
- [ ] `remove_edge_color()` 正确移除颜色，自动清理空边
- [ ] 导入无报错: `from soi.takr_colored_graph import ColoredGraph`

---

### 迭代 10.1: ColoredGraph indel/loss/dup + 结构重排集成

**目标**: 在 `ColoredGraph` 中添加事件检测 + 冲突移除 + 路径覆盖方法。

#### 新增方法

```python
class ColoredGraph:
    # ---------- Indel/Loss ----------
    def find_spanning_edges(self) -> List[Tuple]:
        """找跨越边: unique_edge (a,b) 两端在另一孩子中连通且含中间节点。
        
        对每条 unique_edges:
          用 BFS 检查 a 和 b 之间是否有其他孩子中的路径
          如果连通且路径长度 > 1 → 跨越边
        """
        
    def resolve_indels(self, outgroups=None):
        """删除跨越边，外类群投票确定 gain/loss。
        
        对每条跨越边:
          外类群有 a-b-c-d-e? → 移除跨越边 (loss)
          外类群有 a-e?       → 移除路径边 (gain)
          无外类群 → 简单多数投票
        """

    # ---------- 环检测 + 分类 ----------
    def find_cycles(self) -> List[List]:
        """nx.cycle_basis 包装。"""
        
    def classify_cycle(self, cycle) -> Tuple[Optional[str], List]:
        """分析环颜色模式 → 事件类型。
        
        判断逻辑 (按优先级):
        a. 跨越边模式？ size>4 + 一条边颜色唯一 → indel (已在上一阶段解决)
        b. 环长度=4 + 颜色交替?
           - 端粒参与 → telomere_inversion
           - 否则 → internal_inversion
        c. 颜色交替 + 跨染色体 → RT/URT (端粒参与→URT)
        d. 非交替 + 两端端粒 → EEJ
        e. 非交替 + 内部段为完整染色体 → NCF
        f. 非交替 + 一个邻接可断成两端粒 → fission
        """

    # ---------- 结构重排 ----------
    def resolve_structural_events(self):
        """循环: find_cycles → classify_cycle → 移除冲突边"""
        while True:
            cycles = self.find_cycles()
            if not cycles:
                break
            for cycle in cycles:
                etype, conflict_edges = self.classify_cycle(cycle)
                if etype and conflict_edges:
                    for (h1, h2) in conflict_edges:
                        self.remove_edge_color(h1, h2, ...)
                    self.events.append(...)

    # ---------- 路径覆盖 ----------
    def path_cover(self, telomere_nodes: Set) -> List[List]:
        """端粒约束路径覆盖。
        
        所有节点度 ≤ 2:
          度=0 → 孤立HOG (异常)
          度=1 → 端粒 (启点)
          度=2 → 内部节点
        
        从每个端粒出发走到另一端粒:
          1. 路径上所有节点标记为已访问
          2. 路径长度≥2 → 有效染色体
        
        Returns: [[hog1, hog2, ...], ...]
        """

    # ---------- 转换 ----------
    def to_ancestral_graph(self) -> 'AncestralAdjacencyGraph':
        """path_cover() → 方向化 → AncestralAdjacencyGraph"""

    # ---------- 完整流水线 ----------
    def resolve_all_events(self, outgroups=None) -> List[TAKREvent]:
        """1. duplication resolve → 2. indel/loss → 3. 结构重排 → 4. 路径覆盖"""
```

#### 事件分类核心逻辑

`classify_cycle()` 的核心是**颜色模式分析**:

```python
def classify_cycle(self, cycle):
    nodes = cycle
    edges = [(nodes[i], nodes[(i+1)%len(nodes)]) for i in range(len(nodes))]
    edge_colors = [self.get_colors(u, v) for u, v in edges]
    
    # 特征提取
    n_unique = sum(1 for c in edge_colors if len(c) == 1)    # 单色边数
    n_shared = sum(1 for c in edge_colors if len(c) > 1)     # 共享边数
    
    # 找"冲突边": 单色边中最少的那个颜色的边
    color_counts = defaultdict(int)
    for u, v in edges:
        for col in self.get_colors(u, v):
            color_counts[col] += 1
    
    # 最少出现的颜色 → 该孩子的边是衍生的
    # 移除这些边 = 解决冲突
    conflict_edges = []  # 需要移除的边
    
    # 判断逻辑见设计文档 §4.2
    ...
```

#### 验收标准

- [ ] `find_spanning_edges()` 正确找出 indel/loss 跨越边
- [ ] `classify_cycle()` 正确分类 inversion/RT/EEJ/NCF/fission
- [ ] `path_cover()` 端粒约束正确，每条路径2个端粒
- [ ] `to_ancestral_graph()` 输出格式与 AncestralAdjacencyGraph 兼容

---

### 迭代 10.2: 端到端测试 + 调优

**目标**: 整个流水线跑通，对比 v1 基线 F1。

#### 步骤

1. **单元测试** — ColoredGraph 各方法独立验证

```python
# 测试添加颜色
def test_add_edge_color():
    G = ColoredGraph("N1")
    G.add_edge('A', 'B', 'child1', 0)
    assert G.get_colors('A', 'B') == {('child1', 0)}
    
    G.add_edge('A', 'B', 'child2', 1)
    assert len(G.get_colors('A', 'B')) == 2

# 测试环分类: inversion (颜色交替)
def test_classify_inversion():
    G = ColoredGraph("N1")
    G.add_edge('A', 'B', 'child1', 0)
    G.add_edge('B', 'D', 'child2', 0)
    G.add_edge('D', 'C', 'child1', 0)
    G.add_edge('C', 'A', 'child2', 0)
    # 环 A-B-D-C-A: 颜色 A1,B1,A1,B1 → 交替
    etype, _ = G.classify_cycle(['A', 'B', 'D', 'C'])
    assert etype == 'inversion'

# 测试路径覆盖
def test_path_cover():
    G = ColoredGraph("N1")
    G.graph.add_edge('T1', 'A', colors={('c1', 0)})
    G.graph.add_edge('A', 'B', colors={('c1', 0)})
    G.graph.add_edge('B', 'T2', colors={('c1', 0)})
    paths = G.path_cover({'T1', 'T2', 'T3', 'T4'})
    assert len(paths) == 1
    assert paths[0] == ['T1', 'A', 'B', 'T2']
```

2. **端到端测试** — 流水线完整运行

```bash
cd /media/40T/wlx/zrg/users/zhangrenang/OrthoIndex
export PYTHONPATH="/media/40T/wlx/zrg/users/zhangrenang/OrthoIndex:$PYTHONPATH"

# 生成模拟数据
python3.11 -m soi.evolution_simulator_ak \
  --num-species 4 --num-chroms 6 \
  --inv-rate 5.0 --rt-rate 3.0 --ncf-rate 1.5 --eej-rate 2.0 --fission-rate 0.02 \
  --seed 42 -o tests/sim_data/

# 运行 v2 重建 (通过 AK.py use_v4=True)
python3.11 -c "
from soi.AK import AKR
akr = AKR(
    ogfile='tests/sim_data/ortholog_groups.txt',
    gfffile='tests/sim_data/all_species_gene.gff',
    sptreefile='tests/sim_data/species_tree.nwk',
    outpre='tests/akr_v2/AKR',
    use_v4=True, use_v3=False,
    min_genes=0, timeout=600,
)
akr.run()
print('Done, events:', len(akr.events))
"

# 评估
python3.11 -c "
from soi.eval_ak import load_events, evaluate_branches, print_summary, generate_report
from soi.takr_events import parse_tree
tree, parent_of, children_of = parse_tree('tests/sim_data/species_tree.nwk')
truth = load_events('tests/sim_data/events.tsv', parent_of=parent_of)
detected = load_events('tests/akr_v2/AKR.events.tsv', parent_of=parent_of)
results, global_m = evaluate_branches(truth, detected)
print_summary(results, global_m)
generate_report(results, global_m, 'tests/akr_v2/eval_report.tsv')
"
```

#### 验收标准

- [ ] 单元测试全通过
- [ ] 端到端流水线无异常
- [ ] 染色体数与 truth 一致（N0=8, N1=6 等）
- [ ] Micro F1 > v1 基线 (提升≥0.1)
- [ ] Large-scale F1 > 0.50

---

### 迭代 10.3: AK.py 集成 — reconstruction_algorithm 参数

**目标**: 完善 AKR 入口，支持参数化选择 v2 (ColoredGraph) 或 v1 (当前实现)。

**变更**: `soi/AK.py`

```python
class AKR:
    def __init__(self, reconstruction_algorithm='v3', ...):
        # 'v3' = CP-SAT (原始)
        # 'v4' = v1 event-driven (当前)
        # 'v4_colored' = v2 ColoredGraph (新)
        self.reconstruction_algorithm = reconstruction_algorithm
```

`run()` 方法路由:
- `'v3'` → 原有 CP-SAT 路径
- `'v4'` → 现有 `reconstruct_event_driven()` (v1 共识祖先)
- `'v4_colored'` → 新 `reconstruct_event_driven_v2()` (ColoredGraph)

#### 验收标准

- [ ] `reconstruction_algorithm='v4_colored'` 可执行
- [ ] v3 vs v4 vs v4_colored 三个路径不冲突
- [ ] v4_colored 事件输出格式与 eval_ak.py 兼容

---

### 迭代 10.4: WGD 节点验证

**目标**: 验证 ColoredGraph 在 WGD 节点上的 pre→post 检测。

**关键**: WGD 节点已由 `_collapse_wgd` 拆分为 pre-WGD 图和 post-WGD 图。
ColoredGraph 只需对 pre-WGD 图执行标准流程，不需要特殊 WGD 逻辑。

**变更**: `soi/takr_event_driven.py` 中 `reconstruct_event_driven()` 对 WGD 节点：
1. pre-WGD 图已经通过 `_collapse_wgd` 获得（ploidy × 染色体数）
2. 标准 ColoredGraph 流程 → 祖先图
3. 对比 pre vs post → WGD + fractionation + post-WGD 重排

#### 验收标准

- [ ] WGD 节点 pre→post 事件可检测
- [ ] 染色体数: pre × ploidy + fission - ncf - eej = post
- [ ] pre-WGD 图正确拆分亚基因组

---

## 4. Changelog

### 2026-04-26 v3.2 v4_colored 集成到 AKR 主入口

- [代码] AKR.__init__: `reconstruction_algorithm` 支持 'v4_colored'
- [代码] AKR.run(): 路由到 reconstruct_event_driven_v2
- [代码] 跳过 opt/topdown 阶段（v4_colored 内部处理）
- [测试] `AKR(reconstruction_algorithm='v4_colored').run()` 验证通过
  - N0: 11 chroms, 10 events
  - N1: 9 chroms, 6 events
  - N2: 4 chroms, 3 events

### 2026-04-26 v3.1 ColoredGraph 类实现 + v2入口

- [代码] 新增 `soi/takr_colored_graph.py` — ColoredGraph 完整实现（~600行）
  - 颜色标签管理: add_edge/add_child/shared_edges/unique_edges/remove_edge_color
  - indel/loss: find_spanning_edges + resolve_indels（含外类群投票）
  - 结构重排: find_cycles + classify_cycle（inversion/RT/URT）
  - 路径覆盖: path_cover（端粒约束行走）
  - 转换: to_ancestral_graph() → AncestralAdjacencyGraph
  - 一键执行: resolve_all_events() — 完整 pipeline
  - 单元测试通过（add_edge/get_colors/shared_edges/unique_edges/remove_color/cycle_classification/path_cover）
- [代码] `soi/takr_event_driven.py` 添加 reconstruct_event_driven_v2() 入口

### 2026-04-26 v3.0 ColoredGraph 类设计 + 文档重构

- [设计] 新增 ColoredGraph 类接口设计（§4.9）— 封装颜色标签管理、跨越边检测、环分类、路径覆盖
- [设计] 更新设计文档版本 v4.2→v4.3
- [设计] 设计文档新增 §4.9 ColoredGraph 类设计
- [文档] 重写开发者文档 §3 小计划 — 从旧的 detect→reverse 架构改为 ColoredGraph 架构
- [文档] 重写开发者文档 §7 下一步小计划 — ColoredGraph 实现优先
- [文档] 更新大计划表 (11→14项)
- [文档] 更新已知问题表 (标记已完成项)
- [文档] 移除旧的迭代 6.0-9 内容（已废弃的 detect_events_small/reverse_small 等）
- [文档] 添加 ColoredGraph 代码结构 (soi/takr_colored_graph.py)

### 2026-04-26 v2.1 WGD分支处理完成

- [完成] 迭代 9: WGD 分支处理
  - _handle_wgd_node(): collapse post→pre, 检测pre→post事件
  - 虚拟分支 N2_preWGD-N2 (WGD + fractionation + structural events)
  - 父节点重建时自动使用 pre-WGD 图
  - 事件导出包含 pre_wgd_graphs

- [完成] 迭代 6.0: 实现 takr_event_driven.py
  - 事件驱动重建主入口 reconstruct_event_driven()
  - 共识祖先图构建 _build_consensus_ancestor()
  - 分支级事件检测 _detect_branch_events()
    - inversion (internal_inversion / telomere_inversion)
    - unidirectional translocation
    - EEJ, fission, NCF, RT/URT
  - AK.py 接入: use_v4 标志, 路由到 event-driven 路径
  - 跳过 v3 的 CP-SAT/优化/topdown 检测 (v4 内部处理)

- [完成] 迭代 6.1: takr_events.py 同步 — CSI 移出 CHROM_REARRANGEMENTS, 添加 URT alias
- [完成] 迭代 7: 模拟器 events.tsv 输出改为统一 branch 级格式
  - 列: branch, event_type, genes, chroms, desc, support
  - 过滤 rearrangements 汇总事件, 类型名 canonicalize 化
- [完成] 迭代 7: AKR.events.tsv 输出改为统一 branch 级格式
  - 添加 branch 列 (parent-child), 保留 child_source 列
  - 修复 leaf graph 事件未导出的 bug
- [完成] 迭代 8: 实现 eval_ak.py
  - load_events(): 支持 3 种格式 (unified/old AKR/old simulator)
  - match_events_branch(): 按 event_type + genes Jaccard 匹配
  - evaluate_branches(): 按分支×事件类型输出 TP/FP/FN/F1
- 添加 takr_event_driven.py 模块接口设计 (迭代 6.0)
- 添加事件输出改造方案 — 统一分支级格式 (迭代 7)
- 添加 eval_ak.py 评估模块设计 (迭代 8)
- 添加 WGD 分支处理方案 (迭代 9)
- 更新大计划表 (8→11项)
- 重写文档 §5 测试章节

### 2026-04-26 v1.0 初始版本

- 初始开发者文档，包含大计划、小计划、技术细节

---

## 5. 测试

### 5.1 测试总览

TAKR v4 采用三阶段闭环测试:

```
Phase 1: 生成模拟数据
  evolution_simulator_ak.py
  → events.tsv (truth), ortholog_groups.txt, species_tree.nwk, all_species_gene.gff

Phase 2: 运行祖先重建
  AK.py (--use-v4) 或 takr_event_driven.py
  → AKR.events.tsv (detected), AKR.*.anc.gfa, AKR.tsv

Phase 3: 评估
  eval_ak.py load_events() → evaluate_branches() → generate_report()
  → eval_report.tsv (按分支×事件类型的 TP/FP/FN/Precision/Recall/F1)
```

### 5.2 Phase 1: 生成模拟数据

#### 基本模拟

```bash
cd /media/40T/wlx/zrg/users/zhangrenang/OrthoIndex
export PYTHONPATH="/media/40T/wlx/zrg/users/zhangrenang/OrthoIndex:$PYTHONPATH"

python3.11 -m soi.evolution_simulator_ak \
  --num-species 6 \
  --num-chroms 7 \
  --wgd-rate 0.3 \
  --inv-rate 10.0 \
  --rt-rate 4.0 \
  --ncf-rate 2.0 \
  --eej-rate 2.0 \
  --fission-rate 0.02 \
  --seed 42 \
  -o tests/sim_data/
```

输出文件:
- `species_tree.nwk` — 物种树（含 `[p=N]` 多倍体标注）
- `ortholog_groups.txt` — 直系同源群（HOG输入）
- `all_species_gene.gff` — 所有物种的基因GFF
- `ancestors_gene.gff` — 祖先基因GFF（参考）
- `events.tsv` — **真值事件**（统一branch级格式）
- `ancestors_karyotypes.txt` — 祖先核型参考

#### 带真实物种树的模拟

```bash
python3.11 -m soi.evolution_simulator_ak \
  -t /media/nfs2/wlx/zrg/proj/self/Aethionema_saxatile/orthologs/Species_Tree \
  --num-chroms 8 \
  --wgd-rate 0.2 \
  -o tests/sim_data_real/
```

#### 小规模快速测试（调试用）

```bash
python3.11 -m soi.evolution_simulator_ak \
  --num-species 4 \
  --num-chroms 5 \
  --min-genes 50 \
  --max-genes 200 \
  --inv-rate 3.0 \
  --rt-rate 1.0 \
  --ncf-rate 0.5 \
  --eej-rate 0.5 \
  --fission-rate 0.01 \
  --seed 42 \
  -o tests/sim_data_small/
```

### 5.3 Phase 2: 运行祖先重建

#### 调用 AKR (use_v4)

```bash
cd /media/40T/wlx/zrg/users/zhangrenang/OrthoIndex
export PYTHONPATH="/media/40T/wlx/zrg/users/zhangrenang/OrthoIndex:$PYTHONPATH"

python3.11 -c "
from soi.AK import AKR

akr = AKR(
    ogfile='tests/sim_data/ortholog_groups.txt',
    gfffile='tests/sim_data/all_species_gene.gff',
    sptreefile='tests/sim_data/species_tree.nwk',
    outpre='tests/akr_test_v4/AKR',
    use_v4=True,       # ← 使用事件驱动重建
    use_v3=False,
    min_genes=0,
    timeout=600,
    node_timeout=120,
)
akr.run()
print('Done, events:', len(akr.events))
"
```

#### 调用 takr_event_driven.py (未来)

```bash
python3.11 -c "
from soi.takr_event_driven import reconstruct_event_driven
from soi.AK import AKR

akr = AKR(
    ogfile='tests/sim_data/ortholog_groups.txt',
    gfffile='tests/sim_data/all_species_gene.gff',
    sptreefile='tests/sim_data/species_tree.nwk',
    outpre='tests/akr_v4_test/AKR',
    use_v4=False,
)
akr._build_hogs()
akr._build_leaf_graphs()

anc_graphs = reconstruct_event_driven(akr)
print('Ancestors:', list(anc_graphs.keys()))
"
```

### 5.4 Phase 3: 评估

```bash
cd /media/40T/wlx/zrg/users/zhangrenang/OrthoIndex
export PYTHONPATH="/media/40T/wlx/zrg/users/zhangrenang/OrthoIndex:$PYTHONPATH"

python3.11 -c "
from soi.eval_ak import load_events, evaluate_branches, print_summary, generate_report
from soi.takr_events import parse_tree

tree, parent_of, children_of = parse_tree('tests/sim_data/species_tree.nwk')

# 旧格式: 传 parent_of 供分支推断
truth = load_events('tests/sim_data/events.tsv', parent_of=parent_of)
detected = load_events('tests/akr_test_v4/AKR.events.tsv', parent_of=parent_of)

results, global_m = evaluate_branches(truth, detected)
print_summary(results, global_m)
generate_report(results, global_m, 'tests/akr_test_v4/eval_report.tsv')
"
```

预期输出:
```
=== 评估报告 ===
Micro F1: 0.723  (TP=128 FP=47 FN=52)
  Precision: 0.731  Recall: 0.711

  N1-N2: T=15 D=12 F1=0.733 (P=0.733 R=0.733)
    inversion               T=8  D=6  F1=0.750
    reciprocal_translocation T=2  D=2  F1=1.000
    ncf                     T=3  D=2  F1=0.667
    ...
  N2-Sp_3: T=8 D=6 F1=0.667 ...
```

### 5.5 真实数据测试

```bash
cd /media/40T/wlx/zrg/users/zhangrenang/OrthoIndex

python3.11 -c "
from soi.AK import AKR
akr = AKR(
    ogfile='tests/brassicaceae/ortholog_groups.txt',
    gfffile='tests/brassicaceae/all_species_gene.gff',
    sptreefile='/media/nfs2/wlx/zrg/proj/self/Aethionema_saxatile/orthologs/Species_Tree',
    outpre='tests/brassicaceae/AKR',
    use_v4=True,
    use_v3=False,
    min_genes=0,
)
akr.run()

# 验证: 祖先染色体数应与文献一致
# Crambe_hispanica + Thlaspi_praecox ancestor ≈ 7
# + Arabis_montbretiana + Arabis_sagittata ancestor ≈ 8
for name, aag in akr.anc_graphs.items():
    chrom_count = len(list(aag.chromosomes()))
    print(f'  {name}: {chrom_count} chroms')
"
```

### 5.6 单元测试 (语法/导入检查)

```bash
cd /media/40T/wlx/zrg/users/zhangrenang/OrthoIndex
export PYTHONPATH="/media/40T/wlx/zrg/users/zhangrenang/OrthoIndex:$PYTHONPATH"

python3.11 -c "from soi.AK import RearrangementEvent, AKR; print('AK import OK')"
python3.11 -c "from soi.takr_events import TAKREvent, canonicalize_event_type, branch_id; print('takr_events import OK')"
python3.11 -c "from soi.takr_event_driven import reconstruct_event_driven; print('takr_event_driven import OK')"
python3.11 -c "from soi.eval_ak import load_events, evaluate_branches; print('eval_ak import OK')"

python3.11 -c "
from soi.takr_events import EVENT_TYPES, canonicalize_event_type
assert 'unbalanced_reciprocal_translocation' in EVENT_TYPES
assert canonicalize_event_type('RT') == 'reciprocal_translocation'
assert canonicalize_event_type('translocation') == 'reciprocal_translocation'
print(f'Event type consistency OK ({len(EVENT_TYPES)} types)')
"
```

### 5.7 测试验收标准

| 测试项 | 命令/方法 | 预期结果 |
|--------|-----------|----------|
| 语法检查 | 各模块 import | 无报错 |
| 事件一致性 | canonicalize_event_type 全类型 | 全部可解析 |
| 模拟数据生成 | evolution_simulator_ak + seed | 可复现输出 |
| 模拟输出格式 | events.tsv 含 branch 列 | 列正确, 类型 canonical |
| v4重建 | AKR(use_v4=True).run() | 完成无异常 |
| AKR输出格式 | AKR.events.tsv 含 branch 列 | 列正确, 类型 canonical |
| 评估运行 | eval_ak.py 加载 + 匹配 | 无异常 |
| 事件 F1 | 模拟数据评估 | Micro F1 > 0.85 (目标) |
| WGD分支 | WGD节点的 pre→post 评估 | 分支正确分离 |
| 真实数据 | Brassicaceae | 染色体数与文献一致 |

---

## 6. 已知问题

| 问题 | 状态 | 备注 |
|------|------|------|
| v1 事件检测大规模重排 F1 ≈ 0.025 | 待解决 | v2 ColoredGraph 环检测+冲突移除方案 |
| 共识图碎片化 → 染色体数偏高 (N0: 8→11) | 待解决 | v2 保留所有边+颜色，不移除unique边 |
| 无冲突边检测 + 无路径覆盖 | 待解决 | v2 ColoredGraph 已设计 |
| fractionation 未触发 (无WGD) | 已确认 | 测试数据无WGD时所有HOG拷贝数=1 |
| 模拟器事件输出已统一为 branch 级格式 | 已完成 | 迭代 7 |
| AKR事件输出已添加 branch 列 | 已完成 | 迭代 7 |
| eval_ak.py 评估模块已实现 | 已完成 | 迭代 8 |

---

## 7. 下一步小计划

### 已完成的工作:

- ✅ 迭代 6.0: takr_event_driven.py v1 核心框架 (共识祖先 + 分支级事件检测)
- ✅ 迭代 6.1: takr_events.py 同步 (CSI/URT)
- ✅ 迭代 7: 模拟器 + AKR 事件输出改造 (统一 branch 级格式)
- ✅ 迭代 8: eval_ak.py 评估模块
- ✅ 迭代 9: WGD 分支处理 (模拟器 + AKR)
- ✅ 设计文档 v4.3 ColoredGraph 类接口设计
- ✅ 开发者文档重写 (反映 ColoredGraph 新架构)

### 待完成 (按优先级):

1. **[迭代 10.0] ColoredGraph 类实现** — soi/takr_colored_graph.py
   - `add_edge`, `add_child`, `get_colors`, `shared_edges`, `unique_edges`
   - `remove_edge_color`
   - 验收: 导入正常, 颜色管理正确

2. **[迭代 10.1] ColoredGraph 事件检测 + 路径覆盖**
   - `find_spanning_edges()` + `resolve_indels()` — indel/loss
   - `find_cycles()` + `classify_cycle()` — 结构重排分类
   - `resolve_structural_events()` — 冲突移除
   - `path_cover()` — 端粒约束路径覆盖
   - `to_ancestral_graph()` — 转换
   - `resolve_all_events()` — 完整流水线
   - 验收: 单元测试通过

3. **[迭代 10.2] 端到端测试 + F1对比**
   - 小规模模拟数据 (4种, 6chr, 无WGD)
   - 对比 v1 基线 F1
   - 验证染色体数准确率
   - 验收: 流水线无异常, 染色体数匹配, F1 提升

4. **[迭代 10.3] AK.py reconstruction_algorithm 参数**
   - 支持 'v3' / 'v4' / 'v4_colored'
   - 验收: 三路径不冲突

5. **[迭代 10.4] WGD 节点验证**
   - 验收: pre×ploidy = post 染色体数

---

## 5. Changelog

| 日期 | 版本 | 变更 |
|------|------|------|
| 2026-04-26 | v1.1 | **Outgroup voting for all event types**: `resolve_structural_events()` and `resolve_bridge_events()` now use outgroup adjacency at parent HOG level for polarity. Bridge: outgroup-has-adjacency→fission, else→EEJ/NCF, no-outgroup→bridge_unclassified. Structural: per-edge outgroup count replaces min-count heuristic. `_remove_rare_color_edge()` extracted. `takr_event_driven.py`: `outgroup_leaves_map` precomputation, maps sibling leaves to parent HOG level. |

---

**文档版本**: v1.1
**更新日期**: 2026-04-26
**下次更新**: 迭代 10.3 完成后
