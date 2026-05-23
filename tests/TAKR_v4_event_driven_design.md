# TAKR v4 事件驱动祖先重建设计文档

**文档版本**: v4.3  
**日期**: 2026-04-26  
**状态**: 设计中（待程序员实现）  
**角色**: 产品经理 — 功能规划、生物学背景、验收标准

---

## 1. 项目背景

TAKR（Telomere-Centric Ancestral Karyotype Reconstruction）是基于端粒中心的祖先核型重建工具。v4版本引入**事件驱动重建**范式：先检测重排事件，再逆向恢复祖先状态。

**核心挑战**: 现有v3算法采用"先重建再检测"，导致大量假阳性（unidir_trans FP=94 vs truth=1）和假阴性（fractionation FN=24 vs 0 detected）。v4通过事件优先级和端粒中心模型解决此问题。

---

## 2. 核心原则

### 2.1 事件驱动重建

不是"先重建再检测"，而是：
1. **建初始共识图** — 子节点共享的邻加边加入
2. **检测冲突** — 导致图分支/环的边 = 事件
3. **移除冲突边** — 事件已被"提取"，从共识图中移除冲突边
4. **路径覆盖（端粒约束）** — 干净的线性图重建染色体
5. **外类群投票** — 确定事件极性，分配到正确分支

### 2.2 端粒中心模型（Telomere-Centric）

所有重排检测围绕**端粒**展开：
- 端粒HOG是"锚点"
- EEJ/Fission/NCF 都涉及端粒状态变化
- 端粒倒位必须先处理（改变端粒位置，影响后续检测）
- 路径覆盖阶段必须以端粒为起点，每条染色体恰好两个端粒

### 2.3 事件处理顺序

**原则**: 简单事件先解（简化图），复杂事件后解（在简化图上操作）。每步完成后检查 postcondition，不满足则报错停止。

| Phase | 操作 | Postcondition |
|-------|------|---------------|
| 1 | 每个孩子内部: duplication resolve | 每个孩子图线性 (deg≤2, 无并行边) |
| 2 | 合图 + 单基因 indel/loss/gain | 无跨越边, 无环, 无分支 |
| 3 | 共线性块压缩 | 所有 HOG 已分配到块 |
| 4a | seg_deletion / seg_insertion | 无大段不对称 |
| 4b | unidir_trans | 无单向转移模式 |
| 4c | telomere_inversion (优先!) | 端粒位置一致 |
| 4d | internal_inversion | 块级图无环 |
| 4e | RT / URT | 无跨染色体环 |
| 4f | EEJ / NCF / fission (最后) | 无桥接边 |
| 5 | 端粒约束路径覆盖 | 所有 HOG 覆盖, 每条染色体 2 个端粒 |
| 6 | WGD | — |

### 2.4 无长度限制

任何重排事件都可能涉及**任意数量基因**，不预设长度上限：
- inversion：可能涉及整个染色体臂
- unidir_trans：可能涉及大片段
- NCF：donor染色体可能很长
- 仅设最小阈值（≥3 HOG）避免噪音

---

## 3. 事件类型定义与图结构

线性染色体 = 端粒到端粒的一条路径，边 = HOG 间的邻接。
每条重排可用"断边"和"新边"描述。断边 = 在子代中消失的祖先邻接，新边 = 子代中出现的新邻接。
检测到事件后，断边从共识图中移除，新边不加入共识。

### 3.1 基因级事件（染色体数不变）

| 事件 | 祖先图 | 子代图 | 断边 | 新边 | 图特征 |
|------|--------|--------|------|------|--------|
| gene_loss | T1-A-B-C-D-T2 | T1-A-C-D-T2 | A-B, B-C | A-C | 一个HOG消失，跨越连接 |
| gene_gain | T1-A-B-C-T2 | T1-A-X-B-C-T2 | A-B | A-X, X-B | 一个HOG插入，产生新边 |
| tandem_dup | T1-A-B-C-T2 | T1-A-B-B'-C-T2 | — | B-B' | 相邻重复，并行边 |
| proximal_dup | T1-A-B-C-D-T2 | T1-A-B'-B-C-D-T2 | — | B'-B 或 B-B' | 同染色体附近重复 |
| dispersed_dup | T1-A-B-C-T2, T3-D-E-F-T4 | T1-A-B-C-T2, T3-D-B'-E-F-T4 | D-E | D-B', B'-E | 重复在另一染色体 |
| seg_duplication | 单条或整段 | 整段在另一位置出现 | — | 多条新边 | 整块重复 |
| seg_deletion | T1-A-B-C-D-E-F-T2 | T1-A-E-F-T2 | A-B, D-E | A-E | 整块丢失，跨越连接 |

### 3.2 小规模结构重排（染色体数不变）

| 事件 | 祖先图 | 子代图 | 断边 | 新边 | 图特征 |
|------|--------|--------|------|------|--------|
| internal_inversion | T1-A-B-C-D-E-F-T2 | T1-A-E-D-C-B-F-T2 | A-B, E-F | A-E, B-F | 内部block倒转，两个断点都在内部 |
| telomere_inversion | T1-A-B-C-D-E-F-T2 | T1-D-C-B-A-E-F-T2 | T1-A, D-E | T1-D, A-E | 靠近端粒的block倒转，一个断点在端粒(T1-A) |
| unidir_trans | T1-A-B-C-D-T2, T3-E-F-G-H-T4 | T1-A-T2, T3-E-B-C-D-F-G-H-T4 | A-B, D-T2, E-F | A-T2, E-B, D-F | 单向移动，供体染色体丢失 |

**邻接图特征**:
- **internal_inversion**: 4-cycle A-B-F-E-A，颜色交替，所有边同 chrom_idx
- **telomere_inversion**: 4-cycle T1-D-E-A-T1，颜色交替，端粒节点参与环
- **unidir_trans**: 分叉+汇合，供体染色体消失，受体染色体增长

### 3.3 大尺度重排（染色体数改变）

| 事件 | 祖先图 | 子代图 | 断边 | 新边 | 染色体数 |
|------|--------|--------|------|------|---------|
| EEJ | T1-A-B-T2, T3-C-D-T4 | T1-A-B-C-D-T4 | B-T2, T3-C | B-C | -1 |
| NCF | T1-A-B-C-D-T2, T3-E-F-G-T4 | T1-A-B-E-F-G-C-D-T2 | B-C, T3-E, G-T4 | B-E, G-C | -1 |
| fission | T1-A-B-C-D-T2 | T1-A-B-T3, T4-C-D-T2 | B-C | B-T3, T4-C | +1 |
| RT | T1-A-B-C-D-T2, T3-E-F-G-H-T4 | T1-A-B-G-H-T4, T3-E-F-C-D-T2 | B-C, F-G | B-G, F-C | 0 |
| URT | T1-A-B-C-D-T2, T3-E-F-G-H-T4 | T1-A-B-C-H-T4, T3-E-F-G-D-T2 | D-T2, G-H | D-G, H-T2 | 0 |

**邻接图特征（用于检测和区分事件）**:

| 事件 | 邻接图特征 | 区分方法 |
|------|-----------|---------|
| EEJ | 唯一边连接两个共享分量，两个分量都有端粒 | 端粒参与 + 外类群确认 |
| NCF | 两个唯一边(B-E, G-C)连接 donor 段到 recipient，donor 内部邻接是共享的 | donor 段两端各一个唯一边 |
| fission | 唯一边连接两个共享分量，外类群有该邻接 | 外类群确认：祖先有B-C → fission |
| RT | 4-cycle，边来自不同染色体(chrom_idx不一致) | chrom_idx 不一致 → RT |
| URT | 4-cycle，一端在端粒 | 端粒参与 + chrom_idx 不一致 |

**RT vs inversion 区分**:
- inversion: 4-cycle，所有边同 chrom_idx（同一条染色体内部倒位）
- RT: 4-cycle，边来自不同 chrom_idx（跨染色体交换）

**RT vs URT 区别**:
- RT: 两个断点都在染色体内部 (B-C, F-G)
- URT: 一个断点在端粒 (D-T2, G-H)
- 结果：染色体数都不变，但 URT 改变了端粒归属

### 3.4 特殊事件

- **WGD**: 染色体数 ×ploidy。每条祖先染色体产生 ploidy 条拷贝，HOG 树分裂为 ploidy 个子 HOG。
- **Chromothripsis**: 碎裂后随机重组，无固定图结构。暂不处理。

---

## 4. 关键算法思路

### 4.1 图结构基础

染色体邻接图 = **无向图**（`nx.Graph`）。节点 = HOG，边 = 两个 HOG 在染色体上相邻。
有向边只在路径覆盖后追加（编码链向方向），冲突检测阶段用无向图。

每条边记录颜色标签 **(child_id, chromosome_id)**。完整的颜色由两部分组成——来自哪个孩子，以及在该孩子的哪条染色体上。只有当两个 HOG 在同一孩子的**同一条染色体**上连续时，该边才编码"同色连续"关系。

**颜色定义的层级：**
```
标签: (物种, 染色体ID)
      species_A + chr1  → "A1"
      species_A + chr2  → "A2"
      species_B + chr1  → "B1"
```

**环中颜色模式判断：**
- 两个端点同色 → 环两侧的边来自同一孩子的同一染色体 → inversion
- 两个端点异色 → 环两侧的边来自不同孩子（或同一孩子的不同染色体）→ RT/NCF/EEJ
- 环中存在孤立异色边，其余连续同色 → indel/loss

### 4.2 四种基本图冲突

两种孩子路径在图上产生冲突，表现为**环**或**分叉**。所有重排类型都可以归纳为以下四种冲突之一：

**冲突类型A: 相邻路径分叉（indel/loss）**

可涉及连续多个 HOG 的缺失或插入，环长度可为 4+：

```
孩子1: A-B-C-D-E           (含 B,C,D)
孩子2: A-E                 (不含 B,C,D)
共识: A-B, B-C, C-D, D-E, A-E
环: A-B-C-D-E-A  (5个节点)
颜色: A-B(A1)-B-C(A1)-C-D(A1)-D-E(A1)-E-A(B1)
      → 4条A1 + 1条B1，B1跨越边 A-E 跳过连续 HOGs
解决: 判定为 indel/loss，移除跨越边 A-E
      外类群投票: gain vs loss
```

**冲突类型B: 反向路径（inversion）**
```
孩子1: A-B-C-D
孩子2: A-C-B-D  
断边: A-B, C-D (孩子2没有)
新边: A-C, B-D (孩子2有)
共识: A-B, C-D, A-C, B-D
图: A-B, B-D, D-C, C-A → 环 A-B-D-C-A
颜色: A-B(A1)-B-D(B1)-D-C(A1)-C-A(B1) → 交替
解决: 判定为 inversion，移除交叉边 A-C, B-D
```

**冲突类型C: 交叉路径（RT）**
```
孩子1: B-C, F-G  (不同染色体: B-C chr1, F-G chr2)
孩子2: B-G, F-C  (不同染色体: B-G chr1, F-C chr2)
共识: B-C, F-G, B-G, F-C
图: B-C, C-F, F-G, G-B → 环 B-C-F-G-B
颜色: B-C(A1)-C-F(B1)-F-G(A2)-G-B(B1) → A1, B1, A2, B1
      → 交替 + 跨染色体
解决: 判定为 RT，移除一对交叉边
```

**冲突类型D: 分叉+合并（NCF/EEJ）**
```
NCF: 孩子1有 B-E-G-C (来自另一条染色体)，孩子2有 B-C
     B-E, E-G, G-C, B-C → 分叉 B→(C, E), 汇合 G→C
     环 B-E-G-C-B
颜色: B-E(A2)-E-G(A2)-G-C(A2)-C-B(B1) → 3条A2 + 1条B1
解决: 判定为 NCF，移除跨越边 B-C
```

### 4.3 处理流程

#### 第1步：每个孩子内部 — duplication resolve

**在跨物种比较之前，对每个孩子单独处理**

```
对每个孩子的邻接图:
1. tandem_dup: 相邻拷贝合并（B-B' → B-C）
2. dispersed_dup: 标记事件，保留边（真实重排）
3. proximal_dup: 同染色体附近重复，合并
4. seg_duplication: 整段重复，标记事件

输出: 简化后的孩子图，无多余重复边

Postcondition 检查:
  - 每个节点 degree ≤ 2（无分支）
  - 无并行边（同一对节点之间只有一条边）
  - 失败则报错: "Phase 1 postcondition failed: node X has degree N"
```

#### 第2步：两个孩子合图 + 单基因 indel/loss/gain

```
合并两个孩子为彩色图:
  - 每条边记录颜色 (child_id, chrom_idx)

单基因 indel/loss 检测:
  - 找跨越边（两个端点在另一个孩子中通过中间节点连接）
  - 外类群投票确定 gain vs loss
  - 移除跨越边，记录事件
  - 重复直到无跨越边

gene_gain 检测:
  - 标记孤立 HOG（只在一个孩子中存在）

Postcondition 检查:
  - 无跨越边（no_cross_edges）
  - 无环（no_cycles）
  - 无 degree > 2 节点（无分支）
  - 失败则报错: "Phase 2 postcondition failed: found N crossing edges, M cycles"
```

#### 第3步：共线性块压缩

```
从简化的 HOG 级图构建块:
  - 共享边（≥2 孩子共有）构成块
  - 块 = 共享边的连续路径
  - 单例块 = 不在共享边中的 HOG
构建块级图:
  - 节点 = 块
  - 边 = 块间邻接，颜色 = 孩子来源

Postcondition 检查:
  - 每个 HOG 恰好属于一个块
  - 块内节点 degree ≤ 2
  - 失败则报错: "Phase 3 postcondition failed: HOG X not assigned to any block"
```

#### 第4步：块级事件（按复杂度递增）

```
确保图足够简化后再处理复杂事件
每步完成后检查 postcondition，不满足则报错停止

4a. seg_deletion / seg_insertion
    - 一个孩子有大段，另一个没有
    - 外类群投票区分 deletion vs insertion
    - Postcondition: 无大段不对称 (no_large_asymmetric_segments)

4b. unidir_trans
    - 一条染色体的片段转移到另一条
    - 供体染色体丢失
    - Postcondition: 无单向转移模式 (no_unidirectional_transfers)

4c. telomere_inversion（优先！）
    - 一个断点在端粒，改变端粒位置
    - 用端粒信号决定保留哪种颜色的边
    - Postcondition: 端粒位置一致 (telomere_positions_consistent)

4d. internal_inversion
    - 两个断点都在染色体内部
    - 用 chrom_idx 确认同染色体
    - Postcondition: 块级图无环 (no_cycles)

4e. RT / URT
    - 两条染色体交叉交换片段
    - 用 chrom_idx 区分 RT 和 inversion
    - 用端粒参与区分 RT 和 URT
    - Postcondition: 无跨染色体环 (no_cross_chrom_cycles)

4f. EEJ / NCF / fission（桥接事件，最后处理）
    - 块级共享边图中找连通分量
    - 对每条唯一边: 两端在不同分量 → 桥接候选
    - 外类群投票:
      · 外类群有该邻接 → fission
      · 外类群无该邻接 → NCF（小段）或 EEJ（大段+端粒）
    - 移除/保留边，记录事件
    - Postcondition: 无桥接边 (no_bridge_edges)
```

#### 第5步：端粒约束路径覆盖

```
1. 找共识端粒 HOG（≥2 个孩子中都与端粒相邻）
2. 从共识端粒出发，沿邻接行走
3. 走到另一个端粒 → 一条完整染色体
4. 每条染色体恰好 2 个端粒
5. 覆盖所有 HOG

Postcondition 检查:
  - 所有 HOG 被覆盖 (all_hogs_covered)
  - 每条染色体恰好 2 个端粒 (each_chrom_2_telomeres)
  - 失败则报错: "Phase 5 postcondition failed: N HOGs uncovered, M chromosomes without 2 telomeres"
```

#### 第5步：验证

```
染色体数: anc_chrom + fission - ncf - eej = child_chrom
孤立 HOG: 无
端粒数: 2 × 染色体数
```

### 4.4 端粒HOG识别

端粒HOG = 染色体最左端或最右端的HOG。所有大尺度重排检测围绕端粒展开。

### 4.5 外类群投票

外类群**仅用于确定极性**（哪个状态是祖先），不用于检测：

- EEJ/fission: 外类群的端粒状态决定融合/断裂的方向
- RT: 外类群支持哪组交叉模式
- indel: 外类群有该 HOG? → loss；没有? → gain
- duplication: 外类群有多个拷贝? → ancestral dup；没有? → derived dup

### 4.6 染色体数验证

```
普通节点:   anc_chrom + fission - ncf - eej = child_chrom
pre-WGD节点: anc_chrom * ploidy + fission - ncf - eej = child_chrom
```

### 4.7 事件优先级

1. **基因级事件优先**：移除 loss/dup 造成的环，简化图结构
2. **端粒倒位优先**：改变端粒位置，必须先恢复
3. **小规模先于大尺度**：避免小重排被误判为大重排
4. **EEJ/NCF先断开**（增加染色体数），fission后连接

### 4.8 多倍体叶子 collapse

多倍体叶子（如 Sp_2[p=2]）需要先坍塌到 pre-WGD 状态。

方法（不依赖 CP-SAT）：
1. 映射到父节点 HOG 层级
2. 利用 HOG 树追踪拷贝来源：WGD 后每个父 HOG 产生 ploidy 个子 HOG（.hog0, .hog1, ...）
3. 按子 HOG 分配基因到各亚基因组（不依赖染色体位置，不受亚基因组间重排影响）
4. 亚基因组图 → 标准 event-driven 流程
5. 任意 p 均支持

### 4.9 ColoredGraph 类设计

**定位**: `ColoredGraph` 是事件驱动重建核心的数据结构，封装无向图 + 颜色标签 + 所有图操作。

**类接口**:

```python
class ColoredGraph:
    """彩色邻接图 — 每条边记录一组颜色标签 (child_id, chromosome_id)。
    
    一条边可以有多个颜色标签（多个孩子的同一邻接共享该边）。
    颜色集为空时自动移除边。
    """

    def __init__(self, hog_level: str):
        self.graph = nx.Graph()   # edge attr: 'colors' = set of (child_id, chrom_idx)
        self.hog_level = hog_level
        self.events: List[TAKREvent] = []

    # ---------- 构建 ----------

    def add_edge(self, h1, h2, child_id: str, chrom_idx: int):
        """添加一条边，记录颜色 (child_id, chrom_idx)。
        如果 (h1, h2) 已存在，只把颜色加入已有边的 colors 集。
        """
        ...

    def add_child(self, source_id: str, child_graph: AdjacencyGraph):
        """把一个孩子的所有染色体邻接以该孩子的颜色加入。
        遍历 child_graph 的每条染色体，对染色体上的连续 HOG
        逐一调用 add_edge(self, h1, h2, source_id, chrom_idx)。
        """
        ...

    # ---------- 查询 ----------

    def get_colors(self, h1, h2) -> Set[Tuple[str, int]]:
        """返回 (h1, h2) 上的所有颜色。"""
        ...

    def edge_count(self) -> int:
        """当前总边数（每条边count=1，不计颜色数）。"""
        ...

    def shared_edges(self) -> List[Tuple]:
        """返回有多于一个颜色的边 → 祖先共享邻接。"""
        ...

    def unique_edges(self) -> List[Tuple]:
        """返回只有一个颜色的边 → 可能为衍生边。"""
        ...

    def edges_by_color(self, color: Tuple[str, int]) -> List[Tuple]:
        """返回所有带有指定颜色的边。"""
        ...

    # ---------- 环检测 + 事件分类 ----------

    def find_cycles(self) -> List[List]:
        """nx.cycle_basis 包装。返回图中所有基本环的节点列表。"""
        ...

    def classify_cycle(self, cycle: List) -> Tuple[Optional[str], List]:
        """分析环的颜色模式，判断事件类型。
        
        判断逻辑:
        1. 环中是否有跨越边（两端在同一孩子的另一条路径上）→ indel/loss
        2. 环长度=4 且颜色交替 → 检查边是否端粒:
           是 → telomere_inversion  否 → internal_inversion
        3. 颜色交替 + 跨染色体 → RT/URT（检查端粒参与区分）
        4. 非交替 + 跨越边两端均为端粒 → EEJ
        5. 非交替 + 跨越边串联完整染色体 → NCF
        6. 非交替 + 一个邻接断开成两端粒 → fission
        7. 无法分类 → None
        
        Returns:
            (event_type, conflict_edges) 或 (None, [])
        """
        ...

    # ---------- 边移除 ----------

    def remove_edge_color(self, h1, h2, color: Tuple[str, int]):
        """移除边上的一个颜色。如果该边没有其他颜色，移除整条边。"""
        ...

    def remove_color_from_cycle(self, cycle: List, colors_to_keep: Set = None):
        """从环中移除冲突边（保留共享边时指定 colors_to_keep）。"""
        ...

    # ---------- indel 检测 ----------

    def find_spanning_edges(self) -> List[Tuple]:
        """找跨越边: 边 (a, b) 的端点在另一个孩子的连续路径中相邻。
        
        对每条 unique_edges:
          检查两端 a 和 b 是否在另一孩子的图中通过中间节点连通
          如果连通且路径长度 > 1 → 跨越边
        
        Returns: [(h1, h2, child_id, spanned_hogs), ...]
        """
        ...

    # ---------- 路径覆盖 ----------

    def path_cover(self, telomere_nodes: Set) -> List[List]:
        """端粒约束路径覆盖。
        
        从每个端粒节点出发，沿唯一邻接行走:
        1. 走到另一端粒 → 一条染色体
        2. 所有节点度 ≤ 2 (允许度为1的端粒，度为0的孤立HOG抛出异常)
        3. 每条染色体恰好2个端粒
        4. 覆盖所有HOG
        
        Returns: [[hog1, hog2, ...], ...] 每条染色体一条路径
        """
        ...

    # ---------- 转换 ----------

    def to_ancestral_graph(self) -> 'AncestralAdjacencyGraph':
        """转换为 AncestralAdjacencyGraph。
        
        流程:
        1. path_cover() 获取染色体路径
        2. 每条路径构建方向化邻接
        3. 追加有向边编码链向方向
        4. 返回 AncestralAdjacencyGraph(node_id=self.hog_level)
        """
        ...

    # ---------- 一键执行 ----------

    def resolve_all_events(self, outgroups: Dict = None) -> List[TAKREvent]:
        """按优先级解决所有冲突，返回所有检测到的事件。

        Pipeline: 简单→复杂，每步检查 postcondition，失败则报错停止。
        详见 §2.3 表格和 §4.3 详细流程。
        """
        ...
```

**使用方式**（替代当前 `_build_consensus_ancestor` → `_detect_branch_events` 流程）:

```python
def reconstruct_for_node(node_id, mapped_children, child_source_ids):
    G = ColoredGraph(hog_level=node_id)
    
    # Step 1: duplication resolve (每个孩子单独)
    for mc, sid in zip(mapped_children, child_source_ids):
        G.add_child(sid, mc)
        dup_events = G.resolve_duplications(sid)
    
    # Step 2: indel/loss
    G.resolve_indels(outgroups)
    
    # Step 3: 结构重排
    G.resolve_structural_events()
    
    # Step 4: 路径覆盖 + 转 AncestralAdjacencyGraph
    ancestor = G.to_ancestral_graph()
    ancestor.events = G.events
    return ancestor
```

**优点**:
- 所有颜色操作集中在一处，避免散落函数间
- `resolve_all_events()` 串联完整的处理流水线
- `to_ancestral_graph()` 自动路径覆盖 + 格式转换
- 可单独调用 `find_cycles()` + `classify_cycle()` 调试环模式
- 与 `nx.Graph` 兼容，可利用 networkx 全部图算法

---

## 5. 测试策略

### 5.1 三阶段闭环测试

```
1. 生成模拟数据（evolution_simulator_ak.py）
   → orthogroups.tsv, species.nwk, genes.gff, events.tsv（truth）

2. 运行祖先重建（AK.py --reconstruction-algorithm v4）
   → AKR.events.tsv（detected）, AKR.*.anc.gfa

3. 评估重建质量（python3.11 -m soi.eval_ak）
   → 染色体数对比 + 按规模分类的事件对比 + Micro F1
```

### 5.2 评估指标

| 指标 | 目标值 |
|------|--------|
| 染色体数准确率 | 100% |
| 事件检测Micro F1 | > 0.85 |
| 分支级对比 | 优先（vs 节点级） |
| Large-scale F1 | > 0.80 |
| Small-scale F1 | > 0.80 |
| Gene-level F1 | > 0.70 |

### 5.3 真实数据验证

十字花科数据：
- Crambe_hispanica + Thlaspi_praecox → 祖先 ≈ 7条染色体
- + Arabis_montbretiana + Arabis_sagittata → 祖先 ≈ 8条染色体

---

## 6. 验收标准

### 6.1 功能验收

- [ ] 所有事件类型能被正确检测（II, TI, indel, dup, UT, EEJ, NCF, fission, RT, URT）
- [ ] 染色体数验证100%通过（普通节点 + pre-WGD节点）
- [ ] 外类群投票正确极化gain/loss
- [ ] 分支级事件记录（用于与模拟器truth对比）
- [ ] 亚基因组间重排不破坏pre-WGD collapse

### 6.2 性能验收

- [ ] 模拟数据事件F1 > 0.85
- [ ] Large-scale F1 > 0.80
- [ ] 真实数据染色体数与文献一致
- [ ] 单节点重建时间 < 5分钟

---

## 7. 附录

### 7.1 文件位置

- 本文档：`/media/40T/wlx/zrg/users/zhangrenang/OrthoIndex/TAKR_v4_event_driven_design.md`
- 开发者文档：`/media/40T/wlx/zrg/users/zhangrenang/OrthoIndex/TAKR_v4_developer_doc.md`
- 原始任务书：`/media/40T/wlx/zrg/users/zhangrenang/OrthoIndex/TAKR_v3_develop.md`
- 主代码：`/media/40T/wlx/zrg/users/zhangrenang/OrthoIndex/soi/AK.py`
- 统一事件类型：`/media/40T/wlx/zrg/users/zhangrenang/OrthoIndex/soi/takr_events.py`

### 7.2 版本变更

| 版本 | 日期 | 变更 |
|------|------|------|
| v4.3 | 2026-04-26 | 新增 ColoredGraph 类接口设计（§4.9） — 封装颜色标签管理、环检测+事件分类、indel/loss跨越边检测、路径覆盖、to_ancestral_graph转换 |
| v4.2 | 2026-04-26 | 事件图结构定义（4种冲突类型）、颜色标签(species, chrom)、duplication单物种优先处理、统一冲突检测算法（环识别→判类型→移除冲突边→路径覆盖）、无向图基础、多倍体collapse方案 |
| v4.1 | 2026-04-26 | 移除CSI，添加URT，精简文档 |
| v4.0 | 2026-04-25 | 初始版本 |

---

**文档版本**: v4.2  
**更新日期**: 2026-04-26  
**基于**: TAKR_v3_develop.md v1.1
