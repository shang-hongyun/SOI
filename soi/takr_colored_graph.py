#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
takr_colored_graph.py - ColoredGraph: 彩色邻接图 + 环检测 + 路径覆盖

每条边记录一组颜色标签 (child_id, chromosome_id)，
一条边可以有多个颜色（多个孩子共享该邻接）。
颜色集为空时自动移除边。

使用无向图 (nx.Graph) 存储，方向在路径覆盖后由 AncestralAdjacencyGraph 处理。
"""

import logging
from collections import defaultdict, Counter
from typing import Dict, List, Optional, Set, Tuple

import networkx as nx

from .takr_events import TAKREvent

logger = logging.getLogger(__name__)


def _merge_chromosome_paths(path_a, path_b):
    """Merge two subgenome chromosome paths into one pre-WGD path.
    
    Strategy: union of HOGs, keep path_a order first, append unique from path_b.
    """
    seen = set()
    merged = []
    for h in path_a:
        if h not in seen:
            merged.append(h)
            seen.add(h)
    for h in path_b:
        if h not in seen:
            merged.append(h)
            seen.add(h)
    return merged


class ColoredGraph:
    """彩色邻接图 — 一条边可有多色 (child_id, chrom_idx)。"""

    def __init__(self, hog_level: str):
        self._graph = nx.Graph()  # edge attr: 'colors' = set of (child_id, int)
        self.hog_level = hog_level
        self.events: List[TAKREvent] = []
        # 记录每个 child 的端粒集，用于 indel 检测和路径覆盖
        self._child_telomeres: Dict[str, Set] = {}
        # 记录每个 child 的所有 HOG 集合
        self._child_hogs: Dict[str, set] = {}
        # 记录每个 child 的染色体数，用于验证
        self._child_chrom_counts: Dict[str, int] = {}

    # ==================== 构建 ====================

    def add_edge(self, h1, h2, child_id: str, chrom_idx: int):
        """添加一条边，记录颜色 (child_id, chrom_idx) 和方向。

        direction: +1 = h1 在 h2 左边 (h1→h2), -1 = h2 在 h1 左边。
        如果边已存在，只把颜色和方向加入已有边的集合。
        """
        if self._graph.has_edge(h1, h2):
            colors = self._graph[h1][h2]['colors']
            colors.add((child_id, chrom_idx))
            # 更新方向
            directions = self._graph[h1][h2].get('directions', set())
            directions.add((child_id, chrom_idx, 1))
            self._graph[h1][h2]['directions'] = directions
        else:
            self._graph.add_edge(h1, h2,
                                 colors={(child_id, chrom_idx)},
                                 directions={(child_id, chrom_idx, 1)})

    def add_child(self, source_id: str, child_graph) -> int:
        """把一个孩子的所有染色体邻接以该孩子的颜色加入。
        child_graph: AncestralAdjacencyGraph (含 telomeres 和 chromosomes)
        返回该孩子的染色体数。
        """
        chrom_count = 0
        child_hogs = set()
        child_tels = set()

        for chrom_idx, chrom_nodes in enumerate(child_graph.chromosomes):
            chrom_count += 1
            # 过滤 telomere 节点，只取 HOG 节点
            hogs = [n for n in chrom_nodes if n not in child_graph.telomeres]
            child_hogs.update(hogs)
            for i in range(len(hogs) - 1):
                self.add_edge(hogs[i], hogs[i + 1], source_id, chrom_idx)

        # 收集端粒信息
        for tel in child_graph.telomeres:
            child_tels.add(tel)
        # child graph 中的端粒连接的 HOG
        for h1, h2 in child_graph.get_adjacencies(include_telomere=True):
            if h1 in child_graph.telomeres and h2 in child_graph.gene_nodes:
                child_tels.add(h2)
            elif h2 in child_graph.telomeres and h1 in child_graph.gene_nodes:
                child_tels.add(h1)

        self._child_telomeres[source_id] = child_tels
        self._child_hogs[source_id] = child_hogs
        self._child_chrom_counts[source_id] = chrom_count
        return chrom_count

    # ==================== 查询 ====================

    def get_directions(self, h1, h2) -> set:
        """获取边 (h1, h2) 的方向集合: {(child_id, chrom_idx, direction), ...}"""
        if not self._graph.has_edge(h1, h2):
            return set()
        return self._graph[h1][h2].get('directions', set())

    def edge_has_direction_conflict(self, h1, h2) -> bool:
        """检查边 (h1, h2) 是否有方向冲突 — 不同孩子方向相反。

        如果 child A 认为 h1→h2 (+1)，child B 认为 h2→h1 (-1)，
        则存在方向冲突，说明发生了 inversion。
        """
        directions = self.get_directions(h1, h2)
        if len(directions) <= 1:
            return False
        has_forward = any(d == 1 for _, _, d in directions)
        has_reverse = any(d == -1 for _, _, d in directions)
        return has_forward and has_reverse

    def find_indel_shortcuts(self) -> List[Tuple]:
        """找 indel shortcuts: 跨越边 (spanning edge) 且无方向冲突。

        indel 模式: child A 有 A→B→C, child B 有 A→C。
        在图中: A-C 是 unique edge, A 和 C 在 child A 中通过 B 连通。
        A-C 边的方向: child B 认为 A→C (+1)。
        A-B 和 B-C 边的方向: child A 认为 A→B (+1), B→C (+1)。
        无方向冲突 → 这是 indel shortcut，不是 inversion。

        Returns: [(h1, h2, child_id, spanned_hogs), ...]
        """
        shortcuts = []
        for h1, h2 in self.unique_edges():
            # 检查是否有方向冲突 — 有冲突则不是 indel
            if self.edge_has_direction_conflict(h1, h2):
                continue
            # 检查是否是 spanning edge
            colors = self.get_colors(h1, h2)
            if not colors:
                continue
            child_id, _ = next(iter(colors))
            for other_id in self.children():
                if other_id == child_id:
                    continue
                other_hogs = self._child_hogs.get(other_id, set())
                if h1 not in other_hogs or h2 not in other_hogs:
                    continue
                path = self._shortest_path_through_hogs(h1, h2, other_id)
                if path and len(path) > 2:
                    shortcuts.append((h1, h2, child_id, path[1:-1]))
                    break
        return shortcuts

    def consensus_telomeres(self, min_children: int = 2) -> Set:
        """共识端粒：在 ≥min_children 个孩子中都是端粒邻接的 HOG。

        生物学原理：端粒位置保守，不会无缘无故获得或丢失。
        如果一个 HOG 在多个孩子中都与端粒相邻，它在祖先中也是端粒邻接的。
        """
        tel_counts = defaultdict(int)
        for cid, tels in self._child_telomeres.items():
            for hog in tels:
                tel_counts[hog] += 1
        return {h for h, cnt in tel_counts.items() if cnt >= min_children}

    def child_telomere_set(self) -> Set:
        """所有孩子的端粒邻接 HOG 的并集。"""
        all_tels = set()
        for tels in self._child_telomeres.values():
            all_tels.update(tels)
        return all_tels

    def is_telomere_adjacent(self, hog) -> bool:
        """该 HOG 是否在任意孩子中与端粒相邻。"""
        for tels in self._child_telomeres.values():
            if hog in tels:
                return True
        return False

    def telomere_preserving_color(self, cycle_edges, edge_colors) -> Optional[str]:
        """在交替色环中，找出保留端粒邻接的颜色。

        对每种颜色，检查该颜色的边是否涉及端粒 HOG。
        保留端粒的颜色 = 祖先态（不应被移除）。
        """
        # 统计每种颜色涉及的端粒 HOG 数
        color_tel_count = defaultdict(int)
        color_total = defaultdict(int)
        for (h1, h2), colors in zip(cycle_edges, edge_colors):
            for cid, _ in colors:
                color_total[cid] += 1
                if self.is_telomere_adjacent(h1) or self.is_telomere_adjacent(h2):
                    color_tel_count[cid] += 1

        # 找端粒比例最高的颜色
        best_color = None
        best_ratio = -1
        for cid in color_total:
            ratio = color_tel_count.get(cid, 0) / color_total[cid]
            if ratio > best_ratio:
                best_ratio = ratio
                best_color = cid

        # 只有当端粒比例有明显差异时才返回
        if best_ratio > 0:
            return best_color
        return None

    def get_colors(self, h1, h2) -> Set[Tuple[str, int]]:
        """返回 (h1, h2) 上的颜色集。"""
        if not self._graph.has_edge(h1, h2):
            return set()
        return set(self._graph[h1][h2]['colors'])

    def edge_count(self) -> int:
        return self._graph.number_of_edges()

    def node_count(self) -> int:
        return self._graph.number_of_nodes()

    def all_hogs(self) -> Set:
        """所有 HOG 节点。"""
        return {n for n in self._graph.nodes if not isinstance(n, tuple)
                or len(n) != 2 or n[1] not in ('L', 'R')}

    def shared_edges(self) -> List[Tuple]:
        """返回多于一个颜色的边 → 祖先共享邻接。"""
        result = []
        for h1, h2, data in self._graph.edges(data=True):
            if len(data['colors']) > 1:
                result.append((h1, h2))
        return result

    def unique_edges(self) -> List[Tuple]:
        """返回只有一个颜色的边 → 可能为衍生边。"""
        result = []
        for h1, h2, data in self._graph.edges(data=True):
            if len(data['colors']) == 1:
                result.append((h1, h2))
        return result

    def edges_by_color(self, color: Tuple[str, int]) -> List[Tuple]:
        """返回所有带有指定颜色的边。"""
        result = []
        for h1, h2, data in self._graph.edges(data=True):
            if color in data['colors']:
                result.append((h1, h2))
        return result

    def children(self) -> Set[str]:
        """所有参与颜色的 child_id 集合。"""
        children = set()
        for _, _, data in self._graph.edges(data=True):
            for child_id, _ in data['colors']:
                children.add(child_id)
        return children

    # ==================== 修改 ====================

    def remove_edge_color(self, h1, h2, color: Tuple[str, int]):
        """移除边上的一个颜色。如果该边没有其他颜色，移除整条边。"""
        if not self._graph.has_edge(h1, h2):
            return
        colors = self._graph[h1][h2]['colors']
        colors.discard(color)
        if not colors:
            self._graph.remove_edge(h1, h2)

    def remove_edges_with_color(self, color: Tuple[str, int]):
        """移除所有带有指定颜色的边。"""
        to_remove = list(self.edges_by_color(color))
        for h1, h2 in to_remove:
            self.remove_edge_color(h1, h2, color)

    # ==================== Indel/Loss 检测 ====================

    def find_spanning_edges(self) -> List[Tuple]:
        """找跨越边: unique_edge (a, b) 的端点在另一孩子的路径中连通。

        对每条 unique_edges，检查 a 和 b 是否在同一孩子的另一条路径上:
        - 获取该 unique_edge 的颜色
        - 对于其他孩子的图，检查 a 和 b 之间是否有路径（且路径长度 > 1）

        Returns: [(h1, h2, child_id, spanned_hogs), ...]
        """
        spanning = []
        for h1, h2 in self.unique_edges():
            colors = self.get_colors(h1, h2)
            if not colors:
                continue
            child_id, _ = next(iter(colors))

            # 对每个其他孩子，检查是否有跨越关系
            for other_id in self.children():
                if other_id == child_id:
                    continue
                other_hogs = self._child_hogs.get(other_id, set())
                if h1 not in other_hogs or h2 not in other_hogs:
                    continue

                # 检查 h1 和 h2 是否在 other 的图中通过中间节点连通
                path = self._shortest_path_through_hogs(h1, h2, other_id)
                if path and len(path) > 2:  # 中间至少有一个 HOG
                    spanned = path[1:-1]  # 跳过两端
                    spanning.append((h1, h2, child_id, spanned))
                    break
        return spanning

    def _shortest_path_through_hogs(self, a, b, child_id: str) -> Optional[List]:
        """在指定孩子的子图中找 a→b 的最短路径（仅用该孩子的边）。"""
        G_sub = nx.Graph()
        for h1, h2, data in self._graph.edges(data=True):
            if any(cid == child_id for cid, _ in data['colors']):
                G_sub.add_edge(h1, h2)
        try:
            path = nx.shortest_path(G_sub, a, b)
            return path
        except (nx.NetworkXNoPath, nx.NodeNotFound):
            return None

    def resolve_indels(self, outgroups: Dict = None):
        """检测并解决 indel/loss 冲突。移除所有 indel shortcuts，记录 events。

        使用方向感知: 只移除无方向冲突的跨越边 (indel shortcuts)。
        有方向冲突的跨越边 (inversion) 留给 Phase 4 处理。
        """
        while True:
            shortcuts = self.find_indel_shortcuts()
            if not shortcuts:
                break
            for h1, h2, child_id, spanned_hogs in shortcuts:
                # 确定极性
                if outgroups and child_id in outgroups:
                    og_hogs = outgroups[child_id]
                    # 检查外类群是否包含跨越边覆盖的 HOG
                    any_in_og = any(h in og_hogs for h in spanned_hogs)
                    if any_in_og:
                        # 外类群有这些 HOG → 祖先有 → 当前孩子丢失
                        event_type = 'gene_loss'
                        # 移除跨越边（它是衍生边）
                        color = next(iter(self.get_colors(h1, h2)))
                        self.remove_edge_color(h1, h2, color)
                        logger.debug("  [indel] loss: remove spanning %s-%s (%s)",
                                     h1, h2, child_id)
                    else:
                        # 外类群无这些 HOG → 祖先无 → 当前孩子获得
                        event_type = 'gene_gain'
                        # 移除路径边（它们由插入产生）
                        for ch in spanned_hogs:
                            for u, v in list(self._graph.edges(ch)):
                                colors = self.get_colors(u, v)
                                if any(cid == child_id for cid, _ in colors):
                                    for c in list(colors):
                                        if c[0] == child_id:
                                            self.remove_edge_color(u, v, c)
                        logger.debug("  [indel] gain: remove path edges for %s (%s)",
                                     spanned_hogs, child_id)
                else:
                    event_type = 'gene_indel'
                    # 无外类群 => 移除跨越边（简单方案）
                    color = next(iter(self.get_colors(h1, h2)))
                    self.remove_edge_color(h1, h2, color)

                # 记录事件
                self.events.append(TAKREvent(
                    event_type=event_type,
                    branch=f"{self.hog_level}-{child_id}",
                    genes_involved=list(spanned_hogs),
                    desc=f"{event_type}: {len(spanned_hogs)} HOGs in {child_id}",
                    support=len(spanned_hogs),
                ))

        # 插入节点检测: 度=2 且所有边都是 unique 的节点
        # 模式: h4--[unique]--h5--[unique]--h6, 且 h4-h6 有 shared 边
        # → h5 是插入产物，删除 h5 恢复祖先 h4-h6
        self._resolve_inserted_nodes(outgroups)

    def _resolve_inserted_nodes(self, outgroups: Dict = None):
        """检测并删除插入节点。

        插入节点: 度=2，两条边都是 unique (单孩子)，且两个邻居有 shared 边。
        祖先态: 邻居直接相邻 (shared 边)
        衍生态: 插入了中间节点 (unique 边)
        → 删除插入节点及其边，恢复祖先邻接。
        """
        removed = set()
        for h in list(self._graph.nodes()):
            if h in removed:
                continue
            if h not in self.all_hogs():
                continue
            deg = self._graph.degree(h)
            if deg != 2:
                continue

            # 检查两条边是否都是 unique
            neighbors = list(self._graph.neighbors(h))
            if len(neighbors) != 2:
                continue
            n1, n2 = neighbors
            c1 = self.get_colors(h, n1)
            c2 = self.get_colors(h, n2)
            if len(c1) != 1 or len(c2) != 1:
                continue  # 有 shared 边，不是插入节点

            # 检查两个邻居之间是否有 shared 边
            if not self._graph.has_edge(n1, n2):
                continue
            shared_colors = self.get_colors(n1, n2)
            if len(shared_colors) < 2:
                continue  # 邻居之间也是 unique，不是简单的插入

            # 确认是插入节点: 删除 h 及其边
            child_id = next(iter(c1))[0]
            self._graph.remove_node(h)
            removed.add(h)
            logger.debug("  [indel] inserted node: remove %s between %s-%s (%s)",
                         h, n1, n2, child_id)
            self.events.append(TAKREvent(
                event_type='gene_gain',
                branch=f"{self.hog_level}-{child_id}",
                genes_involved=[h],
                desc=f"gene_gain: {h} inserted between {n1}-{n2} in {child_id}",
                support=1,
            ))

    # ==================== 环检测 + 事件分类 ====================

    def find_cycles(self) -> List[List]:
        """nx.cycle_basis 包装。"""
        try:
            return nx.cycle_basis(self._graph)
        except Exception:
            return []

    def classify_cycle(self, cycle: List) -> Tuple[Optional[str], List, Dict]:
        """分析环的颜色模式，判断事件类型。

        方向感知: 先检查环中是否有方向冲突的边对。
        - 有方向冲突 → 真实结构事件 (inversion/RT)
        - 无方向冲突 → indel shortcut 假环 (不应存在，标记为 gene_indel)

        Returns:
            (event_type, conflict_edges, info_dict)
            event_type = None 表示无法分类
            conflict_edges = 需要移除的边列表 [(h1, h2), ...]
            info_dict = 用于调试的额外信息
        """
        n = len(cycle)
        if n < 3:
            return None, [], {}

        # 构建环的边序列
        edges = [(cycle[i], cycle[(i + 1) % n]) for i in range(n)]
        edge_colors = [self.get_colors(u, v) for u, v in edges]

        info = {
            'cycle': cycle,
            'edges': edges,
            'edge_colors': edge_colors,
            'n_unique': sum(1 for c in edge_colors if len(c) == 1),
            'n_shared': sum(1 for c in edge_colors if len(c) > 1),
        }

        # == 方向冲突检查: 区分真实环 (inversion/RT) 和 indel 假环 ==
        has_dir_conflict = False
        for (u, v) in edges:
            if self.edge_has_direction_conflict(u, v):
                has_dir_conflict = True
                break
        info['has_direction_conflict'] = has_dir_conflict

        if not has_dir_conflict:
            # 无方向冲突 → indel shortcut 假环
            # 这些边不应该形成环，标记为 gene_indel 让 resolve_indels 处理
            info['resolution'] = 'no_direction_conflict_indel'
            return 'gene_indel', list(edges), info

        # 统计每种颜色的出现次数
        color_counts = Counter()
        for colors in edge_colors:
            for c in colors:
                color_counts[c] += 1

        info['color_counts'] = dict(color_counts)

        # 如果颜色数 = 0 (无颜色)，无法分类
        if not color_counts:
            return None, [], info

        # == 情况 A: 4-HOG 环，颜色交替 (A1, B1, A1, B1) ==
        if n >= 4:
            # 检查是否交替: edges[i] 和 edges[i+2] 应有相同的颜色集
            alternating = True
            for i in range(n):
                j = (i + 2) % n
                if edge_colors[i] != edge_colors[j]:
                    alternating = False
                    break

            if alternating and n >= 4:
                # 交替模式 → inversion 或 RT/URT
                # 检查是否涉及不同染色体
                chroms = set()
                for colors in edge_colors:
                    for cid, chrom_idx in colors:
                        chroms.add((cid, chrom_idx))
                # 获取边缘的颜色集
                col_set = set()
                for colors in edge_colors:
                    col_set.update(colors)

                if len(col_set) >= 3:
                    # 3+ 种颜色 → 跨染色体 → RT/URT
                    # 冲突边: 最少出现的颜色的边
                    min_count = min(color_counts.values())
                    rare_colors = {c for c, cnt in color_counts.items()
                                   if cnt == min_count}
                    conflict_edges = []
                    for (h1, h2), colors in zip(edges, edge_colors):
                        if rare_colors & colors:
                            conflict_edges.append((h1, h2))
                    info['conflict_colors'] = rare_colors
                    # 检查是否含端粒 → URT
                    has_tel = any(
                        isinstance(n, tuple) and len(n) == 2 and n[1] in ('L', 'R')
                        for n in cycle
                    )
                    etype = 'unbalanced_reciprocal_translocation' if has_tel else 'reciprocal_translocation'
                    return etype, conflict_edges, info
                else:
                    # 2 种颜色 → inversion
                    # 端粒驱动：保留端粒邻接的颜色，移除另一种
                    tel_color = self.telomere_preserving_color(edges, edge_colors)

                    if tel_color is not None:
                        # 端粒信号明确：移除不保留端粒的颜色
                        remove_colors = {c for c in color_counts if c != tel_color}
                        info['resolution'] = 'telomere-driven'
                        info['kept_color'] = tel_color
                    else:
                        # 无端粒信号：用外类群或回退到出现次数
                        min_count = min(color_counts.values())
                        remove_colors = {c for c, cnt in color_counts.items()
                                         if cnt == min_count}
                        if len(remove_colors) > 1:
                            remove_colors = {next(iter(remove_colors))}
                        info['resolution'] = 'count-based fallback'

                    conflict_edges = []
                    for (h1, h2), colors in zip(edges, edge_colors):
                        if remove_colors & colors:
                            conflict_edges.append((h1, h2))
                    info['conflict_colors'] = remove_colors

                    has_tel = any(
                        isinstance(n, tuple) and len(n) == 2 and n[1] in ('L', 'R')
                        for n in cycle
                    )
                    etype = 'telomere_inversion' if has_tel else 'inversion'
                    return etype, conflict_edges, info

        # == 简单情况: 3-HOG 环 (不可能在染色体图中, 属于错误) ==
        if n == 3:
            return 'gene_indel', list(edges), info

        # 无法分类的环 → 优先用端粒信号，否则移除最少出现的颜色
        tel_color = self.telomere_preserving_color(edges, edge_colors)
        if tel_color is not None:
            remove_colors = {c for c in color_counts if c != tel_color}
            info['resolution'] = 'telomere-driven fallback'
        else:
            min_count = min(color_counts.values())
            remove_colors = {c for c, cnt in color_counts.items() if cnt == min_count}
            info['resolution'] = 'count-based fallback'

        conflict_edges = []
        for (h1, h2), colors in zip(edges, edge_colors):
            if remove_colors & colors:
                conflict_edges.append((h1, h2))
        info['fallback'] = True
        info['fallback_reason'] = 'unclassified_cycle'
        return 'unclassified', conflict_edges, info

    # ==================== 结构重排 resolved ====================

    def resolve_structural_events(self, outgroup_adjacency: Set = None):
        """循环: find_cycles → classify_cycle → 移除外类群确认的衍生方边。

        Outgroup voting 确定极性:
        对环中每条边 (hi, hj), 映射到 parent HOG level 查外类群邻接集。
        - 外类群有该邻接 → 祖先状态 → 保存该边的孩子是祖先方
        - 外类群无该邻接 → 衍生状态 → 有该边的孩子是衍生方
        - 衍生方 = 事件发生方 → 移除该孩子的边

        无外类群 (root): 沿用最少颜色出现次数的启发式。
        """
        iteration = 0
        max_iterations = 100  # 防止无限循环
        while iteration < max_iterations:
            iteration += 1
            cycles = self.find_cycles()
            if not cycles:
                break
            logger.debug("  [structural] iteration %d: %d cycles found",
                         iteration, len(cycles))
            resolved_this_round = 0
            for cycle in cycles:
                etype, conflict_edges, info = self.classify_cycle(cycle)
                if etype is None or not conflict_edges:
                    continue

                # 构建环的边序列（用于外类群投票和回退模式）
                n_cycle = len(cycle)
                cycle_edges = [(cycle[i], cycle[(i + 1) % n_cycle])
                               for i in range(n_cycle)]

                # 用外类群投票确定衍生方（事件发生分支）
                if outgroup_adjacency is not None and etype in (
                    'inversion', 'telomere_inversion',
                    'reciprocal_translocation',
                    'unbalanced_reciprocal_translocation',
                ):
                    cycle_colors = [self.get_colors(u, v) for u, v in cycle_edges]

                    # 对每条边，映射到 parent HOG level 查外类群
                    child_derived_count = defaultdict(int)
                    child_ancestral_count = defaultdict(int)
                    for (u, v), colors in zip(cycle_edges, cycle_colors):
                        u_parent = getattr(u, 'parent',
                                           u.hog_id if hasattr(u, 'hog_id') else str(u))
                        v_parent = getattr(v, 'parent',
                                           v.hog_id if hasattr(v, 'hog_id') else str(v))
                        key = (u_parent, v_parent) if u_parent < v_parent else (v_parent, u_parent)
                        is_ancestral = key in outgroup_adjacency
                        for cid, _ in colors:
                            if is_ancestral:
                                child_ancestral_count[cid] += 1
                            else:
                                child_derived_count[cid] += 1

                    # 出现衍生边的孩子 = 事件发生方
                    if child_derived_count:
                        derived_child_id = max(child_derived_count,
                                               key=child_derived_count.get)
                    else:
                        derived_child_id = None

                    logger.debug("  [structural] cycle %s: og-derived=%s, "
                                 "ancestral=%s, derived=%s",
                                 cycle, derived_child_id,
                                 dict(child_ancestral_count),
                                 dict(child_derived_count))

                    # 只移除衍生方的边
                    if derived_child_id:
                        for h1, h2 in conflict_edges:
                            colors = self.get_colors(h1, h2)
                            for cid, chrom in list(colors):
                                if cid == derived_child_id:
                                    self.remove_edge_color(h1, h2, (cid, chrom))
                    else:
                        # 所有边都是祖先状态？回退到最少颜色出现次数的启发式
                        for h1, h2 in conflict_edges:
                            self._remove_rare_color_edge(h1, h2, cycle, cycle_edges)
                else:
                    # 无外类群 (root) 或 unclassified 环 → 启发式
                    for h1, h2 in conflict_edges:
                        self._remove_rare_color_edge(h1, h2, cycle, cycle_edges)

                # 记录事件
                event = TAKREvent(
                    event_type=etype,
                    branch=f"{self.hog_level}-{derived_child_id}" if derived_child_id else self.hog_level,
                    genes_involved=list(cycle),
                    desc=f"{etype}: {len(cycle)} HOG cycle resolved (iter {iteration})",
                    support=1,
                )
                self.events.append(event)
                resolved_this_round += 1

            if resolved_this_round == 0:
                logger.warning("  [structural] %d cycles found but none resolved",
                               len(cycles))
                break

        if iteration >= max_iterations:
            logger.warning("  [structural] max iterations reached (%d), force-breaking cycles",
                           max_iterations)
            # 强制打破剩余环：对每个环移除最少 child 的所有边
            cycles = self.find_cycles()
            for cycle in cycles:
                etype, conflict_edges, info = self.classify_cycle(cycle)
                for h1, h2 in conflict_edges:
                    colors = self.get_colors(h1, h2)
                    if len(colors) == 1:
                        self.remove_edge_color(h1, h2, next(iter(colors)))
                    elif colors:
                        self.remove_edge_color(h1, h2, next(iter(colors)))
                self.events.append(TAKREvent(
                    event_type='unclassified',
                    branch=self.hog_level,
                    genes_involved=list(cycle),
                    desc=f"force-break: {len(cycle)} HOGs",
                    support=1,
                ))

    # ==================== 边移除辅助 ====================

    def _remove_rare_color_edge(self, h1, h2, cycle, cycle_edges):
        """移除冲突边上出现次数最少的颜色（回退启发式）。"""
        colors = self.get_colors(h1, h2)
        if len(colors) == 1:
            self.remove_edge_color(h1, h2, next(iter(colors)))
        else:
            # 统计各 (cid, chrom) 在环中的出现次数
            child_count = {}
            for (u, v) in cycle_edges:
                for cid, chrom in self.get_colors(u, v):
                    child_count[(cid, chrom)] = child_count.get((cid, chrom), 0) + 1
            if child_count:
                min_child = min(child_count, key=child_count.get)
                if min_child in colors:
                    self.remove_edge_color(h1, h2, min_child)
                else:
                    self.remove_edge_color(h1, h2, next(iter(colors)))
            elif colors:
                self.remove_edge_color(h1, h2, next(iter(colors)))

    # ==================== 路径覆盖 ====================

    def path_cover(self) -> List[List]:
        """端粒约束路径覆盖。

        生物学约束：每条祖先染色体有且仅有 2 个端粒（两端各一个）。
        端粒位置通过多个孩子的共识确定。

        算法：
        1. 找共识端粒 HOG（≥2 个孩子中都与端粒相邻）
        2. 从共识端粒出发，沿邻接行走
        3. 行走时优先选择端粒邻接度高的邻居
        4. 到达另一个端粒时停止 → 一条完整染色体
        5. 处理剩余节点（无端粒的连通分量）

        Returns:
            [[hog1, hog2, ...], ...] — 每条染色体一条路径
        """
        paths = []
        visited = set()

        # 共识端粒：≥2 个孩子中都与端粒相邻的 HOG
        cons_tels = self.consensus_telomeres(min_children=2)
        # 备用：所有孩子的端粒并集
        all_tels = self.child_telomere_set()
        # 度=1 的节点（图结构上的端点）
        degree_ends = {n for n, deg in self._graph.degree()
                       if deg == 1 and n in self.all_hogs()}

        # 选择端粒锚点集：共识端粒 > 度=1 > 所有端粒
        if cons_tels:
            anchor_tels = cons_tels
        elif degree_ends:
            anchor_tels = degree_ends
        else:
            anchor_tels = all_tels

        # 从每个锚点出发行走
        for start in sorted(anchor_tels, key=str):
            if start in visited or not self._graph.has_node(start):
                continue

            path = self._walk_telomere_path(start, visited, anchor_tels)
            if path and len(path) >= 1:
                paths.append(path)
                visited.update(path)

        # 处理剩余未访问节点（无端粒的连通分量 + 孤立节点）
        all_nodes = self.all_hogs()
        unvisited = all_nodes - visited
        if unvisited:
            G_rem = self._graph.subgraph(unvisited)
            for comp_nodes in nx.connected_components(G_rem):
                comp = list(comp_nodes)
                if len(comp) >= 2:
                    start = comp[0]
                    path = self._walk_simple(start, visited)
                    if path:
                        paths.append(path)
                        visited.update(path)
                elif len(comp) == 1:
                    # 孤立节点：不单独成染色体（noise），跳过
                    visited.update(comp)

        return paths

    def _walk_telomere_path(self, start, visited: Set,
                            anchor_tels: Set) -> Optional[List]:
        """从端粒 HOG 出发，走到另一个端粒 HOG。

        行走策略：
        - 度=1 的邻居优先（端点）
        - 端粒锚点 HOG 次之
        - 遇到另一个锚点时停止
        """
        path = [start]
        curr = start
        visited_local = visited | {start}

        while True:
            neighbors = [n for n in self._graph.neighbors(curr)
                         if n not in visited_local and n in self.all_hogs()]
            if not neighbors:
                break

            if len(neighbors) == 1:
                nxt = neighbors[0]
            else:
                # 多个邻居时，按优先级排序
                def neighbor_priority(n):
                    if self._graph.degree(n) == 1:
                        return 0
                    if n in anchor_tels:
                        return 1
                    return 2
                neighbors.sort(key=neighbor_priority)
                nxt = neighbors[0]

            visited_local.add(nxt)
            path.append(nxt)
            curr = nxt

            # 到达另一个锚点 → 停止
            if curr != start and curr in anchor_tels:
                break

        return path if len(path) >= 1 else None

    def _walk_simple(self, start, visited: Set) -> Optional[List]:
        """简单行走（无端粒约束），用于处理剩余连通分量。"""
        path = [start]
        curr = start
        visited_local = visited | {start}

        while True:
            neighbors = [n for n in self._graph.neighbors(curr)
                         if n not in visited_local and n in self.all_hogs()]
            if not neighbors:
                break
            nxt = neighbors[0]
            visited_local.add(nxt)
            path.append(nxt)
            curr = nxt

        return path if len(path) >= 1 else None

    # ==================== 转换 ====================

    def to_ancestral_graph(self):
        """转换为 AncestralAdjacencyGraph。

        复用 resolve_all_events 中缓存的 paths，避免重复计算。
        """
        from .AK import AncestralAdjacencyGraph

        result = AncestralAdjacencyGraph(node_id=self.hog_level)

        # 复用缓存的 paths
        paths = getattr(self, '_cached_paths', None) or self.path_cover()
        if not paths:
            logger.warning("  [colored] no paths found for %s", self.hog_level)
            return result

        for chrom_idx, path in enumerate(paths):
            # 添加 HOG 节点
            for n in path:
                if n not in result.graph:
                    result.graph.add_node(n)
                    result.gene_nodes.add(n)

            # 添加有向边（基于路径顺序）
            for i in range(len(path) - 1):
                result.graph.add_edge(path[i], path[i + 1])

            # 添加端粒
            chrom_name = f"{self.hog_level}_chrom_{chrom_idx}"
            left_tel = (chrom_name, 'L')
            right_tel = (chrom_name, 'R')
            result.graph.add_node(left_tel, telomere=True)
            result.graph.add_node(right_tel, telomere=True)
            result.telomeres.add(left_tel)
            result.telomeres.add(right_tel)
            result.graph.add_edge(left_tel, path[0])
            result.graph.add_edge(path[-1], right_tel)

        result.events = self.events
        return result

    # ==================== 一键执行 ====================

    def collapse_wgd(self, ploidy: int):
        """ColoredGraph-based WGD collapse: post-WGD → pre-WGD.

        post-WGD 图有 ploidy×pre 条染色体，它们是天然分开的（无跨染色体边）。
        WGD collapse = 染色体配对：用 HOG 内容相似度找亚基因组同位染色体。

        Pipeline:
        1. Path_cover → 获取 N 条染色体
        2. 染色体配对: 每对染色体的 HOG Jaccard 相似度 → 找最佳匹配
        3. 每对合并为一条 pre-WGD 染色体（HOG 去重）

        Args:
            ploidy: 倍性 (2=四倍体, 3=六倍体, ...)

        Returns:
            ColoredGraph — pre-WGD 图
        """
        # Step 1: 获取 post-WGD 染色体
        paths = self.path_cover()
        if not paths:
            logger.warning("  [collapse_wgd] no paths found, returning original graph")
            return self

        n_chrom = len(paths)
        logger.debug("  [collapse_wgd] %s: %d post-WGD chromosomes (ploidy=%d)",
                      self.hog_level, n_chrom, ploidy)

        if n_chrom < ploidy:
            logger.warning("  [collapse_wgd] fewer chroms (%d) than ploidy (%d), skipping",
                           n_chrom, ploidy)
            return self

        # Step 2: 构建染色体 × HOG 矩阵
        chrom_hogs = []  # [{hog_id}, ...] per chromosome
        for path in paths:
            chrom_hogs.append(set(path))

        # 计算所有染色体对的 Jaccard 相似度
        n = len(chrom_hogs)
        import itertools
        pairs = []
        for i, j in itertools.combinations(range(n), 2):
            inter = len(chrom_hogs[i] & chrom_hogs[j])
            union = len(chrom_hogs[i] | chrom_hogs[j])
            jaccard = inter / union if union > 0 else 0.0
            if inter > 0:
                pairs.append((jaccard, inter, i, j))

        # Step 3: 贪心匹配 — 从相似度最高的对开始
        pairs.sort(reverse=True)
        paired = set()
        pre_chroms = []  # [(path_i, path_j), ...]

        for jaccard, inter, i, j in pairs:
            if i in paired or j in paired:
                continue
            pre_chroms.append((paths[i], paths[j]))
            paired.add(i)
            paired.add(j)

        # 未配对的染色体作为单条保留
        for i in range(n):
            if i not in paired:
                pre_chroms.append((paths[i],))

        n_pre = len(pre_chroms)
        logger.debug("  [collapse_wgd] %s: %d post-WGD -> %d pre-WGD chroms (paired=%d)",
                      self.hog_level, n_chrom, n_pre, len(pre_chroms))

        # Step 4: 构建 pre-WGD 图（每对染色体合并）
        G2 = ColoredGraph(hog_level=f"{self.hog_level}_preWGD")

        # 记录事件
        paired_count = 0
        for chrom_pair in pre_chroms:
            if len(chrom_pair) == 2:
                # 一对亚基因组同位染色体 → 合并（HOG 去重）
                path_a, path_b = chrom_pair
                merged = _merge_chromosome_paths(path_a, path_b)
                for i in range(len(merged) - 1):
                    G2.add_edge(merged[i], merged[i + 1],
                                f"pre_chr{paired_count}", paired_count)
                paired_count += 1
                # 记录配对事件
                shared_hogs = set(path_a) & set(path_b)
                if shared_hogs:
                    G2.events.append(TAKREvent(
                        event_type='wgd_collapse',
                        branch=f"{self.hog_level}_preWGD-{self.hog_level}",
                        genes_involved=list(shared_hogs),
                        desc=f"wgd_collapse: paired chroms ({len(path_a)}+{len(path_b)} HOGs, {len(shared_hogs)} shared)",
                        support=ploidy,
                    ))
            else:
                # 未配对的单条染色体
                path = chrom_pair[0]
                for i in range(len(path) - 1):
                    G2.add_edge(path[i], path[i + 1],
                                f"pre_chr{paired_count}", paired_count)
                paired_count += 1

        return G2

    # ==================== 共线性块压缩 ====================

    def _build_synteny_blocks(self, min_block_size: int = 2):
        """构建共线性块：最大连续共享边路径。

        块 = 一段 HOG 序列，其中每对相邻 HOG 之间都有共享边（≥2 孩子共有）。
        块间边界 = 唯一边或断裂点。

        self._blocks: {block_id: [hog1, hog2, ...]} — 块内 HOG 顺序
        self._hog_to_block: {hog: block_id} — HOG 到块的映射
        self._block_adj: {bid: set(bid_j)} — 块间邻接关系（用于构建块级图）
        """
        # 共享边图
        shared_G = nx.Graph()
        for h1, h2, data in self._graph.edges(data=True):
            if len(data['colors']) >= 2:
                shared_G.add_edge(h1, h2)

        blocks = {}
        hog_to_block = {}
        assigned = set()

        # 对共享边图的每个连通分量找线性路径
        for comp in nx.connected_components(shared_G):
            sub = shared_G.subgraph(comp)
            deg = dict(sub.degree())

            # 找断面（度 ≠ 2 的节点）作为路径起点
            breakpoints = {n for n, d in deg.items() if d != 2}
            if not breakpoints:
                # 环状分量：任选一个节点打断
                breakpoints = {list(comp)[0]}

            # 从每个断面出发，沿度=2的节点行走
            for bp in list(breakpoints):
                if bp in assigned:
                    continue
                # 向两个方向行走
                for start in (bp,):
                    if start in assigned:
                        continue
                    path = [start]
                    assigned.add(start)
                    # 向前走到下一个断面或末端
                    curr = start
                    prev = None
                    while True:
                        deg2_neighbors = [
                            n for n in sub.neighbors(curr)
                            if n != prev and n not in assigned
                            and deg.get(n, 0) == 2
                        ]
                        if not deg2_neighbors:
                            break
                        nxt = deg2_neighbors[0]
                        assigned.add(nxt)
                        path.append(nxt)
                        prev, curr = curr, nxt

                    if len(path) >= min_block_size:
                        bid = "blk_{}".format(len(blocks))
                        blocks[bid] = path
                        for h in path:
                            hog_to_block[h] = bid
                    else:
                        # 路径太短 (<min_block_size): 不建块，标记为未分配
                        # 留给后续合并阶段处理
                        pass

        # 不在共享图中的 HOG → 可能是 insertion 产物
        # 留给 Phase 2 indel 检测处理，不在此合并

        # 剩余仍未分配的 HOG (无邻居或邻居也未分配) → 单独成块
        for h in self.all_hogs():
            if h not in hog_to_block:
                bid = "blk_{}".format(len(blocks))
                blocks[bid] = [h]
                hog_to_block[h] = bid

        self._blocks = blocks
        self._hog_to_block = hog_to_block

        n_blk = len(blocks)
        n_multi = sum(1 for p in blocks.values() if len(p) >= 2)
        n_single = sum(1 for p in blocks.values() if len(p) == 1)
        logger.debug("  [blocks] %d blocks (%d multi-HOG, %d singleton)",
                     n_blk, n_multi, n_single)
        return blocks

    def _compress_to_block_level(self):
        """将 HOG 级图压缩为块级图。

        block_graph: nx.Graph
          - 节点 = block_id
          - 边 = 块间邻接，colors = 所有跨块 HOG 边颜色的并集

        block_order: {block_id: [block_id, ...]}
          块之间的顺序关系（每个孩子中相邻块的关系）
        """
        if not hasattr(self, '_blocks') or not self._blocks:
            self._build_synteny_blocks()

        block_G = nx.Graph()
        for bid in self._blocks:
            block_G.add_node(bid)

        # 块间边：聚合所有跨块 HOG 边
        for h1, h2, data in self._graph.edges(data=True):
            b1 = self._hog_to_block.get(h1)
            b2 = self._hog_to_block.get(h2)
            if b1 and b2 and b1 != b2:
                if block_G.has_edge(b1, b2):
                    block_G[b1][b2]['colors'].update(data['colors'])
                else:
                    block_G.add_edge(b1, b2, colors=set(data['colors']))

        self._block_graph = block_G

        # 清理无颜色的块间边
        empty_edges = [(b1, b2) for b1, b2, data in block_G.edges(data=True)
                       if not data.get('colors', set())]
        for b1, b2 in empty_edges:
            block_G.remove_edge(b1, b2)
        if empty_edges:
            logger.debug("  [blocks] removed %d empty-colored block edges",
                         len(empty_edges))

        # 块间顺序：对每个孩子，统计相邻块
        child_block_order = defaultdict(list)
        for cid in self.children():
            seen_blocks = set()
            chrom_path = []
            # 收集该孩子所有 HOG → 块 + 顺序
            child_hogs = sorted(self._child_hogs.get(cid, []),
                                key=lambda h: str(h.hog_id if hasattr(h, 'hog_id') else h))
            # 用边顺序重建块序列
            for h1, h2, data in self._graph.edges(data=True):
                if any(c == cid for c, _ in data['colors']):
                    b1 = self._hog_to_block.get(h1)
                    b2 = self._hog_to_block.get(h2)
                    if b1 and b2 and b1 != b2:
                        order = (b1, b2)
                        if order not in child_block_order[cid]:
                            child_block_order[cid].append(order)

        self._child_block_order = dict(child_block_order)

        logger.debug("  [blocks] block graph: %d nodes, %d edges",
                     block_G.number_of_nodes(), block_G.number_of_edges())
        return block_G

    def _find_block_cycles(self) -> List[List]:
        """在块级图上找环。"""
        if not hasattr(self, '_block_graph'):
            return []
        try:
            return nx.cycle_basis(self._block_graph)
        except Exception:
            return []

    def _classify_block_cycle(self, cycle: List) -> Tuple[Optional[str], List, Dict]:
        """分析块级环的颜色模式，判断事件类型。

        块级环 vs HOG 级环：
        - 环节点 = 块（block_id），而不是 HOG
        - 无需考虑块内 HOG 的顺序差异
        - 块间边的颜色 = 两个孩子对该块间邻接的贡献

        Returns:
            (event_type, conflict_block_edges, info)
        """
        n = len(cycle)
        if n < 3:
            return None, [], {}

        bg = self._block_graph

        # 构建环的边序列
        edges = [(cycle[i], cycle[(i + 1) % n]) for i in range(n)]
        edge_colors = []
        for b1, b2 in edges:
            if bg.has_edge(b1, b2):
                edge_colors.append(set(bg[b1][b2].get('colors', set())))
            else:
                edge_colors.append(set())

        info = {
            'cycle': cycle,
            'edges': edges,
            'edge_colors': edge_colors,
            'n_unique': sum(1 for c in edge_colors if len(c) == 1),
            'n_shared': sum(1 for c in edge_colors if len(c) > 1),
        }

        # 统计每种颜色的出现次数
        color_counts = Counter()
        for colors in edge_colors:
            for c in colors:
                color_counts[c] += 1
        info['color_counts'] = dict(color_counts)

        if not color_counts:
            return None, [], info

        # 检查是否交替模式：edges[i] 和 edges[i+2] 应颜色相同
        if n >= 4:
            alternating = all(
                edge_colors[i] == edge_colors[(i + 2) % n]
                for i in range(n)
            )

            if alternating:
                # 交替块边 → inversion 或 RT/URT
                col_set = set()
                for colors in edge_colors:
                    col_set.update(colors)

                if len(col_set) >= 3:
                    # 3+ 种颜色 → 跨染色体 → RT/URT
                    min_count = min(color_counts.values())
                    rare_colors = {c for c, cnt in color_counts.items()
                                   if cnt == min_count}
                    conflict_edges = [(b1, b2) for (b1, b2), colors in zip(edges, edge_colors)
                                      if rare_colors & colors]
                    info['conflict_colors'] = rare_colors
                    etype = 'reciprocal_translocation'
                    return etype, conflict_edges, info
                else:
                    # 2 种颜色 → 染色体内块反转 → inversion
                    min_count = min(color_counts.values())
                    rare_colors = {c for c, cnt in color_counts.items()
                                   if cnt == min_count}
                    if len(rare_colors) > 1 and all(
                            cnt == min_count for cnt in color_counts.values()):
                        rare_colors = {next(iter(color_counts.keys()))}
                    conflict_edges = [(b1, b2) for (b1, b2), colors in zip(edges, edge_colors)
                                      if rare_colors & colors]
                    info['conflict_colors'] = rare_colors
                    etype = 'inversion'
                    return etype, conflict_edges, info

        # 3-cycle 在块级可能是 inversion：两条不同孩子颜色的唯一边 + 一条共享边
        if n == 3:
            # 检查：两条唯一边（颜色不同） + 一条共享边
            unique_group = []
            shared_group = []
            for (b1, b2), colors in zip(edges, edge_colors):
                if len(colors) == 1:
                    unique_group.append((b1, b2, next(iter(colors))[0]))
                else:
                    shared_group.append(colors)

            if len(unique_group) == 2 and len(shared_group) == 1:
                # 两条唯一边来自不同孩子 → 块反转
                c1, c2 = unique_group[0][2], unique_group[1][2]
                if c1 != c2:
                    # 3-cycle 中至少 2 个块是多 HOG 块才算 inversion
                    block_sizes = [len(self._blocks.get(bid, [])) for bid in cycle]
                    n_multi = sum(1 for s in block_sizes if s >= 2)
                    if n_multi < 2:
                        # 多为单例块 → 可能为 indel 边界混淆
                        info['pattern'] = 'block_3cycle_inversion_filtered'
                        return None, [], info
                    # 两条唯一边各来自一个孩子 → 2 色交替 = inversion
                    info['pattern'] = 'block_3cycle_inversion'
                    # 根据唯一边确定冲突方：谁是衍生方由外类群判定
                    # 冲突边 = 两条唯一边
                    conflict_edges = [(b1, b2) for b1, b2, _ in unique_group]
                    conflict_colors = {unique_group[0][2]: None, unique_group[1][2]: None}
                    info['conflict_colors'] = conflict_colors
                    etype = 'inversion'
                    return etype, conflict_edges, info
                elif len(unique_group[0]) == 1 and len(unique_group[1]) == 1:
                    # 两条唯一边来自同一个孩子 → 可能是其他事件
                    pass
            elif len(unique_group) == 3:
                # 三条唯一边 → indel 或错误
                etype = 'gene_indel'
                conflict_edges = [(b1, b2) for b1, b2, _ in unique_group]
                info['pattern'] = 'block_3cycle_all_unique'
                return etype, conflict_edges, info

        # 无法分类的块级环
        min_count = min(color_counts.values())
        rare_colors = {c for c, cnt in color_counts.items() if cnt == min_count}
        conflict_edges = [(b1, b2) for (b1, b2), colors in zip(edges, edge_colors)
                          if rare_colors & colors]
        info['fallback'] = True
        info['fallback_reason'] = 'unclassified_block_cycle'
        return 'unclassified', conflict_edges, info

    def _resolve_block_structural_events(self, outgroup_adjacency: Set = None):
        """在块级图上检测结构重排事件。

        对每个块级环：
        1. classify_cycle → 判断事件类型
        2. outgroup voting → 确定衍生方
        3. 移除衍生方的块间边
        4. 记录事件

        Returns: 检测到的事件数
        """
        if not hasattr(self, '_block_graph'):
            self._compress_to_block_level()

        iteration = 0
        max_iterations = 100
        events_found = 0

        while iteration < max_iterations:
            iteration += 1
            cycles = self._find_block_cycles()
            if not cycles:
                break

            # Step 1: 对所有环分类
            classified = []
            for cycle in cycles:
                etype, conflict_edges, info = self._classify_block_cycle(cycle)
                if etype is not None and conflict_edges:
                    classified.append((etype, conflict_edges, info, cycle))

            if not classified:
                break

            resolved = 0

            # Step 2: 合并相邻 inversion 3-cycle（同一 inversion 的 N-1 个环）
            inv_cycles = [(ce, info, cyc) for et, ce, info, cyc in classified
                          if et == 'inversion']
            other_cycles = [(et, ce, cyc) for et, ce, info, cyc in classified
                            if et != 'inversion']

            # 合并 inversion: 重叠块的环合并为一个事件
            merged_inv_groups = []
            used_inv = set()
            for i, (ce_i, info_i, cyc_i) in enumerate(inv_cycles):
                if i in used_inv:
                    continue
                group = [i]
                blocks_set = set(cyc_i)
                # 找到所有与当前环共享块的相邻环
                changed = True
                while changed:
                    changed = False
                    for j, (ce_j, info_j, cyc_j) in enumerate(inv_cycles):
                        if j in used_inv or j in group:
                            continue
                        if blocks_set & set(cyc_j):
                            group.append(j)
                            blocks_set.update(cyc_j)
                            changed = True
                used_inv.update(group)
                # 收集组内所有冲突边 + 块
                all_conflict = []
                all_blocks = set()
                for g in group:
                    all_conflict.extend(inv_cycles[g][0])
                    all_blocks.update(inv_cycles[g][2])
                merged_inv_groups.append((all_conflict, all_blocks))

            # Step 3: 处理合并后的 inversion 组 + 其他事件
            for conflict_edges, blocks_set in merged_inv_groups:
                # 去重冲突边
                conflict_edges = list(set(conflict_edges))
                # 外类群投票：合并组只有一个极性
                child_derived = defaultdict(int)
                for (b1, b2) in conflict_edges:
                    colors = self._block_graph[b1][b2].get('colors', set())
                    if outgroup_adjacency is not None:
                        hogs1 = self._blocks.get(b1, [])
                        hogs2 = self._blocks.get(b2, [])
                        any_ancestral = False
                        for h1 in hogs1:
                            for h2 in hogs2:
                                if self._graph.has_edge(h1, h2):
                                    u_p = getattr(h1, 'parent',
                                                  h1.hog_id if hasattr(h1, 'hog_id') else str(h1))
                                    v_p = getattr(h2, 'parent',
                                                  h2.hog_id if hasattr(h2, 'hog_id') else str(h2))
                                    key = (u_p, v_p) if u_p < v_p else (v_p, u_p)
                                    if key in outgroup_adjacency:
                                        any_ancestral = True
                                        break
                        for cid, _ in colors:
                            if any_ancestral:
                                pass  # 祖先方不计数
                            else:
                                child_derived[cid] += 1
                    else:
                        for c in colors:
                            child_derived[c[0]] -= 1  # 无外类群回退

                if child_derived:
                    derived_cid = max(child_derived, key=child_derived.get)
                else:
                    # 回退：边上的颜色出现次数较少的
                    cc = Counter()
                    for (b1, b2) in conflict_edges:
                        for c in self._block_graph[b1][b2].get('colors', set()):
                            cc[c] += 1
                    derived_cid = min(cc, key=cc.get)[0] if cc else None

                if derived_cid is None:
                    continue

                # 移除衍生方在冲突边上的颜色
                for (b1, b2) in conflict_edges:
                    colors = self._block_graph[b1][b2].get('colors', set())
                    remaining = set()
                    for cid, chrom in colors:
                        if cid != derived_cid:
                            remaining.add((cid, chrom))
                    if remaining:
                        self._block_graph[b1][b2]['colors'] = remaining
                    else:
                        self._block_graph.remove_edge(b1, b2)

                # 记录合并后的 inversion 事件
                involved_hogs = []
                for bid in blocks_set:
                    involved_hogs.extend(self._blocks.get(bid, []))
                involved_hogs = list(set(involved_hogs))

                event = TAKREvent(
                    event_type='inversion',
                    branch=f"{self.hog_level}-{derived_cid}" if derived_cid else self.hog_level,
                    genes_involved=involved_hogs,
                    desc="inversion (block-level): {} blocks, {} HOGs".format(
                        len(blocks_set), len(involved_hogs)),
                    support=1,
                )
                self.events.append(event)
                resolved += 1
                events_found += 1

                # 同步移除 HOG 级图中衍生方的边
                for (b1, b2) in conflict_edges:
                    for h1 in self._blocks.get(b1, []):
                        for h2 in self._blocks.get(b2, []):
                            if self._graph.has_edge(h1, h2):
                                colors = self.get_colors(h1, h2)
                                for cid, chrom in list(colors):
                                    if cid == derived_cid:
                                        self.remove_edge_color(h1, h2, (cid, chrom))

            # 处理非 inversion 事件（逐环处理）
            for etype, conflict_edges, cycle in other_cycles:
                conflict_edges = list(set(conflict_edges))
                # 检查冲突边是否还在块级图中
                valid_edges = [(b1, b2) for b1, b2 in conflict_edges
                               if self._block_graph.has_edge(b1, b2)]
                if not valid_edges:
                    continue
                # 外类群投票
                child_derived2 = defaultdict(int)
                for (b1, b2) in valid_edges:
                    colors = self._block_graph[b1][b2].get('colors', set())
                    if outgroup_adjacency is not None:
                        hogs1 = self._blocks.get(b1, [])
                        hogs2 = self._blocks.get(b2, [])
                        any_ancestral = False
                        for h1 in hogs1:
                            for h2 in hogs2:
                                if self._graph.has_edge(h1, h2):
                                    u_p = getattr(h1, 'parent',
                                                  h1.hog_id if hasattr(h1, 'hog_id') else str(h1))
                                    v_p = getattr(h2, 'parent',
                                                  h2.hog_id if hasattr(h2, 'hog_id') else str(h2))
                                    key = (u_p, v_p) if u_p < v_p else (v_p, u_p)
                                    if key in outgroup_adjacency:
                                        any_ancestral = True
                                        break
                        for cid, _ in colors:
                            if not any_ancestral:
                                child_derived2[cid] += 1
                if child_derived2:
                    derived_cid2 = max(child_derived2, key=child_derived2.get)
                else:
                    cc2 = Counter()
                    for (b1, b2) in valid_edges:
                        for c in self._block_graph[b1][b2].get('colors', set()):
                            cc2[c] += 1
                    derived_cid2 = min(cc2, key=cc2.get)[0] if cc2 else None

                if derived_cid2 is None:
                    continue

                # 移除衍生方的边
                for (b1, b2) in valid_edges:
                    colors = self._block_graph[b1][b2].get('colors', set())
                    remaining = set()
                    for cid, chrom in colors:
                        if cid != derived_cid2:
                            remaining.add((cid, chrom))
                    if remaining:
                        self._block_graph[b1][b2]['colors'] = remaining
                    else:
                        self._block_graph.remove_edge(b1, b2)

                # 记录事件
                involved_hogs2 = []
                for bid in cycle:
                    involved_hogs2.extend(self._blocks.get(bid, []))
                involved_hogs2 = list(set(involved_hogs2))

                event = TAKREvent(
                    event_type=etype,
                    branch=f"{self.hog_level}-{derived_cid2}" if derived_cid2 else self.hog_level,
                    genes_involved=involved_hogs2,
                    desc="{} (block-level): {} blocks, {} HOGs".format(
                        etype, len(cycle), len(involved_hogs2)),
                    support=1,
                )
                self.events.append(event)
                resolved += 1
                events_found += 1

                # 同步移除 HOG 级边
                for (b1, b2) in valid_edges:
                    for h1 in self._blocks.get(b1, []):
                        for h2 in self._blocks.get(b2, []):
                            if self._graph.has_edge(h1, h2):
                                colors = self.get_colors(h1, h2)
                                for cid, chrom in list(colors):
                                    if cid == derived_cid2:
                                        self.remove_edge_color(h1, h2, (cid, chrom))

            if resolved == 0:
                break

        return events_found

    def _resolve_block_bridge_events(self, outgroup_adjacency: Set = None):
        """在块级图上检测桥接事件。

        块级桥接：块级唯一边连接两个块级共享分量。
        端粒保护：不移除涉及共识端粒块的边。
        """
        if not hasattr(self, '_block_graph'):
            self._compress_to_block_level()

        bg = self._block_graph

        # 块级共享边图
        shared_bg = nx.Graph()
        for b1, b2, data in bg.edges(data=True):
            if len(data.get('colors', set())) > 1:
                shared_bg.add_edge(b1, b2)

        if shared_bg.number_of_edges() == 0:
            return 0

        components = list(nx.connected_components(shared_bg))
        block_to_comp = {}
        for ci, comp in enumerate(components):
            for b in comp:
                block_to_comp[b] = ci

        # 共识端粒块：包含共识端粒 HOG 的块
        cons_tels = self.consensus_telomeres(min_children=2)
        telomere_blocks = set()
        for hog in cons_tels:
            bid = self._hog_to_block.get(hog)
            if bid:
                telomere_blocks.add(bid)

        events_found = 0
        for b1, b2, data in list(bg.edges(data=True)):
            colors = data.get('colors', set())
            if len(colors) != 1:
                continue
            c1 = block_to_comp.get(b1)
            c2 = block_to_comp.get(b2)
            if c1 is not None and c2 is not None and c1 != c2:
                # 端粒保护：不移除涉及端粒块的边
                if b1 in telomere_blocks or b2 in telomere_blocks:
                    continue
                # 块级桥接
                color = next(iter(colors))
                child_id = color[0]

                if outgroup_adjacency is not None:
                    # 检查跨块 HOG 对的 parent 级邻接
                    hogs1 = self._blocks.get(b1, [])
                    hogs2 = self._blocks.get(b2, [])
                    any_ancestral = False
                    for h1 in hogs1:
                        for h2 in hogs2:
                            if self._graph.has_edge(h1, h2):
                                u_p = getattr(h1, 'parent',
                                              h1.hog_id if hasattr(h1, 'hog_id') else str(h1))
                                v_p = getattr(h2, 'parent',
                                              h2.hog_id if hasattr(h2, 'hog_id') else str(h2))
                                key = (u_p, v_p) if u_p < v_p else (v_p, u_p)
                                if key in outgroup_adjacency:
                                    any_ancestral = True
                                    break

                    if any_ancestral:
                        etype = 'fission'
                        suffix = f" (outgroup confirms {child_id} preserved)"
                    else:
                        # 检查 NCF 条件
                        c_size = [len(components[c1]), len(components[c2])]
                        if min(c_size) < 5:
                            etype = 'ncf'
                        else:
                            etype = 'eej'
                        suffix = f" (outgroup: derived in {child_id})"
                else:
                    etype = 'bridge_unclassified'
                    suffix = " (no outgroup, root node)"

                # 展开 HOG
                involved_hogs = list(self._blocks.get(b1, [])) + list(self._blocks.get(b2, []))
                involved_hogs = list(set(involved_hogs))

                logger.debug("  [bridge] %s: %s(%d HOGs,comp%d) <-> %s(%d HOGs,comp%d) via %s",
                            etype, b1, len(self._blocks.get(b1, [])), c1,
                            b2, len(self._blocks.get(b2, [])), c2, child_id)

                self.events.append(TAKREvent(
                    event_type=etype,
                    branch=f"{self.hog_level}-{child_id}",
                    genes_involved=involved_hogs,
                    desc="{} block-bridge: {} + {}{}".format(
                        etype, b1, b2, suffix),
                    support=1,
                ))

                # 移除 HOG 级图中对应边
                for h1 in self._blocks.get(b1, []):
                    for h2 in self._blocks.get(b2, []):
                        if self._graph.has_edge(h1, h2):
                            self.remove_edge_color(h1, h2, color)

                # 移除块级图中的边
                bg.remove_edge(b1, b2)
                events_found += 1

        return events_found

    # ==================== 桥接冲突检测 (EEJ/NCF) ====================

    def resolve_bridge_events(self, outgroup_adjacency: Set = None):
        """检测桥接冲突: unique edge 连接两个共享边连通分量。

        外类群投票区分 fission vs fusion (EEJ/NCF):
        - 外类群 (parent HOG level) 有 (h1, h2) → 祖先有该连接
          → 保存的孩子有边，丢失的孩子没有 → fission on 丢失方的分支
        - 外类群无 (h1, h2) → 祖先无该连接
          → 有边的孩子是融合方 → EEJ/NCF on 有边方的分支
        - 无外类群 (root): 记录为 bridge_unclassified (需上游处理)

        算法:
        1. 建共享边图 (edges with >1 color → ≥2 children)
        2. 找连通分量
        3. 对每条 unique edge (1 color):
           两端在不同分量 → 桥接冲突
           → 用外类群投票分类 → 移除并记录事件
        """
        # Build shared-edge graph
        shared_G = nx.Graph()
        for h1, h2, data in self._graph.edges(data=True):
            if len(data['colors']) > 1:  # shared by ≥2 children
                shared_G.add_edge(h1, h2)

        if shared_G.number_of_edges() == 0:
            return

        # Find connected components in shared graph
        components = list(nx.connected_components(shared_G))
        hog_to_comp = {}
        for ci, comp in enumerate(components):
            for h in comp:
                hog_to_comp[h] = ci

        # Check each unique edge for bridge pattern
        cons_tels = self.consensus_telomeres(min_children=2)
        for h1, h2, data in list(self._graph.edges(data=True)):
            if len(data['colors']) != 1:
                continue  # not a unique edge
            c1 = hog_to_comp.get(h1)
            c2 = hog_to_comp.get(h2)
            if c1 is not None and c2 is not None and c1 != c2:
                # 端粒保护：不移除涉及共识端粒 HOG 的边
                if h1 in cons_tels or h2 in cons_tels:
                    continue
                # Bridge: connects two different shared components
                color = next(iter(data['colors']))
                child_id = color[0]
                comp1_size = len(components[c1])
                comp2_size = len(components[c2])

                # === 外类群投票分类 ===
                if outgroup_adjacency is not None:
                    # 将 h1, h2 映射到 parent HOG level 用于外类群查询
                    h1_parent = getattr(h1, 'parent', h1.hog_id if hasattr(h1, 'hog_id') else str(h1))
                    h2_parent = getattr(h2, 'parent', h2.hog_id if hasattr(h2, 'hog_id') else str(h2))
                    if isinstance(h1_parent, tuple) or isinstance(h1_parent, (int, float)):
                        # HOGrecord 的 parent 是 hog_id string
                        h1_pid = str(h1_parent) if not isinstance(h1_parent, str) else h1_parent
                        h2_pid = str(h2_parent) if not isinstance(h2_parent, str) else h2_parent
                    else:
                        h1_pid = h1_parent
                        h2_pid = h2_parent
                    key = (h1_pid, h2_pid) if h1_pid < h2_pid else (h2_pid, h1_pid)

                    if key in outgroup_adjacency:
                        # 外类群有该邻接 → 祖先有连接 → 另一孩子丢失了它 → fission
                        etype = 'fission'
                        branch_suffix = f" (outgroup confirms {child_id} preserved)"
                    else:
                        # 外类群无该邻接 → 祖先无连接 → 当前孩子融合了它 → EEJ/NCF
                        # NCF: 一端在共享图的端粒位置
                        in_shared_g = h1 in shared_G and h2 in shared_G
                        deg1 = False
                        if in_shared_g:
                            deg1 = shared_G.degree(h1) == 1 or shared_G.degree(h2) == 1
                        if deg1 or min(comp1_size, comp2_size) < 10:
                            etype = 'ncf'
                        else:
                            etype = 'eej'
                        branch_suffix = f" (outgroup: derived in {child_id})"
                else:
                    # 无外类群 (root) → 记录为 unclassified_bridge
                    etype = 'bridge_unclassified'
                    branch_suffix = " (no outgroup, root node)"

                self.events.append(TAKREvent(
                    event_type=etype,
                    branch=f"{self.hog_level}-{child_id}",
                    genes_involved=[h1, h2],
                    desc=f"{etype} bridge: comp{c1} ({comp1_size} HOGs)"
                         f" + comp{c2} ({comp2_size} HOGs)"
                         f" via {child_id} in {color}{branch_suffix}",
                    support=1,
                ))
                self.remove_edge_color(h1, h2, color)

    # ==================== Pipeline postcondition checks ====================

    def _postcondition_linear(self, label: str):
        """Phase 1 postcondition: 每个节点 degree ≤ 2，无并行边。"""
        violations = []
        for n in self._graph.nodes():
            deg = self._graph.degree(n)
            if deg > 2:
                violations.append((n, deg))
        if violations:
            worst = max(violations, key=lambda x: x[1])
            raise RuntimeError(
                f"{label} postcondition failed: {len(violations)} nodes with degree>2, "
                f"worst={worst[0]} deg={worst[1]}")

    def _postcondition_no_cycles(self, label: str):
        """Phase 2/4d postcondition: 图中无环。"""
        cycles = self.find_cycles()
        if cycles:
            raise RuntimeError(
                f"{label} postcondition failed: {len(cycles)} cycles remain, "
                f"first={cycles[0][:5]}...")

    def _postcondition_no_cross_edges(self, label: str):
        """Phase 2 postcondition: 无跨越边（spanning edges）。"""
        spanning = self.find_spanning_edges()
        if spanning:
            raise RuntimeError(
                f"{label} postcondition failed: {len(spanning)} spanning edges remain")

    def _postcondition_all_hogs_assigned(self, label: str):
        """Phase 3 postcondition: 所有 HOG 已分配到块。"""
        if not hasattr(self, '_hog_to_block'):
            raise RuntimeError(f"{label}: _hog_to_block not built")
        unassigned = self.all_hogs() - set(self._hog_to_block.keys())
        if unassigned:
            raise RuntimeError(
                f"{label} postcondition failed: {len(unassigned)} HOGs not in any block")

    def _postcondition_no_bridge_edges(self, label: str):
        """Phase 4f postcondition: 无桥接边。"""
        shared_G = nx.Graph()
        for h1, h2, data in self._graph.edges(data=True):
            if len(data['colors']) > 1:
                shared_G.add_edge(h1, h2)
        if shared_G.number_of_edges() == 0:
            return
        components = list(nx.connected_components(shared_G))
        hog_to_comp = {}
        for ci, comp in enumerate(components):
            for h in comp:
                hog_to_comp[h] = ci
        cons_tels = self.consensus_telomeres(min_children=2)
        bridges = []
        for h1, h2, data in self._graph.edges(data=True):
            if len(data['colors']) != 1:
                continue
            c1 = hog_to_comp.get(h1)
            c2 = hog_to_comp.get(h2)
            if c1 is not None and c2 is not None and c1 != c2:
                if h1 not in cons_tels and h2 not in cons_tels:
                    bridges.append((h1, h2))
        if bridges:
            raise RuntimeError(
                f"{label} postcondition failed: {len(bridges)} bridge edges remain")

    def _postcondition_all_hogs_covered(self, label: str):
        """Phase 5 postcondition: 所有 HOG 被路径覆盖。"""
        # 需要在 path_cover 之后调用，检查 paths 覆盖所有 HOG
        # 由 resolve_all_events 中直接检查

    # ==================== Phase 1: per-child deduplication ====================

    def _deduplicate_children(self, child_graphs: list, child_source_ids: list):
        """Phase 1: 在合图前对每个孩子单独去重。

        处理: tandem_dup (连续重复节点合并), 并行边合并。
        要求: 每个孩子图处理后为线性染色体 (deg<=2, 无并行边)。
        """
        deduped_graphs = []
        for cg, cid in zip(child_graphs, child_source_ids):
            new_cg = self._deduplicate_single_child(cg, cid)
            deduped_graphs.append(new_cg)
        return deduped_graphs

    def _deduplicate_single_child(self, child_graph, source_id: str):
        """对单个孩子图去重：合并连续重复节点，去除并行边。"""
        import copy
        new_graph = copy.deepcopy(child_graph)

        for chrom_idx, chrom_nodes in enumerate(new_graph.chromosomes):
            hogs = [n for n in chrom_nodes if n not in new_graph.telomeres]
            # 合并连续重复节点: A-A-B → A-B
            deduped = []
            for h in hogs:
                if deduped and deduped[-1] == h:
                    logger.debug("  [dedup] %s chrom%d: merge tandem %s",
                                 source_id, chrom_idx, h)
                    continue
                deduped.append(h)

            if len(deduped) != len(hogs):
                # 更新 chromosome 节点列表
                new_chrom = []
                for n in chrom_nodes:
                    if n in new_graph.telomeres:
                        new_chrom.append(n)
                    elif deduped:
                        new_chrom.append(deduped.pop(0))
                new_graph.chromosomes[chrom_idx] = new_chrom

        return new_graph

    def resolve_all_events(self, outgroups: Dict = None,
                           outgroup_adjacency: Set = None,
                           min_hogs: int = 3) -> List[TAKREvent]:
        """Telomere-centric event resolution pipeline.

        Pipeline (简单→复杂，每步检查 postcondition，失败则报错停止):
        ─────────────────────────────────────────────────────────────
        1. 每个孩子内部: duplication resolve (tandem, dispersed, proximal, seg_dup)
           Postcondition: 每个孩子图线性 (deg<=2, no parallel edges)
        2. 合图 + 单基因 indel/loss/gain (HOG 级)
           Postcondition: 无跨越边, 无环, 无分支
        3. 共线性块压缩
           Postcondition: 所有 HOG 已分配到块
        4a. seg_deletion / seg_insertion (多基因 indel)
           Postcondition: 无大段不对称
        4b. unidir_trans (单向转移)
           Postcondition: 无单向转移模式
        4c. telomere_inversion (端粒倒位，优先!)
           Postcondition: 端粒位置一致
        4d. internal_inversion (内部倒位)
           Postcondition: 块级图无环
        4e. RT / URT (相互易位)
           Postcondition: 无跨染色体环
        4f. EEJ / NCF / fission (桥接事件，最后)
           Postcondition: 无桥接边
        5. 端粒约束路径覆盖
           Postcondition: 所有 HOG 覆盖, 每条染色体 2 个端粒
        ─────────────────────────────────────────────────────────────
        """
        logger.info("  [colored] resolve_all_events for %s: %d nodes, %d edges",
                     self.hog_level, self.node_count(), self.edge_count())

        n_events_before = len(self.events)

        # Phase 1 dedup 在 orchestrator 中 add_child 前完成
        logger.info("  [colored] Phase 1 (dedup): done in orchestrator, %d children",
                     len(self._child_chrom_counts))

        # ====== Phase 2: 单基因 indel/loss/gain ======
        n0, e0 = self.node_count(), self.edge_count()
        self.resolve_indels(outgroups)
        n1, e1 = self.node_count(), self.edge_count()
        logger.info("  [colored] Phase 2 (indels): %d events, nodes %d→%d (-%d), edges %d→%d (-%d)",
                     len(self.events) - n_events_before,
                     n0, n1, n0 - n1, e0, e1, e0 - e1)
        try:
            self._postcondition_no_cross_edges("Phase 2")
        except RuntimeError as e:
            logger.warning("  [colored] %s", e)

        # ====== Phase 3: 共线性块压缩 ======
        self._build_synteny_blocks()
        self._compress_to_block_level()
        logger.info("  [colored] Phase 3 (blocks): %d blocks (%.1f HOG/block avg)",
                     len(self._blocks),
                     self.node_count() / max(len(self._blocks), 1))
        try:
            self._postcondition_all_hogs_assigned("Phase 3")
        except RuntimeError as e:
            logger.warning("  [colored] %s", e)

        # ====== Phase 4: 块级事件 (按复杂度递增) ======

        # Phase 4b: unidir_trans
        n_ut = self._resolve_unidir_trans(outgroup_adjacency=outgroup_adjacency)
        if n_ut:
            logger.info("  [colored] Phase 4b (unidir_trans): %d events", n_ut)

        # ====== Phase 4c-4e: 块级结构重排 (inversion, RT) ======
        n_before_struct = len(self.events)
        n2, e2 = self.node_count(), self.edge_count()
        n_struct = self._resolve_block_structural_events(
            outgroup_adjacency=outgroup_adjacency)
        n3, e3 = self.node_count(), self.edge_count()
        logger.info("  [colored] Phase 4c-4e (structural): %d events, nodes %d→%d (-%d), edges %d→%d (-%d)",
                     len(self.events) - n_before_struct,
                     n2, n3, n2 - n3, e2, e3, e2 - e3)

        # ====== Phase 4f: 块级桥接 (EEJ/NCF/fission) ======
        n_before_bridge = len(self.events)
        n4, e4 = self.node_count(), self.edge_count()
        n_bridge = self._resolve_block_bridge_events(
            outgroup_adjacency=outgroup_adjacency)
        n5, e5 = self.node_count(), self.edge_count()
        logger.info("  [colored] Phase 4f (bridge): %d events, nodes %d→%d (-%d), edges %d→%d (-%d)",
                     len(self.events) - n_before_bridge,
                     n4, n5, n4 - n5, e4, e5, e4 - e5)
        try:
            self._postcondition_no_bridge_edges("Phase 4f")
        except RuntimeError as e:
            logger.warning("  [colored] %s", e)

        # ====== Phase 5: 端粒约束路径覆盖 ======
        paths = self.path_cover()
        self._cached_paths = paths  # 缓存给 to_ancestral_graph 复用
        all_hogs = self.all_hogs()
        covered = set()
        for p in paths:
            covered.update(p)
        uncovered = all_hogs - covered
        if uncovered:
            logger.warning("  [colored] Phase 5: %d HOGs not covered by path_cover",
                           len(uncovered))
        n_chrom = len(paths)
        logger.info("  [colored] Phase 5 (path cover): %d chromosomes", n_chrom)

        # 过滤小事件
        kept = []
        for e in self.events:
            if e.event_type in ('eej', 'ncf', 'fission', 'bridge_unclassified',
                                'reciprocal_translocation',
                                'unbalanced_reciprocal_translocation',
                                'inversion', 'telomere_inversion'):
                kept.append(e)
            elif len(e.genes_involved) >= min_hogs:
                kept.append(e)
        self.events = kept

        from collections import Counter
        type_counts = Counter(e.event_type for e in self.events)
        type_str = ', '.join(f"{t}={c}" for t, c in sorted(type_counts.items()))
        logger.info("  [colored] done: %d events %s (after min_hogs=%d), %d chroms",
                     len(self.events), type_str, min_hogs, n_chrom)

        # === 染色体数一致性检查 (按孩子分支) ===
        # 每个孩子独立推算祖先染色体数:
        #   祖先_chroms = child_chroms - fissions + fusions
        # 所有孩子推算结果应一致
        new_events = self.events[n_events_before:]
        # 按分支分组事件
        from collections import defaultdict
        branch_events = defaultdict(lambda: {'fission': 0, 'fusion': 0})
        for e in new_events:
            # branch 格式: "N1-Sp_1"
            branch = e.branch
            if e.event_type == 'fission':
                branch_events[branch]['fission'] += 1
            elif e.event_type in ('eej', 'ncf'):
                branch_events[branch]['fusion'] += 1

        inferred_ancestor_chroms = {}
        for cid, child_chroms in self._child_chrom_counts.items():
            branch_key = f"{self.hog_level}-{cid}"
            ev = branch_events.get(branch_key, {'fission': 0, 'fusion': 0})
            # 祖先 = 孩子 - fissions + fusions
            anc_chroms = child_chroms - ev['fission'] + ev['fusion']
            inferred_ancestor_chroms[cid] = anc_chroms
            logger.info("  [chrom] %s: %d chroms, fission=%d, fusion=%d → ancestor=%d",
                        cid, child_chroms, ev['fission'], ev['fusion'], anc_chroms)

        # 检查所有孩子推算的祖先染色体数是否一致
        values = list(inferred_ancestor_chroms.values())
        if len(set(values)) > 1:
            logger.error("  [chrom] INCONSISTENT ancestor chroms: %s", inferred_ancestor_chroms)
        else:
            inferred = values[0]
            logger.info("  [chrom] consistent: ancestor inferred %d chromosomes", inferred)
            # 对比 path_cover 实际结果
            if inferred != n_chrom:
                logger.error("  [chrom] MISMATCH: inferred=%d, actual path_cover=%d",
                             inferred, n_chrom)
            else:
                logger.info("  [chrom] OK: inferred=%d == path_cover=%d",
                            inferred, n_chrom)

        return self.events

    # ==================== Phase 4a/4b: seg events & unidir_trans ====================

    def _resolve_seg_events(self, outgroup_adjacency: Set = None) -> int:
        """Phase 4a: 检测 seg_deletion / seg_insertion。

        模式: 一个孩子有连续 N 个块 (N>=3)，另一个孩子完全缺失。
        在块级图中检测: 共享边路径中，某个孩子缺失一段连续块。
        """
        if not hasattr(self, '_blocks'):
            return 0

        events_found = 0
        # 对每个孩子，检查其块序列中是否有大段缺失
        for cid in self.children():
            child_blocks = set()
            for bid, hogs in self._blocks.items():
                for h in hogs:
                    if h in self._child_hogs.get(cid, set()):
                        child_blocks.add(bid)
                        break

            # 检查共享边路径中缺失的连续块
            # 共享边路径 = 所有孩子共有的块序列
            all_blocks = set(self._blocks.keys())
            missing = all_blocks - child_blocks
            if not missing:
                continue

            # 检查缺失块是否在共享边路径中连续
            # 简化: 标记缺失块为 seg_deletion (需要外类群确认)
            if len(missing) >= 3:
                # 检查缺失块是否形成连续路径
                bg = self._block_graph
                missing_sub = bg.subgraph(missing) if missing else None
                if missing_sub and nx.is_connected(missing_sub):
                    involved_hogs = []
                    for bid in missing:
                        involved_hogs.extend(self._blocks.get(bid, []))

                    etype = 'seg_deletion'
                    self.events.append(TAKREvent(
                        event_type=etype,
                        branch=self.hog_level,
                        genes_involved=involved_hogs,
                        desc=f"seg_deletion: {len(missing)} blocks missing in {cid}",
                        support=1,
                    ))
                    events_found += 1

        return events_found

    def _resolve_unidir_trans(self, outgroup_adjacency: Set = None) -> int:
        """Phase 4b: 检测 unidirectional translocation。

        模式: A 有 a-b-c-d，B 有 a-x-d (b,c 被 x 替换，x 只在 B 中)。
        在块级图中: 两个孩子共享 a 和 d，但中间路径不同。
        """
        # 简化实现: 检查块级图中的非循环结构差异
        # 完整实现需要对每个孩子的块序列做比对
        # 当前返回 0，后续补充
        return 0

    # ==================== 可视化 ====================

    def draw_block_graph(self, outpath: str, dpi: int = 200,
                         title: str = None, show_hogs: bool = False):
        """Draw block-level colored adjacency graph using Graphviz.

        Nodes = synteny blocks (sized by HOG count).
        Edges:
          - Shared (multi-child): thick black solid
          - Unique (single-child): thin colored solid, one color per child
        Telomere-adjacent blocks: double-octagon shape, orange fill.

        Args:
            outpath: Output PNG path
            dpi: Resolution
            title: Plot title
            show_hogs: If True, show individual HOG names in block labels
        """
        try:
            import graphviz
            # Check if dot binary is available
            import shutil
            if not shutil.which('dot'):
                raise ImportError('dot binary not found')
        except (ImportError, Exception):
            logger.debug("graphviz dot not available, using matplotlib fallback")
            self._draw_block_graph_mpl(outpath, dpi, title)
            return

        if not hasattr(self, '_block_graph') or not self._block_graph:
            self._build_synteny_blocks()
            self._compress_to_block_level()

        bg = self._block_graph
        if bg.number_of_nodes() == 0:
            return

        # Build telomere block set
        telomere_blocks = set()
        for cid, tels in self._child_telomeres.items():
            for tel_hog in tels:
                bid = self._hog_to_block.get(tel_hog)
                if bid:
                    telomere_blocks.add(bid)

        # Cycle nodes
        cycles = []
        try:
            cycles = nx.cycle_basis(bg)
        except Exception:
            pass
        cycle_nodes = set()
        for cyc in cycles:
            cycle_nodes.update(cyc)

        # Child color palette (saturated, distinct)
        children = sorted(self.children())
        palette = ['#c0392b', '#2980b9', '#27ae60', '#f39c12',
                   '#8e44ad', '#16a085', '#d35400', '#2c3e50']
        child_colors = {c: palette[i % len(palette)]
                        for i, c in enumerate(children)}

        # Build graphviz DOT
        dot = graphviz.Digraph(format='png')
        dot.attr(rankdir='LR', label=title or f'Block Graph: {self.hog_level}',
                 labelloc='t', fontsize='18', fontname='Helvetica-Bold',
                 fontcolor='#1a1a2e', bgcolor='white', pad='0.5',
                 nodesep='0.4', ranksep='0.6')
        dot.attr('node', fontname='Helvetica-Bold', fontsize='11',
                 fontcolor='#1a1a2e', margin='0.15,0.08')
        dot.attr('edge', fontname='Helvetica', fontsize='8')

        # Nodes
        for n in bg.nodes():
            hog_list = self._blocks.get(n, [])
            hog_count = len(hog_list)
            label = f'{n}\n({hog_count})'
            if show_hogs and hog_count <= 5:
                label += '\n' + ', '.join(str(h) for h in hog_list[:5])

            attrs = {'label': label,
                     'width': str(max(0.4, hog_count * 0.03)),
                     'height': str(max(0.3, hog_count * 0.02)),
                     'fontcolor': '#1a1a2e', 'fontsize': '11'}

            if n in telomere_blocks:
                attrs['shape'] = 'doubleoctagon'
                attrs['style'] = 'filled,bold'
                attrs['fillcolor'] = '#fff3cd'  # light yellow
                attrs['color'] = '#e67e22'       # orange border
                attrs['penwidth'] = '2.5'
                attrs['fontcolor'] = '#7d5a00'   # dark gold text
            elif n in cycle_nodes:
                attrs['shape'] = 'octagon'
                attrs['style'] = 'filled,bold'
                attrs['fillcolor'] = '#fce4ec'   # light pink
                attrs['color'] = '#c62828'       # dark red border
                attrs['penwidth'] = '2.5'
                attrs['fontcolor'] = '#b71c1c'   # dark red text
            else:
                attrs['shape'] = 'box'
                attrs['style'] = 'filled,rounded'
                attrs['fillcolor'] = '#e3f2fd'   # light blue
                attrs['color'] = '#1565c0'       # blue border
                attrs['penwidth'] = '1.5'

            dot.node(str(n), **attrs)

        # Edges
        for h1, h2, data in bg.edges(data=True):
            colors = data.get('colors', set())
            child_ids = sorted(set(c for c, _ in colors))

            if len(child_ids) > 1:
                # Shared edge: dark navy, very thick
                dot.edge(str(h1), str(h2),
                         color='#1a1a2e', penwidth='3.5',
                         style='solid')
            elif len(child_ids) == 1:
                cid = child_ids[0]
                color = child_colors.get(cid, '#7f8c8d')
                dot.edge(str(h1), str(h2),
                         color=color, penwidth='2.0',
                         style='solid')

        # Legend as subgraph
        with dot.subgraph(name='cluster_legend') as legend:
            legend.attr(label='Legend', style='dashed', fontsize='12',
                        fontcolor='#1a1a2e', bgcolor='#f5f5f5')
            legend.node('_lg_shared', 'shared\n(multi-child)',
                        shape='box', style='filled,rounded',
                        fillcolor='#e3f2fd', color='#1565c0',
                        fontcolor='#1a1a2e')
            legend.node('_lg_cycle', 'in cycle',
                        shape='octagon', style='filled,bold',
                        fillcolor='#fce4ec', color='#c62828',
                        fontcolor='#b71c1c')
            legend.node('_lg_tel', 'telomere',
                        shape='doubleoctagon', style='filled,bold',
                        fillcolor='#fff3cd', color='#e67e22',
                        fontcolor='#7d5a00')
            for cid in children:
                color = child_colors.get(cid, '#7f8c8d')
                legend.node(f'_lg_{cid}', f'{cid}\n(unique)',
                            shape='box', style='filled,rounded',
                            fillcolor=color, fontcolor='white')
            legend.edge('_lg_shared', '_lg_cycle', style='invis')

        # Render
        base = outpath.rsplit('.', 1)[0]
        dot.render(base, cleanup=True)
        logger.info("  [viz] block graph saved to %s", base + '.png')

    def _draw_block_graph_mpl(self, outpath: str, dpi: int = 200,
                              title: str = None):
        """Fallback matplotlib block graph (original implementation)."""
        try:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
        except ImportError:
            logger.warning("matplotlib not available, skipping block graph")
            return

        if not hasattr(self, '_block_graph') or not self._block_graph:
            self._build_synteny_blocks()
            self._compress_to_block_level()

        bg = self._block_graph
        if bg.number_of_nodes() == 0:
            return

        fig, ax = plt.subplots(1, 1, figsize=(12, 10))
        try:
            pos = nx.kamada_kawai_layout(bg)
        except Exception:
            pos = nx.spring_layout(bg, k=2.0, seed=42)

        node_sizes = []
        for n in bg.nodes():
            hog_count = len(self._blocks.get(n, []))
            node_sizes.append(max(200, hog_count * 15))

        cycles = []
        try:
            cycles = nx.cycle_basis(bg)
        except Exception:
            pass
        cycle_nodes = set()
        for cyc in cycles:
            cycle_nodes.update(cyc)

        # Telomere blocks
        telomere_blocks = set()
        for cid, tels in self._child_telomeres.items():
            for tel_hog in tels:
                bid = self._hog_to_block.get(tel_hog)
                if bid:
                    telomere_blocks.add(bid)

        node_colors = []
        node_edge_colors = []
        for n in bg.nodes():
            if n in telomere_blocks:
                node_colors.append('#fff3cd')   # light yellow
                node_edge_colors.append('#e67e22')
            elif n in cycle_nodes:
                node_colors.append('#fce4ec')   # light pink
                node_edge_colors.append('#c62828')
            else:
                node_colors.append('#e3f2fd')   # light blue
                node_edge_colors.append('#1565c0')

        nx.draw_networkx_nodes(bg, pos, ax=ax, node_size=node_sizes,
                               node_color=node_colors, alpha=0.95,
                               edgecolors=node_edge_colors, linewidths=2.0)

        labels = {n: f'{n}\n({len(self._blocks.get(n, []))})' for n in bg.nodes()}
        nx.draw_networkx_labels(bg, pos, labels, ax=ax, font_size=8,
                                font_color='#1a1a2e', font_weight='bold')

        children = sorted(self.children())
        palette = ['#c0392b', '#2980b9', '#27ae60', '#f39c12',
                   '#8e44ad', '#16a085', '#d35400', '#2c3e50']
        child_colors = {c: palette[i % len(palette)]
                        for i, c in enumerate(children)}

        for h1, h2, data in bg.edges(data=True):
            colors = data.get('colors', set())
            child_ids = sorted(set(c for c, _ in colors))
            if len(child_ids) > 1:
                ax.plot([pos[h1][0], pos[h2][0]], [pos[h1][1], pos[h2][1]],
                        color='#1a1a2e', linewidth=3.5, alpha=0.9, zorder=1)
            elif len(child_ids) == 1:
                cid = child_ids[0]
                color = child_colors.get(cid, '#7f8c8d')
                ax.plot([pos[h1][0], pos[h2][0]], [pos[h1][1], pos[h2][1]],
                        color=color, linewidth=2.0, alpha=0.8, zorder=0)

        legend_handles = [
            plt.Line2D([0], [0], marker='o', color='w',
                       markerfacecolor='#fff3cd', markersize=12,
                       markeredgecolor='#e67e22', markeredgewidth=2,
                       label='telomere'),
            plt.Line2D([0], [0], marker='o', color='w',
                       markerfacecolor='#fce4ec', markersize=12,
                       markeredgecolor='#c62828', markeredgewidth=2,
                       label='in cycle'),
            plt.Line2D([0], [0], color='#1a1a2e', linewidth=3.5,
                       label='shared (multi-child)'),
        ]
        for cid in children:
            color = child_colors.get(cid, '#7f8c8d')
            legend_handles.append(
                plt.Line2D([0], [0], color=color, linewidth=2.5,
                           label=f'{cid} (unique)'))

        ax.legend(handles=legend_handles, loc='upper left', fontsize=8)
        ax.set_title(title or f'Block Graph: {self.hog_level}',
                     fontsize=12, fontweight='bold')
        ax.axis('off')
        fig.tight_layout()
        fig.savefig(outpath, dpi=dpi, bbox_inches='tight')
        plt.close(fig)
        logger.info("  [viz] block graph saved to %s", outpath)

    def draw_adjacency_heatmap(self, outpath: str, dpi: int = 200,
                               title: str = None):
        """Draw adjacency matrix heatmap: rows/cols = HOGs, cells = edge colors.

        Useful for seeing the global structure of shared vs unique adjacencies.
        """
        try:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
            import numpy as np
        except ImportError:
            logger.warning("matplotlib not available, skipping heatmap")
            return

        # Order HOGs by block, then by position within block
        if not hasattr(self, '_blocks'):
            self._build_synteny_blocks()

        hog_order = []
        for bid in sorted(self._blocks.keys()):
            hog_order.extend(self._blocks[bid])

        # Add HOGs not in any block
        all_hogs = set()
        for n in self._graph.nodes():
            all_hogs.add(n)
        for h in sorted(all_hogs, key=str):
            if h not in set(hog_order):
                hog_order.append(h)

        n = len(hog_order)
        if n == 0 or n > 500:
            logger.info("  [viz] skipping heatmap: %d HOGs (too many or zero)", n)
            return

        hog_idx = {h: i for i, h in enumerate(hog_order)}

        # Build color matrix: 0=none, 1=unique, 2=shared
        mat = np.zeros((n, n), dtype=int)
        children = sorted(self.children())
        child_list = list(children)

        for h1, h2, data in self._graph.edges(data=True):
            if h1 in hog_idx and h2 in hog_idx:
                i, j = hog_idx[h1], hog_idx[h2]
                colors = data.get('colors', set())
                if len(colors) > 1:
                    mat[i][j] = 2
                    mat[j][i] = 2
                else:
                    mat[i][j] = 1
                    mat[j][i] = 1

        fig, ax = plt.subplots(1, 1, figsize=(10, 10))
        cmap = matplotlib.colors.ListedColormap(['#ffffff', '#4ecdc4', '#ff6b6b'])
        ax.imshow(mat, cmap=cmap, interpolation='nearest', aspect='equal')

        # Block boundaries
        offset = 0
        for bid in sorted(self._blocks.keys()):
            blen = len(self._blocks[bid])
            ax.axhline(y=offset + blen - 0.5, color='gray', linewidth=0.5)
            ax.axvline(x=offset + blen - 0.5, color='gray', linewidth=0.5)
            offset += blen

        ax.set_title(title or f'Adjacency Matrix: {self.hog_level}',
                     fontsize=12, fontweight='bold')
        ax.set_xlabel('HOG index')
        ax.set_ylabel('HOG index')

        legend_handles = [
            plt.Line2D([0], [0], marker='s', color='w',
                       markerfacecolor='#ffffff', markersize=10, label='no edge'),
            plt.Line2D([0], [0], marker='s', color='w',
                       markerfacecolor='#4ecdc4', markersize=10, label='unique'),
            plt.Line2D([0], [0], marker='s', color='w',
                       markerfacecolor='#ff6b6b', markersize=10, label='shared'),
        ]
        ax.legend(handles=legend_handles, loc='upper right', fontsize=8)

        fig.tight_layout()
        fig.savefig(outpath, dpi=dpi, bbox_inches='tight')
        plt.close(fig)
        logger.info("  [viz] adjacency heatmap saved to %s", outpath)

    def __repr__(self):
        return (f"ColoredGraph(hog_level={self.hog_level!r}, "
                f"nodes={self.node_count()}, edges={self.edge_count()}, "
                f"events={len(self.events)})")
