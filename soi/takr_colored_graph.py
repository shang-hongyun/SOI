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
        """添加一条边，记录颜色 (child_id, chrom_idx)。
        如果边已存在，只把颜色加入已有边的 colors 集。
        """
        if self._graph.has_edge(h1, h2):
            colors = self._graph[h1][h2]['colors']
            colors.add((child_id, chrom_idx))
        else:
            self._graph.add_edge(h1, h2, colors={(child_id, chrom_idx)})

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
        """检测并解决 indel/loss 冲突。

        外类群投票确定方向:
        - 外类群有 HOG(spanning) → 祖先有这些 HOG → loss in this child
        - 外类群无 HOG(spanning) → 祖先无这些 HOG → gain in this child
        """
        while True:
            spanning = self.find_spanning_edges()
            if not spanning:
                break
            for h1, h2, child_id, spanned_hogs in spanning:
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

    # ==================== 环检测 + 事件分类 ====================

    def find_cycles(self) -> List[List]:
        """nx.cycle_basis 包装。"""
        try:
            return nx.cycle_basis(self._graph)
        except Exception:
            return []

    def classify_cycle(self, cycle: List) -> Tuple[Optional[str], List, Dict]:
        """分析环的颜色模式，判断事件类型。

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

            if alternating and n == 4:
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
                    # 2 种颜色 → 同一 child 的染色体内某段倒置 → inversion
                    # 只需移除一个 child 的冲突边（另一组的边保留）
                    # 选 color_counts 中出现次数较少的颜色组移除
                    min_count = min(color_counts.values())
                    rare_colors = {c for c, cnt in color_counts.items()
                                   if cnt == min_count}
                    # 如果两种颜色出现次数相同（完美交替），选第一个
                    if len(rare_colors) > 1 and list(color_counts.values()).count(min_count) == len(color_counts):
                        # 所有颜色出现次数相同 → 选任一种，这里选第一个
                        rare_colors = {next(iter(color_counts.keys()))}
                    conflict_edges = []
                    for (h1, h2), colors in zip(edges, edge_colors):
                        if rare_colors & colors:
                            conflict_edges.append((h1, h2))
                    info['conflict_colors'] = rare_colors
                    # 检查是否含端粒
                    has_tel = any(
                        isinstance(n, tuple) and len(n) == 2 and n[1] in ('L', 'R')
                        for n in cycle
                    )
                    etype = 'telomere_inversion' if has_tel else 'inversion'
                    return etype, conflict_edges, info

        # == 简单情况: 3-HOG 环 (不可能在染色体图中, 属于错误) ==
        if n == 3:
            return 'gene_indel', list(edges), info

        # 无法分类的环 → 移除最少出现的颜色的边
        min_count = min(color_counts.values())
        rare_colors = {c for c, cnt in color_counts.items() if cnt == min_count}
        conflict_edges = []
        for (h1, h2), colors in zip(edges, edge_colors):
            if rare_colors & colors:
                conflict_edges.append((h1, h2))
        info['fallback'] = True
        info['fallback_reason'] = 'unclassified_cycle'
        return 'unclassified', conflict_edges, info

    # ==================== 结构重排 resolved ====================

    def resolve_structural_events(self):
        """循环: find_cycles → classify_cycle → 移除冲突边。"""
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

                # 找到冲突边的颜色 (来自最少出现的 child)
                # 对于 inversion/RT: 移除单色边
                for h1, h2 in conflict_edges:
                    colors = self.get_colors(h1, h2)
                    if len(colors) == 1:
                        self.remove_edge_color(h1, h2, next(iter(colors)))
                    else:
                        # 共享边 → 移除最少出现的 child 的颜色
                        # 统计各 child 在环中的出现次数
                        child_count = {}
                        for i_c in range(len(cycle)):
                            u, v = cycle[i_c], cycle[(i_c + 1) % len(cycle)]
                            for cid, chrom in self.get_colors(u, v):
                                child_count[(cid, chrom)] = child_count.get((cid, chrom), 0) + 1
                        # 找到出现最少的颜色
                        min_child = min(child_count, key=child_count.get)
                        if min_child in colors:
                            self.remove_edge_color(h1, h2, min_child)
                        else:
                            # 该边没有该 child 的颜色 → 移除任何一个
                            self.remove_edge_color(h1, h2, next(iter(colors)))

                # 记录事件
                event = TAKREvent(
                    event_type=etype,
                    branch=self.hog_level,
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

    # ==================== 路径覆盖 ====================

    def path_cover(self) -> List[List]:
        """端粒约束路径覆盖。

        使用当前图中的所有剩余边（所有颜色），
        从端粒节点出发沿唯一邻接行走。

        Returns:
            [[hog1, hog2, ...], ...] — 每条染色体一条路径
            路径中不包含 telomere 节点
        """
        paths = []
        visited = set()

        # 找端粒节点: 度为 1 的节点
        degree_counts = dict(self._graph.degree())
        telomeres = {n for n, deg in degree_counts.items()
                     if deg == 1 and n in self.all_hogs()}

        # 从每个端粒出发行走
        for start in telomeres:
            if start in visited:
                continue

            path = self._walk_path(start, visited | {start})
            if path:
                paths.append(path)
                visited.update(path)

        # 检查是否有未访问的连通分量（环）
        all_nodes = self.all_hogs()
        unvisited = all_nodes - visited
        if unvisited:
            # 对每个未访问的连通分量，选一个节点并切断最低支撑边
            G_rem = self._graph.subgraph(unvisited)
            for comp_nodes in nx.connected_components(G_rem):
                comp = list(comp_nodes)
                if len(comp) >= 2:
                    # 找到该分量中第一个有边可走的节点
                    start = comp[0]
                    path = self._walk_path_no_telomere(start)
                    if path:
                        paths.append(path)
                        visited.update(path)

        return paths

    def _walk_path(self, start, visited: Set) -> Optional[List]:
        """从 start 出发，沿唯一邻接行走直到另一个端粒或没有邻居。"""
        path = [start]
        curr = start
        while True:
            neighbors = [n for n in self._graph.neighbors(curr)
                         if n not in visited]
            if not neighbors:
                break
            if len(neighbors) > 1:
                # 分叉 → 选择边数最少的邻居（度最小的）
                neighbors.sort(key=lambda n: self._graph.degree(n))
            nxt = neighbors[0]
            visited.add(nxt)
            path.append(nxt)
            curr = nxt
            # 如果遇到端粒（度=1），停止
            if self._graph.degree(curr) == 1 and curr != start:
                break

        if len(path) < 2:
            return None
        # 过滤 telomere 节点
        path = [n for n in path if not (isinstance(n, tuple)
                                        and len(n) == 2 and n[1] in ('L', 'R'))]
        return path if len(path) >= 2 else None

    def _walk_path_no_telomere(self, start) -> Optional[List]:
        """从一个非端粒节点出发行走，用于处理环。"""
        path = [start]
        visited = {start}
        curr = start
        prev = None
        while True:
            neighbors = [n for n in self._graph.neighbors(curr)
                         if n != prev]
            if not neighbors:
                break
            nxt = neighbors[0]
            if nxt in visited:
                break
            visited.add(nxt)
            path.append(nxt)
            prev, curr = curr, nxt
        return path if len(path) >= 2 else None

    # ==================== 转换 ====================

    def to_ancestral_graph(self):
        """转换为 AncestralAdjacencyGraph。

        调用 path_cover() 获取染色体路径，然后构建有向图。
        """
        from .AK import AncestralAdjacencyGraph

        result = AncestralAdjacencyGraph(node_id=self.hog_level)

        paths = self.path_cover()
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

        # 添加所有孤立的 HOG
        for n in self._graph.nodes:
            if n in self.all_hogs() and n not in result.gene_nodes:
                result.graph.add_node(n)
                result.gene_nodes.add(n)

        result.events = self.events
        return result

    # ==================== 一键执行 ====================

    def collapse_wgd(self, ploidy: int):
        """ColoredGraph-based WGD collapse: post-WGD → pre-WGD.

        利用边颜色中的染色体来源 (child_id) 作为亚基因组标识：
        - 跨孩子共享的边 = pre-WGD 保守邻接（保留）
        - 单个孩子独有的边 = post-WGD 衍生（移除）
        然后路径覆盖 → pre-WGD 图。

        Args:
            ploidy: 倍性 (2=四倍体, 3=六倍体, ...)

        Returns:
            ColoredGraph — pre-WGD 图（已 resolve 所有冲突）
        """
        # Step 1: 统计每条边的总颜色数（跨全部孩子的 (child_id, chrom_idx) 对）
        # p=2 → 每条 pre-WGD 保守边应出现至少 2 次（两个亚基因组）
        # p=3 → 至少 3 次
        threshold = max(ploidy, 1)
        edge_color_counts = {}  # (h1, h2) -> total color count across all children
        for h1, h2, data in self._graph.edges(data=True):
            edge_color_counts[(h1, h2)] = len(data['colors'])

        # Step 2: 构建 pre-WGD 图（仅保留颜色数 >= ploidy 的边）
        pre_G = ColoredGraph(hog_level=self.hog_level + "_preWGD")

        # 先添加所有节点
        for n in self._graph.nodes:
            pre_G._graph.add_node(n)

        # 保留颜色数 >= ploidy 的边 → pre-WGD 保守邻接
        kept = 0
        removed = 0
        for (h1, h2), color_count in edge_color_counts.items():
            if h1 not in pre_G._graph or h2 not in pre_G._graph:
                continue
            if color_count >= threshold:
                # 保留
                for child_id, chrom_idx in self._graph[h1][h2]['colors']:
                    pre_G.add_edge(h1, h2, child_id, chrom_idx)
                kept += 1
            else:
                # 颜色数不足 ploidy → post-WGD 衍生
                removed += 1
                event_type = 'post_wgd_rearrangement'
                if len(pre_G.events) < 100:
                    pre_G.events.append(TAKREvent(
                        event_type=event_type,
                        branch=self.hog_level + "_preWGD-" + self.hog_level,
                        genes_involved=[h1, h2],
                        desc=f"{event_type}: {color_count} colors < {threshold}",
                        support=1,
                    ))

        logger.debug("  [collapse_wgd] %s: %d shared kept, %d unique removed (threshold=%d)",
                      self.hog_level, kept, removed, threshold)

        # Step 3: 路径覆盖 → 还原染色体
        # 如果 pre_G 为空图，回退到原图
        if pre_G.edge_count() == 0:
            logger.warning("  [collapse_wgd] no shared edges, using original graph")
            return self

        return pre_G

    def resolve_all_events(self, outgroups: Dict = None,
                           min_hogs: int = 3) -> List[TAKREvent]:
        """按优先级解决所有冲突。

        Pipeline:
        1. indel/loss resolve
        2. 结构重排 (循环: find_cycles → classify_cycle → 移除冲突边)
        3. 路径覆盖

        Returns: 所有检测到的事件
        """
        logger.info("  [colored] resolve_all_events for %s: %d nodes, %d edges",
                     self.hog_level, self.node_count(), self.edge_count())

        # Step 1: indel/loss
        self.resolve_indels(outgroups)
        logger.info("  [colored] after indel: %d nodes, %d edges, %d events",
                     self.node_count(), self.edge_count(), len(self.events))

        # Step 2: 结构重排
        self.resolve_structural_events()
        logger.info("  [colored] after structural: %d nodes, %d edges, %d events",
                     self.node_count(), self.edge_count(), len(self.events))

        # 过滤小事件 (min_hogs)
        self.events = [e for e in self.events
                       if len(e.genes_involved) >= min_hogs]
        logger.info("  [colored] done: %d events (after min_hogs=%d)",
                     len(self.events), min_hogs)

        return self.events

    def __repr__(self):
        return (f"ColoredGraph(hog_level={self.hog_level!r}, "
                f"nodes={self.node_count()}, edges={self.edge_count()}, "
                f"events={len(self.events)})")
