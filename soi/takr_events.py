#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""takr_events.py - Unified event types and tree parsing for TAKR framework.

This module provides:
  1. Canonical event type definitions shared by simulator, AKR, and eval
  2. Newick tree parsing utilities (wraps tree.py)
  3. Event serialization/deserialization helpers

Usage:
    from soi.takr_events import (
        EVENT_TYPES, CHROM_REARRANGEMENTS, CHROM_COUNT_EVENTS,
        canonicalize_event_type, is_chrom_rearrangement, chrom_count_delta,
        parse_tree, TAKREvent
    )

Event type taxonomy (重建需要精确区分的重排类型):
================================================================================

一、染色体重排 (Chromosomal Rearrangements) - 重建必须精确区分

1. Inversion (INV, 倒位)
   缩写: inv
   定义: 染色体内部片段翻转180度
   染色体数变化: 0
   
   子类型 (重建时必须区分，影响端粒-中心模型):
   a) Internal Inversion (内部倒位)
      - 两个断点都在染色体内部
      - 不影响端粒位置
      - 重建: 直接翻转片段即可
   
   b) Telomere Inversion (端粒倒位, TI)
      - 至少一个断点在端粒
      - 端粒位置会改变！影响后续的EEJ/NCF检测
      - 重建: 必须先恢复端粒位置，再检测大尺度重排
      - 检测方法: 倒位块包含端粒HOG (bstart==0 或 bend==len)
   
   c) Centromere-Spanning Inversion (着丝粒跨越倒位, CSI)
      - 倒位跨越着丝粒
      - 可能改变染色体臂的比例
      - 重建: 需要特殊处理着丝粒位置
   
   重建优先级: 端粒倒位 > 内部倒位 > 大尺度重排
   原因: 端粒倒位改变端粒位置，必须先恢复才能正确检测EEJ/NCF

2. Reciprocal Translocation (RT, 相互易位)
   缩写: rt
   定义: 两条染色体交换片段，双方都保留片段
   染色体数变化: 0
   重建注意: 与URT的关键区别是交换平衡

3. Unbalanced Reciprocal Translocation (URT, 不平衡相互易位)
   缩写: urt
   定义: 相互易位的子类型，其中一条断点发生在端粒，导致交换极不平衡
   染色体数变化: 0
   重建注意: 图结构可能与RT不同，需单独检测
   与RT区别: URT的一个断点在端粒，另一个在染色体内部

4. Unidirectional Translocation (UT, 单向易位)
   缩写: utrans
   定义: 片段从一条染色体移到另一条，原位置丢失
   染色体数变化: 0
   重建注意: 与RT/URT的区别是只有单向移动

5. Nested Chromosome Fusion (NCF, 嵌套染色体融合)
   缩写: ncf
   定义: 一条染色体插入到另一条染色体的内部
   染色体数变化: -1 (两条变一条)
   重建注意: donor两端都不是端粒，与EEJ区分

6. End-End Join (EEJ, 端粒-端粒连接)
   缩写: eej
   定义: 两条染色体的端粒直接连接
   染色体数变化: -1 (两条变一条)
   重建注意: 两个断点都是端粒，与NCF区分

7. Fission (FIS, 染色体断裂)
   缩写: fis
   定义: 一条染色体断裂为两条
   染色体数变化: +1 (一条变两条)
   重建注意: EEJ的逆事件

二、基因水平事件 (Gene-level Events) - 重建可统一处理

8. Gene Loss / Fractionation (基因丢失/去冗余)
   缩写: loss, frac
   定义: 基因从基因组中丢失
   染色体数变化: 0
   重建注意: 重建时不区分loss和fractionation，统一视为缺失
   模拟器区分: loss=随机丢失, fractionation=WGD后的去冗余

9. Gene Gain (基因获得)
   缩写: gain
   定义: 新基因加入基因组
   染色体数变化: 0
   重建注意: 与外类群投票结合确定

10. Tandem Duplication (TDUP, 串联重复)
    缩写: tdup
    定义: 相邻位置产生基因拷贝
    染色体数变化: 0

11. Dispersed Duplication (DDUP, 分散重复)
    缩写: ddup
    定义: 远距离位置产生基因拷贝
    染色体数变化: 0

12. Segmental Duplication (SDUP, 片段重复)
    缩写: segdup
    定义: 大片段重复
    染色体数变化: 0

13. Segmental Deletion (SDEL, 片段删除)
    缩写: segdel
    定义: 大片段删除
    染色体数变化: 0

三、特殊事件

14. Whole Genome Duplication (WGD, 全基因组加倍)
    缩写: wgd
    定义: 整个基因组复制
    染色体数变化: *ploidy (乘以倍性)

15. Chromothripsis (CHT, 染色体碎裂)
    缩写: cht
    定义: 染色体灾难性碎裂后重组
    染色体数变化: 不确定

重建策略说明:
================================================================================
- 精确区分: EEJ vs NCF vs Fission vs RT vs URT vs UT vs Inversion
  这些改变染色体结构，必须准确检测才能正确重建祖先
  
- 统一处理: Gene Loss + Fractionation + Segmental Deletion
  都表现为基因缺失，重建时统一处理，不需要区分具体机制
  
- 统一处理: Tandem Dup + Dispersed Dup + Segmental Dup
  都表现为基因重复，重建时统一处理
  
- 外类群投票: 所有gain/loss类事件都需要外类群确定祖先状态
"""

from collections import defaultdict
import re


# ============================================================
# Canonical event type definitions
# ============================================================

# All canonical event types used across simulator, AKR, and eval
EVENT_TYPES = {
    # === Large-scale rearrangements (must be precisely distinguished) ===
    'inversion',                           # INV (generic)
    'internal_inversion',                  # II: both breakpoints internal
    'telomere_inversion',                  # TI: one breakpoint at telomere
    'centromere_spanning_inversion',       # CSI: spans centromere
    'reciprocal_translocation',            # RT
    'unbalanced_reciprocal_translocation', # URT
    'unidir_trans',                        # UT
    'ncf',                                 # NCF
    'eej',                                 # EEJ
    'fission',                             # FIS

    # === Gene-level events (can be unified in reconstruction) ===
    'gene_gain',                           # gain
    'gene_loss',                           # loss
    'fractionation',                       # frac
    'tandem_dup',                          # tdup
    'dispersed_dup',                       # ddup
    'proximal_dup',
    'seg_duplication',                     # segdup
    'seg_deletion',                        # segdel
    'indel',

    # === Special ===
    'chromothripsis',                      # cht
    'WGD',
}

# Chromosomal rearrangements that change chromosome structure
# These MUST be precisely distinguished during reconstruction
CHROM_REARRANGEMENTS = {
    # Inversions (all subtypes - must be detected before large-scale events)
    'inversion',
    'internal_inversion',
    'telomere_inversion',
    # Note: centromere_spanning_inversion (CSI) removed in v4.1
    # centromere position unknown, cannot reconstruct
    # Other large-scale rearrangements
    'reciprocal_translocation',
    'unbalanced_reciprocal_translocation',
    'unidir_trans',
    'ncf',
    'eej',
    'fission',
}

# Events that change chromosome count
CHROM_COUNT_EVENTS = {
    'ncf': -1,      # nested fusion: 2 chroms -> 1
    'eej': -1,      # end-end join: 2 chroms -> 1
    'fission': +1,  # fission: 1 chrom -> 2
    'WGD': None,    # depends on ploidy
}

# Events that can be unified during reconstruction
# (reconstruction doesn't need to distinguish these subtypes)
UNIFIED_EVENTS = {
    # All deletion/loss types -> unified as 'deletion'
    'gene_loss': 'deletion',
    'fractionation': 'deletion',
    'seg_deletion': 'deletion',
    'indel': 'deletion',
    # All duplication types -> unified as 'duplication'
    'tandem_dup': 'duplication',
    'dispersed_dup': 'duplication',
    'proximal_dup': 'duplication',
    'seg_duplication': 'duplication',
}

# Simulator short names -> canonical names
SIM_SHORT_NAMES = {
    'inv': 'inversion',
    'ii': 'internal_inversion',
    'ti': 'telomere_inversion',
    'csi': 'centromere_spanning_inversion',
    'rt': 'reciprocal_translocation',
    'urt': 'unbalanced_reciprocal_translocation',
    'ncf': 'ncf',
    'eej': 'eej',
    'fis': 'fission',
    'utrans': 'unidir_trans',
    'tdup': 'tandem_dup',
    'ddup': 'dispersed_dup',
    'segdup': 'seg_duplication',
    'segdel': 'seg_deletion',
    'frac': 'fractionation',
    'gain': 'gene_gain',
    'loss': 'gene_loss',
    'cht': 'chromothripsis',
}

# AKR detector event type aliases -> canonical names
AKR_ALIASES = {
    'internal_inversion': 'inversion',
    'telomere_inversion': 'inversion',
    'translocation': 'reciprocal_translocation',
    'unbalanced_reciprocal_translocation': 'reciprocal_translocation',
    'chrom_gain': 'gene_gain',
    'chrom_loss': 'gene_loss',
}

# Full name -> abbreviation mapping (for display)
EVENT_ABBREVIATIONS = {
    'inversion': 'INV',
    'internal_inversion': 'II',
    'telomere_inversion': 'TI',
    'centromere_spanning_inversion': 'CSI',
    'reciprocal_translocation': 'RT',
    'unbalanced_reciprocal_translocation': 'URT',
    'unidir_trans': 'UT',
    'ncf': 'NCF',
    'eej': 'EEJ',
    'fission': 'FIS',
    'gene_gain': 'gain',
    'gene_loss': 'loss',
    'fractionation': 'frac',
    'tandem_dup': 'TDUP',
    'dispersed_dup': 'DDUP',
    'proximal_dup': 'PDUP',
    'seg_duplication': 'SDUP',
    'seg_deletion': 'SDEL',
    'indel': 'indel',
    'chromothripsis': 'CHT',
    'WGD': 'WGD',
}

# Reverse mapping: abbreviation -> full name
ABBREV_TO_EVENT = {v.lower(): k for k, v in EVENT_ABBREVIATIONS.items()}


def canonicalize_event_type(event_type, source='auto'):
    """Convert any event type name to canonical form.
    
    Args:
        event_type: Input event type string
        source: 'sim', 'akr', or 'auto' (try both)
    
    Returns:
        Canonical event type name, or original if unknown
    """
    et = event_type.lower().strip()
    
    # Direct match
    if et in EVENT_TYPES:
        return et
    
    # Try simulator short names
    if source in ('auto', 'sim'):
        if et in SIM_SHORT_NAMES:
            return SIM_SHORT_NAMES[et]
    
    # Try AKR aliases
    if source in ('auto', 'akr'):
        if et in AKR_ALIASES:
            return AKR_ALIASES[et]
    
    # Try abbreviation
    if et in ABBREV_TO_EVENT:
        return ABBREV_TO_EVENT[et]
    
    # Try case-insensitive match
    for canonical in EVENT_TYPES:
        if et == canonical.lower():
            return canonical
    
    return event_type


def unify_event_type(event_type):
    """Unify event types for reconstruction.
    
    Reconstruction doesn't need to distinguish:
    - loss / fractionation / seg_deletion / indel -> 'deletion'
    - tandem_dup / dispersed_dup / proximal_dup / seg_duplication -> 'duplication'
    
    But MUST keep distinct:
    - inversion, RT, URT, UT, NCF, EEJ, Fission
    """
    canonical = canonicalize_event_type(event_type)
    return UNIFIED_EVENTS.get(canonical, canonical)


def is_chrom_rearrangement(event_type):
    """Check if event type is a chromosomal rearrangement.
    
    These are the events that MUST be precisely distinguished
    during ancestral reconstruction.
    """
    return canonicalize_event_type(event_type) in CHROM_REARRANGEMENTS


def chrom_count_delta(event_type):
    """Get chromosome count change for an event type.
    
    Returns:
        int: chromosome count change
        None: for WGD (depends on ploidy)
    """
    canonical = canonicalize_event_type(event_type)
    return CHROM_COUNT_EVENTS.get(canonical, 0)


def get_abbreviation(event_type):
    """Get abbreviation for display."""
    canonical = canonicalize_event_type(event_type)
    return EVENT_ABBREVIATIONS.get(canonical, canonical.upper()[:4])


# ============================================================
# Tree parsing utilities (wraps tree.py for unified interface)
# ============================================================

def parse_tree(tree_file):
    """Parse a tree file using the standard number_nodes function.
    
    This wraps tree.number_nodes to ensure consistent tree labeling
    across simulator, AKR, and eval.
    
    Args:
        tree_file: Path to Newick tree file
    
    Returns:
        tree: ete3 Tree object with labeled internal nodes
        parent_of: dict {child_node: parent_node}
        children_of: dict {parent_node: [child1, child2]}
    """
    from .tree import number_nodes
    
    tree = number_nodes(tree_file)
    
    # Build parent-child mappings from ete3 tree
    parent_of = {}
    children_of = defaultdict(list)
    
    for node in tree.traverse():
        if not node.is_root():
            parent_of[node.name] = node.up.name
        for child in node.children:
            children_of[node.name].append(child.name)
    
    return tree, parent_of, dict(children_of)


def get_ancestor_nodes(parent_of):
    """Get all ancestor (internal) nodes from parent_of mapping.
    
    Returns:
        set of ancestor node names
    """
    ancestors = set()
    for child, parent in parent_of.items():
        ancestors.add(parent)
    return ancestors


def get_leaf_nodes(parent_of):
    """Get all leaf (species) nodes from parent_of mapping.
    
    Returns:
        set of leaf node names
    """
    all_nodes = set(parent_of.keys()) | set(parent_of.values())
    ancestors = get_ancestor_nodes(parent_of)
    return all_nodes - ancestors


def get_branches(parent_of):
    """Get all branches as (parent, child) tuples.
    
    Returns:
        list of (parent, child) tuples
    """
    return [(parent, child) for child, parent in parent_of.items()]


def branch_id(parent, child):
    """Create a branch identifier string."""
    return '{}-{}'.format(parent, child)


def parse_branch_id(branch_str):
    """Parse a branch identifier into (parent, child).
    
    Returns:
        (parent, child) tuple, or (branch_str, None) if not a branch id
    """
    parts = branch_str.split('-', 1)
    if len(parts) == 2:
        return parts[0], parts[1]
    return branch_str, None


# ============================================================
# Event record helpers
# ============================================================

class TAKREvent:
    """Unified event record for TAKR framework.
    
    Used by both simulator (forward) and AKR (backward).
    
    Attributes:
        event_type: Canonical event type name
        branch: Branch identifier (e.g., "N2-Sp_1" or just "N2")
        genes_involved: List of gene/HOG identifiers
        desc: Human-readable description
        support: Support count/confidence
        details: Additional type-specific details dict
    """
    
    def __init__(self, event_type, branch, genes_involved=None,
                 desc='', support=None, details=None):
        """
        Args:
            event_type: Canonical event type name
            branch: Branch identifier (e.g., "N2-Sp_1" or just "N2")
            genes_involved: List of gene/HOG identifiers
            desc: Human-readable description
            support: Support count/confidence
            details: Additional type-specific details dict
        """
        self.event_type = canonicalize_event_type(event_type)
        self.branch = branch
        self.genes_involved = genes_involved or []
        self.desc = desc
        self.support = support
        self.details = details or {}
    
    def to_dict(self):
        d = {
            'event_type': self.event_type,
            'branch': self.branch,
            'genes_involved': self.genes_involved,
            'desc': self.desc,
            'support': self.support,
        }
        if self.details:
            d['details'] = self.details
        return d
    
    @classmethod
    def from_dict(cls, d):
        return cls(
            event_type=d.get('event_type', 'unknown'),
            branch=d.get('branch', ''),
            genes_involved=d.get('genes_involved', []),
            desc=d.get('desc', ''),
            support=d.get('support'),
            details=d.get('details', {})
        )
    
    def __repr__(self):
        return '<{} on {}: {}>'.format(self.event_type, self.branch, self.desc)


# ============================================================
# Validation helpers
# ============================================================

def validate_chrom_count(node, parent_of, events_by_branch, ploidy=1):
    """Validate chromosome count using event counts.
    
    Formula: child_chrom = parent_chrom * ploidy + fission - ncf - eej
    
    Args:
        node: Node name
        parent_of: parent-child mapping
        events_by_branch: {branch_id: [TAKREvent, ...]}
        ploidy: Ploidy level (1 for normal, 2+ for post-WGD)
    
    Returns:
        (expected_chrom, actual_chrom) or None if can't compute
    """
    # This is a placeholder - actual implementation needs chromosome counts
    pass


def get_events_for_node(node, events_by_branch, parent_of):
    """Get all events associated with a node (all child branches).
    
    Args:
        node: Node name
        events_by_branch: {branch_id: [TAKREvent, ...]}
        parent_of: {child: parent} mapping
    
    Returns:
        list of TAKREvent
    """
    events = []
    for child, parent in parent_of.items():
        if parent == node:
            branch = branch_id(node, child)
            events.extend(events_by_branch.get(branch, []))
    return events


def count_events_by_type(events):
    """Count events by canonical type.
    
    Args:
        events: list of TAKREvent or event type strings
    
    Returns:
        Counter of canonical event types
    """
    from collections import Counter
    if events and isinstance(events[0], TAKREvent):
        return Counter(e.event_type for e in events)
    return Counter(canonicalize_event_type(e) for e in events)


# Backward compatibility: keep old names
RearrangementEvent = TAKREvent
