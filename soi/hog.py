import sys
import re
import networkx as nx
from collections import defaultdict, Counter

from .OrthoFinder import OrthoMCLGroup
from .mcscan import ColinearGroups
from .tree import number_nodes
from .RunCmdsMP import logger

class HOG:
	def __init__(self, ogfile=None, orthfiles=None, sptreefile=None, outfile = "HOGs.tsv",
		  paralog=False):
		self.ogfile = ogfile
		self.orthfiles = orthfiles
		self.sptreefile = sptreefile
		self.outfile = outfile
		self.noparalog = not paralog
	def pipe(self, write_tsv=True):
		logger.info(f'Reading and Numbering species tree from {self.sptreefile}')
		self.tree = sptree = number_nodes(self.sptreefile)
		self.species = sptree.get_leaf_names()
		
		logger.info(f'Loading orthologs from {self.orthfiles}')
		graph = ColinearGroups(self.orthfiles, spsd=self.species, 
					noparalog=self.noparalog).graph
		
		logger.info(f'Splitting HOGs from {self.ogfile}')
		prefix = 'hog'
		# 存储所有HOG：key = hog_id, value = HOG对象/dict
		self.all_hogs = {}
		self.all_genes = []
		for og in OrthoMCLGroup(self.ogfile):
			og_genes = og.genes
			self.all_genes  += og_genes
			og_species = set(og.species)   # 这个OG包含的物种
			og_id = og.ogid
			og_spdict = og.spdict		 # species -> [genes]
			subgraph = graph.subgraph(og_genes).copy()

			# 临时存储：node_id → 该节点下的所有HOGs
			node_to_hogs = defaultdict(list)

			# ------------------------------
			# 遍历物种树（HOG标准：后序）
			# ------------------------------
			for node in sptree.traverse(strategy="postorder"):
				if self.noparalog and node.is_leaf():
					continue
				node_id = node.name
				node_species = set(node.get_leaf_names())

				# --------------------------
				# 核心：无交集 → 跳过
				# --------------------------
				if not node_species & og_species:
					continue

				# --------------------------
				# 从 og_spdict 直接取 subset_genes
				# --------------------------
				subset_genes = []
				for sp in node_species & og_species:
					subset_genes.extend(og_spdict[sp])

				if not subset_genes:
					continue

				# --------------------------
				# 共线性连通分量 = 单个HOG
				# --------------------------
				hog_subgraph = subgraph.subgraph(subset_genes)
				connected_components = list(nx.connected_components(hog_subgraph))

				# --------------------------
				# 为每个连通分量创建 HOG
				# --------------------------
				for idx, cc_genes in enumerate(connected_components):
					hog_id = f"{og_id}.{node_id}.{prefix}{idx}"

					# 构建 HOG 条目（含层级、父子关系）
					# hog = {
						# "hog_id": hog_id,
						# "og_id": og_id,
						# "node_id": node_id,
						# "genes": list(cc_genes),
						# "species": [sp for sp in og_spdict if sp in node_species],
						# "parent": None,		  # 关键：父HOG
						# "children": [],		  # 关键：子HOG
					# }
					hog = HOGrecord(
						hog_id=hog_id,
						og_id=og_id,
						node_id=node_id,
						genes=list(cc_genes),
						species=[sp for sp in og_spdict if sp in node_species],
						parent=None,
						children=[]
					)
					self.all_hogs[hog_id] = hog
					node_to_hogs[node_id].append(hog)

			# ------------------------------
			# 最关键步骤：建立 HOG 层级父子关系
			# ------------------------------
			for node in sptree.traverse(strategy="postorder"):
				node_id = node.name
				current_hogs = node_to_hogs.get(node_id, [])

				if not current_hogs:
					continue

				# 找父节点
				parent_node = node.up
				if not parent_node:
					continue  # 根节点无parent

				p_node_id = parent_node.name
				parent_hogs = node_to_hogs.get(p_node_id, [])

				# 每个子HOG必须找到包含它的父HOG
				for hog in current_hogs:
					hog_genes = set(hog["genes"])
					for phog in parent_hogs:
						phog_genes = set(phog["genes"])
						if hog_genes.issubset(phog_genes):
							hog["parent"] = phog["hog_id"]
							phog["children"].append(hog["hog_id"])
							break

		logger.info(f"Processed {len(self.all_hogs)} HOGs")
		print(Counter(v['node_id'] for v in self.all_hogs.values()))
		logger.info("All HOGs with hierarchy built successfully!")
		
		if write_tsv:
			with open(self.outfile, 'w', encoding='utf-8') as f:
				self.write_all_hogs_in_one_file(f)
			logger.info(f"HOGs written to {self.outfile}")
		return self.all_hogs
		
	def write_all_hogs_in_one_file(self, fout=sys.stdout):
		"""
		最终输出格式：HOG  OG  Tree_Node  Parent  Genes
		Genes: 空格分隔
		Parent: 父HOG，无父节点则为 Root
		"""
		# 表头
		header = ["HOG", "OG", "Tree_Node", "Parent", "Genes"]
		fout.write("\t".join(header) + "\n")

		for hog in self.all_hogs.values():
			hog_id = hog["hog_id"]
			og_id = hog["og_id"]
			node_id = hog["node_id"]
			parent_id = hog["parent"] if hog["parent"] is not None else "Root"
			gene_str = " ".join(sorted(hog["genes"]))  # 空格分隔

			# 写入一行
			row = [hog_id, og_id, node_id, parent_id, gene_str]
			fout.write("\t".join(row) + "\n")

class HOGrecord:
    __slots__ = ["hog_id", "og_id", "node_id", "genes", "species", "parent", "children"]
    
    def __init__(self, hog_id, og_id, node_id, genes, species, parent=None, children=None):
        # 核心字段
        self.hog_id = hog_id
        self.og_id = og_id
        self.node_id = node_id
        self.genes = genes
        self.species = species
        # 层级关系
        self.parent = parent
        self.children = children if children is not None else []

    # 🔥 支持字典访问：hog["hog_id"]
    def __getitem__(self, key):
        return getattr(self, key)

    # 🔥 支持字典赋值：hog["parent"] = xxx
    def __setitem__(self, key, value):
        setattr(self, key, value)

    # 🔥 可哈希：唯一标识为 hog_id
    def __hash__(self):
        return hash(self.hog_id)

    # 🔥 相等判断：两个HOG的hog_id相同则视为同一个对象
    def __eq__(self, other):
        return isinstance(other, HOG) and self.hog_id == other.hog_id

    # 打印美化（调试用）
    def __repr__(self):
        return f"{self.hog_id}"

def main():
	HOG(ogfile=sys.argv[1], orthfiles=sys.argv[2], sptreefile=sys.argv[3]).pipe()
if __name__ == '__main__':
	main()
