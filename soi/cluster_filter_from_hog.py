#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
从HOG信息创建单拷贝OG文件
根据HOG（Hierarchical Orthologous Groups）信息，对OG文件进行处理，
确保每个物种在每个OG中只保留一个基因。
"""

import sys
import argparse
from collections import defaultdict
from .hog import HOG, number_nodes
from ete3 import Tree


def load_hog_info_from_class(ogfile, orthfiles, sptreefile, paralog=False):
	"""
	从HOG类加载HOG信息
	:param ogfile: OG文件路径
	:param orthfiles: 正交文件路径列表
	:param sptreefile: 物种树文件路径
	:param paralog: 是否包含旁系同源，默认False
	:return: all_hogs字典；tree对象
	"""
	# 创建HOG实例并处理
	hog_instance = HOG(ogfile=ogfile, orthfiles=orthfiles, sptreefile=sptreefile, paralog=paralog)
	all_hogs = hog_instance.pipe(write_tsv=False)  # 不写入文件，只获取HOG
	
	# 建立基因到HOG的映射（虽然目前未使用，但保持接口兼容性）
	# 移除未使用的gene_to_hog_dict变量，直接返回all_hogs和tree
	return all_hogs, hog_instance.tree


def get_all_descendants(hog_id, all_hogs):
	"""
	递归获取某个HOG的所有后代HOG ID
	"""
	descendants = set()
	hog_record = all_hogs.get(hog_id)
	if hog_record:
		for child_id in hog_record['children']:
			descendants.add(child_id)
			descendants.update(get_all_descendants(child_id, all_hogs))
	return descendants


def build_deletion_map(all_hogs, tree):
	"""
	遍历树一次，构建要删除的基因集合
	"""
	# 预构建节点到HOG的映射，避免在循环中遍历所有HOG
	node_to_hogs = defaultdict(list)
	# 创建HOG记录对象到其ID的映射，用于排序
	record_to_hog_id = {}
	for hog_id, hog_record in all_hogs.items():
		node_to_hogs[hog_record['node_id']].append(hog_record)
		record_to_hog_id[hog_record] = hog_id
	
	# 确保每次遍历时节点顺序一致 - 按层级顺序遍历（从根到叶）
	traversal_nodes = tree.traverse(strategy="levelorder")
	
	# 初始化要删除的基因集合
	leaf_deleted_genes = set()  # 记录在叶节点处理中被删除的基因
	internal_deleted_genes = set()  # 记录在内部节点处理中被删除的基因
	treeorder = open('tree_order.txt', 'w')
	# 遍历树的节点（从根到叶）
	for node in traversal_nodes:  # levelorder = 广度优先，确保从根到叶的处理顺序
		node_name = node.name
		treeorder.write(f"{node_name}\n")
		# 获取当前节点对应的HOGs (使用预构建的映射)
		node_hogs = node_to_hogs.get(node_name, [])
		
		# 确保每次处理HOG的顺序一致 - 直接使用HOG对象的hog_id属性进行排序
		sorted_node_hogs = sorted(node_hogs, key=lambda x: record_to_hog_id.get(x, ""))
		
		for hog_record in sorted_node_hogs:
			# hog_genes = set(hog_record['genes'])  # 这个变量未被使用，保留注释说明
			
			if node.is_leaf():
				# 叶子节点：获取该节点的第一个父节点的HOG，对这些HOG中该叶节点物种的基因进行处理
				# 即在每个HOG中，对该物种的多拷贝基因只保留一个
				parent_node = node.up
				if parent_node is not None:
					parent_node_name = parent_node.name
					parent_node_hogs = node_to_hogs.get(parent_node_name, [])
					
					# 遍历父节点的每个HOG，对属于当前叶节点物种的基因进行处理
					for parent_hog_record in parent_node_hogs:
						parent_hog_genes = set(parent_hog_record['genes'])
						
						# 获取这个父HOG中的当前叶节点物种的基因
						leaf_species_genes = set()
						for gene in parent_hog_genes:
							if gene.split('|')[0] == node_name and gene not in leaf_deleted_genes:
								leaf_species_genes.add(gene)
						
						# 对这个HOG中的该物种的基因，如果多于一个，则只保留一个
						leaf_species_list = sorted(list(leaf_species_genes))  # 排序确保一致性
						if len(leaf_species_list) > 1:
							# 保留第一个基因，将其余的标记为删除
							for gene in leaf_species_list[1:]:
								leaf_deleted_genes.add(gene)
			else:
				# 中间节点：对于当前HOG的每一个子HOG，按您描述的新逻辑处理
				# 直接使用HOG记录中的children属性获取当前HOG的子HOG
				children_hogs = []
				for child_hog_id in sorted(hog_record['children']):
					child_hog = all_hogs.get(child_hog_id)
					if child_hog:
						children_hogs.append(child_hog)
				
				# 过滤掉已经被删除的HOG的后代
				# 检查方式：如果子HOG中的基因在internal_deleted_genes中，则该子HOG应被跳过
				valid_children_hogs = []
				for child_hog in children_hogs:
					# 获取这个child_hog中所有的基因
					child_all_genes = set(child_hog['genes'])
					
					# 检查这个child_hog是否包含任何已被删除的基因
					# 如果这个child_hog中有任何基因在internal_deleted_genes中，
					# 说明这个child_hog是之前已被删除的HOG
					has_deleted_genes = any(g in internal_deleted_genes for g in child_all_genes)
					
					if has_deleted_genes:
						continue  # 跳过已被删除的HOG
					
					# 传递整个child_hog给后续处理，不在此处过滤基因
					valid_children_hogs.append((child_hog, child_all_genes))
				
				if valid_children_hogs:
					# 按节点对子HOG进行分组
					children_by_node = defaultdict(list)
					for child_hog, child_target_genes in valid_children_hogs:
						children_by_node[child_hog['node_id']].append((child_hog, child_target_genes))
					
					# 检查是否有两个以上的子节点有HOG
					nodes_with_hogs = list(children_by_node.keys())  # 不再排序，保持原始顺序
					if len(nodes_with_hogs) >= 2:
						# 如果有多个子节点都有HOG，每个节点保留物种数最多的那个子HOG
						for node_id in nodes_with_hogs:
							node_children_data = children_by_node[node_id]
							if len(node_children_data) > 1:
								# 找出这个节点中物种数最多的子HOG
								max_species_count = 0
								selected_child_hog = None
								
								for child_hog, child_target_genes in node_children_data:
									# 计算这个子HOG中的物种数
									child_species = {gene.split('|')[0] for gene in child_target_genes}
									if len(child_species) > max_species_count:
										max_species_count = len(child_species)
										selected_child_hog = child_hog
								# 删除其他子HOG中的基因（同一节点的其他HOG）
								for child_hog, child_target_genes in node_children_data:
									if child_hog != selected_child_hog:
										for gene in child_target_genes:
											internal_deleted_genes.add(gene)
	treeorder.close()
	return leaf_deleted_genes, internal_deleted_genes


def process_og_with_hog(og_file, hog_args, output_file, restore_gene=False, restore_log=False):
	"""
	根据HOG信息处理OG文件
	:param og_file: 输入OG文件路径 (MCL格式)
	:param hog_args: 包含HOG参数的字典 (ogfile, orthfiles, sptreefile, paralog)
	:param output_file: 输出文件路径
	:param restore_gene: 当一个物种没有任何基因保留时，是否恢复该物种的第一个基因
	:param restore_log: 记录恢复基因的日志文件路径
	"""
	# 从HOG类加载信息
	all_hogs, tree = load_hog_info_from_class(
		ogfile=hog_args['ogfile'], 
		orthfiles=hog_args['orthfiles'], 
		sptreefile=hog_args['sptreefile'],
		paralog=hog_args.get('paralog', False)
	)
	
	# 构建删除基因的映射，只遍历一次树
	leaf_deleted_genes, internal_deleted_genes = build_deletion_map(all_hogs, tree)
	all_deleted_genes = leaf_deleted_genes.union(internal_deleted_genes)	
	
	# 打开日志文件（如果指定了路径）
	log_file = None
	if restore_log:
		log_file = open(restore_log, 'w')

	# 处理OG文件
	with open(og_file, 'r') as fin, open(output_file, 'w') as fout:
		for line in fin:
			line = line.strip()
			if not line:
				continue
				
			# 解析OG行 (MCL格式: "OG_id: gene1 gene2 gene3...")
			if ':' in line:
				og_id, genes_part = line.split(':', 1)
				og_id = og_id.strip()
				genes = genes_part.strip().split()
			else:
				# 如果没有OG ID，生成一个默认ID
				og_id = "OG_" + str(hash(line) % 10000)
				genes = line.split()
			
			# 按物种分组基因并预先计算基因集合
			genes_by_species = defaultdict(list)
			gene_set = set(genes)  # 预先创建基因集合，加快查找速度
			for gene in genes:
				if '|' in gene:
					# 提取第一个|前面的部分作为物种名称
					species = gene.split('|')[0]
					genes_by_species[species].append(gene)
			
			# 按照预计算的删除集合来处理基因
			kept_genes = set()
			
			# 遍历当前OG中的基因，根据预计算的删除集合来决定哪些基因保留
			restored_species_for_og = []  # 记录当前OG中恢复的物种和基因
			for species, species_genes in genes_by_species.items():
				# 从当前物种的基因中找出未被删除的基因
				active_genes = [g for g in species_genes if g not in all_deleted_genes]
				
				# 按字母顺序排序，确保一致性
				active_genes.sort()
				
				if active_genes:
					# 如果有保留的基因，保留第一个（按字母顺序）
					kept_genes.update(active_genes)
				elif species_genes:  # 如果没有保留的基因，记录应该恢复的基因（不管restore_gene是否为True）
					# 如果没有保留的基因，但启用了恢复功能，则恢复该物种的第一个基因
					original_gene = species_genes[0]
					
					# 即使不实际恢复基因，也要记录本应恢复的基因信息
					restored_species_for_og.append(f"{species}:{original_gene}")
					
					# 只有在restore_gene为True时才真正添加到kept_genes
					if restore_gene:
						kept_genes.add(original_gene)
			
			# 如果有恢复的物种，将它们写入日志文件（一个OG一行）
			if len(restored_species_for_og) > 0 and log_file:
				log_file.write(f"{og_id}\t{';'.join(restored_species_for_og)}\n")
			
			# 最终保留的基因
			final_kept = []
			for gene in kept_genes:
				if gene in gene_set:  # 确保基因属于当前OG
					final_kept.append(gene)
			
			# 输出处理后的OG行
			if final_kept:
				output_line = "{}: {}".format(og_id, ' '.join(sorted(final_kept)))
				fout.write(output_line + '\n')
	
	# 关闭日志文件（如果打开）
	if log_file:
		log_file.close()


def main():
	parser = argparse.ArgumentParser(description='根据HOG信息创建单拷贝OG文件')
	parser.add_argument('-og', '--og-file', required=True, type=str,
						help='输入OG文件路径 (MCL格式)')
	parser.add_argument('-orthfiles', '--ortholog-files', required=True, type=str, nargs='+',
						help='正交文件路径列表')
	parser.add_argument('-sptree', '--species-tree-file', required=True, type=str,
						help='物种树文件路径')
	parser.add_argument('-paralog', '--paralog', action='store_true',
						help='包含旁系同源 (默认: False)')
	parser.add_argument('-o', '--output', required=True, type=str,
						help='输出文件路径')
	
	args = parser.parse_args()
	
	hog_args = {
		'ogfile': args.og_file,
		'orthfiles': args.ortholog_files,
		'sptreefile': args.species_tree_file,
		'paralog': args.paralog
	}
	
	print("处理OG文件: {} 使用HOG信息".format(args.og_file))
	process_og_with_hog(args.og_file, hog_args, args.output)
	print("结果已保存至: {}".format(args.output))


if __name__ == "__main__":
	main()