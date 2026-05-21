#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
从HOG信息创建单拷贝OG文件
根据HOG（Hierarchical Orthologous Groups）信息，对OG文件进行处理，
确保每个物种在每个OG中只保留一个基因。
"""

import time
from collections import defaultdict
from .hog import HOG
from .OrthoFinder import OrthoMCLGroup, gene_format_o
from .RunCmdsMP import logger
from ete3 import Tree


def load_hog_info_from_class(ogfile, orthfiles, sptreefile, paralog=False):
	"""
	从HOG类加载HOG信息
	:param ogfile: OG文件路径
	:param orthfiles: 正交文件路径列表
	:param sptreefile: 物种树文件路径
	:param paralog: 是否包含旁系同源，默认False
	:return: all_hogs字典, tree对象
	"""
	# 创建HOG实例并处理
	hog_instance = HOG(ogfile=ogfile, orthfiles=orthfiles, sptreefile=sptreefile, paralog=paralog)
	all_hogs = hog_instance.pipe(write_tsv=False)	# 不写入文件，只获取HOG
	
	# 建立基因到HOG的映射
			
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

def remove_hog_and_descendants(hog, node_to_hogs, all_hogs):
	"""递归删除 hog 及其所有后代 HOG 的基因，并从 node_to_hogs 中移除"""
	# 使用栈迭代，避免递归深度过大
	stack = [hog]
	while stack:
		cur = stack.pop()
		# 从它所在节点的列表中移除该 HOG
		node_list = node_to_hogs.get(cur.node_id)
		if node_list:
			try:
				node_list.remove(cur)
			except ValueError:
				pass  # 可能已经被移除了
		# 将它的孩子加入栈，以便处理后代
		for child_id in cur.children:
			child_hog = all_hogs.get(child_id)
			if child_hog:
				stack.append(child_hog)

def build_deletion_map(all_hogs, tree, gene_degree=None):
	"""
	遍历树一次，构建要删除的基因集合
	:param all_hogs: HOG字典
	:param tree: 物种树对象
	"""
	
	# 预构建节点到HOG的映射，避免在循环中遍历所有HOG
	node_to_hogs = defaultdict(list)
	for hog_id, hog_record in all_hogs.items():
		node_to_hogs[hog_record['node_id']].append(hog_record)
	
	# 确保每次遍历时节点顺序一致 - 按层级顺序遍历（从根到叶）
	traversal_nodes = tree.traverse(strategy="levelorder")
	
	# 初始化要删除的基因集合
	deleted_genes = set()
	
	# 遍历树的节点（从根到叶）
	for node in traversal_nodes:	# levelorder = 广度优先，确保从根到叶的处理顺序
		node_name = node.name
		node_hogs = node_to_hogs.get(node_name, [])

		if node_name == 'N0':
			# 按 og_id 分组（同 SOG 的 HOG 具有相同的 og_id）
			hogs_by_og = defaultdict(list)
			for hog in node_hogs:
				hogs_by_og[hog['og_id']].append(hog)	
			# 对每组sog，只保留物种数最多的 HOG（已按物种数降序），其余删除
			for og_id, hogs in hogs_by_og.items():
				if len(hogs) > 1:
					hogs.sort(key=lambda x: len({gene_format_o(g)[0] for g in x['genes']}), reverse=True)
					# hogs[0] 已为物种数最多，删除其余
					for hog in hogs[1:]:
						deleted_genes.update(hog['genes'])
						node_to_hogs[node_name].remove(hog)

		if node.is_leaf():
			# 对于叶子节点，获取其父节点的HOGs
			if node.up:
				node_hogs = node_to_hogs.get(node.up.name, [])

		logger.debug('Processing node %s (%d HOGs)', node_name, len(node_hogs))
		if not node_hogs:
			logger.debug('No HOGs for node %s (may have been removed by previous step)', node_name)

		# 一次性对当前节点的HOG进行排序，避免多次排序
		sorted_node_hogs = sorted(node_hogs, key=lambda x: x['hog_id'])

		for hog_record in sorted_node_hogs:
			# 检查是否为需要调试的OG			
			if node.is_leaf():
				# 叶子节点：获取该节点的第一个父节点的HOG，对这些HOG中该叶节点物种的基因进行处理
				leaf_hog_genes = set(hog_record['genes'])
				# 获取这个父HOG中的当前叶节点物种的基因
				leaf_species_genes = []
				for gene in leaf_hog_genes:
					sp, _ = gene_format_o(gene)
					if sp == node_name and gene not in deleted_genes:
						leaf_species_genes.append(gene)
						
			# 对这个HOG中的该物种的基因，如果多于一个，则只保留一个
			if len(leaf_species_genes) > 1:
				if gene_degree:
					# 保留 degree 最高的基因
					leaf_species_genes.sort(key=lambda g: gene_degree.get(g, 0), reverse=True)
				else:
					leaf_species_genes.sort()  # fallback: 字母序
				# 保留第一个基因，将其余的标记为删除
				for gene in leaf_species_genes[1:]:
					deleted_genes.add(gene)
			else:
				# 中间节点：对于当前HOG的每一个子HOG，按您描述的新逻辑处理
				# 直接使用HOG记录中的children属性获取当前HOG的子HOG
				children_by_node = defaultdict(list)
				for child_hog_id in sorted(hog_record['children']):  # 预先排序子HOG
					child_hog = all_hogs.get(child_hog_id)
					if child_hog:
						children_by_node[child_hog.node_id].append((child_hog, set(child_hog.genes)))
				#if len(children_by_node.keys()) < 2:
					#print(f"Warning: now Node {node_name}; Hog: {hog_record['hog_id']} has 0 child hog in an tree child node....")
				# 过滤掉已经被删除的HOG的后代
				for node_id, node_children_data in children_by_node.items():
					if len(node_children_data) > 1:
						# 初始化：先把第一个当作“当前最大”
						selected_item = node_children_data[0]
						# 遍历其余元素
						for item in node_children_data[1:]:
							child_hog, child_genes = item
							species_count = len({gene_format_o(g)[0] for g in child_genes})
							selected_species_count = len({gene_format_o(g)[0] for g in selected_item[1]})
							if species_count > selected_species_count:
								# 新的更大 → 把旧的（原最大）整个删除
								deleted_genes.update(selected_item[1])
								remove_hog_and_descendants(selected_item[0], node_to_hogs, all_hogs)
								selected_item = item
							else:
								# 不是更大 → 直接删除当前这个
								deleted_genes.update(child_genes)
								remove_hog_and_descendants(child_hog, node_to_hogs, all_hogs)	
	return deleted_genes


def process_og_with_hog(og_file, hog_args, output_file, restore_gene=False, restore_log=None, debug_output=None):
	"""
	根据HOG信息处理OG文件
	:param og_file: 输入OG文件路径 (MCL格式)
	:param hog_args: 包含HOG参数的字典 (ogfile, orthfiles, sptreefile, paralog)
	:param output_file: 输出文件路径
	:param restore_gene: 当一个物种没有任何基因保留时，是否恢复该物种的第一个基因
	:param restore_log: 记录恢复基因的日志文件路径
	:param debug_output: 调试输出文件路径，用于输出all_deleted_genes集合内容
	"""
	start_time = time.time()  # 添加开始时间记录
	# 从HOG类加载信息
	all_hogs, tree = load_hog_info_from_class(
		ogfile=hog_args['ogfile'], 
		orthfiles=hog_args['orthfiles'], 
		sptreefile=hog_args['sptreefile'],
		paralog=hog_args.get('paralog', False)
	)
	
	# 从 orthfiles 构建基因 degree 字典（复用 Collinearity 解析器兼容各种格式）
	from .mcscan import Collinearity
	gene_degree = defaultdict(int)
	for f in hog_args['orthfiles']:
		for rc in Collinearity(f):
			for g1, g2 in zip(rc.genes1, rc.genes2):
				gene_degree[g1] += 1
				gene_degree[g2] += 1
		
	# 构建删除基因的映射，只遍历一次树
	all_deleted_genes = build_deletion_map(all_hogs, tree, gene_degree)
	end_time = time.time()
	logger.info('build_deletion_map took %.2f s', end_time - start_time)
	
	# 如果需要调试输出，将all_deleted_genes写入文件
	if debug_output:
		with open(debug_output, 'w') as f_debug:
			f_debug.write(f"Total deleted genes count: {len(all_deleted_genes)}\n")
			for gene in sorted(all_deleted_genes):	# 排序以确保一致性
				f_debug.write(f"{gene}\n")
	
	# 打开日志文件（如果指定了路径）
	log_file = None
	if restore_log:
		log_file = open(restore_log, 'w')

	# 处理OG文件，使用 OrthoMCLGroup 类解析
	if isinstance(output_file, str):
		fout = open(output_file, 'w')
	else:
		fout = output_file
	try:
		for og in OrthoMCLGroup(og_file):
			og_id = og.ogid
			genes = og.genes
			
			# 按物种分组基因
			genes_by_species = defaultdict(list)
			for gene in genes:
				sp, _ = gene_format_o(gene)
				if sp is not None:
					genes_by_species[sp].append(gene)
			
			# 按照预计算的删除集合来处理基因
			kept_genes = set()
			
			# 遍历当前OG中的基因，根据预计算的删除集合来决定哪些基因保留
			restored_species_for_og = []	# 记录当前OG中恢复的物种和基因
			for species, species_genes in genes_by_species.items():
				# 从当前物种的基因中找出未被删除的基因
				active_genes = [g for g in species_genes if g not in all_deleted_genes]
				
				# 按字母顺序排序，确保一致性
				active_genes.sort()
				
				if active_genes:
					# 如果有保留的基因
					kept_genes.update(active_genes)
				elif species_genes:	# 如果没有保留的基因，记录应该恢复的基因（不管restore_gene是否为True）
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
			
		# 最终保留的基因（restore时保留所有，否则排除被删基因）
		if restore_gene:
			final_kept = sorted(kept_genes)
		else:
			final_kept = sorted(kept_genes - all_deleted_genes)
			
			# 输出处理后的OG行
		if final_kept:
			output_line = "{}: {}".format(og_id, ' '.join(sorted(final_kept)))
			fout.write(output_line + '\n')
	finally:
		if isinstance(output_file, str):
			fout.close()

	# 关闭日志文件（如果打开）
	if log_file:
		log_file.close()