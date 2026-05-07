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
from .hog import HOG


def load_hog_info_from_class(ogfile, orthfiles, sptreefile, paralog=False):
	"""
	从HOG类加载HOG信息
	:param ogfile: OG文件路径
	:param orthfiles: 正交文件路径列表
	:param sptreefile: 物种树文件路径
	:param paralog: 是否包含旁系同源，默认False
	:return: gene_to_hog_dict 字典，映射基因到HOG
	"""
	# 创建HOG实例并处理
	hog_instance = HOG(ogfile=ogfile, orthfiles=orthfiles, sptreefile=sptreefile, paralog=paralog)
	all_hogs = hog_instance.pipe(write_tsv=False)  # 不写入文件，只获取HOG
	
	gene_to_hog_dict = {}
	
	# 遍历所有HOG，建立基因到HOG的映射
	for hog_id, hog_record in all_hogs.items():
		for gene in hog_record['genes']:
			gene_to_hog_dict[gene] = hog_id
	
	return gene_to_hog_dict


def process_og_with_hog(og_file, hog_args, output_file):
	"""
	根据HOG信息处理OG文件
	:param og_file: 输入OG文件路径 (MCL格式)
	:param hog_args: 包含HOG参数的字典 (ogfile, orthfiles, sptreefile, paralog)
	:param output_file: 输出文件路径
	"""
	# 从HOG类加载信息
	gene_to_hog = load_hog_info_from_class(
		ogfile=hog_args['ogfile'], 
		orthfiles=hog_args['orthfiles'], 
		sptreefile=hog_args['sptreefile'],
		paralog=hog_args.get('paralog', False)
	)
	
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
			
			# 按物种分组基因
			genes_by_species = defaultdict(list)
			for gene in genes:
				if '|' in gene:
					# 假设基因名格式为 "gene_id|species" 或 "species|gene_id"
					if gene.count('|') == 1:
						parts = gene.split('|')
						if len(parts) == 2:
							# 判断哪一部分是物种（通常第一部分或第二部分是物种）
							# 根据常见的命名约定，通常是 species|gene 或 gene|species
							# 这里我们假设格式为 species|gene
							species = parts[0]
							genes_by_species[species].append(gene)
			
			# 为每个物种选择一个基因
			selected_genes = []
			for species, species_genes in genes_by_species.items():
				if not species_genes:
					continue
				
				# 检查此物种是否有基因存在于HOG中
				species_genes_in_hog = [g for g in species_genes if g in gene_to_hog]
				species_genes_not_in_hog = [g for g in species_genes if g not in gene_to_hog]
				
				# 将不在HOG中的基因直接加入结果
				selected_genes.extend(species_genes_not_in_hog)
				
				if not species_genes_in_hog:
					# 如果此物种没有基因在HOG中，跳过后续HOG处理
					continue
				
				# 按HOG分组，然后从每组中选择一个基因
				hog_groups = defaultdict(list)
				for gene in species_genes_in_hog:
					hog_id = gene_to_hog.get(gene, "unknown_HOG")
					hog_groups[hog_id].append(gene)
				
				# 按物种收集所有基因（不管属于哪个HOG），每个物种只保留一个基因
				species_to_gene = {}  # 映射：物种 -> 该物种的一个基因
				
				for hog_id, genes_in_hog in hog_groups.items():
					if genes_in_hog:
						# 对HOG中的基因进行排序，选择第一个
						selected_gene = sorted(genes_in_hog)[0]
						
						# 提取基因所属的物种
						if '|' in selected_gene:
							species = selected_gene.split('|')[0]  # 假设基因名格式为 species|gene
						else:
							species = selected_gene  # 如果没有物种信息，则使用基因名本身
							
						# 如果这个物种还没有被记录，则添加这个基因
						if species not in species_to_gene:
							species_to_gene[species] = selected_gene
				
				# 将每个物种选择的基因添加到结果中
				selected_genes.extend(species_to_gene.values())
			
			# 输出处理后的OG行
			if selected_genes:
				output_line = "{}: {}".format(og_id, ' '.join(sorted(selected_genes)))
				fout.write(output_line + '\n')


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