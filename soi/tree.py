import re
from ete3 import Tree

class SpeciesTree:
	pass

def number_nodes(sptreefile):
	treestr = convertNHX(sptreefile)
#	print(treestr)
	tree = load_tree_smart(treestr)
#	try: tree = Tree(treestr, format=0)
#	except :
#		tree = Tree(treestr, format=1)
	i = 0
	for node in tree.traverse():
		node.show = True
		if node.is_leaf():
			continue
		if node.name:
			continue
		name = 'N{}'.format(i)
		node.name = name
		i += 1
	# 压缩单子节点：删掉只有一个孩子的内部节点，祖孙直连
	for node in list(tree.traverse()):
		if not node.is_leaf() and not node.is_root() and len(node.children) == 1:
			node.delete(prevent_nondicotomic=False)
	tree.write(outfile=sptreefile + ".labeled.nwk", format=1)
	return tree 


def load_tree_smart(treestr):
    """
    智能判断：检测内部节点标签是纯数字(support)还是文本(name)
    """
    # 提取所有 ")LABEL:" 或 ")LABEL;" 模式的内部节点标签
    # 匹配 ) 后到 : 或 ; 或 ) 前的内容
    internal_labels = re.findall(r'\)([^:;,\(\)]+)(?=:)', treestr)
    
    # 判断是否有非数字标签（允许空，允许小数）
    has_text_label = any(
        label and not re.match(r'^\d*\.?\d+$', label) 
        for label in internal_labels
    )
    
    if has_text_label:
        # 有文本标签（如 N1），用 format=1
        return Tree(treestr, format=1)
    else:
        # 只有数字或空，用 format=0（保留 support 值）
        return Tree(treestr, format=0)

def convert_newick(line: str) -> str:
    '''(A,B[p=2]); (A,B[p=2]:0.1); (A,B:0.1[p=2]); (A,B[&&NHX:p=2])'''
    # 1. [tag]:length → :length[&&NHX:tag]（带安全检查）
    def patch_tag_length(m):
        tag, length = m.groups()
        if not tag.startswith("&&NHX:"):
            tag = f"&&NHX:{tag}"
        return f":{length}[{tag}]"
    
    line = re.sub(r'\[([^]]+)\]:([\d.eE-]+)', patch_tag_length, line)

    # 2. :length[tag] → :length[&&NHX:tag]
    def patch_length_tag(m):
        length, tag = m.groups()
        if not tag.startswith("&&NHX:"):
            tag = f"&&NHX:{tag}"
        return f":{length}[{tag}]"
    
    line = re.sub(r':([\d.eE-]+)\[([^]]+)\]', patch_length_tag, line)

    # 3. [tag]（孤立无长度）→ [&&NHX:tag]
    def patch_isolated(m):
        tag = m.group(1)
        if not tag.startswith("&&NHX:"):
            tag = f"&&NHX:{tag}"
        return f"[{tag}]"
    
    line = re.sub(r'\[([^]]+)\](?![\d:\[])', patch_isolated, line)
#    print(line)
    return line
	
def convertNHX(inNwk, ):
    nwk = []
    for line in open(inNwk):
        nwk += [convert_newick(line)]
    return ''.join(nwk)