import re
from ete3 import Tree

class SpeciesTree:
	pass

def number_nodes(sptreefile):
	treestr = convertNHX(sptreefile)
	tree = Tree(treestr)
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
	tree.write(outfile=sptreefile + ".labeled.nwk", format=1)
	return tree 


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


