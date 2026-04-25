import sys,os
import argparse
from collections import OrderedDict
import numpy as np
import matplotlib as mpl
mpl.use("Agg")
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42

def Args():
	parser = argparse.ArgumentParser()
	parser.add_argument('--kaks', metavar='kaks file', type=str, default=None, help='')
	parser.add_argument('--pair', metavar='species pair control file', type=str, default=None, help='')
	parser.add_argument('--figure', metavar='output figure file', type=str, default=None, help='[default=--pair + .pdf]')
	parser.add_argument('--max-ks', metavar='max ks', type=float, default=2, help='[default=%(default)s]')
	parser.add_argument('--bins-per-ks', metavar='bins per ks unit', type=int, default=50, help='[default=%(default)s]')
	parser.add_argument('--normed', action='store_true', default=False, help='[default=%(default)s]')
	parser.add_argument('--yn00', action='store_true', default=False, help='turn to YN00[default=%(default)s]')
	parser.add_argument('--method', metavar='meth', type=str, default='NG86', help='[default=%(default)s]]')
	parser.add_argument('--fdtv', action='store_true', default=False, help='turn to 4DTV[default=%(default)s]')
	parser.add_argument('--xlabel', metavar='x label', type=str, default='Ks', help='[default=%(default)s]]')
	parser.add_argument('--homology-class', metavar='class of homology', type=str, default=None, help='allow repeat[default=%(default)s]]')
	parser.add_argument('--only-class', action='store_true', default=False, help='only plot homology in --homology-class. [default=%(default)s]]')
	parser.add_argument('--plot-others', action='store_true', default=False, help='plot homology not in --homology-class but in --kaks as a single line. [default=%(default)s]')
	parser.add_argument('--out-data', action='store_true', default=False, help='output plot data? [default=%(default)s]')
	parser.add_argument('--remove_same_chr', action='store_true', default=False, 
				help='remove data that on the same chromsome [default=%(default)s]')
	parser.add_argument('--gff', type=str, default=None, help='mcscan gff file[default=%(default)s]]')

	args = parser.parse_args()
	if args.fdtv and args.xlabel == 'Ks':
		args.xlabel = '4dTV'
	suffix = args.xlabel + '.yn00' if args.yn00  else args.xlabel
	if args.figure is None:
		if args.pair is not None:
			args.figure = args.pair + '.' + suffix + '.pdf'
		else:
			args.figure = args.kaks + '.' + suffix + '.pdf'
	return args
def main(args=Args()):
##	if args.pair is not None:
	pairs = parse_pair(args.pair)
##	else:
##		pairs = None
##	if args.homology_class is not None:
	homology_class = parse_homology_class(args.homology_class)
	d_gff = parse_gff(args.gff) if args.remove_same_chr else None
##	else:
##		homology_class = None
##	print len(homology_class)
	d_ks = parse_kaks(kaks=args.kaks, pairs=pairs, max_ks=args.max_ks, fdtv=args.fdtv, 
		yn00=args.yn00, method=args.method, d_gff=d_gff,
		homology_class=homology_class, only_class=args.only_class, plot_others=args.plot_others)
	if homology_class is not None and pairs is not None:
		pairs = sorted(sorted(d_ks.keys()), key=lambda x: pairs.index((x[0], x[1])))
#	print 'keys:', d_ks.keys()
	if args.out_data:
		args.out_data = args.figure + '.data.xls'
		with open(args.out_data, 'w') as f:
			print('#parameters:', vars(args), file=f)
	plot_hist(d_ks, args.figure, order=pairs, normed=args.normed, nbins=int(args.max_ks*args.bins_per_ks), xlabel=args.xlabel, out_data=args.out_data, max_ks=args.max_ks)
#	test_diff(d_ks)
	
def parse_homology_class(homology_class):
	if homology_class is None:
		return None
	d = {}
	for line in open(homology_class):
		temp = line.strip().split()
		pair = tuple(temp[:2])
		try: CLASS = temp[2]
		except IndexError: CLASS = None
		if CLASS is not None:
			try: d[pair] += [CLASS]
			except KeyError: d[pair] = [CLASS]
		else:
			d[pair] = []
	return d
def test_diff(d_ks):
	from itertools import combinations
	from scipy import stats
	for p1, p2 in combinations(list(d_ks.keys()), 2):
		v1, v2 = d_ks[p1], d_ks[p2]
		test = stats.kruskal(v1, v2)
		test = stats.ttest_ind(v1, v2)
		print(p1, p2, test.pvalue)

def plot_hist(d_data, outFig, order=None, normed=False, xlabel='Ks', nbins=300, ylabel=' of syntenic gene pairs', out_data=False, max_ks=3, density=True):
	import matplotlib.pyplot as plt
	
	if order is None:
		order = sorted(d_data.keys())
	den_ylabel = 'Percent' + ylabel
	if normed:
		ylabel = 'Percent' + ylabel
	else:
		ylabel = 'Number' + ylabel
	fig, ax = plt.subplots()
	x_list = []
	y_list = []
	pair_list = []
#	for pair, kss in d_data.items():
#		print pair, len(kss)
	print('\t'.join(['species1', 'species2', 'category', 'genePairNumber', 'mean', 'median', 'peak', '95% CI']))
	for pair in order:
		data = d_data[pair]
#		print pair[0], pair[1], len(data), np.mean(data), np.median(data), np.min(data), np.max(data)
		pair_list += [generate_label(pair)]
		data2 = [0, max_ks] + data
		n,bins,patches=ax.hist(data2, bins=nbins, normed=normed,label=pair, )
		y_hist = [(bins[i] + bins[i+1])/2 for i in range(len(bins)-1)]
#		print bins[0], bins[1], bins[-1]
		x_list.append(y_hist)
		y_list.append(n)
		y_hist = format_ks(y_hist, nbins)
		peak_ks = max(list(zip(y_hist,n)), key=lambda x:x[1])[0]
		if out_data:
			with open(out_data, 'a') as f:
				print('\t'.join(list(pair)+['x']+ list(map(str, y_hist)) ), file=f)
				print('\t'.join(list(pair)+['y']+ list(map(str, n)) ), file=f)
		_tile2_5 = np.percentile(data, 2.5)
		_tile97_5 = np.percentile(data, 97.5)
		_ci95 = '{:.3f}-{:.3f}'.format(abs(_tile2_5), _tile97_5)
		category = '-' if len(pair) == 2 else pair[2]
		#print pair
		print('\t'.join(map(str, [ pair[0], pair[1], category, len(data), np.mean(data), np.median(data), peak_ks, _ci95])))
	ax.cla()
	for i in range(len(pair_list)):
		ax.plot(x_list[i],y_list[i],label=pair_list[i],linewidth=1.5)
		legend = ax.legend(fontsize='x-small')
	ax.minorticks_on()
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)
	plt.savefig(outFig)
	
	# density
	prefix = os.path.splitext(outFig)[0]
	datafile = prefix + '.density.data'
	outfig = prefix + '.density.pdf'
	with open(datafile, 'w') as f:
		print('{}\t{}\t{}'.format('pair', 'value', 'homology'), file=f)
		for pair in order:
			label = generate_label(pair)
			if len(set(pair)) == 1:
				logs = 'inparalogs'
			else:
				logs = 'orthologs'
			for value in d_data[pair]:
				print('{}\t{}\t{}'.format(label, value, logs), file=f)
	rsrc = prefix + '.density.r'
	with open(rsrc, 'w') as f:
		print('''datafile = '{datafile}'
data = read.table(datafile, head=T)
library(ggplot2)
p <- ggplot(data, aes(x=value, color=pair, lty=homology)) + geom_line(stat="density", size=1) + xlab('{xlabel}') + ylab('{ylabel}') + xlim(0, {max_ks}) + scale_colour_hue(l=45) + theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.position = c(0.7, 0.6))

ggsave('{outfig}', p, width=7, height=5)

'''.format(datafile=datafile, outfig=outfig, xlabel=xlabel, ylabel=den_ylabel, max_ks=max_ks, ), file=f)
	cmd = 'Rscript {}'.format(rsrc)
	os.system(cmd)
def generate_label(pair):
	if len(pair)> 2:
		return pair[2]
	return '-'.join(pair)
def format_ks(values, bins):
	return [int(round(v*bins,0))*1.0/bins for v in values]
def parse_kaks(kaks, pairs=None, max_ks=None, fdtv=False, yn00=False, method='NG86', 
		homology_class=None, only_class=False, plot_others=False, d_gff=None):
	from .mcscan import KaKs, KaKsParser
#	print pairs
	if pairs is not None:	
		pairs = set(pairs)
	d_ks = {}
#	for line in open(kaks):
	for rc in KaKsParser(kaks):
#		temp = line.strip().split()
#		if temp[0] == 'Sequence':
#			if temp[1] == 'dS-YN00':
#				yn00 = True
#			elif temp[1] == '4D_Sites':
#				fdtv = True
#			continue
#		rc = KaKs(temp, fdtv=fdtv, yn00=yn00, method=method)
		sp1, sp2 = rc.species
		g1, g2 = rc.pair
		if d_gff:
			if d_gff[g1] == d_gff[g2]:
				continue
		if pairs is None:
			pair = (sp1, sp2)
		else:
			if (sp1, sp2) in pairs:
				pair = (sp1, sp2)
			elif (sp2, sp1) in pairs:
				pair = (sp2, sp1)
			else:
				continue
#		print pair
		if rc.ks <= 0:
			continue
		if max_ks is not None and rc.ks > max_ks:
			continue

		Xpairs = [pair] # legend
		if homology_class is not None:
			if only_class:
				Xpairs = []
			else:
				Xpairs = [tuple(list(pair)+['ALL'])]
##			print (g1,g2)
			if (g1, g2) in homology_class:
				classes = homology_class[(g1, g2)]
			elif (g2, g1) in homology_class:
				classes = homology_class[(g2, g1)]
			else:
				if only_class:
					continue
				elif plot_others:
					classes = ['OTHERS']
				else:
					classes = []
			for CLASS in classes:
				Xpairs += [tuple(list(pair)+[CLASS])]
#				Xpairs += [tuple([CLASS])]
#		print Xpairs
		for pair in Xpairs:
			try:
				d_ks[pair] += [rc.ks]
			except KeyError: d_ks[pair] = [rc.ks]
		
	return d_ks
def parse_pair(pair):
	if pair is None:
		return None
	pairs = []
	for line in open(pair):
		if line.startswith('#'):
			continue
		temp = line.strip().split()
		p1, p2 = temp[:2]
		pairs += [(p1, p2)]
	return pairs
def parse_gff(gff):
	d = {}
	for line in open(gff):
		line = line.strip().split()
		chr, g = line[:2]
		d[g] = chr
	return d
if __name__ == '__main__':
	main()
