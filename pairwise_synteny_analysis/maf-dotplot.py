#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
-----
Program to make a dot-plot from a pairwise MAF alignment.

Author: Pedro Almeida
-----
'''

import sys, gzip, argparse
from operator import itemgetter
from collections import defaultdict
import matplotlib
import matplotlib.collections as mcoll
import matplotlib.patches as patches


def sortQaligns(m, maxgap=10000):
	'''
	return a list of query names sorted by reference start position
	this assumes there's only one reference sequence

	m = [name1, name2, [start1, end1], [start2, end2]]
	'''
	refStart = {}
	qAlignLen = {}
	start2, end2 = {}, {}
	bestQ = defaultdict(int)
	prevEnd2 = defaultdict(lambda: float('inf'))
	prevAlignLen = defaultdict(int)

	for row in m:
		qAlignLen[row[1]] = abs(row[3][1] - row[3][0] + 1)

		start2[row[1]], end2[row[1]] = sorted((row[2][0], row[2][1]))

		if row[1] not in refStart:
			# if the first block, keep this refStart
			refStart[row[1]] = start2[row[1]]

		if start2[row[1]] - maxgap <= prevEnd2[row[1]]:
			qAlignLen[row[1]] += prevAlignLen[row[1]]

		if bestQ[row[1]] < qAlignLen[row[1]]:
			bestQ[row[1]] = qAlignLen[row[1]]
			# when blocks are "broken", add the new block refStart
			refStart[row[1]] = start2[row[1]]

		prevAlignLen[row[1]] = qAlignLen[row[1]]
		prevEnd2[row[1]] = end2[row[1]]

	q_view = [(v, k) for k, v in refStart.items()]
	q_view.sort()  # natively sort tuples by first element
	orderedNames = []
	for v, k in q_view:
		orderedNames.append(k)

	return orderedNames


def updateCoords(m, dR, dQ):
	'''
	increment the coordinates in matrix
	'''
	for i, row in enumerate(m):
		if dR:
			m[i][2] = [x + dR[row[0]] for x in row[2] if row[0] in dR]
		elif dQ:
			m[i][3] = [x + dQ[row[1]] for x in row[3] if row[1] in dQ]

	return m


def updateChrSizes(d, lst):
	''' uses order in lst '''
	upD = {}
	add = 0
	last = lst[-1]
	for i in lst:
		if i == last:
			upD[i] = add
		else:
			upD[i] = add
		try:
			add += d[i]
		except KeyError:
			pass

	return upD


def uniqElements(lst):
	''' returns unique elements from list preserving order '''
	uniq = []
	[uniq.append(x) for x in lst if x not in uniq]
	return uniq


def chrSizes(d, name, size):
	if name not in d:
		d[name] = size
	else:
		if d[name] != size:
			raise Exception('Sequences have different lengths: {}'.format(name))
	return d


def formatAlignments(aligns):
	'''
	returns a coordinate-like format
	startR, endR, startQ, endQ, nameR, lenR, nameQ, lenQ
	'''
	for aln in aligns:
		b1 = aln[4] + 1  # assumes seq1 is always forward
		e1 = aln[4] + aln[6]
		if aln[5] < 0:  # seq2 in reverse orientation
			b2 = -(aln[5] + aln[6])
			e2 = b2 + aln[6] - 1
			yield b1, e1, e2, b2, aln[0], aln[1], aln[2], aln[3]
		else:
			b2 = aln[5] + 1
			e2 = aln[5] + aln[6]
			yield b1, e1, b2, e2, aln[0], aln[1], aln[2], aln[3]


def sortAlignments(alignments):
	''' sort all alignments by ref name and start '''
	alignments.sort(key=itemgetter(0, 4))
	return alignments


def destructBlocks(block):
	'''
	Get alignments and sequence lengths

	Return:
	targetName targetLength queryName queryLength [(alignBlocks)]

	where alignBlocks is a list of tuples consisting of:
	targetStart queryStart gaplessLength

	Code from "last-dotplot" : http://last.cbrc.jp
	'''
	mafCount = 0
	for b in block:
		w = b.split()
		if w[0] == "s":  # MAF sequence line
			if mafCount == 0:
				chr1, beg1, seqlen1, seq1 = w[1], int(w[2]), int(w[5]), w[6]
				if w[4] == '-':
					beg1 -= seqlen1
				mafCount = 1
			else:
				chr2, beg2, seqlen2, seq2 = w[1], int(w[2]), int(w[5]), w[6]
				if w[4] == '-':
					beg2 -= seqlen2
				blocks = mafBlocks(beg1, beg2, seq1, seq2)
				yield chr1, seqlen1, chr2, seqlen2, blocks
				mafCount = 0


def mafBlocks(beg1, beg2, seq1, seq2):
	'''
	Get the gapless blocks of an alignment, from MAF format.
	Code from "last-dotplot" : http://last.cbrc.jp
	'''
	size = 0  # length of each alignment, excluding gaps
	for x, y in zip(seq1, seq2):
		if x == '-':  # a gap
			if size:
				yield beg1, beg2, size
				beg1 += size
				beg2 += size
				size = 0
			beg2 += 1
		elif y == '-':  # a gap
			if size:
				yield beg1, beg2, size
				beg1 += size
				beg2 += size
				size = 0
			beg1 += 1
		else:
			size += 1

	# the last alignment block
	if size:
		yield beg1, beg2, size


def bucketAlignments(alignments, blocks):
	n1, l1, n2, l2, bcs = blocks
	for b in bcs:
		alignments.append((n1, l1, n2, l2, *b))
	return alignments


def readMAF(fp):
	'''
	Read the contents of a MAF alignment file.
	Returns a generator with blocks (a list).
	The block has size 1 if it's a comment otherwise size is always > 1

	# with open(filename, 'r') as fp:
	# 	for block in readMAF(fp):
		...
	'''

	maf_words = {'a', 's', 'q', 'p', 'i', 'e'}
	l_return = []
	for line in fp:
		line = line.strip()
		try:
			line = line.decode('utf-8')
		except AttributeError:
			pass
		if line:
			c = line[0]
			if c in maf_words:
				if c == 'a':
					if l_return:
						yield l_return
						l_return = []  # reset list
				l_return.append(line)
			elif c == '#':
				yield [line]
			else:
				continue
	# return the last block
	# if there are comment lines below this, those are returned before
	yield l_return


def readGaplessMAF(fname, namesR=None, namesQ=None, bestStrand=False):
	'''
	Returns a matrix of gapless alignments, chromosome sizes in Ref, chromosome sizes in Qry
	'''
	aligns = []
	with fopen(fname, 'r') as fp:
		for block in readMAF(fp):
			if len(block) == 1:  # skip comment lines
				continue
			else:
				for alignBlocks in destructBlocks(block):
					aligns = bucketAlignments(aligns, alignBlocks)

	sizesR = {}
	sizesQ = {}
	fwQaligns = defaultdict(int)
	rvQaligns = defaultdict(int)
	mat = []     # R | Q | Y | X (note that X-axis is related to seq2/query)
	for startR, endR, startQ, endQ, nameR, sizeR, nameQ, sizeQ in formatAlignments(
																		sortAlignments(aligns)):
		if namesQ and nameQ not in namesQ:
			continue
		if namesR and nameR not in namesR:
			continue

		sizesR = chrSizes(sizesR, nameR, sizeR)
		sizesQ = chrSizes(sizesQ, nameQ, sizeQ)
		mat.append([nameR, nameQ, [startR, endR], [startQ, endQ]])

		if bestStrand:
			# count the total length of forward/reverse alignments for each query
			if startQ > endQ:
				rvQaligns[nameQ] += startQ - endQ + 1
			else:
				fwQaligns[nameQ] += endQ - startQ + 1

	if bestStrand:
		for i, row in enumerate(mat):
			if rvQaligns[row[1]] > fwQaligns[row[1]]:
				# if length of alignments in reverse strand is higher than in forward
				# just assume the query is in reverse compliment
				mat[i][3] = [sizesQ[row[1]] - row[3][1] + 1, sizesQ[row[1]] - row[3][0] + 1][::-1]

	return mat, sizesR, sizesQ


def readBED(fname, sizesR, sizesQ):
	bedR = []
	bedQ = []
	with fopen(fname, 'r') as fp:
		for line in fp:
			line = line.strip('\r\n').split()
			# n, s, e, c = line[0], int(line[1]), int(line[2]), line[3]
			n, s, e = line[0], int(line[1]), int(line[2])
			try:
				c = line[3]
			except IndexError:
				c = 'lightgrey'
			if line[0] in sizesR:
				bedR.append((n, s + sizesR[n], e + sizesR[n], c))
			if line[0] in sizesQ:
				bedQ.append((n, s + sizesQ[n], e + sizesQ[n], c))
	return bedR, bedQ


def fopen(file, mode):
	''' Opens a (gzip) file in open mode '''
	if file == '-':
		fobj = sys.stdin
	elif file.lower().endswith('.gz'):
		fobj = gzip.open(file, mode)
	else:
		fobj = open(file, mode)
	return fobj


def drawVHlines(coords, ax, axis=0, **kwargs):
	'''
	axis=0, vertical lines
	axis=1, horizontal lines
	'''
	assert axis in {0, 1}, 'axis must be 0 or 1'

	if axis == 0:
		[ax.axvline(c, **kwargs) for c in coords]
	else:
		[ax.axhline(c, **kwargs) for c in coords]


def setPlotLims(ax, xlim=None, ylim=None):
	if xlim:
		ax.set(xlim=(xlim[0], xlim[1]))
	if ylim:
		ax.set(ylim=(ylim[0], ylim[1]))


def setTicks(addSizes, chrSizes, names, ax, axis=0, **kwargs):
	'''
	axis=0, x-axis
	axis=1, y-axis
	'''
	_ticks, labs = list(addSizes.values()), list(addSizes.keys())
	_ticks = _ticks + [_ticks[-1] + chrSizes[names[-1]]]  # required to make last tick
	ticks = [(_ticks[i] + _ticks[i + 1]) / 2 for i in range(0, len(_ticks) - 1)]  # make ticks in middle

	if axis == 0:
		ax.set_xticks(ticks)
		ax.set_xticklabels(labs, **kwargs)
	elif axis == 1:
		ax.set_yticks(ticks)
		ax.set_yticklabels(labs, **kwargs)
	ax.tick_params(axis='both', length=0)


def addRectangle(ax, bed, axis, rectmax):
	'''
	axis=0, x-axis (vertical)
	axis=1, y-axis (horizontal)
	rectmax is the height/width depending on the axis (x/y respect.)
	'''
	try:
		kwargs = {'linewidth': 0, 'facecolor': bed[3], 'alpha': bed[4]}
	except IndexError:
		kwargs = {'linewidth': 0, 'facecolor': bed[3], 'alpha': 0.4}
	if axis == 0:
		rect = patches.Rectangle((bed[1], 0), bed[2] - bed[1], rectmax, **kwargs)
	elif axis == 1:
		rect = patches.Rectangle((0, bed[1]), rectmax, bed[2] - bed[1], **kwargs)

	ax.add_patch(rect)


def get_parser():
	''' get parser object for script '''

	parser = argparse.ArgumentParser(
			description=__doc__,
			formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument(
		'maf', metavar='MAF',
		help='input maf file (reads from pipe with -)')
	parser.add_argument(
		'out', metavar='FILE',
		help='output file name')
	parser.add_argument(
		'--png', action='store_true',
		help='generate plots in png format instead of pdf')

	data = parser.add_argument_group(title='Data options')
	data.add_argument(
			'-q', '--qry', metavar='STR', action='append',
			help='name of query sequences to plot (can be given multiple times)')
	data.add_argument(
			'-r', '--ref', metavar='STR', action='append',
			help='name of reference sequences to plot (can be given multiple times)')
	data.add_argument(
			'-b', '--bed', metavar='FILE',
			help='BED file with regions to shade')
	data.add_argument(
			'-g', '--maxgap', metavar='INT', type=int, default=10000,
			help='cluster neighbour alignments with less than INT gap size [%(default)s]')
	data.add_argument(
			'-o', '--orient', action='store_true',
			help='reverse strands if length of reverse alignments is longer [%(default)s]')

	aesthetics = parser.add_argument_group(title='Aesthetics options')
	aesthetics.add_argument(
			'-f', '--fontsize', metavar='INT', type=int, default=6,
			help='font size for chromosome/contig names [%(default)s]')
	aesthetics.add_argument(
			'-l', '--linewidth', metavar='FLOAT', type=float, default=0.75,
			help='line width for connecting alignments [%(default)s]')
	aesthetics.add_argument(
			'-c', '--fcol', metavar='STR', default='b',
			help='colour for forward alignments [%(default)s]')
	aesthetics.add_argument(
			'-C', '--rcol', metavar='STR', default='r',
			help='colour for reverse alignments [%(default)s]')
	aesthetics.add_argument(
			'-V', '--no_vlines', action='store_true',
			help='do not draw vertical lines for query')
	aesthetics.add_argument(
			'-t', '--ticks', metavar='INT', type=float, default=0,
			help='convert to INT and draw length ticks on axis instead of contig names [%(default)s]')
	aesthetics.add_argument(
			'-T', '--title', metavar='STR', help='title for the plot')
	aesthetics.add_argument(
			'-F', '--figsize', metavar='LIST', default='8,6',
			help='size of figure in inches (width,height) [%(default)s]')

	args = parser.parse_args()

	return args


if __name__ == '__main__':

	args = get_parser()

	if args.png:
		matplotlib.use('agg')
		plot_kwargs = {'format': 'png', 'dpi': 400}
	else:
		matplotlib.use('pdf')
		plot_kwargs = {'format': 'pdf'}

	import matplotlib.pyplot as plt

	figsize = list(map(float, args.figsize.split(',')))
	namesQ = args.qry  # read from args
	namesR = args.ref  # read from args

	# ----- Read MAF into matrix
	alignMat, sizesR, sizesQ = readGaplessMAF(args.maf, namesR, namesQ, args.orient)

	# ----- use R|Q names/order from file, otherwise use input order
	if not namesR:
		namesR = uniqElements(list(zip(*alignMat))[0])
	addSizesR = updateChrSizes(sizesR, namesR)
	alignMat = updateCoords(alignMat, addSizesR, None)

	if not namesQ:
		# if names is not given sort by ref start
		namesQ = sortQaligns(alignMat, args.maxgap)
	addSizesQ = updateChrSizes(sizesQ, namesQ)
	alignMat = updateCoords(alignMat, None, addSizesQ)

	# ----- Read BED file
	if args.bed:
		bedR, bedQ = readBED(args.bed, addSizesR, addSizesQ)

	# ----- dot-plot -----
	fig, ax = plt.subplots(figsize=figsize)

	# thanks to unutbu @ https://stackoverflow.com/questions/53035858/plot-multiple-values-with-matplotlib-without-loop/53036005#53036005
	# LineCollection wants pairs of points of the form ((x0, y0), (x1, y1))
	segments = [[[r[3][0], r[2][0]], [r[3][1], r[2][1]]] for r in alignMat]
	colors = [args.fcol if r[3][0] <= r[3][1] else args.rcol for r in alignMat]
	line_segments = mcoll.LineCollection(
								segments, colors=colors,
								linewidth=args.linewidth, capstyle='round')

	ax.add_collection(line_segments)

	# ----- add chr/scaffold boundaries
	vhl_kwards = {'color': 'k', 'linestyle': '-', 'linewidth': 0.25, 'alpha': 0.5}
	drawVHlines(list(addSizesR.values()), ax, axis=1, **vhl_kwards)
	if not args.no_vlines:
		drawVHlines(list(addSizesQ.values()), ax, axis=0, **vhl_kwards)

	# ----- set x/y lims
	xMax = sum(sizesQ.values())
	yMax = sum(sizesR.values())
	setPlotLims(ax, xlim=(1, xMax), ylim=(1, yMax))

	# ----- ticks
	if args.ticks > 1:
		# instead of drawing chromosome/contig names, draws ticks and tick labels
		# most useful when plotting a single chromosome
		_ = ax.set_xticklabels(
					['{}'.format(int(x / args.ticks)) for x in ax.get_xticks().tolist()],
					fontsize=args.fontsize)
		_ = ax.set_yticklabels(
					['{}'.format(int(x / args.ticks)) for x in ax.get_xticks().tolist()],
					fontsize=args.fontsize)
	else:
		xticks_kwargs = {'rotation': 45, 'fontsize': args.fontsize, 'horizontalalignment': 'right'}
		yticks_kwargs = {'rotation': 0, 'fontsize': args.fontsize}
		setTicks(addSizesQ, sizesQ, namesQ, ax, axis=0, **xticks_kwargs)
		setTicks(addSizesR, sizesR, namesR, ax, axis=1, **yticks_kwargs)

	# ----- add rectangles
	if args.bed:
		for b in bedQ:
			addRectangle(ax, b, 0, yMax)
		for b in bedR:
			addRectangle(ax, b, 1, xMax)

	# ----- Finalize figure
	if args.title:
		ax.set_title(args.title)

	tightlayout = dict(pad=1, w_pad=0, h_pad=0.0)
	fig.set_tight_layout(tightlayout)
	fig.savefig(args.out, **plot_kwargs)

