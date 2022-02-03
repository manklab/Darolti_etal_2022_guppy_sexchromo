#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
-----
Program to filter maf alignment blocks for E values.
See: http://last.cbrc.jp/doc/last-evalues.html

Note that as for now, the smallest e-value is 1e-300.

Author: Pedro Almeida
-----
'''

import sys, gzip, argparse


def maf_mismapFilter(block, threshold=1e-5):
	mm = block[0].split()[-1]
	if mm.startswith('E='):
		mm = float(mm.replace('E=', ''))
		if mm <= threshold:
			return True
		else:
			return False
	else:
		print('ERROR: could not find e-value in MAF block!')
		sys.exit(1)


def maf_qCovFilter(block, cov=1, perc=None):
	# this just works for 2 sequence lines in maf block
	alen = int(block[2].split()[3])
	sSize = int(block[2].split()[5])
	if perc is None:
		if alen < cov:
			return False
	else:
		pqCov = (alen / sSize) * 100
		if pqCov < perc:
			return False

	return True


def write_maf_block(block):
	[sys.stdout.write('{}\n'.format(line)) for line in block]
	sys.stdout.write('\n')


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


def fopen(file, mode):
	''' Opens a (gzip) file in open mode '''
	if file == '-':
		fobj = sys.stdin
	elif file.lower().endswith('.gz'):
		fobj = gzip.open(file, mode)
	else:
		fobj = open(file, mode)
	return fobj


def get_parser():
	''' get parser object for script '''

	epilog = '''
Notes:
* filters are applied either to query alignment length (-c or -p) OR mismap probability (-m)
* in the former, by default filters are based on aligned query bp unless -p/--perc is specified
	'''

	parser = argparse.ArgumentParser(
			description=__doc__,
			epilog=epilog,
			formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument(
		'maf', metavar='FILE',
		help='input maf file')
	parser.add_argument(
		'-e', '--eval', metavar='FLOAT', type=float,
		help='maximum e-value threshold for block alignments')
	parser.add_argument(
		'-c', '--cov', metavar='INT', type=int,
		help='minimum aligned query bp in each block')
	parser.add_argument(
		'-p', '--perc', metavar='INT', type=int, choices=range(0, 101),
		help='minimum percentage of aligned query in each block')

	args = parser.parse_args()
	if not args.cov and not args.perc and args.eval is None:
		parser.error('need to choose one of -c/--cov OR -p/--perc OR -e/--eval')
	return args


if __name__ == '__main__':

	args = get_parser()

	reject = 0
	with fopen(args.maf, 'r') as fp:
		for block in readMAF(fp):
			if len(block) == 1:
				sys.stdout.write('{}\n'.format(*block))
			else:
				# if block[1].split()[1] == block[2].split()[1]:
				# 	# skip self alignments, if present
				# 	continue

				if args.cov or args.perc:
					if maf_qCovFilter(block, args.cov, args.perc):
						write_maf_block(block)
					else:
						reject += 1
				elif args.eval is not None:
					if maf_mismapFilter(block, args.eval):
						write_maf_block(block)
					else:
						reject += 1

	print('# rejected blocks: {}'.format(reject), file=sys.stderr)

