import re

from sys_cmd import *

# merge bcf files and convert to vcf
def bcf_2_vcf(BCFTOOLS, infile_list, outfile):
	if not exists_and_is_nonZero(outfile):
		if BCFTOOLS == '':
			print '\nError: No BCFtools exe specified.\n'
			exit(1)
		if len(infile_list) == 1:
			cmd = BCFTOOLS + ' view -O v -o ' + outfile + ' ' + infile_list[0]
			exe(cmd)
		elif len(infile_list) > 1:
			cmd = BCFTOOLS + ' merge -O v -o ' + outfile + ' ' + ' '.join(infile_list) + ' --force-samples'
			exe(cmd)

###############################################################################
#                                                                             #
#	PARSE SV FILES --> PUT INTO A STANDARD FORMAT:                            #
#                                                                             #
### DEL:                                                                      #
#       ['DEL', 'chr1', s_pos, e_pos, s_buff, e_buff, meta]                   #
#                                                                             #
#       sequence [s_pos+1,e_pos] has been removed with respect to reference   #
#       s_pos is right-most coordinate involved in event                      #
#                                                                             #
### DUP:                                                                      #
#       ['DUP', 'chr1', s_pos, e_pos, s_buff, e_buff, meta]                   #
#                                                                             #
#       sequence [s_pos+1,e_pos] has been duplicated and inserted EITHER at   #
#       s_pos or e_pos+1                                                      #
#                                                                             #
### INS:                                                                      #
#       ['INS', 'chr1', s_pos, len, s_buff, l_buff, meta]                     #
#                                                                             #
#       novel sequence of length "len" is inserted, beginning at s_pos+1      #
#                                                                             #
### INV:                                                                      #
#       ['INV', 'chr1', s_pos, e_pos, s_buff, e_buff, meta]                   #
#                                                                             #
#       sequence [s_pos+1,e_pos] has been inverted with respect to reference  #
#                                                                             #
### BND:                                                                      #
#       ['BND', 'chr1', s_pos, 'chr2', e_pos, meta]                           #
#                                                                             #
#       (chr1, s_pos+1) is connected to (chr2, e_pos)                         #
#                                                                             #
###############################################################################
#                                                                             #
#	ALL COORDINATES HERE ARE 1-INDEXED, EXACTLY LIKE THE VCF FILES            #
#                                                                             #
###############################################################################
#                                                                             #
#	meta = (sv_caller, normalized_delta_read_depth)                           #
#                                                                             #
###############################################################################

#
#
#
def parse_vcf(fn,VARIANT_FILTERS,sv_caller=None):
	f = open(fn,'r')
	sv_list = []
	for line in f:
		if line[0] != '#':
			splt  = line.strip().split('\t')
			myChr = splt[0]
			myPos = int(splt[1])
			sv_type = re.findall(r";SVTYPE=.*?(?=;)",';'+splt[7]+';')[0][8:]
			sv_filt = splt[6]

			if sv_type == 'BND':
				pass
			else:
				myEnd = int(re.findall(r";END=.*?(?=;)",';'+splt[7]+';')[0][5:])
				ciPos = [int(n) for n in re.findall(r";CIPOS=.*?(?=;)",';'+splt[7]+';')[0][7:].split(',')]
				ciEnd = [int(n) for n in re.findall(r";CIEND=.*?(?=;)",';'+splt[7]+';')[0][7:].split(',')]
				sortC = sorted([[myPos,ciPos],[myEnd,ciEnd]])
				sv_size = sortC[1][0] - sortC[0][0]

				myMeta = [sv_caller,None]

				# apply filters
				passed_all_filters = True
				if VARIANT_FILTERS['QUAL']:
					if (sv_filt != '.' and sv_filt != 'PASS'):
						passed_all_filters = False
				if sv_size < VARIANT_FILTERS['SIZE_MIN']:
					passed_all_filters = False
				if sv_size > VARIANT_FILTERS['SIZE_MAX']:
					passed_all_filters = False
				if not passed_all_filters:
					continue

				if sv_type == 'DEL':
					sv_list.append(['DEL',myChr,sortC[0][0],sortC[1][0],tuple(sortC[0][1]),tuple(sortC[1][1]),tuple(myMeta)])
				elif sv_type == 'DUP':
					sv_list.append(['DUP',myChr,sortC[0][0],sortC[1][0],tuple(sortC[0][1]),tuple(sortC[1][1]),tuple(myMeta)])
				elif sv_type == 'INS':
					# I haven't seen any of these thus far, so lets skip it for now
					continue
				elif sv_type == 'INV':
					sv_list.append(['INV',myChr,sortC[0][0],sortC[1][0],tuple(sortC[0][1]),tuple(sortC[1][1]),tuple(myMeta)])
	f.close()
	return sv_list

#
#
#
def parse_cnvnator(fn,VARIANT_FILTERS):
	f = open(fn,'r')
	sv_list = []
	for line in f:
		splt  = line.strip().split('\t')
		splt2 = splt[1].split(':')
		myChr = splt2[0]
		[myPos,myEnd] = [int(n) for n in splt2[1].split('-')]
		[ciPos,ciEnd] = [[0,0],[0,0]]
		sortC   = sorted([[myPos,ciPos],[myEnd,ciEnd]])
		sv_type = splt[0]
		sv_size = myEnd - myPos + 1

		my_norm_readDepth = float(splt[3])

		# apply filters
		passed_all_filters = True
		if sv_size < VARIANT_FILTERS['SIZE_MIN']:
			passed_all_filters = False
		if sv_size > VARIANT_FILTERS['SIZE_MAX']:
			passed_all_filters = False
		if not passed_all_filters:
			continue

		myMeta = ['cnvnator',my_norm_readDepth]

		if sv_type == 'deletion':
			sv_list.append(['DEL',myChr,sortC[0][0],sortC[1][0],tuple(sortC[0][1]),tuple(sortC[1][1]),tuple(myMeta)])
		elif sv_type == 'duplication':
			sv_list.append(['DUP',myChr,sortC[0][0],sortC[1][0],tuple(sortC[0][1]),tuple(sortC[1][1]),tuple(myMeta)])

	f.close()
	return sv_list

#
#
#
def parse_wandy(fn,VARIANT_FILTERS):
	f = open(fn,'r')
	sv_list = []
	isFirstLine = True
	for line in f:
		# skip first line which contains header information
		if isFirstLine:
			isFirstLine = False
			continue

		splt  = line.strip().split('\t')
		myChr = splt[1]
		# need to do int(float()) here because format is stupid and allows ints to be printed in scientific notation
		[myPos,myEnd] = [int(float(splt[2])),int(float(splt[3]))]
		[ciPos,ciEnd] = [[0,0],[0,0]]
		sortC   = sorted([[myPos,ciPos],[myEnd,ciEnd]])
		sv_size = myEnd - myPos

		log2_norm_readDepth = float(splt[5])
		my_norm_readDepth   = 2.**log2_norm_readDepth

		# apply filters
		passed_all_filters = True
		if sv_size < VARIANT_FILTERS['SIZE_MIN']:
			passed_all_filters = False
		if sv_size > VARIANT_FILTERS['SIZE_MAX']:
			passed_all_filters = False
		if log2_norm_readDepth > VARIANT_FILTERS['WANDY_DEL'] and log2_norm_readDepth < VARIANT_FILTERS['WANDY_DUP']:
			passed_all_filters = False
		if not passed_all_filters:
			continue

		myMeta = ['wandy',my_norm_readDepth]

		if log2_norm_readDepth <= VARIANT_FILTERS['WANDY_DEL']:
			sv_list.append(['DEL',myChr,sortC[0][0],sortC[1][0],tuple(sortC[0][1]),tuple(sortC[1][1]),tuple(myMeta)])
		elif log2_norm_readDepth >= VARIANT_FILTERS['WANDY_DUP']:
			sv_list.append(['DUP',myChr,sortC[0][0],sortC[1][0],tuple(sortC[0][1]),tuple(sortC[1][1]),tuple(myMeta)])

	f.close()
	return sv_list


