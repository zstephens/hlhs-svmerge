import re
import copy

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

# write output files
def write_output_file(filtered_variants,VAR_FILT,fn):
	f = open(fn,'w')
	# write variant filter parameters in header
	for k in sorted(VAR_FILT.keys()):
		f.write('##\t'+k+'\t'+str(VAR_FILT[k])+'\n')
	# write all filtered SVs
	for k in sorted(filtered_variants.keys()):
		for svcs in filtered_variants[k]:
			f.write('#\t'+k+'\t'+'\t'.join([str(n) for n in svcs[0][0:3]])+'\n')
			for svc in svcs[1]:
				f.write('\t'.join([str(n) for n in svc])+'\n')
	f.close()

# read input file
def read_input_file(IN_FILE):
	FILTERED_SVS = {}
	f = open(IN_FILE,'r')
	current_list = []
	for line in f:
		splt = line.strip().split('\t')
		if splt[0] == '##':
			continue
		elif splt[0] == '#':
			if len(current_list):
				FILTERED_SVS[samp_name].append([(current_cc,current_lb,current_ub),copy.deepcopy(current_list)])
			samp_name = splt[1]
			if samp_name not in FILTERED_SVS:
				FILTERED_SVS[samp_name] = []
			current_cc = splt[2]
			current_lb = int(splt[3])
			current_ub = int(splt[4])
			current_list = []
		else:
			tup1 = splt[4][1:-1].split(',')
			tup2 = splt[5][1:-1].split(',')
			tup3 = splt[6][1:-1].split(',')
			if 'None' in tup3[1]:
				tup3_parse = [tup3[0][1:-1],None]
			else:
				tup3_parse = [tup3[0][1:-1],float(tup3[1])]
			current_list.append([splt[0],splt[1],int(splt[2]),int(splt[3]),(int(tup1[0]),int(tup1[1])),(int(tup2[0]),int(tup2[1])),(tup3_parse[0],tup3_parse[1])])
	if len(current_list):
			FILTERED_SVS[samp_name].append([(current_cc,current_lb,current_ub),copy.deepcopy(current_list)])
	f.close()
	return FILTERED_SVS


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


###############################################################################
#                                                                             #
#	SV FILTERS                                                                #
#                                                                             #
###############################################################################

# remove exact duplicates
def svfilter_unique(inList):
	countDict = {}
	outList = []
	for i in xrange(len(inList)):
		myKey = (inList[i][0],inList[i][1],inList[i][2],inList[i][3],inList[i][6][0])
		if myKey not in countDict:
			outList.append(inList[i])
		countDict[myKey] = True
	return outList

# remove events on called on reference sequences we're not interested in
def svfilter_chromWhitelist(inList,whitelist):
	outList = []
	for i in xrange(len(inList)):
		ccoord = inList[i][1]
		if ccoord in whitelist:
			outList.append(inList[i])
	return outList

# separate events that intersect bed file (e.g. mappability, gene regions)
#
# intersect=True --> return events that intersect bed file
# intersect=False -> return events that are outside bed regions
#
# input_is_cluster=True --> inList is a list of lists of events (clustered svs) instead of just list of events.
#
# mode="endpoint" --> only consider regions buffered around endpoint coordinates
# mode="span" ------> consider entire span of event
#
def svfilter_bedIntersect(inList,mapTrack,intersect,mode,mapBuff=None,mapMax=None,input_is_cluster=False):
	if mode == "endpoint":
		if mapBuff == None or mapMax == None:
			print "\nError: bedIntersect 'endpoint' method requires specifying mapBuff and mapMax.\n"
			exit(1)
	if mode == "span":
		if mapMax == None:
			print "\nError: bedIntersect 'span' method requires specifying mapMax.\n"
			exit(1)

	outList = []
	for i in xrange(len(inList)):
		if inList[i][0] in ['DUP','DEL','INV']:
			ccoord = inList[i][1]
			lcoord = inList[i][2]# + inList[i][4][0]
			ucoord = inList[i][3]# + inList[i][5][1]

			isInside = False

			if mode == "endpoint":
				lrange = (max([0,lcoord-mapBuff]),lcoord+mapBuff)
				urange = (max([0,ucoord-mapBuff]),ucoord+mapBuff)
				unmap_count_l   = mapTrack.query_range_faster(ccoord,lrange[0],lrange[1],query_endPointInclusive=True)
				unmap_count_u   = mapTrack.query_range_faster(ccoord,urange[0],urange[1],query_endPointInclusive=True)
				unmap_percent_l = unmap_count_l/float(2.*mapBuff)
				unmap_percent_u = unmap_count_u/float(2.*mapBuff)
				if unmap_percent_l > mapMax or unmap_percent_u > mapMax:
					isInside = True

			elif mode == "span":
				unmap_count   = mapTrack.query_range_faster(ccoord,lcoord,ucoord,query_endPointInclusive=True)
				unmap_percent = unmap_count/float(ucoord-lcoord)
				if unmap_percent > mapMax:
					isInside = True

			else:
				print "\nError: Unknown bedIntersect method.\n"
				exit(1)
			
			if (isInside and intersect) or (not(isInside) and not(intersect)):
				outList.append(inList[i])
	return outList

# remove events not supported by enough separate callers
def svfilter_concordance(inCluster,concVal):
	outCluster = []
	for i in xrange(len(inCluster)):
		callersFound = {}
		for n in inCluster[i][1]:
			callersFound[n[6][0]] = True
		if len(callersFound) >= concVal:
			outCluster.append(copy.deepcopy(inCluster[i]))
	return outCluster



