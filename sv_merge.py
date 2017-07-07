#!/usr/bin/env python
# encoding: utf-8
""" ////////////////////////////////////////////////////////////////////////////////
   ///                                                                          ///
  ///       sv_merge.py                                                        ///
 ///                                                                          ///
///////     Merging and filtering of SV calls for HLHS datasets              //////
   ///                                                                         ///
  ///       Written by:     Zach Stephens                                     ///
 ///        For:            Mayo Clinic (Dr. Timothy Olson et al.)           ///
////////    Date:           July 5, 2017                                    ////////
    ///     Contact:        zstephe2@illinois.edu                               ///
   ///                      zstephens2718@gmail.edu                            ///
  ///                                                                         ///       
 ///                                                                         ///
/////////////////////////////////////////////////////////////////////////////// """

import os
import sys
import copy
import time
import argparse
import cPickle as pickle

# absolute path to this script
SIM_PATH = '/'.join(os.path.realpath(__file__).split('/')[:-1])
sys.path.append(SIM_PATH+'/py/')

from sys_cmd import *
from parse import parse_vcf, parse_cnvnator, parse_wandy, bcf_2_vcf
from mappability import MappabilityTrack
from sv_comparison import sv_filter

DELLY_FILE_NAMES = ['delly.BND.bcf','delly.DEL.bcf','delly.DUP.bcf','delly.INV.bcf']

CHR_WHITELIST = [str(n) for n in xrange(1,23)] + ['X','Y','M','Mt']
CHR_WHITELIST += ['chr'+n for n in CHR_WHITELIST]

# default values
VAR_FILT = {'WHITELIST':True,	# only use variants from desired chromosomes (see CHR_WHITELIST)
            'QUAL':True,		# only use variants with PASS
            'SIZE_MAX':100000,	# maximum size of SV
            'SIZE_MIN':50,		# minimum size of SV
            'WANDY_DEL':-0.6,	# maximum log_2(depth) for deletions in WANDY output
            'WANDY_DUP':0.4,	# minimum log_2(depth) for dups in WANDY output
            'MAP_BUFF':150,		# how much to buffer around breakpoint regions when computing mappability
            'MAP_MAX':0.20,		# maximum percentage of sv allowed to overlap unmappable regions
            'CONC':2}			# how many SV callers must support the same event for it to be considered?

parser = argparse.ArgumentParser(description='sv_merge.py')
parser.add_argument('-i',             type=str,   required=True,  metavar='<str>',                                    help="* base_directory/")
parser.add_argument('--sv-size-min',  type=int,   required=False, metavar='<int>',     default=VAR_FILT['SIZE_MIN'],  help="minimum allowable SV span [%(default)s]")
parser.add_argument('--sv-size-max',  type=int,   required=False, metavar='<int>',     default=VAR_FILT['SIZE_MAX'],  help="maximum allowable SV span [%(default)s]")
parser.add_argument('--wandy-del',    type=float, required=False, metavar='<float>',   default=VAR_FILT['WANDY_DEL'], help="maximum log_2(depth) for WANDY DEL calls [%(default)s]")
parser.add_argument('--wandy-dup',    type=float, required=False, metavar='<float>',   default=VAR_FILT['WANDY_DUP'], help="minimum log_2(depth) for WANDY DUP calls [%(default)s]")
parser.add_argument('--map-track',    type=str,   required=False, metavar='<str>',     default='',                    help="mappability track bed file")
parser.add_argument('--map-buff',     type=int,   required=False, metavar='<int>',     default=VAR_FILT['MAP_BUFF'],  help="mappability buffer around breakpoints [%(default)s]")
parser.add_argument('--map-max',      type=float, required=False, metavar='<float>',   default=VAR_FILT['MAP_MAX'],   help="maximum percentage of unmappable buffer [%(default)s]")
parser.add_argument('--concordance',  type=int,   required=False, metavar='<int>',     default=VAR_FILT['CONC'],      help="minimum concordance (#SV callers, <2 = OFF) [%(default)s]")
parser.add_argument('--vcf-nofilter',             required=False, action='store_true', default=False,                 help='include VCF entries that fail QUAL filter')
parser.add_argument('--no-whitelist',             required=False, action='store_true', default=False,                 help='include SVs from all chr (not just whitelist)')
parser.add_argument('--bcftools',     type=str,   required=False, metavar='<str>',     default='',                    help="/path/to/bcftools")
parser.add_argument('--file-in',      type=str,   required=False, metavar='<str>',     default='',                    help="(skip parsing and read SVs from here instead)")
parser.add_argument('--file-out',     type=str,   required=False, metavar='<str>',     default='',                    help="(save parsed SVs here for speedier future runs)")
args = parser.parse_args()

# TODO: input value sanity-checking code

BASE_INPUT_DIRECTORY = dir_format(args.i)

VAR_FILT['SIZE_MIN']  = args.sv_size_min
VAR_FILT['SIZE_MAX']  = args.sv_size_max
VAR_FILT['WANDY_DEL'] = args.wandy_del
VAR_FILT['WANDY_DUP'] = args.wandy_dup
VAR_FILT['MAP_BUFF']  = args.map_buff
VAR_FILT['MAP_MAX']   = args.map_max
VAR_FILT['CONC']      = args.concordance

if args.vcf_nofilter:
	VAR_FILT['QUAL'] = False

if args.no_whitelist:
	VAR_FILT['WHITELIST'] = False

(IN_FILE, OUT_FILE) = (args.file_in, args.file_out)

BCFTOOLS = args.bcftools
#BCFTOOLS = '/Users/zach/Desktop/bioinformatics/tools/bcftools'

MAPPABILITY_TRACK = args.map_track
if MAPPABILITY_TRACK == '':
	tryDefault = SIM_PATH+'/data/hg19_e2_l400_mappabilityTrack.bed'
	if exists_and_is_nonZero(tryDefault):
		print 'mappability track not specified, using default:'
		print '--',tryDefault
	MAPPABILITY_TRACK = tryDefault
#MAPPABILITY_TRACK = '/Users/zach/Desktop/hlhs_summer_2017/hg19_e2_l300_mappabilityTrack.bed'



################################################################################
#																			   #
#	MAIN																	   #
#																			   #
################################################################################


def main():

	myMappabilityTrack = MappabilityTrack(MAPPABILITY_TRACK)

	FILTERED_SVS = {}
	if IN_FILE == '':
		#
		#	READ INPUT FOLDER
		#
		print 'finding input files...'
		tt = time.time()
		ALL_FILE_NAMES = {}
		listing = [n for n in os.listdir(BASE_INPUT_DIRECTORY) if '_family_cnv_sv_files' in n]
		for family_dir in listing:
			fd = BASE_INPUT_DIRECTORY+family_dir+'/'
			family_prefix = family_dir.split('_')[0]
			listing2 = [n for n in os.listdir(fd) if os.path.isdir(fd+n)]
			all_samp_names = {}
			for sv_caller in listing2:
				svd = fd+sv_caller+'/'
				if sv_caller == 'cnvnator':
					sublisting = [n for n in os.listdir(svd) if os.path.isdir(svd+n)]
					for subdir in sublisting:
						svd_sd = svd+subdir+'/'
						sample_suffix = subdir[-3:]			# HARD-CODED FILENAMES. HOPE THIS DOESN'T BREAK!!
						sample_id     = family_prefix+'-'+sample_suffix
						if sample_suffix not in all_samp_names:
							all_samp_names[sample_suffix] = []
							ALL_FILE_NAMES[sample_id] = {}
						ALL_FILE_NAMES[sample_id]['cnvnator'] = svd_sd+'all_chrs.cnv.bin100.out'
						all_samp_names[sample_suffix].append('cnvnator')
				elif sv_caller == 'delly':
					sublisting = [n for n in os.listdir(svd) if os.path.isdir(svd+n)]
					for subdir in sublisting:
						svd_sd = svd+subdir+'/'
						delly_list = [svd_sd+n for n in os.listdir(svd_sd) if n in DELLY_FILE_NAMES]
						bcf_2_vcf(BCFTOOLS, delly_list, svd_sd+'delly.MERGED.vcf')
						sample_suffix = subdir[-3:]			# HARD-CODED FILENAMES. HOPE THIS DOESN'T BREAK!!
						sample_id     = family_prefix+'-'+sample_suffix
						if sample_suffix not in all_samp_names:
							all_samp_names[sample_suffix] = []
							ALL_FILE_NAMES[sample_id] = {}
						ALL_FILE_NAMES[sample_id]['delly'] = svd_sd+'delly.MERGED.vcf'
						all_samp_names[sample_suffix].append('delly')
				elif sv_caller == 'lumpy':
					vcf_files = [n for n in os.listdir(svd) if n[-4:] == '.vcf']
					for vcf_file in vcf_files:
						sample_suffix = vcf_file[-7:-4]		# HARD-CODED FILENAMES. HOPE THIS DOESN'T BREAK!!
						sample_id     = family_prefix+'-'+sample_suffix
						if sample_suffix not in all_samp_names:
							all_samp_names[sample_suffix] = []
							ALL_FILE_NAMES[sample_id] = {}
						ALL_FILE_NAMES[sample_id]['lumpy'] = svd+vcf_file
						all_samp_names[sample_suffix].append('lumpy')
				elif sv_caller == 'wandy':
					wandy_files = [n for n in os.listdir(svd) if '_CNV_segmentation.txt' in n]
					for wandy_file in wandy_files:
						sample_suffix = wandy_file[10:13]	# HARD-CODED FILENAMES. HOPE THIS DOESN'T BREAK!!
						sample_id     = family_prefix+'-'+sample_suffix
						if sample_suffix not in all_samp_names:
							all_samp_names[sample_suffix] = []
							ALL_FILE_NAMES[sample_id] = {}
						ALL_FILE_NAMES[sample_id]['wandy'] = svd+wandy_file
						all_samp_names[sample_suffix].append('wandy')
			#for k in sorted(all_samp_names.keys()):
			#	print k, sorted(all_samp_names[k])
			#for k in sorted(ALL_FILE_NAMES.keys()):
			#	print '---',k
			#	for k2 in sorted(ALL_FILE_NAMES[k].keys()):
			#		print k2, ALL_FILE_NAMES[k][k2]
			#break

		#
		# CONVERT ALL SVS TO SAME INTERNAL FORMAT
		#
		for samp_name in sorted(ALL_FILE_NAMES.keys()):

			all_svs_for_sample = []
			for sv_method in sorted(ALL_FILE_NAMES[samp_name].keys()):
				
				print 'parsing sample: ['+samp_name+'] - sv_method: ['+sv_method+']'
				if sv_method == 'cnvnator':
					all_svs_for_sample.extend(parse_cnvnator(ALL_FILE_NAMES[samp_name][sv_method],VAR_FILT))

				elif sv_method == 'delly':
					all_svs_for_sample.extend(parse_vcf(ALL_FILE_NAMES[samp_name][sv_method],VAR_FILT,sv_caller='delly'))

				elif sv_method == 'lumpy':
					all_svs_for_sample.extend(parse_vcf(ALL_FILE_NAMES[samp_name][sv_method],VAR_FILT,sv_caller='lumpy'))

				elif sv_method == 'wandy':
					all_svs_for_sample.extend(parse_wandy(ALL_FILE_NAMES[samp_name][sv_method],VAR_FILT))

			print '--',len(all_svs_for_sample),'total SVs...'

			#
			# PRUNE IDENTICAL VARIANTS
			#
			countDict = {}
			all_svs_for_sample_unique = []
			for i in xrange(len(all_svs_for_sample)):
				myKey = (all_svs_for_sample[i][0],all_svs_for_sample[i][1],all_svs_for_sample[i][2],all_svs_for_sample[i][3],all_svs_for_sample[i][6][0])
				if myKey not in countDict:
					all_svs_for_sample_unique.append(all_svs_for_sample[i])
				countDict[myKey] = True
			del all_svs_for_sample
			print '--',len(all_svs_for_sample_unique),'unique SVs...'

			#
			# FILTER BY CHROMOSOME
			#
			print '* applying chromosome whitelist filter...'
			all_svs_for_sample_whitelist = []
			for i in xrange(len(all_svs_for_sample_unique)):
				ccoord = all_svs_for_sample_unique[i][1]
				if ccoord in CHR_WHITELIST:
					all_svs_for_sample_whitelist.append(all_svs_for_sample_unique[i])
			del all_svs_for_sample_unique
			print '--',len(all_svs_for_sample_whitelist),'SVs...'

			#
			# FILTER BY MAPPABILITY
			#
			print '* applying mappability filter...'
			all_svs_for_sample_mappability_filtered = []
			for i in xrange(len(all_svs_for_sample_whitelist)):
				if all_svs_for_sample_whitelist[i][0] in ['DUP','DEL','INV']:
					ccoord = all_svs_for_sample_whitelist[i][1]
					lcoord = all_svs_for_sample_whitelist[i][2]# + all_svs_for_sample_whitelist[i][4][0]
					ucoord = all_svs_for_sample_whitelist[i][3]# + all_svs_for_sample_whitelist[i][5][1]
					lrange = (max([0,lcoord-VAR_FILT['MAP_BUFF']]),lcoord+VAR_FILT['MAP_BUFF'])
					urange = (max([0,ucoord-VAR_FILT['MAP_BUFF']]),ucoord+VAR_FILT['MAP_BUFF'])
					#unmap_count = myMappabilityTrack.query_range_faster(ccoord,lcoord,ucoord,query_endPointInclusive=True)
					#unmap_perct = 100.*unmap_count/float(ucoord-lcoord+1)

					unmap_count_l = myMappabilityTrack.query_range_faster(ccoord,lrange[0],lrange[1],query_endPointInclusive=True)
					unmap_count_u = myMappabilityTrack.query_range_faster(ccoord,urange[0],urange[1],query_endPointInclusive=True)
					unmap_percent_l = unmap_count_l/float(2.*VAR_FILT['MAP_BUFF'])
					unmap_percent_u = unmap_count_u/float(2.*VAR_FILT['MAP_BUFF'])

					if unmap_percent_l <= VAR_FILT['MAP_MAX'] and unmap_percent_u <= VAR_FILT['MAP_MAX']:
						all_svs_for_sample_mappability_filtered.append(all_svs_for_sample_whitelist[i])
			del all_svs_for_sample_whitelist
			print '--',len(all_svs_for_sample_mappability_filtered),'SVs...'

			#
			# SORT AND CLUSTER EVENTS
			#
			print '* sorting and clustering by position...'
			sorted_svs  = sorted(all_svs_for_sample_mappability_filtered)
			sv_clusters = [[sorted_svs[0]]]
			for i in xrange(1,len(sorted_svs)):
				inPrevCluster = False
				for n in sv_clusters[-1]:
					x_chr = n[1]
					y_chr = sorted_svs[i][1]
					x1    = n[2] + n[4][0]
					x2    = n[3] + n[5][1]
					y1    = sorted_svs[i][2] + sorted_svs[i][4][0]
					y2    = sorted_svs[i][3] + sorted_svs[i][5][1]
					if x_chr == y_chr and x1 <= y2 and y1 <= x2:
						inPrevCluster = True
				if inPrevCluster:
					sv_clusters[-1].append(sorted_svs[i])
				else:
					sv_clusters.append([sorted_svs[i]])
			del sorted_svs
			print '--',len(sv_clusters),'SVs...'

			####sorted_clusters = []
			####for i in xrange(len(sv_clusters)):
			####	lb = min([n[2]+n[4][0] for n in sv_clusters[i]])
			####	ub = min([n[3]+n[5][1] for n in sv_clusters[i]])
			####	sorted_clusters.append([(sv_clusters[i][0][1],lb,ub),copy.deepcopy(sv_clusters[i])])
			####sorted_clusters = sorted(sorted_clusters)
			####del sv_clusters

			#
			# APPLY CONCORDANCE FILTER
			#
			print '* applying concordance filter...'
			FILTERED_SVS[samp_name] = []
			for i in xrange(len(sv_clusters)):
				callersFound = {}
				for n in sv_clusters[i]:
					callersFound[n[6][0]] = True
				if len(callersFound) >= VAR_FILT['CONC']:
					lb = min([n[2]+n[4][0] for n in sv_clusters[i]])
					ub = min([n[3]+n[5][1] for n in sv_clusters[i]])
					FILTERED_SVS[samp_name].append([(sv_clusters[i][0][1],lb,ub),copy.deepcopy(sv_clusters[i])])
			del sv_clusters
			FILTERED_SVS[samp_name] = sorted(FILTERED_SVS[samp_name])
			print '--',len(FILTERED_SVS[samp_name]),'SVs...'

		#print time.time()-tt,'(sec)'

	#
	# OR SKIP ALL THAT AND READ IN FROM FILE, IF SPECIFIED.
	#
	else:
		print 'reading input file...'
		tt = time.time()
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
		print time.time()-tt,'(sec)'

	####for k in sorted(FILTERED_SVS.keys()):
	####	print '----',k
	####	for n in FILTERED_SVS[k]:
	####		print n
	####exit(1)

	#
	# SAVE SORTED/FILTERED SVS TO OUTPUT FILE, IF DESIRED.
	#
	if OUT_FILE != '':
		print 'saving to output file...'
		tt = time.time()
		f = open(OUT_FILE,'w')
		# write variant filter parameters in header
		for k in sorted(VAR_FILT.keys()):
			f.write('##\t'+k+'\t'+str(VAR_FILT[k])+'\n')
		# write all filtered SVs
		for k in sorted(FILTERED_SVS.keys()):
			for svcs in FILTERED_SVS[k]:
				f.write('#\t'+k+'\t'+'\t'.join([str(n) for n in svcs[0][0:3]])+'\n')
				for svc in svcs[1]:
					f.write('\t'.join([str(n) for n in svc])+'\n')
		f.close()
		print time.time()-tt,'(sec)'


	for k in sorted(FILTERED_SVS.keys()):
		print k, len(FILTERED_SVS[k])
	print '\n********************\n'

	#
	#
	#
	lab = [n[1][0][0] for n in FILTERED_SVS['011H-001']]
	print 'DEL', lab.count('DEL')
	print 'DUP', lab.count('DUP')
	print 'INV', lab.count('INV')

	print ''
	lens = [len(n[1]) for n in FILTERED_SVS['011H-001']]
	for i in xrange(1,5):
		print i, lens.count(i)

	one_count = {}
	for n in FILTERED_SVS['011H-001']:
		nCallers = {}
		for m in n[1]:
			nCallers[m[6][0]] = True
		if len(nCallers) == 1:
			k = nCallers.keys()[0]
			if k not in one_count:
				one_count[k] = 0
			one_count[k] += 1
	for k in one_count.keys():
		print k, one_count[k],'{0:.2f}%'.format(100.*float(one_count[k])/8946.)

	# remove father
	print '\n*', len(FILTERED_SVS['011H-001']), '-->',
	de_novo = sv_filter(FILTERED_SVS['011H-001'],FILTERED_SVS['011H-003'])
	print len(de_novo),'-->',
	# remove mother
	de_novo = sv_filter(de_novo,FILTERED_SVS['011H-004'])
	print len(de_novo),'-->',
	# remove sister
	de_novo = sv_filter(de_novo,FILTERED_SVS['011H-100'])
	print len(de_novo),'-->',
	# remove brother
	de_novo = sv_filter(de_novo,FILTERED_SVS['011H-101'])
	print len(de_novo),'de novo SVs'


	lab = [n[1][0][0] for n in de_novo]
	print 'DEL', lab.count('DEL')
	print 'DUP', lab.count('DUP')
	print 'INV', lab.count('INV')

	print ''
	lens = [len(n[1]) for n in de_novo]
	for i in xrange(1,5):
		print i, lens.count(i)

	for i in xrange(len(de_novo)):
		if len(de_novo[i][1]) == 4:
			print ''
			print de_novo[i][1]

	one_count = {}
	for n in de_novo:
		nCallers = {}
		for m in n[1]:
			nCallers[m[6][0]] = True
		if len(nCallers) == 1:
			k = nCallers.keys()[0]
			if k not in one_count:
				one_count[k] = 0
			one_count[k] += 1
	for k in one_count.keys():
		print k, one_count[k],'{0:.2f}%'.format(100.*float(one_count[k])/3890.)


if __name__ == '__main__':

	#myTrack = MappabilityTrack('/Users/zach/Desktop/hlhs_summer_2017/test.bed')
	#print (80,100), 1, myTrack.query_range('chr1',80,100)
	#print (80,120), 20, myTrack.query_range('chr1',80,120)
	#print (120,140), 11, myTrack.query_range('chr1',120,140)
	#print (130,140), 1, myTrack.query_range('chr1',130,140)
	#print (120,220), 30, myTrack.query_range('chr1',120,220)
	#exit(1)

	main()


