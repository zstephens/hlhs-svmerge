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
from sv_comparison import subtractive_filter, intersection_filter

DELLY_FILE_NAMES = ['delly.BND.bcf','delly.DEL.bcf','delly.DUP.bcf','delly.INV.bcf']

# default values
VAR_FILT = {'QUAL':True,		# only use variants with PASS
            'SIZE_MAX':100000,	# maximum size of SV
            'SIZE_MIN':50,		# minimum size of SV
            'WANDY_DEL':-0.6,	# maximum log_2(depth) for deletions in WANDY output
            'WANDY_DUP':0.4,	# minimum log_2(depth) for dups in WANDY output
            'MAP_BUFF':150,		# how much to buffer around breakpoint regions when computing mappability
            'MAP_MAX':0.20}	# maximum percentage of sv allowed to overlap unmappable regions

parser = argparse.ArgumentParser(description='sv_merge.py')
parser.add_argument('-i',             type=str,   required=True,  metavar='<str>',                                    help="* base_directory/")
parser.add_argument('--sv-size-min',  type=int,   required=False, metavar='<int>',     default=VAR_FILT['SIZE_MIN'],  help="minimum allowable SV span [%(default)s]")
parser.add_argument('--sv-size-max',  type=int,   required=False, metavar='<int>',     default=VAR_FILT['SIZE_MAX'],  help="maximum allowable SV span [%(default)s]")
parser.add_argument('--wandy-del',    type=float, required=False, metavar='<float>',   default=VAR_FILT['WANDY_DEL'], help="maximum log_2(depth) for WANDY DEL calls [%(default)s]")
parser.add_argument('--wandy-dup',    type=float, required=False, metavar='<float>',   default=VAR_FILT['WANDY_DUP'], help="minimum log_2(depth) for WANDY DUP calls [%(default)s]")
parser.add_argument('--map-track',    type=str,   required=False, metavar='<str>',     default='',                    help="mappability track bed file")
parser.add_argument('--map-buff',     type=int,   required=False, metavar='<int>',     default=VAR_FILT['MAP_BUFF'],  help="mappability buffer around breakpoints [%(default)s]")
parser.add_argument('--map-max',      type=float, required=False, metavar='<float>',   default=VAR_FILT['MAP_MAX'],   help="maximum percentage of unmappable buffer [%(default)s]")
parser.add_argument('--vcf-nofilter',             required=False, action='store_true', default=False,                 help='include VCF entries that fail QUAL filter')
parser.add_argument('--bcftools',     type=str,   required=False, metavar='<str>',     default='',                    help="/path/to/bcftools")
parser.add_argument('--pickle-in',    type=str,   required=False, metavar='<str>',     default='',                    help="(skip parsing and read SVs from here instead)")
parser.add_argument('--pickle-out',   type=str,   required=False, metavar='<str>',     default='',                    help="(save parsed SVs here for speedier future runs)")
args = parser.parse_args()

# TODO: input value sanity-checking code

# this is the worst variable name I could possibly use, but I'm doing it anyway for brevity!!
D = dir_format(args.i)

if args.vcf_nofilter:
	VAR_FILT['QUAL'] = False

(IN_PICKLE, OUT_PICKLE) = (args.pickle_in, args.pickle_out)

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

	if IN_PICKLE == '':
		#
		#	READ INPUT FOLDER
		#
		print 'finding input files...'
		tt = time.time()
		IN_FILES = {}
		listing = [n for n in os.listdir(D) if '_family_cnv_sv_files' in n]
		for family_dir in listing:
			fd = D+family_dir+'/'
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
							IN_FILES[sample_id] = {}
						IN_FILES[sample_id]['cnvnator'] = svd_sd+'all_chrs.cnv.bin100.out'
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
							IN_FILES[sample_id] = {}
						IN_FILES[sample_id]['delly'] = svd_sd+'delly.MERGED.vcf'
						all_samp_names[sample_suffix].append('delly')
				elif sv_caller == 'lumpy':
					vcf_files = [n for n in os.listdir(svd) if n[-4:] == '.vcf']
					for vcf_file in vcf_files:
						sample_suffix = vcf_file[-7:-4]		# HARD-CODED FILENAMES. HOPE THIS DOESN'T BREAK!!
						sample_id     = family_prefix+'-'+sample_suffix
						if sample_suffix not in all_samp_names:
							all_samp_names[sample_suffix] = []
							IN_FILES[sample_id] = {}
						IN_FILES[sample_id]['lumpy'] = svd+vcf_file
						all_samp_names[sample_suffix].append('lumpy')
				elif sv_caller == 'wandy':
					wandy_files = [n for n in os.listdir(svd) if '_CNV_segmentation.txt' in n]
					for wandy_file in wandy_files:
						sample_suffix = wandy_file[10:13]	# HARD-CODED FILENAMES. HOPE THIS DOESN'T BREAK!!
						sample_id     = family_prefix+'-'+sample_suffix
						if sample_suffix not in all_samp_names:
							all_samp_names[sample_suffix] = []
							IN_FILES[sample_id] = {}
						IN_FILES[sample_id]['wandy'] = svd+wandy_file
						all_samp_names[sample_suffix].append('wandy')
			#for k in sorted(all_samp_names.keys()):
			#	print k, sorted(all_samp_names[k])
			#for k in sorted(IN_FILES.keys()):
			#	print '---',k
			#	for k2 in sorted(IN_FILES[k].keys()):
			#		print k2, IN_FILES[k][k2]
			#break

		#
		# CONVERT ALL SVS TO SAME INTERNAL FORMAT
		#
		SORTED_SVS_BY_SAMPLE = {}
		for samp_name in sorted(IN_FILES.keys()):

			all_svs_for_sample = []
			for sv_method in sorted(IN_FILES[samp_name].keys()):
				
				print 'parsing sample: ['+samp_name+'] - sv_method: ['+sv_method+']'
				if sv_method == 'cnvnator':
					all_svs_for_sample.extend(parse_cnvnator(IN_FILES[samp_name][sv_method],VAR_FILT))

				elif sv_method == 'delly':
					all_svs_for_sample.extend(parse_vcf(IN_FILES[samp_name][sv_method],VAR_FILT,sv_caller='delly'))

				elif sv_method == 'lumpy':
					all_svs_for_sample.extend(parse_vcf(IN_FILES[samp_name][sv_method],VAR_FILT,sv_caller='lumpy'))

				elif sv_method == 'wandy':
					all_svs_for_sample.extend(parse_wandy(IN_FILES[samp_name][sv_method],VAR_FILT))

			print '--',len(all_svs_for_sample),'total SVs...'

			#
			# FILTER BY MAPPABILITY
			#
			print '* applying mappability filter...'
			all_svs_for_sample_mappability_filtered = []
			for i in xrange(len(all_svs_for_sample)):
				if all_svs_for_sample[i][0] in ['DUP','DEL','INV']:
					ccoord = all_svs_for_sample[i][1]
					lcoord = all_svs_for_sample[i][2]# + all_svs_for_sample[i][4][0]
					ucoord = all_svs_for_sample[i][3]# + all_svs_for_sample[i][5][1]
					lrange = (max([0,lcoord-VAR_FILT['MAP_BUFF']]),lcoord+VAR_FILT['MAP_BUFF'])
					urange = (max([0,ucoord-VAR_FILT['MAP_BUFF']]),ucoord+VAR_FILT['MAP_BUFF'])
					#unmap_count = myMappabilityTrack.query_range_faster(ccoord,lcoord,ucoord,query_endPointInclusive=True)
					#unmap_perct = 100.*unmap_count/float(ucoord-lcoord+1)

					unmap_count_l = myMappabilityTrack.query_range_faster(ccoord,lrange[0],lrange[1],query_endPointInclusive=True)
					unmap_count_u = myMappabilityTrack.query_range_faster(ccoord,urange[0],urange[1],query_endPointInclusive=True)
					unmap_percent_l = unmap_count_l/float(2.*VAR_FILT['MAP_BUFF'])
					unmap_percent_u = unmap_count_u/float(2.*VAR_FILT['MAP_BUFF'])

					if unmap_percent_l <= VAR_FILT['MAP_MAX'] and unmap_percent_u <= VAR_FILT['MAP_MAX']:
						all_svs_for_sample_mappability_filtered.append(all_svs_for_sample[i])
			del all_svs_for_sample
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
					x1 = n[2] + n[4][0]
					x2 = n[3] + n[5][1]
					y1 = sorted_svs[i][2] + sorted_svs[i][4][0]
					y2 = sorted_svs[i][3] + sorted_svs[i][5][1]
					if x1 <= y2 and y1 <= x2:
						inPrevCluster = True
				if inPrevCluster:
					sv_clusters[-1].append(sorted_svs[i])
				else:
					sv_clusters.append([sorted_svs[i]])
			del sorted_svs
			print '--',len(sv_clusters),'SVs...'

			SORTED_SVS_BY_SAMPLE[samp_name] = []
			for i in xrange(len(sv_clusters)):
				lb = min([n[2]+n[4][0] for n in sv_clusters[i]])
				ub = min([n[3]+n[5][1] for n in sv_clusters[i]])
				SORTED_SVS_BY_SAMPLE[samp_name].append([(lb,ub),copy.deepcopy(sv_clusters[i])])
			SORTED_SVS_BY_SAMPLE[samp_name] = sorted(SORTED_SVS_BY_SAMPLE[samp_name])
			del sv_clusters

		print time.time()-tt,'(sec)'

	#
	# OR SKIP ALL THAT AND READ IN FROM PICKLE, IF SPECIFIED.
	#
	else:
		print 'reading input pickle...'
		tt = time.time()
		SORTED_SVS_BY_SAMPLE = pickle.load(open(IN_PICKLE,"rb"))
		print time.time()-tt,'(sec)'

	#
	# SAVE SORTED/FILTERED SVS TO OUTPUT PICKLE, IF DESIRED.
	#
	if OUT_PICKLE != '':
		pickle.dump(SORTED_SVS_BY_SAMPLE,open(OUT_PICKLE,"wb"))

	for k in sorted(SORTED_SVS_BY_SAMPLE.keys()):
		print k, len(SORTED_SVS_BY_SAMPLE[k])
	print '\n********************\n'

	#
	#
	#
	lab = [n[1][0][0] for n in SORTED_SVS_BY_SAMPLE['011H-001']]
	print 'DEL', lab.count('DEL')
	print 'DUP', lab.count('DUP')
	print 'INV', lab.count('INV')

	print ''
	lens = [len(n[1]) for n in SORTED_SVS_BY_SAMPLE['011H-001']]
	for i in xrange(1,5):
		print i, lens.count(i)

	one_count = {}
	for n in SORTED_SVS_BY_SAMPLE['011H-001']:
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
	print '\n*', len(SORTED_SVS_BY_SAMPLE['011H-001']), '-->',
	de_novo = subtractive_filter(SORTED_SVS_BY_SAMPLE['011H-001'],SORTED_SVS_BY_SAMPLE['011H-003'])
	print len(de_novo),'-->',
	# remove mother
	de_novo = subtractive_filter(de_novo,SORTED_SVS_BY_SAMPLE['011H-004'])
	print len(de_novo),'-->',
	# remove sister
	de_novo = subtractive_filter(de_novo,SORTED_SVS_BY_SAMPLE['011H-100'])
	print len(de_novo),'-->',
	# remove brother
	de_novo = subtractive_filter(de_novo,SORTED_SVS_BY_SAMPLE['011H-101'])
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


