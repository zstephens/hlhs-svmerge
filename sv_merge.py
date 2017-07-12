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
from parse import *
from mappability import MappabilityTrack
from sv_comparison import get_denovo, pairwise_all_samples

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
            'GAP_BUFF':150,		# how much to buffer around breakpoint regions when intersecting with gap track
            'GAP_MAX':0.20,		# maximum percentage of sv allowed to overlap gap regions
            'CONC':2}			# how many SV callers must support the same event for it to be considered?

parser = argparse.ArgumentParser(description='sv_merge.py')
parser.add_argument('-i',             type=str,   required=True,  metavar='<str>',                                    help="* base_directory/")
parser.add_argument('-o',             type=str,   required=True,  metavar='<str>',                                    help="* output_dir/")
parser.add_argument('--sv-size-min',  type=int,   required=False, metavar='<int>',     default=VAR_FILT['SIZE_MIN'],  help="minimum allowable SV span [%(default)s]")
parser.add_argument('--sv-size-max',  type=int,   required=False, metavar='<int>',     default=VAR_FILT['SIZE_MAX'],  help="maximum allowable SV span [%(default)s]")
parser.add_argument('--wandy-del',    type=float, required=False, metavar='<float>',   default=VAR_FILT['WANDY_DEL'], help="maximum log_2(depth) for WANDY DEL calls [%(default)s]")
parser.add_argument('--wandy-dup',    type=float, required=False, metavar='<float>',   default=VAR_FILT['WANDY_DUP'], help="minimum log_2(depth) for WANDY DUP calls [%(default)s]")
parser.add_argument('--alt-track',    type=str,   required=False, metavar='<str>',     default='',                    help="alt track bed file")
parser.add_argument('--gap-track',    type=str,   required=False, metavar='<str>',     default='',                    help="gap track bed file")
parser.add_argument('--map-track',    type=str,   required=False, metavar='<str>',     default='',                    help="mappability track bed file")
parser.add_argument('--track-buff',   type=int,   required=False, metavar='<int>',     default=VAR_FILT['MAP_BUFF'],  help="up/downstream buffer for intersecting bed [%(default)s]")
parser.add_argument('--track-thresh', type=float, required=False, metavar='<float>',   default=VAR_FILT['MAP_MAX'],   help="percent threshold for intersecting bed [%(default)s]")
parser.add_argument('--concordance',  type=int,   required=False, metavar='<int>',     default=VAR_FILT['CONC'],      help="minimum concordance (#SV callers, <2 = OFF) [%(default)s]")
parser.add_argument('--vcf-nofilter',             required=False, action='store_true', default=False,                 help='include VCF entries that fail QUAL filter')
parser.add_argument('--no-whitelist',             required=False, action='store_true', default=False,                 help='include SVs from all chr (not just whitelist)')
parser.add_argument('--bcftools',     type=str,   required=False, metavar='<str>',     default='',                    help="/path/to/bcftools")
parser.add_argument('--file-in',      type=str,   required=False, metavar='<str>',     default='',                    help="(skip parsing and read SVs from here instead)")
parser.add_argument('--file-out',     type=str,   required=False, metavar='<str>',     default='',                    help="(save parsed SVs here for speedier future runs)")
args = parser.parse_args()

# TODO: input value sanity-checking code

BASE_INPUT_DIRECTORY  = dir_format(args.i)
RESULTS_DIR           = dir_format(args.o)

VAR_FILT['SIZE_MIN']  = args.sv_size_min
VAR_FILT['SIZE_MAX']  = args.sv_size_max
VAR_FILT['WANDY_DEL'] = args.wandy_del
VAR_FILT['WANDY_DUP'] = args.wandy_dup
VAR_FILT['MAP_BUFF']  = args.track_buff
VAR_FILT['MAP_MAX']   = args.track_thresh
VAR_FILT['CONC']      = args.concordance

if args.vcf_nofilter:
	VAR_FILT['QUAL'] = False

if args.no_whitelist:
	VAR_FILT['WHITELIST'] = False

(IN_FILE, OUT_FILE) = (args.file_in, args.file_out)

BCFTOOLS = args.bcftools

ALT_TRACK = args.alt_track
if ALT_TRACK == '':
	tryDefault = SIM_PATH+'/data/hg38_alt.bed'
	if exists_and_is_nonZero(tryDefault):
		print 'alt track not specified, using default:'
		print '--',tryDefault
	ALT_TRACK = tryDefault

GAP_TRACK = args.gap_track
if GAP_TRACK == '':
	tryDefault = SIM_PATH+'/data/hg38_gap.bed'
	if exists_and_is_nonZero(tryDefault):
		print 'gap track not specified, using default:'
		print '--',tryDefault
	GAP_TRACK = tryDefault

MAPPABILITY_TRACK = args.map_track
if MAPPABILITY_TRACK == '':
	tryDefault = SIM_PATH+'/data/hg19_e2_l400_mappabilityTrack.bed'
	if exists_and_is_nonZero(tryDefault):
		print 'mappability track not specified, using default:'
		print '--',tryDefault
	MAPPABILITY_TRACK = tryDefault


################################################################################
#																			   #
#	MAIN																	   #
#																			   #
################################################################################


def main():

	myAltTrack = MappabilityTrack(ALT_TRACK)
	myGapTrack = MappabilityTrack(GAP_TRACK)
	myMapTrack = MappabilityTrack(MAPPABILITY_TRACK)

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
			all_svs_for_sample = svfilter_unique(all_svs_for_sample)
			print '--',len(all_svs_for_sample),'unique SVs...'

			#
			# FILTER BY CHROMOSOME
			#
			all_svs_for_sample = svfilter_chromWhitelist(all_svs_for_sample,CHR_WHITELIST)
			print '--',len(all_svs_for_sample),'SVs on whitelisted chroms...'

			#
			# FILTER BY MAPPABILITY AND GAP TRACK
			#
			all_svs_for_sample = svfilter_bedIntersect(all_svs_for_sample,myAltTrack,False,"endpoint",mapBuff=VAR_FILT['MAP_BUFF'],mapMax=VAR_FILT['MAP_MAX'])
			print '--',len(all_svs_for_sample),'SVs in non-alt regions...'
			all_svs_for_sample = svfilter_bedIntersect(all_svs_for_sample,myGapTrack,False,"endpoint",mapBuff=VAR_FILT['MAP_BUFF'],mapMax=VAR_FILT['MAP_MAX'])
			print '--',len(all_svs_for_sample),'SVs in non-gap regions...'
			all_svs_for_sample = svfilter_bedIntersect(all_svs_for_sample,myMapTrack,False,"endpoint",mapBuff=VAR_FILT['MAP_BUFF'],mapMax=VAR_FILT['MAP_MAX'])
			print '--',len(all_svs_for_sample),'SVs in mappable regions...'

			#
			# SORT AND CLUSTER EVENTS
			#
			sorted_svs  = sorted(all_svs_for_sample)
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
			FILTERED_SVS[samp_name] = []
			for i in xrange(len(sv_clusters)):
				lb = min([n[2]+n[4][0] for n in sv_clusters[i]])
				ub = min([n[3]+n[5][1] for n in sv_clusters[i]])
				FILTERED_SVS[samp_name].append([(sv_clusters[i][0][1],lb,ub),copy.deepcopy(sv_clusters[i])])
			FILTERED_SVS[samp_name] = sorted(FILTERED_SVS[samp_name])
			del sv_clusters
			print '--',len(FILTERED_SVS[samp_name]),'SV clusters...'

			#
			# APPLY CONCORDANCE FILTER
			#
			FILTERED_SVS[samp_name] = svfilter_concordance(FILTERED_SVS[samp_name],VAR_FILT['CONC'])
			print '--',len(FILTERED_SVS[samp_name]),'SV clusters that pass concordance filter...'


	#
	# OR SKIP ALL THAT AND READ IN FROM FILE, IF SPECIFIED.
	#
	else:
		print 'reading input file...'
		tt = time.time()
		FILTERED_SVS = read_input_file(IN_FILE)
		print time.time()-tt,'(sec)'

	#
	# SAVE SORTED/FILTERED SVS TO OUTPUT FILE, IF DESIRED.
	#
	if OUT_FILE != '':
		print 'saving to output file...'
		tt = time.time()
		write_output_file(FILTERED_SVS,VAR_FILT,OUT_FILE)
		print time.time()-tt,'(sec)'

	#
	# OUTPUT DENOVO SVS
	#
	print 'finding denovo SVs...'
	for k in sorted(FILTERED_SVS.keys()):
		suffix = k.split('-')[1]
		rdir   = RESULTS_DIR+k+'_denovo/'
		if suffix == '001':
			print '--',k
			makedir(rdir)
			myDenovo = {}
			myDenovo[k] = get_denovo(FILTERED_SVS,k)
			if len(myDenovo[k]):
				rout = rdir+k+'_denovo.txt'
				write_output_file(myDenovo,VAR_FILT,rout)
	
	#pairwise_all_samples(FILTERED_SVS)


if __name__ == '__main__':

	#myTrack = MappabilityTrack('/Users/zach/Desktop/hlhs_summer_2017/test.bed')
	#print (80,100), 1, myTrack.query_range('chr1',80,100)
	#print (80,120), 20, myTrack.query_range('chr1',80,120)
	#print (120,140), 11, myTrack.query_range('chr1',120,140)
	#print (130,140), 1, myTrack.query_range('chr1',130,140)
	#print (120,220), 30, myTrack.query_range('chr1',120,220)
	#exit(1)

	main()


