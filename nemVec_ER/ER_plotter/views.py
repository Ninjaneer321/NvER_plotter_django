from django.http import HttpResponse
from django.shortcuts import render
from django.forms import formset_factory
from django.core.paginator import Paginator, EmptyPage, PageNotAnInteger
from ER_plotter.models import Fasta, Regen_cpm, Embryo_cpm, Annotation, Regen_SE, Mfuzz, Regen_log_SE, Embryo_SE, de_table
from .forms import Gene_searchForm, NvERTxForm, ConvertForm, comparisonForm
import math, re, diggPaginator
import os
import sys
import tempfile
from blastplus import utils
from blastplus.forms import BlastForm, TBlastnForm, BlastpForm, BlastxForm
from blastplus.settings import BLAST_CORRECT_PARAMS
from blastplus.settings import EVALUE_BLAST_DEFAULT, BLAST_MAX_NUMBER_SEQ_IN_INPUT
from blastplus.settings import EXAMPLE_FASTA_NUCL_FILE_PATH, EXAMPLE_FASTA_PROT_FILE_PATH
from blastplus.settings import BLAST_DB_NUCL_LIST

from django.db.models import Transform

class AbsoluteValue(Transform):
    lookup_name = 'abs'
    function = 'ABS'
    
from django.db.models import IntegerField
IntegerField.register_lookup(AbsoluteValue)

def results(request):
	regen_point = [0,2,4,8,12,16,20,24,36,48,60,72,96,120,144]
	gene_search_form = Gene_searchForm(request.POST or None)
	nvertx_form = NvERTxForm(request.POST or None)
	convert_form = ConvertForm(request.POST or None)
	
	#get id
	nvid = request.GET.getlist('Nvid', '')

	nvertx_1 = False
	nvertx_2 = False
	nvertx_3 = False
	nvertx_4 = False
	nvertx_5 = False
	log2 = False

	if nvid :
		nvertx_1 = nvid[0]
		if len(nvid) >= 2 :
			nvertx_2 = nvid[1]
		if len(nvid) >= 3 :
			nvertx_3 = nvid[2]
		if len(nvid) >= 4 :
			nvertx_4 = nvid[3]
		if len(nvid) == 5 :
			nvertx_5 = nvid[4]
		nvertx_search = True


	if nvertx_form.is_valid():
		nvertx_1 = nvertx_form.cleaned_data['nvertx_1']
		if nvertx_1[0] != 'N' :
			nvertx_1 = 'NvERTx.4.' + nvertx_1
		if not re.match('NvERTx.4.[0-9]+$', nvertx_1):
			nvertx_1 = 'NA'
		nvertx_2 = nvertx_form.cleaned_data['nvertx_2']
		if nvertx_2 and nvertx_2[0] != 'N' :
			nvertx_2 = 'NvERTx.4.' + nvertx_2
		if nvertx_2 and not re.match('NvERTx.4.[0-9]+$', nvertx_2):
			nvertx_2 = 'NA'
		nvertx_3 = nvertx_form.cleaned_data['nvertx_3']
		if nvertx_3 and nvertx_3[0] != 'N' :
			nvertx_3 = 'NvERTx.4.' + nvertx_3
		if nvertx_3 and not re.match('NvERTx.4.[0-9]+$', nvertx_3):
			nvertx_3 = 'NA'
		nvertx_4 = nvertx_form.cleaned_data['nvertx_4']
		if nvertx_4 and nvertx_4[0] != 'N' :
			nvertx_4 = 'NvERTx.4.' + nvertx_4
		if nvertx_4 and not re.match('NvERTx.4.[0-9]+$', nvertx_4):
			nvertx_4 = 'NA'
		nvertx_5 = nvertx_form.cleaned_data['nvertx_5']
		if nvertx_5 and nvertx_5[0] != 'N' :
			nvertx_5 = 'NvERTx.4.' + nvertx_5
		if nvertx_5 and not re.match('NvERTx.4.[0-9]+$', nvertx_5):
			nvertx_5 = 'NA'
		log2 = nvertx_form.cleaned_data['log2']
		nvertx_search = True
		
	if nvertx_1 :
		nvertx_1_embryo_warner_invalid = False
		nvertx_1_embryo_fischer_invalid = False
		nvertx_1_embryo_helm_invalid = False
		try :
			sequence_fasta_1 = Fasta.objects.get(nvertx_id=nvertx_1).fasta_sequence
		except :
			sequence_fasta_1 = 'No fasta sequence'
		try :
			if not log2 :
				regen_UC_1 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_1).regen_anc_UC+1,2),2)
				regen_0_1 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_1).regen_anc_0HPA+1,2),2)
				regen_2_1 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_1).regen_anc_2HPA+1,2),2)
				regen_4_1 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_1).regen_anc_4HPA+1,2),2)
				regen_8_1 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_1).regen_anc_8HPA+1,2),2)
				regen_12_1 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_1).regen_anc_12HPA+1,2),2)
				regen_16_1 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_1).regen_anc_16HPA+1,2),2)
				regen_20_1 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_1).regen_anc_20HPA+1,2),2)
				regen_24_1 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_1).regen_anc_24HPA+1,2),2)
				regen_36_1 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_1).regen_anc_36HPA+1,2),2)
				regen_48_1 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_1).regen_anc_48HPA+1,2),2)
				regen_60_1 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_1).regen_anc_60HPA+1,2),2)
				regen_72_1 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_1).regen_anc_72HPA+1,2),2)
				regen_96_1 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_1).regen_anc_96HPA+1,2),2)
				regen_120_1 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_1).regen_anc_120HPA+1,2),2)
				regen_144_1 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_1).regen_anc_144HPA+1,2),2)
				regen_se_UC_1 = round(Regen_log_SE.objects.get(nvertx_id=nvertx_1).regen_log_se_UC,2)
				regen_se_0_1 = round(Regen_log_SE.objects.get(nvertx_id=nvertx_1).regen_log_se_0HPA,2)
				regen_se_2_1 = round(Regen_log_SE.objects.get(nvertx_id=nvertx_1).regen_log_se_2HPA,2)
				regen_se_4_1 = round(Regen_log_SE.objects.get(nvertx_id=nvertx_1).regen_log_se_4HPA,2)
				regen_se_8_1 = round(Regen_log_SE.objects.get(nvertx_id=nvertx_1).regen_log_se_8HPA,2)
				regen_se_12_1 = round(Regen_log_SE.objects.get(nvertx_id=nvertx_1).regen_log_se_12HPA,2)
				regen_se_16_1 = round(Regen_log_SE.objects.get(nvertx_id=nvertx_1).regen_log_se_16HPA,2)
				regen_se_20_1 = round(Regen_log_SE.objects.get(nvertx_id=nvertx_1).regen_log_se_20HPA,2)
				regen_se_24_1 = round(Regen_log_SE.objects.get(nvertx_id=nvertx_1).regen_log_se_24HPA,2)
				regen_se_36_1 = round(Regen_log_SE.objects.get(nvertx_id=nvertx_1).regen_log_se_36HPA,2)
				regen_se_48_1 = round(Regen_log_SE.objects.get(nvertx_id=nvertx_1).regen_log_se_48HPA,2)
				regen_se_60_1 = round(Regen_log_SE.objects.get(nvertx_id=nvertx_1).regen_log_se_60HPA,2)
				regen_se_72_1 = round(Regen_log_SE.objects.get(nvertx_id=nvertx_1).regen_log_se_72HPA,2)
				regen_se_96_1 = round(Regen_log_SE.objects.get(nvertx_id=nvertx_1).regen_log_se_96HPA,2)
				regen_se_120_1 = round(Regen_log_SE.objects.get(nvertx_id=nvertx_1).regen_log_se_120HPA,2)
				regen_se_144_1 = round(Regen_log_SE.objects.get(nvertx_id=nvertx_1).regen_log_se_144HPA,2)
			else :
				regen_UC_1 = round(Regen_cpm.objects.get(nvertx_id=nvertx_1).regen_anc_UC,2)
				regen_0_1 = round(Regen_cpm.objects.get(nvertx_id=nvertx_1).regen_anc_0HPA,2)
				regen_2_1 = round(Regen_cpm.objects.get(nvertx_id=nvertx_1).regen_anc_2HPA,2)
				regen_4_1 = round(Regen_cpm.objects.get(nvertx_id=nvertx_1).regen_anc_4HPA,2)
				regen_8_1 = round(Regen_cpm.objects.get(nvertx_id=nvertx_1).regen_anc_8HPA,2)
				regen_12_1 = round(Regen_cpm.objects.get(nvertx_id=nvertx_1).regen_anc_12HPA,2)
				regen_16_1 = round(Regen_cpm.objects.get(nvertx_id=nvertx_1).regen_anc_16HPA,2)
				regen_20_1 = round(Regen_cpm.objects.get(nvertx_id=nvertx_1).regen_anc_20HPA,2)
				regen_24_1 = round(Regen_cpm.objects.get(nvertx_id=nvertx_1).regen_anc_24HPA,2)
				regen_36_1 = round(Regen_cpm.objects.get(nvertx_id=nvertx_1).regen_anc_36HPA,2)
				regen_48_1 = round(Regen_cpm.objects.get(nvertx_id=nvertx_1).regen_anc_48HPA,2)
				regen_60_1 = round(Regen_cpm.objects.get(nvertx_id=nvertx_1).regen_anc_60HPA,2)
				regen_72_1 = round(Regen_cpm.objects.get(nvertx_id=nvertx_1).regen_anc_72HPA,2)
				regen_96_1 = round(Regen_cpm.objects.get(nvertx_id=nvertx_1).regen_anc_96HPA,2)
				regen_120_1 = round(Regen_cpm.objects.get(nvertx_id=nvertx_1).regen_anc_120HPA,2)
				regen_144_1 = round(Regen_cpm.objects.get(nvertx_id=nvertx_1).regen_anc_144HPA,2)
				regen_se_UC_1 = round(Regen_SE.objects.get(nvertx_id=nvertx_1).regen_se_UC,2)
				regen_se_0_1 = round(Regen_SE.objects.get(nvertx_id=nvertx_1).regen_se_0HPA,2)
				regen_se_2_1 = round(Regen_SE.objects.get(nvertx_id=nvertx_1).regen_se_2HPA,2)
				regen_se_4_1 = round(Regen_SE.objects.get(nvertx_id=nvertx_1).regen_se_4HPA,2)
				regen_se_8_1 = round(Regen_SE.objects.get(nvertx_id=nvertx_1).regen_se_8HPA,2)
				regen_se_12_1 = round(Regen_SE.objects.get(nvertx_id=nvertx_1).regen_se_12HPA,2)
				regen_se_16_1 = round(Regen_SE.objects.get(nvertx_id=nvertx_1).regen_se_16HPA,2)
				regen_se_20_1 = round(Regen_SE.objects.get(nvertx_id=nvertx_1).regen_se_20HPA,2)
				regen_se_24_1 = round(Regen_SE.objects.get(nvertx_id=nvertx_1).regen_se_24HPA,2)
				regen_se_36_1 = round(Regen_SE.objects.get(nvertx_id=nvertx_1).regen_se_36HPA,2)
				regen_se_48_1 = round(Regen_SE.objects.get(nvertx_id=nvertx_1).regen_se_48HPA,2)
				regen_se_60_1 = round(Regen_SE.objects.get(nvertx_id=nvertx_1).regen_se_60HPA,2)
				regen_se_72_1 = round(Regen_SE.objects.get(nvertx_id=nvertx_1).regen_se_72HPA,2)
				regen_se_96_1 = round(Regen_SE.objects.get(nvertx_id=nvertx_1).regen_se_96HPA,2)
				regen_se_120_1 = round(Regen_SE.objects.get(nvertx_id=nvertx_1).regen_se_120HPA,2)
				regen_se_144_1 = round(Regen_SE.objects.get(nvertx_id=nvertx_1).regen_se_144HPA,2)
		except :
			nvertx_1_regen_invalid = True
		try :
			embryo_warner_24_1 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_1).warner_anc_24HPF,2)
			embryo_warner_48_1 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_1).warner_anc_48HPF,2)
			embryo_warner_72_1 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_1).warner_anc_72HPF,2)
			embryo_warner_96_1 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_1).warner_anc_96HPF,2)
			embryo_warner_120_1 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_1).warner_anc_120HPF,2)
			embryo_warner_144_1 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_1).warner_anc_144HPF,2)
			embryo_warner_168_1 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_1).warner_anc_168HPF,2)
			embryo_warner_192_1 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_1).warner_anc_192HPF,2)
			embryo_warner_se_24_1 = round(Embryo_SE.objects.get(nvertx_id=nvertx_1).warner_se_24HPF,2)
			embryo_warner_se_48_1 = round(Embryo_SE.objects.get(nvertx_id=nvertx_1).warner_se_48HPF,2)
			embryo_warner_se_72_1 = round(Embryo_SE.objects.get(nvertx_id=nvertx_1).warner_se_72HPF,2)
			embryo_warner_se_96_1 = round(Embryo_SE.objects.get(nvertx_id=nvertx_1).warner_se_96HPF,2)
			embryo_warner_se_120_1 = round(Embryo_SE.objects.get(nvertx_id=nvertx_1).warner_se_120HPF,2)
			embryo_warner_se_144_1 = round(Embryo_SE.objects.get(nvertx_id=nvertx_1).warner_se_144HPF,2)
			embryo_warner_se_168_1 = round(Embryo_SE.objects.get(nvertx_id=nvertx_1).warner_se_168HPF,2)
			embryo_warner_se_192_1 = round(Embryo_SE.objects.get(nvertx_id=nvertx_1).warner_se_192HPF,2)
		except :
			nvertx_1_embryo_warner_invalid = True
		try :
			embryo_fischer_0_1 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_1).fischer_anc_0HPF,2)
			embryo_fischer_1_1 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_1).fischer_anc_1HPF,2)
			embryo_fischer_2_1 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_1).fischer_anc_2HPF,2)
			embryo_fischer_3_1 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_1).fischer_anc_3HPF,2)
			embryo_fischer_4_1 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_1).fischer_anc_4HPF,2)
			embryo_fischer_5_1 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_1).fischer_anc_5HPF,2)
			embryo_fischer_6_1 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_1).fischer_anc_6HPF,2)
			embryo_fischer_7_1 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_1).fischer_anc_7HPF,2)
			embryo_fischer_8_1 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_1).fischer_anc_8HPF,2)
			embryo_fischer_9_1 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_1).fischer_anc_9HPF,2)	
			embryo_fischer_10_1 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_1).fischer_anc_10HPF,2)
			embryo_fischer_11_1 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_1).fischer_anc_11HPF,2)
			embryo_fischer_12_1 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_1).fischer_anc_12HPF,2)
			embryo_fischer_13_1 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_1).fischer_anc_13HPF,2)
			embryo_fischer_14_1 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_1).fischer_anc_14HPF,2)
			embryo_fischer_15_1 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_1).fischer_anc_15HPF,2)
			embryo_fischer_16_1 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_1).fischer_anc_16HPF,2)		
			embryo_fischer_17_1 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_1).fischer_anc_17HPF,2)
			embryo_fischer_18_1 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_1).fischer_anc_18HPF,2)
			embryo_fischer_19_1 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_1).fischer_anc_19HPF,2)
			embryo_fischer_se_0_1 = round(Embryo_SE.objects.get(nvertx_id=nvertx_1).fischer_se_0HPF,2)
			embryo_fischer_se_1_1 = round(Embryo_SE.objects.get(nvertx_id=nvertx_1).fischer_se_1HPF,2)
			embryo_fischer_se_2_1 = round(Embryo_SE.objects.get(nvertx_id=nvertx_1).fischer_se_2HPF,2)
			embryo_fischer_se_3_1 = round(Embryo_SE.objects.get(nvertx_id=nvertx_1).fischer_se_3HPF,2)
			embryo_fischer_se_4_1 = round(Embryo_SE.objects.get(nvertx_id=nvertx_1).fischer_se_4HPF,2)
			embryo_fischer_se_5_1 = round(Embryo_SE.objects.get(nvertx_id=nvertx_1).fischer_se_5HPF,2)
			embryo_fischer_se_6_1 = round(Embryo_SE.objects.get(nvertx_id=nvertx_1).fischer_se_6HPF,2)
			embryo_fischer_se_7_1 = round(Embryo_SE.objects.get(nvertx_id=nvertx_1).fischer_se_7HPF,2)
			embryo_fischer_se_8_1 = round(Embryo_SE.objects.get(nvertx_id=nvertx_1).fischer_se_8HPF,2)
			embryo_fischer_se_9_1 = round(Embryo_SE.objects.get(nvertx_id=nvertx_1).fischer_se_9HPF,2)
			#embryo_fischer_se_10_1 = [round(embryo_fischer_10_1 - Embryo_SE.objects.get(nvertx_id=nvertx_1).fischer_se_10HPF,2),round(embryo_fischer_10_1 + Embryo_SE.objects.get(nvertx_id=nvertx_1).fischer_se_10HPF,2)]
			embryo_fischer_se_11_1 = round(Embryo_SE.objects.get(nvertx_id=nvertx_1).fischer_se_11HPF,2)
			embryo_fischer_se_12_1 = round(Embryo_SE.objects.get(nvertx_id=nvertx_1).fischer_se_12HPF,2)
			embryo_fischer_se_13_1 = round(Embryo_SE.objects.get(nvertx_id=nvertx_1).fischer_se_13HPF,2)
			embryo_fischer_se_14_1 = round(Embryo_SE.objects.get(nvertx_id=nvertx_1).fischer_se_14HPF,2)
			embryo_fischer_se_15_1 = round(Embryo_SE.objects.get(nvertx_id=nvertx_1).fischer_se_15HPF,2)
			embryo_fischer_se_16_1 = round(Embryo_SE.objects.get(nvertx_id=nvertx_1).fischer_se_16HPF,2)
			embryo_fischer_se_17_1 = round(Embryo_SE.objects.get(nvertx_id=nvertx_1).fischer_se_17HPF,2)
			embryo_fischer_se_18_1 = round(Embryo_SE.objects.get(nvertx_id=nvertx_1).fischer_se_18HPF,2)
			embryo_fischer_se_19_1 = round(Embryo_SE.objects.get(nvertx_id=nvertx_1).fischer_se_19HPF,2)
		except :
			nvertx_1_embryo_fischer_invalid = True
		try :
			embryo_helm_2_1 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_1).helm_anc_2HPF,2)
			embryo_helm_7_1 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_1).helm_anc_7HPF,2)
			embryo_helm_12_1 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_1).helm_anc_12HPF,2)
			embryo_helm_24_1 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_1).helm_anc_24HPF,2)
			embryo_helm_120_1 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_1).helm_anc_120HPF,2)
			embryo_helm_240_1 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_1).helm_anc_240HPF,2)
			embryo_helm_se_2_1 = round(Embryo_SE.objects.get(nvertx_id=nvertx_1).helm_se_2HPF,2)
			embryo_helm_se_7_1 = round(Embryo_SE.objects.get(nvertx_id=nvertx_1).helm_se_7HPF,2)
			embryo_helm_se_12_1 = round(Embryo_SE.objects.get(nvertx_id=nvertx_1).helm_se_12HPF,2)
			embryo_helm_se_24_1 = round(Embryo_SE.objects.get(nvertx_id=nvertx_1).helm_se_24HPF,2)
			embryo_helm_se_120_1 = round(Embryo_SE.objects.get(nvertx_id=nvertx_1).helm_se_120HPF,2)
			embryo_helm_se_240_1 = round(Embryo_SE.objects.get(nvertx_id=nvertx_1).helm_se_240HPF,2)
		except :
			nvertx_1_embryo_helm_invalid = True
		if not nvertx_1_embryo_warner_invalid and not nvertx_1_embryo_fischer_invalid and not nvertx_1_embryo_helm_invalid :
			try :
				embryo_mean_2_1 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_1).mean_2HPF,2)
				embryo_mean_7_1 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_1).mean_7HPF,2)
				embryo_mean_12_1 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_1).mean_12HPF,2)
				embryo_mean_24_1 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_1).mean_24HPF,2)
				embryo_mean_120_1 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_1).mean_120HPF,2)
				embryo_mean_0_1 = embryo_fischer_0_1
				embryo_mean_1_1 = embryo_fischer_1_1
				embryo_mean_3_1 = embryo_fischer_3_1
				embryo_mean_4_1 = embryo_fischer_4_1
				embryo_mean_5_1 = embryo_fischer_5_1
				embryo_mean_6_1 = embryo_fischer_6_1
				embryo_mean_8_1 = embryo_fischer_8_1
				embryo_mean_9_1 = embryo_fischer_9_1
				embryo_mean_10_1 = embryo_fischer_10_1
				embryo_mean_11_1 = embryo_fischer_11_1
				embryo_mean_13_1 = embryo_fischer_13_1
				embryo_mean_14_1 = embryo_fischer_14_1
				embryo_mean_15_1 = embryo_fischer_15_1
				embryo_mean_16_1 = embryo_fischer_16_1
				embryo_mean_17_1 = embryo_fischer_17_1
				embryo_mean_18_1 = embryo_fischer_18_1
				embryo_mean_19_1 = embryo_fischer_19_1
				embryo_mean_48_1 = embryo_warner_48_1
				embryo_mean_72_1 = embryo_warner_72_1
				embryo_mean_96_1 = embryo_warner_96_1
				embryo_mean_144_1 = embryo_warner_144_1
				embryo_mean_168_1 = embryo_warner_168_1
				embryo_mean_192_1 = embryo_warner_192_1
				embryo_mean_240_1 = embryo_helm_240_1
			except :
				nvertx_1_embryo_mean_invalid = True
		else :
			nvertx_1_embryo_mean_invalid = True
		try :
			annot_nemve1_tophit_1 = Annotation.objects.get(nvertx_id=nvertx_1).Nemve1_tophit
			annot_nemve1_e_val_1 = round(Annotation.objects.get(nvertx_id=nvertx_1).Nemve1_e_val,150)
			if annot_nemve1_e_val_1 == 99 :
				annot_nemve1_e_val_1 = "NA"
			annot_mfuzz_r_clust_1 = Annotation.objects.get(nvertx_id=nvertx_1).Mfuzz_R_Clust
			annot_mfuzz_r_score_1 = round(Annotation.objects.get(nvertx_id=nvertx_1).Mfuzz_R_Score,2)
			if annot_mfuzz_r_score_1 == -1 :
				annot_mfuzz_r_score_1 = "NA"
			annot_mfuzz_e_clust_1 = Annotation.objects.get(nvertx_id=nvertx_1).Mfuzz_E_Clust
			annot_mfuzz_e_score_1 = round(Annotation.objects.get(nvertx_id=nvertx_1).Mfuzz_E_Score,2)
			if annot_mfuzz_e_score_1 == -1 :
				annot_mfuzz_e_score_1 = "NA"
			annot_uniprot_id_1 = Annotation.objects.get(nvertx_id=nvertx_1).Uniprot_ID
			annot_uniprot_description_1 = Annotation.objects.get(nvertx_id=nvertx_1).Uniprot_Description
			annot_top_nr_hit_eval_1 = Annotation.objects.get(nvertx_id=nvertx_1).Top_nr_hit_eval
			if annot_top_nr_hit_eval_1 != "NA" :
				annot_top_nr_hit_eval_1_split = annot_top_nr_hit_eval_1.split('|',4)
				annot_nr_beg_1 = annot_top_nr_hit_eval_1_split[0] + "|" + annot_top_nr_hit_eval_1_split[1] + "|" + annot_top_nr_hit_eval_1_split[2]
				annot_nr_link_1 = annot_top_nr_hit_eval_1_split[3]
				annot_nr_end_1 = annot_top_nr_hit_eval_1_split[4]
				annot_other_nr_hits_1 = Annotation.objects.get(nvertx_id=nvertx_1).Other_nr_hits
				nr_hit_graph_1 = re.search('[\[\- \w]+\]', annot_top_nr_hit_eval_1).group(0)
		except :
			nvertx_1_annot_invalid = True
		try :
			ncbi_1 = annot_top_nr_hit_eval_1.split('|')[1]
			prot_1 = annot_top_nr_hit_eval_1.split('|')[3]
		except :
			nvertx_1_links_invalid = True

		if nvertx_2 :
			nvertx_2_embryo_warner_invalid = False
			nvertx_2_embryo_fischer_invalid = False
			nvertx_2_embryo_helm_invalid = False
			try :
				sequence_fasta_2 = Fasta.objects.get(nvertx_id=nvertx_2).fasta_sequence
			except :
				sequence_fasta_2 = 'No fasta sequence'
			try :
				if not log2 :
					regen_UC_2 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_2).regen_anc_UC+1,2),2)
					regen_0_2 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_2).regen_anc_0HPA+1,2),2)
					regen_2_2 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_2).regen_anc_2HPA+1,2),2)
					regen_4_2 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_2).regen_anc_4HPA+1,2),2)
					regen_8_2 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_2).regen_anc_8HPA+1,2),2)
					regen_12_2 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_2).regen_anc_12HPA+1,2),2)
					regen_16_2 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_2).regen_anc_16HPA+1,2),2)
					regen_20_2 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_2).regen_anc_20HPA+1,2),2)
					regen_24_2 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_2).regen_anc_24HPA+1,2),2)
					regen_36_2 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_2).regen_anc_36HPA+1,2),2)
					regen_48_2 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_2).regen_anc_48HPA+1,2),2)
					regen_60_2 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_2).regen_anc_60HPA+1,2),2)
					regen_72_2 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_2).regen_anc_72HPA+1,2),2)
					regen_96_2 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_2).regen_anc_96HPA+1,2),2)
					regen_120_2 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_2).regen_anc_120HPA+1,2),2)
					regen_144_2 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_2).regen_anc_144HPA+1,2),2)
					regen_se_UC_2 = round(Regen_log_SE.objects.get(nvertx_id=nvertx_2).regen_log_se_UC,2)
					regen_se_0_2 = round(Regen_log_SE.objects.get(nvertx_id=nvertx_2).regen_log_se_0HPA,2)
					regen_se_2_2 = round(Regen_log_SE.objects.get(nvertx_id=nvertx_2).regen_log_se_2HPA,2)
					regen_se_4_2 = round(Regen_log_SE.objects.get(nvertx_id=nvertx_2).regen_log_se_4HPA,2)
					regen_se_8_2 = round(Regen_log_SE.objects.get(nvertx_id=nvertx_2).regen_log_se_8HPA,2)
					regen_se_12_2 = round(Regen_log_SE.objects.get(nvertx_id=nvertx_2).regen_log_se_12HPA,2)
					regen_se_16_2 = round(Regen_log_SE.objects.get(nvertx_id=nvertx_2).regen_log_se_16HPA,2)
					regen_se_20_2 = round(Regen_log_SE.objects.get(nvertx_id=nvertx_2).regen_log_se_20HPA,2)
					regen_se_24_2 = round(Regen_log_SE.objects.get(nvertx_id=nvertx_2).regen_log_se_24HPA,2)
					regen_se_36_2 = round(Regen_log_SE.objects.get(nvertx_id=nvertx_2).regen_log_se_36HPA,2)
					regen_se_48_2 = round(Regen_log_SE.objects.get(nvertx_id=nvertx_2).regen_log_se_48HPA,2)
					regen_se_60_2 = round(Regen_log_SE.objects.get(nvertx_id=nvertx_2).regen_log_se_60HPA,2)
					regen_se_72_2 = round(Regen_log_SE.objects.get(nvertx_id=nvertx_2).regen_log_se_72HPA,2)
					regen_se_96_2 = round(Regen_log_SE.objects.get(nvertx_id=nvertx_2).regen_log_se_96HPA,2)
					regen_se_120_2 = round(Regen_log_SE.objects.get(nvertx_id=nvertx_2).regen_log_se_120HPA,2)
					regen_se_144_2 = round(Regen_log_SE.objects.get(nvertx_id=nvertx_2).regen_log_se_144HPA,2)
				else :
					regen_UC_2 = round(Regen_cpm.objects.get(nvertx_id=nvertx_2).regen_anc_UC,2)
					regen_0_2 = round(Regen_cpm.objects.get(nvertx_id=nvertx_2).regen_anc_0HPA,2)
					regen_2_2 = round(Regen_cpm.objects.get(nvertx_id=nvertx_2).regen_anc_2HPA,2)
					regen_4_2 = round(Regen_cpm.objects.get(nvertx_id=nvertx_2).regen_anc_4HPA,2)
					regen_8_2 = round(Regen_cpm.objects.get(nvertx_id=nvertx_2).regen_anc_8HPA,2)
					regen_12_2 = round(Regen_cpm.objects.get(nvertx_id=nvertx_2).regen_anc_12HPA,2)
					regen_16_2 = round(Regen_cpm.objects.get(nvertx_id=nvertx_2).regen_anc_16HPA,2)
					regen_20_2 = round(Regen_cpm.objects.get(nvertx_id=nvertx_2).regen_anc_20HPA,2)
					regen_24_2 = round(Regen_cpm.objects.get(nvertx_id=nvertx_2).regen_anc_24HPA,2)
					regen_36_2 = round(Regen_cpm.objects.get(nvertx_id=nvertx_2).regen_anc_36HPA,2)
					regen_48_2 = round(Regen_cpm.objects.get(nvertx_id=nvertx_2).regen_anc_48HPA,2)
					regen_60_2 = round(Regen_cpm.objects.get(nvertx_id=nvertx_2).regen_anc_60HPA,2)
					regen_72_2 = round(Regen_cpm.objects.get(nvertx_id=nvertx_2).regen_anc_72HPA,2)
					regen_96_2 = round(Regen_cpm.objects.get(nvertx_id=nvertx_2).regen_anc_96HPA,2)
					regen_120_2 = round(Regen_cpm.objects.get(nvertx_id=nvertx_2).regen_anc_120HPA,2)
					regen_144_2 = round(Regen_cpm.objects.get(nvertx_id=nvertx_2).regen_anc_144HPA,2)
					regen_se_UC_2 = round(Regen_SE.objects.get(nvertx_id=nvertx_2).regen_se_UC,2)
					regen_se_0_2 = round(Regen_SE.objects.get(nvertx_id=nvertx_2).regen_se_0HPA,2)
					regen_se_2_2 = round(Regen_SE.objects.get(nvertx_id=nvertx_2).regen_se_2HPA,2)
					regen_se_4_2 = round(Regen_SE.objects.get(nvertx_id=nvertx_2).regen_se_4HPA,2)
					regen_se_8_2 = round(Regen_SE.objects.get(nvertx_id=nvertx_2).regen_se_8HPA,2)
					regen_se_12_2 = round(Regen_SE.objects.get(nvertx_id=nvertx_2).regen_se_12HPA,2)
					regen_se_16_2 = round(Regen_SE.objects.get(nvertx_id=nvertx_2).regen_se_16HPA,2)
					regen_se_20_2 = round(Regen_SE.objects.get(nvertx_id=nvertx_2).regen_se_20HPA,2)
					regen_se_24_2 = round(Regen_SE.objects.get(nvertx_id=nvertx_2).regen_se_24HPA,2)
					regen_se_36_2 = round(Regen_SE.objects.get(nvertx_id=nvertx_2).regen_se_36HPA,2)
					regen_se_48_2 = round(Regen_SE.objects.get(nvertx_id=nvertx_2).regen_se_48HPA,2)
					regen_se_60_2 = round(Regen_SE.objects.get(nvertx_id=nvertx_2).regen_se_60HPA,2)
					regen_se_72_2 = round(Regen_SE.objects.get(nvertx_id=nvertx_2).regen_se_72HPA,2)
					regen_se_96_2 = round(Regen_SE.objects.get(nvertx_id=nvertx_2).regen_se_96HPA,2)
					regen_se_120_2 = round(Regen_SE.objects.get(nvertx_id=nvertx_2).regen_se_120HPA,2)
					regen_se_144_2 = round(Regen_SE.objects.get(nvertx_id=nvertx_2).regen_se_144HPA,2)
			except :
				nvertx_2_regen_invalid = True
			try :
				embryo_warner_24_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).warner_anc_24HPF,2)
				embryo_warner_48_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).warner_anc_48HPF,2)
				embryo_warner_72_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).warner_anc_72HPF,2)
				embryo_warner_96_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).warner_anc_96HPF,2)
				embryo_warner_120_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).warner_anc_120HPF,2)
				embryo_warner_144_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).warner_anc_144HPF,2)
				embryo_warner_168_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).warner_anc_168HPF,2)
				embryo_warner_192_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).warner_anc_192HPF,2)
				embryo_warner_se_24_2 = round(Embryo_SE.objects.get(nvertx_id=nvertx_2).warner_se_24HPF,2)
				embryo_warner_se_48_2 = round(Embryo_SE.objects.get(nvertx_id=nvertx_2).warner_se_48HPF,2)
				embryo_warner_se_72_2 = round(Embryo_SE.objects.get(nvertx_id=nvertx_2).warner_se_72HPF,2)
				embryo_warner_se_96_2 = round(Embryo_SE.objects.get(nvertx_id=nvertx_2).warner_se_96HPF,2)
				embryo_warner_se_120_2 = round(Embryo_SE.objects.get(nvertx_id=nvertx_2).warner_se_120HPF,2)
				embryo_warner_se_144_2 = round(Embryo_SE.objects.get(nvertx_id=nvertx_2).warner_se_144HPF,2)
				embryo_warner_se_168_2 = round(Embryo_SE.objects.get(nvertx_id=nvertx_2).warner_se_168HPF,2)
				embryo_warner_se_192_2 = round(Embryo_SE.objects.get(nvertx_id=nvertx_2).warner_se_192HPF,2)
			except :
				nvertx_2_embryo_warner_invalid = True
			try :
				embryo_fischer_0_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).fischer_anc_0HPF,2)
				embryo_fischer_1_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).fischer_anc_1HPF,2)
				embryo_fischer_2_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).fischer_anc_2HPF,2)
				embryo_fischer_3_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).fischer_anc_3HPF,2)
				embryo_fischer_4_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).fischer_anc_4HPF,2)
				embryo_fischer_5_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).fischer_anc_5HPF,2)
				embryo_fischer_6_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).fischer_anc_6HPF,2)
				embryo_fischer_7_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).fischer_anc_7HPF,2)
				embryo_fischer_8_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).fischer_anc_8HPF,2)
				embryo_fischer_9_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).fischer_anc_9HPF,2)	
				embryo_fischer_10_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).fischer_anc_10HPF,2)
				embryo_fischer_11_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).fischer_anc_11HPF,2)
				embryo_fischer_12_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).fischer_anc_12HPF,2)
				embryo_fischer_13_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).fischer_anc_13HPF,2)
				embryo_fischer_14_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).fischer_anc_14HPF,2)
				embryo_fischer_15_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).fischer_anc_15HPF,2)
				embryo_fischer_16_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).fischer_anc_16HPF,2)		
				embryo_fischer_17_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).fischer_anc_17HPF,2)
				embryo_fischer_18_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).fischer_anc_18HPF,2)
				embryo_fischer_19_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).fischer_anc_19HPF,2)
				embryo_fischer_se_0_2 = round(Embryo_SE.objects.get(nvertx_id=nvertx_2).fischer_se_0HPF,2)
				embryo_fischer_se_1_2 = round(Embryo_SE.objects.get(nvertx_id=nvertx_2).fischer_se_1HPF,2)
				embryo_fischer_se_2_2 = round(Embryo_SE.objects.get(nvertx_id=nvertx_2).fischer_se_2HPF,2)
				embryo_fischer_se_3_2 = round(Embryo_SE.objects.get(nvertx_id=nvertx_2).fischer_se_3HPF,2)
				embryo_fischer_se_4_2 = round(Embryo_SE.objects.get(nvertx_id=nvertx_2).fischer_se_4HPF,2)
				embryo_fischer_se_5_2 = round(Embryo_SE.objects.get(nvertx_id=nvertx_2).fischer_se_5HPF,2)
				embryo_fischer_se_6_2 = round(Embryo_SE.objects.get(nvertx_id=nvertx_2).fischer_se_6HPF,2)
				embryo_fischer_se_7_2 = round(Embryo_SE.objects.get(nvertx_id=nvertx_2).fischer_se_7HPF,2)
				embryo_fischer_se_8_2 = round(Embryo_SE.objects.get(nvertx_id=nvertx_2).fischer_se_8HPF,2)
				embryo_fischer_se_9_2 = round(Embryo_SE.objects.get(nvertx_id=nvertx_2).fischer_se_9HPF,2)
				#embryo_fischer_se_10_1 = [round(embryo_fischer_10_1 - Embryo_SE.objects.get(nvertx_id=nvertx_1).fischer_se_10HPF,2),round(embryo_fischer_10_1 + Embryo_SE.objects.get(nvertx_id=nvertx_1).fischer_se_10HPF,2)]
				embryo_fischer_se_11_2 = round(Embryo_SE.objects.get(nvertx_id=nvertx_2).fischer_se_11HPF,2)
				embryo_fischer_se_12_2 = round(Embryo_SE.objects.get(nvertx_id=nvertx_2).fischer_se_12HPF,2)
				embryo_fischer_se_13_2 = round(Embryo_SE.objects.get(nvertx_id=nvertx_2).fischer_se_13HPF,2)
				embryo_fischer_se_14_2 = round(Embryo_SE.objects.get(nvertx_id=nvertx_2).fischer_se_14HPF,2)
				embryo_fischer_se_15_2 = round(Embryo_SE.objects.get(nvertx_id=nvertx_2).fischer_se_15HPF,2)
				embryo_fischer_se_16_2 = round(Embryo_SE.objects.get(nvertx_id=nvertx_2).fischer_se_16HPF,2)
				embryo_fischer_se_17_2 = round(Embryo_SE.objects.get(nvertx_id=nvertx_2).fischer_se_17HPF,2)
				embryo_fischer_se_18_2 = round(Embryo_SE.objects.get(nvertx_id=nvertx_2).fischer_se_18HPF,2)
				embryo_fischer_se_19_2 = round(Embryo_SE.objects.get(nvertx_id=nvertx_2).fischer_se_19HPF,2)
			except :
				nvertx_2_embryo_fischer_invalid = True
			try :
				embryo_helm_2_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).helm_anc_2HPF,2)
				embryo_helm_7_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).helm_anc_7HPF,2)
				embryo_helm_12_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).helm_anc_12HPF,2)
				embryo_helm_24_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).helm_anc_24HPF,2)
				embryo_helm_120_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).helm_anc_120HPF,2)
				embryo_helm_240_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).helm_anc_240HPF,2)
				embryo_helm_se_2_2 = round(Embryo_SE.objects.get(nvertx_id=nvertx_2).helm_se_2HPF,2)
				embryo_helm_se_7_2 = round(Embryo_SE.objects.get(nvertx_id=nvertx_2).helm_se_7HPF,2)
				embryo_helm_se_12_2 = round(Embryo_SE.objects.get(nvertx_id=nvertx_2).helm_se_12HPF,2)
				embryo_helm_se_24_2 = round(Embryo_SE.objects.get(nvertx_id=nvertx_2).helm_se_24HPF,2)
				embryo_helm_se_120_2 = round(Embryo_SE.objects.get(nvertx_id=nvertx_2).helm_se_120HPF,2)
				embryo_helm_se_240_2 = round(Embryo_SE.objects.get(nvertx_id=nvertx_2).helm_se_240HPF,2)
			except :
				nvertx_2_embryo_helm_invalid = True
			if not nvertx_2_embryo_warner_invalid and not nvertx_2_embryo_fischer_invalid and not nvertx_2_embryo_helm_invalid :
				try :
					embryo_mean_2_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).mean_2HPF,2)
					embryo_mean_7_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).mean_7HPF,2)
					embryo_mean_12_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).mean_12HPF,2)
					embryo_mean_24_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).mean_24HPF,2)
					embryo_mean_120_2 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_2).mean_120HPF,2)
					embryo_mean_0_2 = embryo_fischer_0_2
					embryo_mean_1_2 = embryo_fischer_1_2
					embryo_mean_3_2 = embryo_fischer_3_2
					embryo_mean_4_2 = embryo_fischer_4_2
					embryo_mean_5_2 = embryo_fischer_5_2
					embryo_mean_6_2 = embryo_fischer_6_2
					embryo_mean_8_2 = embryo_fischer_8_2
					embryo_mean_9_2 = embryo_fischer_9_2
					embryo_mean_10_2 = embryo_fischer_10_2
					embryo_mean_11_2 = embryo_fischer_11_2
					embryo_mean_13_2 = embryo_fischer_13_2
					embryo_mean_14_2 = embryo_fischer_14_2
					embryo_mean_15_2 = embryo_fischer_15_2
					embryo_mean_16_2 = embryo_fischer_16_2
					embryo_mean_17_2 = embryo_fischer_17_2
					embryo_mean_18_2 = embryo_fischer_18_2
					embryo_mean_19_2 = embryo_fischer_19_2
					embryo_mean_48_2 = embryo_warner_48_2
					embryo_mean_72_2 = embryo_warner_72_2
					embryo_mean_96_2 = embryo_warner_96_2
					embryo_mean_144_2 = embryo_warner_144_2
					embryo_mean_168_2 = embryo_warner_168_2
					embryo_mean_192_2 = embryo_warner_192_2
					embryo_mean_240_2 = embryo_helm_240_2					
				except :
					nvertx_2_embryo_mean_invalid = True
			else :
				nvertx_2_embryo_mean_invalid = True
			try :
				annot_nemve1_tophit_2 = Annotation.objects.get(nvertx_id=nvertx_2).Nemve1_tophit
				annot_nemve1_e_val_2 = round(Annotation.objects.get(nvertx_id=nvertx_2).Nemve1_e_val,150)
				if annot_nemve1_e_val_2 == 99 :
					annot_nemve1_e_val_2 = "NA"
				annot_mfuzz_r_clust_2 = Annotation.objects.get(nvertx_id=nvertx_2).Mfuzz_R_Clust
				annot_mfuzz_r_score_2 = round(Annotation.objects.get(nvertx_id=nvertx_2).Mfuzz_R_Score,2)
				if annot_mfuzz_r_score_2 == -1 :
					annot_mfuzz_r_score_2 = "NA"
				annot_mfuzz_e_clust_2 = Annotation.objects.get(nvertx_id=nvertx_2).Mfuzz_E_Clust
				annot_mfuzz_e_score_2 = round(Annotation.objects.get(nvertx_id=nvertx_2).Mfuzz_E_Score,2)
				if annot_mfuzz_e_score_2 == -1 :
					annot_mfuzz_e_score_2 = "NA"
				annot_uniprot_id_2 = Annotation.objects.get(nvertx_id=nvertx_2).Uniprot_ID
				annot_uniprot_description_2 = Annotation.objects.get(nvertx_id=nvertx_2).Uniprot_Description
				annot_top_nr_hit_eval_2 = Annotation.objects.get(nvertx_id=nvertx_2).Top_nr_hit_eval
				if annot_top_nr_hit_eval_2 != "NA" :
					annot_top_nr_hit_eval_2_split = annot_top_nr_hit_eval_2.split('|',4)
					annot_nr_beg_2 = annot_top_nr_hit_eval_2_split[0] + "|" + annot_top_nr_hit_eval_2_split[1] + "|" + annot_top_nr_hit_eval_2_split[2]
					annot_nr_link_2 = annot_top_nr_hit_eval_2_split[3]
					annot_nr_end_2 = annot_top_nr_hit_eval_2_split[4]
					annot_other_nr_hits_2 = Annotation.objects.get(nvertx_id=nvertx_2).Other_nr_hits
					nr_hit_graph_2 = re.search('[\[\- \w]+\]', annot_top_nr_hit_eval_2).group(0)
			except :
				nvertx_2_annot_invalid = True
			try :
				ncbi_2 = annot_top_nr_hit_eval_2.split('|')[1]
				prot_2 = annot_top_nr_hit_eval_2.split('|')[3]
			except :
				nvertx_2_links_invalid = True

		if nvertx_3 :
			nvertx_3_embryo_warner_invalid = False
			nvertx_3_embryo_fischer_invalid = False
			nvertx_3_embryo_helm_invalid = False
			try :
				sequence_fasta_3 = Fasta.objects.get(nvertx_id=nvertx_3).fasta_sequence
			except :
				sequence_fasta_3 = 'No fasta sequence'
			try :
				if not log2 :
					regen_UC_3 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_3).regen_anc_UC+1,2),2)
					regen_0_3 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_3).regen_anc_0HPA+1,2),2)
					regen_2_3 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_3).regen_anc_2HPA+1,2),2)
					regen_4_3 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_3).regen_anc_4HPA+1,2),2)
					regen_8_3 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_3).regen_anc_8HPA+1,2),2)
					regen_12_3 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_3).regen_anc_12HPA+1,2),2)
					regen_16_3 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_3).regen_anc_16HPA+1,2),2)
					regen_20_3 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_3).regen_anc_20HPA+1,2),2)
					regen_24_3 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_3).regen_anc_24HPA+1,2),2)
					regen_36_3 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_3).regen_anc_36HPA+1,2),2)
					regen_48_3 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_3).regen_anc_48HPA+1,2),2)
					regen_60_3 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_3).regen_anc_60HPA+1,2),2)
					regen_72_3 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_3).regen_anc_72HPA+1,2),2)
					regen_96_3 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_3).regen_anc_96HPA+1,2),2)
					regen_120_3 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_3).regen_anc_120HPA+1,2),2)
					regen_144_3 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_3).regen_anc_144HPA+1,2),2)
					regen_se_UC_3 = round(Regen_log_SE.objects.get(nvertx_id=nvertx_3).regen_log_se_UC,2)
					regen_se_0_3 = round(Regen_log_SE.objects.get(nvertx_id=nvertx_3).regen_log_se_0HPA,2)
					regen_se_2_3 = round(Regen_log_SE.objects.get(nvertx_id=nvertx_3).regen_log_se_2HPA,2)
					regen_se_4_3 = round(Regen_log_SE.objects.get(nvertx_id=nvertx_3).regen_log_se_4HPA,2)
					regen_se_8_3 = round(Regen_log_SE.objects.get(nvertx_id=nvertx_3).regen_log_se_8HPA,2)
					regen_se_12_3 = round(Regen_log_SE.objects.get(nvertx_id=nvertx_3).regen_log_se_12HPA,2)
					regen_se_16_3 = round(Regen_log_SE.objects.get(nvertx_id=nvertx_3).regen_log_se_16HPA,2)
					regen_se_20_3 = round(Regen_log_SE.objects.get(nvertx_id=nvertx_3).regen_log_se_20HPA,2)
					regen_se_24_3 = round(Regen_log_SE.objects.get(nvertx_id=nvertx_3).regen_log_se_24HPA,2)
					regen_se_36_3 = round(Regen_log_SE.objects.get(nvertx_id=nvertx_3).regen_log_se_36HPA,2)
					regen_se_48_3 = round(Regen_log_SE.objects.get(nvertx_id=nvertx_3).regen_log_se_48HPA,2)
					regen_se_60_3 = round(Regen_log_SE.objects.get(nvertx_id=nvertx_3).regen_log_se_60HPA,2)
					regen_se_72_3 = round(Regen_log_SE.objects.get(nvertx_id=nvertx_3).regen_log_se_72HPA,2)
					regen_se_96_3 = round(Regen_log_SE.objects.get(nvertx_id=nvertx_3).regen_log_se_96HPA,2)
					regen_se_120_3 = round(Regen_log_SE.objects.get(nvertx_id=nvertx_3).regen_log_se_120HPA,2)
					regen_se_144_3 = round(Regen_log_SE.objects.get(nvertx_id=nvertx_3).regen_log_se_144HPA,2)
				else :
					regen_UC_3 = round(Regen_cpm.objects.get(nvertx_id=nvertx_3).regen_anc_UC,2)
					regen_0_3 = round(Regen_cpm.objects.get(nvertx_id=nvertx_3).regen_anc_0HPA,2)
					regen_2_3 = round(Regen_cpm.objects.get(nvertx_id=nvertx_3).regen_anc_2HPA,2)
					regen_4_3 = round(Regen_cpm.objects.get(nvertx_id=nvertx_3).regen_anc_4HPA,2)
					regen_8_3 = round(Regen_cpm.objects.get(nvertx_id=nvertx_3).regen_anc_8HPA,2)
					regen_12_3 = round(Regen_cpm.objects.get(nvertx_id=nvertx_3).regen_anc_12HPA,2)
					regen_16_3 = round(Regen_cpm.objects.get(nvertx_id=nvertx_3).regen_anc_16HPA,2)
					regen_20_3 = round(Regen_cpm.objects.get(nvertx_id=nvertx_3).regen_anc_20HPA,2)
					regen_24_3 = round(Regen_cpm.objects.get(nvertx_id=nvertx_3).regen_anc_24HPA,2)
					regen_36_3 = round(Regen_cpm.objects.get(nvertx_id=nvertx_3).regen_anc_36HPA,2)
					regen_48_3 = round(Regen_cpm.objects.get(nvertx_id=nvertx_3).regen_anc_48HPA,2)
					regen_60_3 = round(Regen_cpm.objects.get(nvertx_id=nvertx_3).regen_anc_60HPA,2)
					regen_72_3 = round(Regen_cpm.objects.get(nvertx_id=nvertx_3).regen_anc_72HPA,2)
					regen_96_3 = round(Regen_cpm.objects.get(nvertx_id=nvertx_3).regen_anc_96HPA,2)
					regen_120_3 = round(Regen_cpm.objects.get(nvertx_id=nvertx_3).regen_anc_120HPA,2)
					regen_144_3 = round(Regen_cpm.objects.get(nvertx_id=nvertx_3).regen_anc_144HPA,2)
					regen_se_UC_3 = round(Regen_SE.objects.get(nvertx_id=nvertx_3).regen_se_UC,2)
					regen_se_0_3 = round(Regen_SE.objects.get(nvertx_id=nvertx_3).regen_se_0HPA,2)
					regen_se_2_3 = round(Regen_SE.objects.get(nvertx_id=nvertx_3).regen_se_2HPA,2)
					regen_se_4_3 = round(Regen_SE.objects.get(nvertx_id=nvertx_3).regen_se_4HPA,2)
					regen_se_8_3 = round(Regen_SE.objects.get(nvertx_id=nvertx_3).regen_se_8HPA,2)
					regen_se_12_3 = round(Regen_SE.objects.get(nvertx_id=nvertx_3).regen_se_12HPA,2)
					regen_se_16_3 = round(Regen_SE.objects.get(nvertx_id=nvertx_3).regen_se_16HPA,2)
					regen_se_20_3 = round(Regen_SE.objects.get(nvertx_id=nvertx_3).regen_se_20HPA,2)
					regen_se_24_3 = round(Regen_SE.objects.get(nvertx_id=nvertx_3).regen_se_24HPA,2)
					regen_se_36_3 = round(Regen_SE.objects.get(nvertx_id=nvertx_3).regen_se_36HPA,2)
					regen_se_48_3 = round(Regen_SE.objects.get(nvertx_id=nvertx_3).regen_se_48HPA,2)
					regen_se_60_3 = round(Regen_SE.objects.get(nvertx_id=nvertx_3).regen_se_60HPA,2)
					regen_se_72_3 = round(Regen_SE.objects.get(nvertx_id=nvertx_3).regen_se_72HPA,2)
					regen_se_96_3 = round(Regen_SE.objects.get(nvertx_id=nvertx_3).regen_se_96HPA,2)
					regen_se_120_3 = round(Regen_SE.objects.get(nvertx_id=nvertx_3).regen_se_120HPA,2)
					regen_se_144_3 = round(Regen_SE.objects.get(nvertx_id=nvertx_3).regen_se_144HPA,2)
			except :
				nvertx_3_regen_invalid = True
			try :
				embryo_warner_24_3 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_3).warner_anc_24HPF,2)
				embryo_warner_48_3 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_3).warner_anc_48HPF,2)
				embryo_warner_72_3 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_3).warner_anc_72HPF,2)
				embryo_warner_96_3 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_3).warner_anc_96HPF,2)
				embryo_warner_120_3 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_3).warner_anc_120HPF,2)
				embryo_warner_144_3 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_3).warner_anc_144HPF,2)
				embryo_warner_168_3 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_3).warner_anc_168HPF,2)
				embryo_warner_192_3 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_3).warner_anc_192HPF,2)
				embryo_warner_se_24_3 = round(Embryo_SE.objects.get(nvertx_id=nvertx_3).warner_se_24HPF,2)
				embryo_warner_se_48_3 = round(Embryo_SE.objects.get(nvertx_id=nvertx_3).warner_se_48HPF,2)
				embryo_warner_se_72_3 = round(Embryo_SE.objects.get(nvertx_id=nvertx_3).warner_se_72HPF,2)
				embryo_warner_se_96_3 = round(Embryo_SE.objects.get(nvertx_id=nvertx_3).warner_se_96HPF,2)
				embryo_warner_se_120_3 = round(Embryo_SE.objects.get(nvertx_id=nvertx_3).warner_se_120HPF,2)
				embryo_warner_se_144_3 = round(Embryo_SE.objects.get(nvertx_id=nvertx_3).warner_se_144HPF,2)
				embryo_warner_se_168_3 = round(Embryo_SE.objects.get(nvertx_id=nvertx_3).warner_se_168HPF,2)
				embryo_warner_se_192_3 = round(Embryo_SE.objects.get(nvertx_id=nvertx_3).warner_se_192HPF,2)
			except :
				nvertx_3_embryo_warner_invalid = True
			try :
				embryo_fischer_0_3 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_3).fischer_anc_0HPF,2)
				embryo_fischer_1_3 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_3).fischer_anc_1HPF,2)
				embryo_fischer_2_3 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_3).fischer_anc_2HPF,2)
				embryo_fischer_3_3 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_3).fischer_anc_3HPF,2)
				embryo_fischer_4_3 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_3).fischer_anc_4HPF,2)
				embryo_fischer_5_3 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_3).fischer_anc_5HPF,2)
				embryo_fischer_6_3 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_3).fischer_anc_6HPF,2)
				embryo_fischer_7_3 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_3).fischer_anc_7HPF,2)
				embryo_fischer_8_3 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_3).fischer_anc_8HPF,2)
				embryo_fischer_9_3 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_3).fischer_anc_9HPF,2)	
				embryo_fischer_10_3 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_3).fischer_anc_10HPF,2)
				embryo_fischer_11_3 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_3).fischer_anc_11HPF,2)
				embryo_fischer_12_3 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_3).fischer_anc_12HPF,2)
				embryo_fischer_13_3 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_3).fischer_anc_13HPF,2)
				embryo_fischer_14_3 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_3).fischer_anc_14HPF,2)
				embryo_fischer_15_3 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_3).fischer_anc_15HPF,2)
				embryo_fischer_16_3 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_3).fischer_anc_16HPF,2)		
				embryo_fischer_17_3 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_3).fischer_anc_17HPF,2)
				embryo_fischer_18_3 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_3).fischer_anc_18HPF,2)
				embryo_fischer_19_3 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_3).fischer_anc_19HPF,2)
				embryo_fischer_se_0_3 = round(Embryo_SE.objects.get(nvertx_id=nvertx_3).fischer_se_0HPF,2)
				embryo_fischer_se_1_3 = round(Embryo_SE.objects.get(nvertx_id=nvertx_3).fischer_se_1HPF,2)
				embryo_fischer_se_2_3 = round(Embryo_SE.objects.get(nvertx_id=nvertx_3).fischer_se_2HPF,2)
				embryo_fischer_se_3_3 = round(Embryo_SE.objects.get(nvertx_id=nvertx_3).fischer_se_3HPF,2)
				embryo_fischer_se_4_3 = round(Embryo_SE.objects.get(nvertx_id=nvertx_3).fischer_se_4HPF,2)
				embryo_fischer_se_5_3 = round(Embryo_SE.objects.get(nvertx_id=nvertx_3).fischer_se_5HPF,2)
				embryo_fischer_se_6_3 = round(Embryo_SE.objects.get(nvertx_id=nvertx_3).fischer_se_6HPF,2)
				embryo_fischer_se_7_3 = round(Embryo_SE.objects.get(nvertx_id=nvertx_3).fischer_se_7HPF,2)
				embryo_fischer_se_8_3 = round(Embryo_SE.objects.get(nvertx_id=nvertx_3).fischer_se_8HPF,2)
				embryo_fischer_se_9_3 = round(Embryo_SE.objects.get(nvertx_id=nvertx_3).fischer_se_9HPF,2)
				#embryo_fischer_se_10_1 = [round(embryo_fischer_10_1 - Embryo_SE.objects.get(nvertx_id=nvertx_1).fischer_se_10HPF,2),round(embryo_fischer_10_1 + Embryo_SE.objects.get(nvertx_id=nvertx_1).fischer_se_10HPF,2)]
				embryo_fischer_se_11_3 = round(Embryo_SE.objects.get(nvertx_id=nvertx_3).fischer_se_11HPF,2)
				embryo_fischer_se_12_3 = round(Embryo_SE.objects.get(nvertx_id=nvertx_3).fischer_se_12HPF,2)
				embryo_fischer_se_13_3 = round(Embryo_SE.objects.get(nvertx_id=nvertx_3).fischer_se_13HPF,2)
				embryo_fischer_se_14_3 = round(Embryo_SE.objects.get(nvertx_id=nvertx_3).fischer_se_14HPF,2)
				embryo_fischer_se_15_3 = round(Embryo_SE.objects.get(nvertx_id=nvertx_3).fischer_se_15HPF,2)
				embryo_fischer_se_16_3 = round(Embryo_SE.objects.get(nvertx_id=nvertx_3).fischer_se_16HPF,2)
				embryo_fischer_se_17_3 = round(Embryo_SE.objects.get(nvertx_id=nvertx_3).fischer_se_17HPF,2)
				embryo_fischer_se_18_3 = round(Embryo_SE.objects.get(nvertx_id=nvertx_3).fischer_se_18HPF,2)
				embryo_fischer_se_19_3 = round(Embryo_SE.objects.get(nvertx_id=nvertx_3).fischer_se_19HPF,2)
			except :
				nvertx_3_embryo_fischer_invalid = True
			try :
				embryo_helm_2_3 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_3).helm_anc_2HPF,2)
				embryo_helm_7_3 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_3).helm_anc_7HPF,2)
				embryo_helm_12_3 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_3).helm_anc_12HPF,2)
				embryo_helm_24_3 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_3).helm_anc_24HPF,2)
				embryo_helm_120_3 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_3).helm_anc_120HPF,2)
				embryo_helm_240_3 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_3).helm_anc_240HPF,2)
				embryo_helm_se_2_3 = round(Embryo_SE.objects.get(nvertx_id=nvertx_3).helm_se_2HPF,2)
				embryo_helm_se_7_3 = round(Embryo_SE.objects.get(nvertx_id=nvertx_3).helm_se_7HPF,2)
				embryo_helm_se_12_3 = round(Embryo_SE.objects.get(nvertx_id=nvertx_3).helm_se_12HPF,2)
				embryo_helm_se_24_3 = round(Embryo_SE.objects.get(nvertx_id=nvertx_3).helm_se_24HPF,2)
				embryo_helm_se_120_3 = round(Embryo_SE.objects.get(nvertx_id=nvertx_3).helm_se_120HPF,2)
				embryo_helm_se_240_3 = round(Embryo_SE.objects.get(nvertx_id=nvertx_3).helm_se_240HPF,2)
			except :
				nvertx_3_embryo_helm_invalid = True
			if not nvertx_3_embryo_warner_invalid and not nvertx_3_embryo_fischer_invalid and not nvertx_3_embryo_helm_invalid :
				try :
					embryo_mean_2_3 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_3).mean_2HPF,2)
					embryo_mean_7_3 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_3).mean_7HPF,2)
					embryo_mean_12_3 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_3).mean_12HPF,2)
					embryo_mean_24_3 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_3).mean_24HPF,2)
					embryo_mean_120_3 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_3).mean_120HPF,2)
					embryo_mean_0_3 = embryo_fischer_0_3
					embryo_mean_1_3 = embryo_fischer_1_3
					embryo_mean_3_3 = embryo_fischer_3_3
					embryo_mean_4_3 = embryo_fischer_4_3
					embryo_mean_5_3 = embryo_fischer_5_3
					embryo_mean_6_3 = embryo_fischer_6_3
					embryo_mean_8_3 = embryo_fischer_8_3
					embryo_mean_9_3 = embryo_fischer_9_3
					embryo_mean_10_3 = embryo_fischer_10_3
					embryo_mean_11_3 = embryo_fischer_11_3
					embryo_mean_13_3 = embryo_fischer_13_3
					embryo_mean_14_3 = embryo_fischer_14_3
					embryo_mean_15_3 = embryo_fischer_15_3
					embryo_mean_16_3 = embryo_fischer_16_3
					embryo_mean_17_3 = embryo_fischer_17_3
					embryo_mean_18_3 = embryo_fischer_18_3
					embryo_mean_19_3 = embryo_fischer_19_3
					embryo_mean_48_3 = embryo_warner_48_3
					embryo_mean_72_3 = embryo_warner_72_3
					embryo_mean_96_3 = embryo_warner_96_3
					embryo_mean_144_3 = embryo_warner_144_3
					embryo_mean_168_3 = embryo_warner_168_3
					embryo_mean_192_3 = embryo_warner_192_3
					embryo_mean_240_3 = embryo_helm_240_3
				except :
					nvertx_3_embryo_mean_invalid = True
			else :
				nvertx_3_embryo_mean_invalid = True
			try :
				annot_nemve1_tophit_3 = Annotation.objects.get(nvertx_id=nvertx_3).Nemve1_tophit
				annot_nemve1_e_val_3 = round(Annotation.objects.get(nvertx_id=nvertx_3).Nemve1_e_val,150)
				if annot_nemve1_e_val_3 == 99 :
					annot_nemve1_e_val_3 = "NA"
				annot_mfuzz_r_clust_3 = Annotation.objects.get(nvertx_id=nvertx_3).Mfuzz_R_Clust
				annot_mfuzz_r_score_3 = round(Annotation.objects.get(nvertx_id=nvertx_3).Mfuzz_R_Score,2)
				if annot_mfuzz_r_score_3 == -1 :
					annot_mfuzz_r_score_3 = "NA"
				annot_mfuzz_e_clust_3 = Annotation.objects.get(nvertx_id=nvertx_3).Mfuzz_E_Clust
				annot_mfuzz_e_score_3 = round(Annotation.objects.get(nvertx_id=nvertx_3).Mfuzz_E_Score,2)
				if annot_mfuzz_e_score_3 == -1 :
					annot_mfuzz_e_score_3 = "NA"
				annot_uniprot_id_3 = Annotation.objects.get(nvertx_id=nvertx_3).Uniprot_ID
				annot_uniprot_description_3 = Annotation.objects.get(nvertx_id=nvertx_3).Uniprot_Description
				annot_top_nr_hit_eval_3 = Annotation.objects.get(nvertx_id=nvertx_3).Top_nr_hit_eval
				if annot_top_nr_hit_eval_3 != "NA" :
					annot_top_nr_hit_eval_3_split = annot_top_nr_hit_eval_3.split('|',4)
					annot_nr_beg_3 = annot_top_nr_hit_eval_3_split[0] + "|" + annot_top_nr_hit_eval_3_split[1] + "|" + annot_top_nr_hit_eval_3_split[2]
					annot_nr_link_3 = annot_top_nr_hit_eval_3_split[3]
					annot_nr_end_3 = annot_top_nr_hit_eval_3_split[4]
					annot_other_nr_hits_3 = Annotation.objects.get(nvertx_id=nvertx_3).Other_nr_hits
					nr_hit_graph_3 = re.search('[\[\- \w]+\]', annot_top_nr_hit_eval_3).group(0)
			except :
				nvertx_3_annot_invalid = True
			try :
				ncbi_3 = annot_top_nr_hit_eval_3.split('|')[1]
				prot_3 = annot_top_nr_hit_eval_3.split('|')[3]
			except :
				nvertx_3_links_invalid = True

		if nvertx_4 :
			nvertx_4_embryo_warner_invalid = False
			nvertx_4_embryo_fischer_invalid = False
			nvertx_4_embryo_helm_invalid = False
			try :
				sequence_fasta_4 = Fasta.objects.get(nvertx_id=nvertx_4).fasta_sequence
			except :
				sequence_fasta_4 = 'No fasta sequence'
			try :
				if not log2 :
					regen_UC_4 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_4).regen_anc_UC+1,2),2)
					regen_0_4 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_4).regen_anc_0HPA+1,2),2)
					regen_2_4 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_4).regen_anc_2HPA+1,2),2)
					regen_4_4 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_4).regen_anc_4HPA+1,2),2)
					regen_8_4 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_4).regen_anc_8HPA+1,2),2)
					regen_12_4 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_4).regen_anc_12HPA+1,2),2)
					regen_16_4 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_4).regen_anc_16HPA+1,2),2)
					regen_20_4 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_4).regen_anc_20HPA+1,2),2)
					regen_24_4 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_4).regen_anc_24HPA+1,2),2)
					regen_36_4 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_4).regen_anc_36HPA+1,2),2)
					regen_48_4 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_4).regen_anc_48HPA+1,2),2)
					regen_60_4 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_4).regen_anc_60HPA+1,2),2)
					regen_72_4 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_4).regen_anc_72HPA+1,2),2)
					regen_96_4 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_4).regen_anc_96HPA+1,2),2)
					regen_120_4 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_4).regen_anc_120HPA+1,2),2)
					regen_144_4 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_4).regen_anc_144HPA+1,2),2)
					regen_se_UC_4 = round(Regen_log_SE.objects.get(nvertx_id=nvertx_4).regen_log_se_UC,2)
					regen_se_0_4 = round(Regen_log_SE.objects.get(nvertx_id=nvertx_4).regen_log_se_0HPA,2)
					regen_se_2_4 = round(Regen_log_SE.objects.get(nvertx_id=nvertx_4).regen_log_se_2HPA,2)
					regen_se_4_4 = round(Regen_log_SE.objects.get(nvertx_id=nvertx_4).regen_log_se_4HPA,2)
					regen_se_8_4 = round(Regen_log_SE.objects.get(nvertx_id=nvertx_4).regen_log_se_8HPA,2)
					regen_se_12_4 = round(Regen_log_SE.objects.get(nvertx_id=nvertx_4).regen_log_se_12HPA,2)
					regen_se_16_4 = round(Regen_log_SE.objects.get(nvertx_id=nvertx_4).regen_log_se_16HPA,2)
					regen_se_20_4 = round(Regen_log_SE.objects.get(nvertx_id=nvertx_4).regen_log_se_20HPA,2)
					regen_se_24_4 = round(Regen_log_SE.objects.get(nvertx_id=nvertx_4).regen_log_se_24HPA,2)
					regen_se_36_4 = round(Regen_log_SE.objects.get(nvertx_id=nvertx_4).regen_log_se_36HPA,2)
					regen_se_48_4 = round(Regen_log_SE.objects.get(nvertx_id=nvertx_4).regen_log_se_48HPA,2)
					regen_se_60_4 = round(Regen_log_SE.objects.get(nvertx_id=nvertx_4).regen_log_se_60HPA,2)
					regen_se_72_4 = round(Regen_log_SE.objects.get(nvertx_id=nvertx_4).regen_log_se_72HPA,2)
					regen_se_96_4 = round(Regen_log_SE.objects.get(nvertx_id=nvertx_4).regen_log_se_96HPA,2)
					regen_se_120_4 = round(Regen_log_SE.objects.get(nvertx_id=nvertx_4).regen_log_se_120HPA,2)
					regen_se_144_4 = round(Regen_log_SE.objects.get(nvertx_id=nvertx_4).regen_log_se_144HPA,2)
				else :
					regen_UC_4 = round(Regen_cpm.objects.get(nvertx_id=nvertx_4).regen_anc_UC,2)
					regen_0_4 = round(Regen_cpm.objects.get(nvertx_id=nvertx_4).regen_anc_0HPA,2)
					regen_2_4 = round(Regen_cpm.objects.get(nvertx_id=nvertx_4).regen_anc_2HPA,2)
					regen_4_4 = round(Regen_cpm.objects.get(nvertx_id=nvertx_4).regen_anc_4HPA,2)
					regen_8_4 = round(Regen_cpm.objects.get(nvertx_id=nvertx_4).regen_anc_8HPA,2)
					regen_12_4 = round(Regen_cpm.objects.get(nvertx_id=nvertx_4).regen_anc_12HPA,2)
					regen_16_4 = round(Regen_cpm.objects.get(nvertx_id=nvertx_4).regen_anc_16HPA,2)
					regen_20_4 = round(Regen_cpm.objects.get(nvertx_id=nvertx_4).regen_anc_20HPA,2)
					regen_24_4 = round(Regen_cpm.objects.get(nvertx_id=nvertx_4).regen_anc_24HPA,2)
					regen_36_4 = round(Regen_cpm.objects.get(nvertx_id=nvertx_4).regen_anc_36HPA,2)
					regen_48_4 = round(Regen_cpm.objects.get(nvertx_id=nvertx_4).regen_anc_48HPA,2)
					regen_60_4 = round(Regen_cpm.objects.get(nvertx_id=nvertx_4).regen_anc_60HPA,2)
					regen_72_4 = round(Regen_cpm.objects.get(nvertx_id=nvertx_4).regen_anc_72HPA,2)
					regen_96_4 = round(Regen_cpm.objects.get(nvertx_id=nvertx_4).regen_anc_96HPA,2)
					regen_120_4 = round(Regen_cpm.objects.get(nvertx_id=nvertx_4).regen_anc_120HPA,2)
					regen_144_4 = round(Regen_cpm.objects.get(nvertx_id=nvertx_4).regen_anc_144HPA,2)
					regen_se_UC_4 = round(Regen_SE.objects.get(nvertx_id=nvertx_4).regen_se_UC,2)
					regen_se_0_4 = round(Regen_SE.objects.get(nvertx_id=nvertx_4).regen_se_0HPA,2)
					regen_se_2_4 = round(Regen_SE.objects.get(nvertx_id=nvertx_4).regen_se_2HPA,2)
					regen_se_4_4 = round(Regen_SE.objects.get(nvertx_id=nvertx_4).regen_se_4HPA,2)
					regen_se_8_4 = round(Regen_SE.objects.get(nvertx_id=nvertx_4).regen_se_8HPA,2)
					regen_se_12_4 = round(Regen_SE.objects.get(nvertx_id=nvertx_4).regen_se_12HPA,2)
					regen_se_16_4 = round(Regen_SE.objects.get(nvertx_id=nvertx_4).regen_se_16HPA,2)
					regen_se_20_4 = round(Regen_SE.objects.get(nvertx_id=nvertx_4).regen_se_20HPA,2)
					regen_se_24_4 = round(Regen_SE.objects.get(nvertx_id=nvertx_4).regen_se_24HPA,2)
					regen_se_36_4 = round(Regen_SE.objects.get(nvertx_id=nvertx_4).regen_se_36HPA,2)
					regen_se_48_4 = round(Regen_SE.objects.get(nvertx_id=nvertx_4).regen_se_48HPA,2)
					regen_se_60_4 = round(Regen_SE.objects.get(nvertx_id=nvertx_4).regen_se_60HPA,2)
					regen_se_72_4 = round(Regen_SE.objects.get(nvertx_id=nvertx_4).regen_se_72HPA,2)
					regen_se_96_4 = round(Regen_SE.objects.get(nvertx_id=nvertx_4).regen_se_96HPA,2)
					regen_se_120_4 = round(Regen_SE.objects.get(nvertx_id=nvertx_4).regen_se_120HPA,2)
					regen_se_144_4 = round(Regen_SE.objects.get(nvertx_id=nvertx_4).regen_se_144HPA,2)
			except :
				nvertx_4_regen_invalid = True
			try :
				embryo_warner_24_4 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_4).warner_anc_24HPF,2)
				embryo_warner_48_4 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_4).warner_anc_48HPF,2)
				embryo_warner_72_4 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_4).warner_anc_72HPF,2)
				embryo_warner_96_4 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_4).warner_anc_96HPF,2)
				embryo_warner_120_4 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_4).warner_anc_120HPF,2)
				embryo_warner_144_4 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_4).warner_anc_144HPF,2)
				embryo_warner_168_4 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_4).warner_anc_168HPF,2)
				embryo_warner_192_4 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_4).warner_anc_192HPF,2)
				embryo_warner_se_24_4 = round(Embryo_SE.objects.get(nvertx_id=nvertx_4).warner_se_24HPF,2)
				embryo_warner_se_48_4 = round(Embryo_SE.objects.get(nvertx_id=nvertx_4).warner_se_48HPF,2)
				embryo_warner_se_72_4 = round(Embryo_SE.objects.get(nvertx_id=nvertx_4).warner_se_72HPF,2)
				embryo_warner_se_96_4 = round(Embryo_SE.objects.get(nvertx_id=nvertx_4).warner_se_96HPF,2)
				embryo_warner_se_120_4 = round(Embryo_SE.objects.get(nvertx_id=nvertx_4).warner_se_120HPF,2)
				embryo_warner_se_144_4 = round(Embryo_SE.objects.get(nvertx_id=nvertx_4).warner_se_144HPF,2)
				embryo_warner_se_168_4 = round(Embryo_SE.objects.get(nvertx_id=nvertx_4).warner_se_168HPF,2)
				embryo_warner_se_192_4 = round(Embryo_SE.objects.get(nvertx_id=nvertx_4).warner_se_192HPF,2)
			except :
				nvertx_4_embryo_warner_invalid = True
			try :
				embryo_fischer_0_4 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_4).fischer_anc_0HPF,2)
				embryo_fischer_1_4 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_4).fischer_anc_1HPF,2)
				embryo_fischer_2_4 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_4).fischer_anc_2HPF,2)
				embryo_fischer_3_4 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_4).fischer_anc_3HPF,2)
				embryo_fischer_4_4 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_4).fischer_anc_4HPF,2)
				embryo_fischer_5_4 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_4).fischer_anc_5HPF,2)
				embryo_fischer_6_4 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_4).fischer_anc_6HPF,2)
				embryo_fischer_7_4 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_4).fischer_anc_7HPF,2)
				embryo_fischer_8_4 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_4).fischer_anc_8HPF,2)
				embryo_fischer_9_4 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_4).fischer_anc_9HPF,2)	
				embryo_fischer_10_4 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_4).fischer_anc_10HPF,2)
				embryo_fischer_11_4 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_4).fischer_anc_11HPF,2)
				embryo_fischer_12_4 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_4).fischer_anc_12HPF,2)
				embryo_fischer_13_4 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_4).fischer_anc_13HPF,2)
				embryo_fischer_14_4 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_4).fischer_anc_14HPF,2)
				embryo_fischer_15_4 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_4).fischer_anc_15HPF,2)
				embryo_fischer_16_4 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_4).fischer_anc_16HPF,2)		
				embryo_fischer_17_4 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_4).fischer_anc_17HPF,2)
				embryo_fischer_18_4 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_4).fischer_anc_18HPF,2)
				embryo_fischer_19_4 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_4).fischer_anc_19HPF,2)
				embryo_fischer_se_0_4 = round(Embryo_SE.objects.get(nvertx_id=nvertx_4).fischer_se_0HPF,2)
				embryo_fischer_se_1_4 = round(Embryo_SE.objects.get(nvertx_id=nvertx_4).fischer_se_1HPF,2)
				embryo_fischer_se_2_4 = round(Embryo_SE.objects.get(nvertx_id=nvertx_4).fischer_se_2HPF,2)
				embryo_fischer_se_3_4 = round(Embryo_SE.objects.get(nvertx_id=nvertx_4).fischer_se_3HPF,2)
				embryo_fischer_se_4_4 = round(Embryo_SE.objects.get(nvertx_id=nvertx_4).fischer_se_4HPF,2)
				embryo_fischer_se_5_4 = round(Embryo_SE.objects.get(nvertx_id=nvertx_4).fischer_se_5HPF,2)
				embryo_fischer_se_6_4 = round(Embryo_SE.objects.get(nvertx_id=nvertx_4).fischer_se_6HPF,2)
				embryo_fischer_se_7_4 = round(Embryo_SE.objects.get(nvertx_id=nvertx_4).fischer_se_7HPF,2)
				embryo_fischer_se_8_4 = round(Embryo_SE.objects.get(nvertx_id=nvertx_4).fischer_se_8HPF,2)
				embryo_fischer_se_9_4 = round(Embryo_SE.objects.get(nvertx_id=nvertx_4).fischer_se_9HPF,2)
				#embryo_fischer_se_10_1 = [round(embryo_fischer_10_1 - Embryo_SE.objects.get(nvertx_id=nvertx_1).fischer_se_10HPF,2),round(embryo_fischer_10_1 + Embryo_SE.objects.get(nvertx_id=nvertx_1).fischer_se_10HPF,2)]
				embryo_fischer_se_11_4 = round(Embryo_SE.objects.get(nvertx_id=nvertx_4).fischer_se_11HPF,2)
				embryo_fischer_se_12_4 = round(Embryo_SE.objects.get(nvertx_id=nvertx_4).fischer_se_12HPF,2)
				embryo_fischer_se_13_4 = round(Embryo_SE.objects.get(nvertx_id=nvertx_4).fischer_se_13HPF,2)
				embryo_fischer_se_14_4 = round(Embryo_SE.objects.get(nvertx_id=nvertx_4).fischer_se_14HPF,2)
				embryo_fischer_se_15_4 = round(Embryo_SE.objects.get(nvertx_id=nvertx_4).fischer_se_15HPF,2)
				embryo_fischer_se_16_4 = round(Embryo_SE.objects.get(nvertx_id=nvertx_4).fischer_se_16HPF,2)
				embryo_fischer_se_17_4 = round(Embryo_SE.objects.get(nvertx_id=nvertx_4).fischer_se_17HPF,2)
				embryo_fischer_se_18_4 = round(Embryo_SE.objects.get(nvertx_id=nvertx_4).fischer_se_18HPF,2)
				embryo_fischer_se_19_4 = round(Embryo_SE.objects.get(nvertx_id=nvertx_4).fischer_se_19HPF,2)
			except :
				nvertx_4_embryo_fischer_invalid = True
			try :
				embryo_helm_2_4 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_4).helm_anc_2HPF,2)
				embryo_helm_7_4 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_4).helm_anc_7HPF,2)
				embryo_helm_12_4 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_4).helm_anc_12HPF,2)
				embryo_helm_24_4 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_4).helm_anc_24HPF,2)
				embryo_helm_120_4 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_4).helm_anc_120HPF,2)
				embryo_helm_240_4 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_4).helm_anc_240HPF,2)
				embryo_helm_se_2_4 = round(Embryo_SE.objects.get(nvertx_id=nvertx_4).helm_se_2HPF,2)
				embryo_helm_se_7_4 = round(Embryo_SE.objects.get(nvertx_id=nvertx_4).helm_se_7HPF,2)
				embryo_helm_se_12_4 = round(Embryo_SE.objects.get(nvertx_id=nvertx_4).helm_se_12HPF,2)
				embryo_helm_se_24_4 = round(Embryo_SE.objects.get(nvertx_id=nvertx_4).helm_se_24HPF,2)
				embryo_helm_se_120_4 = round(Embryo_SE.objects.get(nvertx_id=nvertx_4).helm_se_120HPF,2)
				embryo_helm_se_240_4 = round(Embryo_SE.objects.get(nvertx_id=nvertx_4).helm_se_240HPF,2)
			except :
				nvertx_4_embryo_helm_invalid = True
			if not nvertx_4_embryo_warner_invalid and not nvertx_4_embryo_fischer_invalid and not nvertx_4_embryo_helm_invalid :
				try :
					embryo_mean_2_4 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_4).mean_2HPF,2)
					embryo_mean_7_4 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_4).mean_7HPF,2)
					embryo_mean_12_4 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_4).mean_12HPF,2)
					embryo_mean_24_4 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_4).mean_24HPF,2)
					embryo_mean_120_4 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_4).mean_120HPF,2)
					embryo_mean_0_4 = embryo_fischer_0_4
					embryo_mean_1_4 = embryo_fischer_1_4
					embryo_mean_3_4 = embryo_fischer_3_4
					embryo_mean_4_4 = embryo_fischer_4_4
					embryo_mean_5_4 = embryo_fischer_5_4
					embryo_mean_6_4 = embryo_fischer_6_4
					embryo_mean_8_4 = embryo_fischer_8_4
					embryo_mean_9_4 = embryo_fischer_9_4
					embryo_mean_10_4 = embryo_fischer_10_4
					embryo_mean_11_4 = embryo_fischer_11_4
					embryo_mean_13_4 = embryo_fischer_13_4
					embryo_mean_14_4 = embryo_fischer_14_4
					embryo_mean_15_4 = embryo_fischer_15_4
					embryo_mean_16_4 = embryo_fischer_16_4
					embryo_mean_17_4 = embryo_fischer_17_4
					embryo_mean_18_4 = embryo_fischer_18_4
					embryo_mean_19_4 = embryo_fischer_19_4
					embryo_mean_48_4 = embryo_warner_48_4
					embryo_mean_72_4 = embryo_warner_72_4
					embryo_mean_96_4 = embryo_warner_96_4
					embryo_mean_144_4 = embryo_warner_144_4
					embryo_mean_168_4 = embryo_warner_168_4
					embryo_mean_192_4 = embryo_warner_192_4
					embryo_mean_240_4 = embryo_helm_240_4
				except :
					nvertx_4_embryo_mean_invalid = True
			else :
				nvertx_4_embryo_mean_invalid = True
			try :
				annot_nemve1_tophit_4 = Annotation.objects.get(nvertx_id=nvertx_4).Nemve1_tophit
				annot_nemve1_e_val_4 = round(Annotation.objects.get(nvertx_id=nvertx_4).Nemve1_e_val,150)
				if annot_nemve1_e_val_4 == 99 :
					annot_nemve1_e_val_4 = "NA"
				annot_mfuzz_r_clust_4 = Annotation.objects.get(nvertx_id=nvertx_4).Mfuzz_R_Clust
				annot_mfuzz_r_score_4 = round(Annotation.objects.get(nvertx_id=nvertx_4).Mfuzz_R_Score,2)
				if annot_mfuzz_r_score_4 == -1 :
					annot_mfuzz_r_score_4 = "NA"
				annot_mfuzz_e_clust_4 = Annotation.objects.get(nvertx_id=nvertx_4).Mfuzz_E_Clust
				annot_mfuzz_e_score_4 = round(Annotation.objects.get(nvertx_id=nvertx_4).Mfuzz_E_Score,2)
				if annot_mfuzz_e_score_4 == -1 :
					annot_mfuzz_e_score_4 = "NA"
				annot_uniprot_id_4 = Annotation.objects.get(nvertx_id=nvertx_4).Uniprot_ID
				annot_uniprot_description_4 = Annotation.objects.get(nvertx_id=nvertx_4).Uniprot_Description
				annot_top_nr_hit_eval_4 = Annotation.objects.get(nvertx_id=nvertx_4).Top_nr_hit_eval
				if annot_top_nr_hit_eval_4 != "NA" :
					annot_top_nr_hit_eval_4_split = annot_top_nr_hit_eval_4.split('|',4)
					annot_nr_beg_4 = annot_top_nr_hit_eval_4_split[0] + "|" + annot_top_nr_hit_eval_4_split[1] + "|" + annot_top_nr_hit_eval_4_split[2]
					annot_nr_link_4 = annot_top_nr_hit_eval_4_split[3]
					annot_nr_end_4 = annot_top_nr_hit_eval_4_split[4]
					annot_other_nr_hits_4 = Annotation.objects.get(nvertx_id=nvertx_4).Other_nr_hits
					nr_hit_graph_4 = re.search('[\[\- \w]+\]', annot_top_nr_hit_eval_4).group(0)
			except :
				nvertx_4_annot_invalid = True
			try :
				ncbi_4 = annot_top_nr_hit_eval_4.split('|')[1]
				prot_4 = annot_top_nr_hit_eval_4.split('|')[3]
			except :
				nvertx_4_links_invalid = True

		if nvertx_5 :
			nvertx_5_embryo_warner_invalid = False
			nvertx_5_embryo_fischer_invalid = False
			nvertx_5_embryo_helm_invalid = False
			try :
				sequence_fasta_5 = Fasta.objects.get(nvertx_id=nvertx_5).fasta_sequence
			except :
				sequence_fasta_5 = 'No fasta sequence'
			try :
				if not log2 :
					regen_UC_5 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_5).regen_anc_UC+1,2),2)
					regen_0_5 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_5).regen_anc_0HPA+1,2),2)
					regen_2_5 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_5).regen_anc_2HPA+1,2),2)
					regen_4_5 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_5).regen_anc_4HPA+1,2),2)
					regen_8_5 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_5).regen_anc_8HPA+1,2),2)
					regen_12_5 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_5).regen_anc_12HPA+1,2),2)
					regen_16_5 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_5).regen_anc_16HPA+1,2),2)
					regen_20_5 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_5).regen_anc_20HPA+1,2),2)
					regen_24_5 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_5).regen_anc_24HPA+1,2),2)
					regen_36_5 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_5).regen_anc_36HPA+1,2),2)
					regen_48_5 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_5).regen_anc_48HPA+1,2),2)
					regen_60_5 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_5).regen_anc_60HPA+1,2),2)
					regen_72_5 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_5).regen_anc_72HPA+1,2),2)
					regen_96_5 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_5).regen_anc_96HPA+1,2),2)
					regen_120_5 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_5).regen_anc_120HPA+1,2),2)
					regen_144_5 = round(math.log(Regen_cpm.objects.get(nvertx_id=nvertx_5).regen_anc_144HPA+1,2),2)
					regen_se_UC_5 = round(Regen_log_SE.objects.get(nvertx_id=nvertx_5).regen_log_se_UC,2)
					regen_se_0_5 = round(Regen_log_SE.objects.get(nvertx_id=nvertx_5).regen_log_se_0HPA,2)
					regen_se_2_5 = round(Regen_log_SE.objects.get(nvertx_id=nvertx_5).regen_log_se_2HPA,2)
					regen_se_4_5 = round(Regen_log_SE.objects.get(nvertx_id=nvertx_5).regen_log_se_4HPA,2)
					regen_se_8_5 = round(Regen_log_SE.objects.get(nvertx_id=nvertx_5).regen_log_se_8HPA,2)
					regen_se_12_5 = round(Regen_log_SE.objects.get(nvertx_id=nvertx_5).regen_log_se_12HPA,2)
					regen_se_16_5 = round(Regen_log_SE.objects.get(nvertx_id=nvertx_5).regen_log_se_16HPA,2)
					regen_se_20_5 = round(Regen_log_SE.objects.get(nvertx_id=nvertx_5).regen_log_se_20HPA,2)
					regen_se_24_5 = round(Regen_log_SE.objects.get(nvertx_id=nvertx_5).regen_log_se_24HPA,2)
					regen_se_36_5 = round(Regen_log_SE.objects.get(nvertx_id=nvertx_5).regen_log_se_36HPA,2)
					regen_se_48_5 = round(Regen_log_SE.objects.get(nvertx_id=nvertx_5).regen_log_se_48HPA,2)
					regen_se_60_5 = round(Regen_log_SE.objects.get(nvertx_id=nvertx_5).regen_log_se_60HPA,2)
					regen_se_72_5 = round(Regen_log_SE.objects.get(nvertx_id=nvertx_5).regen_log_se_72HPA,2)
					regen_se_96_5 = round(Regen_log_SE.objects.get(nvertx_id=nvertx_5).regen_log_se_96HPA,2)
					regen_se_120_5 = round(Regen_log_SE.objects.get(nvertx_id=nvertx_5).regen_log_se_120HPA,2)
					regen_se_144_5 = round(Regen_log_SE.objects.get(nvertx_id=nvertx_5).regen_log_se_144HPA,2)
				else :
					regen_UC_5 = round(Regen_cpm.objects.get(nvertx_id=nvertx_5).regen_anc_UC,2)
					regen_0_5 = round(Regen_cpm.objects.get(nvertx_id=nvertx_5).regen_anc_0HPA,2)
					regen_2_5 = round(Regen_cpm.objects.get(nvertx_id=nvertx_5).regen_anc_2HPA,2)
					regen_4_5 = round(Regen_cpm.objects.get(nvertx_id=nvertx_5).regen_anc_4HPA,2)
					regen_8_5 = round(Regen_cpm.objects.get(nvertx_id=nvertx_5).regen_anc_8HPA,2)
					regen_12_5 = round(Regen_cpm.objects.get(nvertx_id=nvertx_5).regen_anc_12HPA,2)
					regen_16_5 = round(Regen_cpm.objects.get(nvertx_id=nvertx_5).regen_anc_16HPA,2)
					regen_20_5 = round(Regen_cpm.objects.get(nvertx_id=nvertx_5).regen_anc_20HPA,2)
					regen_24_5 = round(Regen_cpm.objects.get(nvertx_id=nvertx_5).regen_anc_24HPA,2)
					regen_36_5 = round(Regen_cpm.objects.get(nvertx_id=nvertx_5).regen_anc_36HPA,2)
					regen_48_5 = round(Regen_cpm.objects.get(nvertx_id=nvertx_5).regen_anc_48HPA,2)
					regen_60_5 = round(Regen_cpm.objects.get(nvertx_id=nvertx_5).regen_anc_60HPA,2)
					regen_72_5 = round(Regen_cpm.objects.get(nvertx_id=nvertx_5).regen_anc_72HPA,2)
					regen_96_5 = round(Regen_cpm.objects.get(nvertx_id=nvertx_5).regen_anc_96HPA,2)
					regen_120_5 = round(Regen_cpm.objects.get(nvertx_id=nvertx_5).regen_anc_120HPA,2)
					regen_144_5 = round(Regen_cpm.objects.get(nvertx_id=nvertx_5).regen_anc_144HPA,2)
					regen_se_UC_5 = round(Regen_SE.objects.get(nvertx_id=nvertx_5).regen_se_UC,2)
					regen_se_0_5 = round(Regen_SE.objects.get(nvertx_id=nvertx_5).regen_se_0HPA,2)
					regen_se_2_5 = round(Regen_SE.objects.get(nvertx_id=nvertx_5).regen_se_2HPA,2)
					regen_se_4_5 = round(Regen_SE.objects.get(nvertx_id=nvertx_5).regen_se_4HPA,2)
					regen_se_8_5 = round(Regen_SE.objects.get(nvertx_id=nvertx_5).regen_se_8HPA,2)
					regen_se_12_5 = round(Regen_SE.objects.get(nvertx_id=nvertx_5).regen_se_12HPA,2)
					regen_se_16_5 = round(Regen_SE.objects.get(nvertx_id=nvertx_5).regen_se_16HPA,2)
					regen_se_20_5 = round(Regen_SE.objects.get(nvertx_id=nvertx_5).regen_se_20HPA,2)
					regen_se_24_5 = round(Regen_SE.objects.get(nvertx_id=nvertx_5).regen_se_24HPA,2)
					regen_se_36_5 = round(Regen_SE.objects.get(nvertx_id=nvertx_5).regen_se_36HPA,2)
					regen_se_48_5 = round(Regen_SE.objects.get(nvertx_id=nvertx_5).regen_se_48HPA,2)
					regen_se_60_5 = round(Regen_SE.objects.get(nvertx_id=nvertx_5).regen_se_60HPA,2)
					regen_se_72_5 = round(Regen_SE.objects.get(nvertx_id=nvertx_5).regen_se_72HPA,2)
					regen_se_96_5 = round(Regen_SE.objects.get(nvertx_id=nvertx_5).regen_se_96HPA,2)
					regen_se_120_5 = round(Regen_SE.objects.get(nvertx_id=nvertx_5).regen_se_120HPA,2)
					regen_se_144_5 = round(Regen_SE.objects.get(nvertx_id=nvertx_5).regen_se_144HPA,2)
			except :
				nvertx_5_regen_invalid = True
			try :
				embryo_warner_24_5 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_5).warner_anc_24HPF,2)
				embryo_warner_48_5 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_5).warner_anc_48HPF,2)
				embryo_warner_72_5 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_5).warner_anc_72HPF,2)
				embryo_warner_96_5 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_5).warner_anc_96HPF,2)
				embryo_warner_120_5 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_5).warner_anc_120HPF,2)
				embryo_warner_144_5 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_5).warner_anc_144HPF,2)
				embryo_warner_168_5 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_5).warner_anc_168HPF,2)
				embryo_warner_192_5 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_5).warner_anc_192HPF,2)
				embryo_warner_se_24_5 = round(Embryo_SE.objects.get(nvertx_id=nvertx_5).warner_se_24HPF,2)
				embryo_warner_se_48_5 = round(Embryo_SE.objects.get(nvertx_id=nvertx_5).warner_se_48HPF,2)
				embryo_warner_se_72_5 = round(Embryo_SE.objects.get(nvertx_id=nvertx_5).warner_se_72HPF,2)
				embryo_warner_se_96_5 = round(Embryo_SE.objects.get(nvertx_id=nvertx_5).warner_se_96HPF,2)
				embryo_warner_se_120_5 = round(Embryo_SE.objects.get(nvertx_id=nvertx_5).warner_se_120HPF,2)
				embryo_warner_se_144_5 = round(Embryo_SE.objects.get(nvertx_id=nvertx_5).warner_se_144HPF,2)
				embryo_warner_se_168_5 = round(Embryo_SE.objects.get(nvertx_id=nvertx_5).warner_se_168HPF,2)
				embryo_warner_se_192_5 = round(Embryo_SE.objects.get(nvertx_id=nvertx_5).warner_se_192HPF,2)
			except :
				nvertx_5_embryo_warner_invalid = True
			try :
				embryo_fischer_0_5 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_5).fischer_anc_0HPF,2)
				embryo_fischer_1_5 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_5).fischer_anc_1HPF,2)
				embryo_fischer_2_5 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_5).fischer_anc_2HPF,2)
				embryo_fischer_3_5 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_5).fischer_anc_3HPF,2)
				embryo_fischer_4_5 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_5).fischer_anc_4HPF,2)
				embryo_fischer_5_5 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_5).fischer_anc_5HPF,2)
				embryo_fischer_6_5 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_5).fischer_anc_6HPF,2)
				embryo_fischer_7_5 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_5).fischer_anc_7HPF,2)
				embryo_fischer_8_5 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_5).fischer_anc_8HPF,2)
				embryo_fischer_9_5 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_5).fischer_anc_9HPF,2)	
				embryo_fischer_10_5 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_5).fischer_anc_10HPF,2)
				embryo_fischer_11_5 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_5).fischer_anc_11HPF,2)
				embryo_fischer_12_5 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_5).fischer_anc_12HPF,2)
				embryo_fischer_13_5 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_5).fischer_anc_13HPF,2)
				embryo_fischer_14_5 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_5).fischer_anc_14HPF,2)
				embryo_fischer_15_5 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_5).fischer_anc_15HPF,2)
				embryo_fischer_16_5 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_5).fischer_anc_16HPF,2)		
				embryo_fischer_17_5 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_5).fischer_anc_17HPF,2)
				embryo_fischer_18_5 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_5).fischer_anc_18HPF,2)
				embryo_fischer_19_5 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_5).fischer_anc_19HPF,2)
				embryo_fischer_se_0_5 = round(Embryo_SE.objects.get(nvertx_id=nvertx_5).fischer_se_0HPF,2)
				embryo_fischer_se_1_5 = round(Embryo_SE.objects.get(nvertx_id=nvertx_5).fischer_se_1HPF,2)
				embryo_fischer_se_2_5 = round(Embryo_SE.objects.get(nvertx_id=nvertx_5).fischer_se_2HPF,2)
				embryo_fischer_se_3_5 = round(Embryo_SE.objects.get(nvertx_id=nvertx_5).fischer_se_3HPF,2)
				embryo_fischer_se_4_5 = round(Embryo_SE.objects.get(nvertx_id=nvertx_5).fischer_se_4HPF,2)
				embryo_fischer_se_5_5 = round(Embryo_SE.objects.get(nvertx_id=nvertx_5).fischer_se_5HPF,2)
				embryo_fischer_se_6_5 = round(Embryo_SE.objects.get(nvertx_id=nvertx_5).fischer_se_6HPF,2)
				embryo_fischer_se_7_5 = round(Embryo_SE.objects.get(nvertx_id=nvertx_5).fischer_se_7HPF,2)
				embryo_fischer_se_8_5 = round(Embryo_SE.objects.get(nvertx_id=nvertx_5).fischer_se_8HPF,2)
				embryo_fischer_se_9_5 = round(Embryo_SE.objects.get(nvertx_id=nvertx_5).fischer_se_9HPF,2)
				#embryo_fischer_se_10_1 = [round(embryo_fischer_10_1 - Embryo_SE.objects.get(nvertx_id=nvertx_1).fischer_se_10HPF,2),round(embryo_fischer_10_1 + Embryo_SE.objects.get(nvertx_id=nvertx_1).fischer_se_10HPF,2)]
				embryo_fischer_se_11_5 = round(Embryo_SE.objects.get(nvertx_id=nvertx_5).fischer_se_11HPF,2)
				embryo_fischer_se_12_5 = round(Embryo_SE.objects.get(nvertx_id=nvertx_5).fischer_se_12HPF,2)
				embryo_fischer_se_13_5 = round(Embryo_SE.objects.get(nvertx_id=nvertx_5).fischer_se_13HPF,2)
				embryo_fischer_se_14_5 = round(Embryo_SE.objects.get(nvertx_id=nvertx_5).fischer_se_14HPF,2)
				embryo_fischer_se_15_5 = round(Embryo_SE.objects.get(nvertx_id=nvertx_5).fischer_se_15HPF,2)
				embryo_fischer_se_16_5 = round(Embryo_SE.objects.get(nvertx_id=nvertx_5).fischer_se_16HPF,2)
				embryo_fischer_se_17_5 = round(Embryo_SE.objects.get(nvertx_id=nvertx_5).fischer_se_17HPF,2)
				embryo_fischer_se_18_5 = round(Embryo_SE.objects.get(nvertx_id=nvertx_5).fischer_se_18HPF,2)
				embryo_fischer_se_19_5 = round(Embryo_SE.objects.get(nvertx_id=nvertx_5).fischer_se_19HPF,2)
			except :
				nvertx_5_embryo_fischer_invalid = True
			try :
				embryo_helm_2_5 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_5).helm_anc_2HPF,2)
				embryo_helm_7_5 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_5).helm_anc_7HPF,2)
				embryo_helm_12_5 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_5).helm_anc_12HPF,2)
				embryo_helm_24_5 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_5).helm_anc_24HPF,2)
				embryo_helm_120_5 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_5).helm_anc_120HPF,2)
				embryo_helm_240_5 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_5).helm_anc_240HPF,2)
				embryo_helm_se_2_5 = round(Embryo_SE.objects.get(nvertx_id=nvertx_5).helm_se_2HPF,2)
				embryo_helm_se_7_5 = round(Embryo_SE.objects.get(nvertx_id=nvertx_5).helm_se_7HPF,2)
				embryo_helm_se_12_5 = round(Embryo_SE.objects.get(nvertx_id=nvertx_5).helm_se_12HPF,2)
				embryo_helm_se_24_5 = round(Embryo_SE.objects.get(nvertx_id=nvertx_5).helm_se_24HPF,2)
				embryo_helm_se_120_5 = round(Embryo_SE.objects.get(nvertx_id=nvertx_5).helm_se_120HPF,2)
				embryo_helm_se_240_5 = round(Embryo_SE.objects.get(nvertx_id=nvertx_5).helm_se_240HPF,2)
			except :
				nvertx_5_embryo_helm_invalid = True
			if not nvertx_5_embryo_warner_invalid and not nvertx_5_embryo_fischer_invalid and not nvertx_5_embryo_helm_invalid :
				try :
					embryo_mean_2_5 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_5).mean_2HPF,2)
					embryo_mean_7_5 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_5).mean_7HPF,2)
					embryo_mean_12_5 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_5).mean_12HPF,2)
					embryo_mean_24_5 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_5).mean_24HPF,2)
					embryo_mean_120_5 = round(Embryo_cpm.objects.get(nvertx_id=nvertx_5).mean_120HPF,2)
					embryo_mean_0_5 = embryo_fischer_0_5
					embryo_mean_1_5 = embryo_fischer_1_5
					embryo_mean_3_5 = embryo_fischer_3_5
					embryo_mean_4_5 = embryo_fischer_4_5
					embryo_mean_5_5 = embryo_fischer_5_5
					embryo_mean_6_5 = embryo_fischer_6_5
					embryo_mean_8_5 = embryo_fischer_8_5
					embryo_mean_9_5 = embryo_fischer_9_5
					embryo_mean_10_5 = embryo_fischer_10_5
					embryo_mean_11_5 = embryo_fischer_11_5
					embryo_mean_13_5 = embryo_fischer_13_5
					embryo_mean_14_5 = embryo_fischer_14_5
					embryo_mean_15_5 = embryo_fischer_15_5
					embryo_mean_16_5 = embryo_fischer_16_5
					embryo_mean_17_5 = embryo_fischer_17_5
					embryo_mean_18_5 = embryo_fischer_18_5
					embryo_mean_19_5 = embryo_fischer_19_5
					embryo_mean_48_5 = embryo_warner_48_5
					embryo_mean_72_5 = embryo_warner_72_5
					embryo_mean_96_5 = embryo_warner_96_5
					embryo_mean_144_5 = embryo_warner_144_5
					embryo_mean_168_5 = embryo_warner_168_5
					embryo_mean_192_5 = embryo_warner_192_5
					embryo_mean_240_5 = embryo_helm_240_5
				except :
					nvertx_5_embryo_mean_invalid = True
			else :
				nvertx_5_embryo_mean_invalid = True
			try :
				annot_nemve1_tophit_5 = Annotation.objects.get(nvertx_id=nvertx_5).Nemve1_tophit
				annot_nemve1_e_val_5 = round(Annotation.objects.get(nvertx_id=nvertx_5).Nemve1_e_val,150)
				if annot_nemve1_e_val_5 == 99 :
					annot_nemve1_e_val_5 = "NA"
				annot_mfuzz_r_clust_5 = Annotation.objects.get(nvertx_id=nvertx_5).Mfuzz_R_Clust
				annot_mfuzz_r_score_5 = round(Annotation.objects.get(nvertx_id=nvertx_5).Mfuzz_R_Score,2)
				if annot_mfuzz_r_score_5 == -1 :
					annot_mfuzz_r_score_5 = "NA"
				annot_mfuzz_e_clust_5 = Annotation.objects.get(nvertx_id=nvertx_5).Mfuzz_E_Clust
				annot_mfuzz_e_score_5 = round(Annotation.objects.get(nvertx_id=nvertx_5).Mfuzz_E_Score,2)
				if annot_mfuzz_e_score_5 == -1 :
					annot_mfuzz_e_score_5 = "NA"
				annot_uniprot_id_5 = Annotation.objects.get(nvertx_id=nvertx_5).Uniprot_ID
				annot_uniprot_description_5 = Annotation.objects.get(nvertx_id=nvertx_5).Uniprot_Description
				annot_top_nr_hit_eval_5 = Annotation.objects.get(nvertx_id=nvertx_5).Top_nr_hit_eval
				if annot_top_nr_hit_eval_5 != "NA" :
					annot_top_nr_hit_eval_5_split = annot_top_nr_hit_eval_5.split('|',4)
					annot_nr_beg_5 = annot_top_nr_hit_eval_5_split[0] + "|" + annot_top_nr_hit_eval_5_split[1] + "|" + annot_top_nr_hit_eval_5_split[2]
					annot_nr_link_5 = annot_top_nr_hit_eval_5_split[3]
					annot_nr_end_5 = annot_top_nr_hit_eval_5_split[4]
					annot_other_nr_hits_5 = Annotation.objects.get(nvertx_id=nvertx_5).Other_nr_hits
					nr_hit_graph_5 = re.search('[\[\- \w]+\]', annot_top_nr_hit_eval_5).group(0)
			except :
				nvertx_5_annot_invalid = True
			try :
				ncbi_5 = annot_top_nr_hit_eval_5.split('|')[1]
				prot_5 = annot_top_nr_hit_eval_5.split('|')[3]
			except :
				nvertx_5_links_invalid = True
		
		#This chunk cats the multifastas
		multifasta = ''
		if nvertx_1 :
			multifasta = '>' + nvertx_1 + '\n' + sequence_fasta_1
		if nvertx_2 :
			multifasta = multifasta + '\n' + '>' + nvertx_2 + '\n' + sequence_fasta_2
		if nvertx_3 :
			multifasta = multifasta + '\n' + '>' + nvertx_3 + '\n' + sequence_fasta_3
		if nvertx_4 :
			multifasta = multifasta + '\n' + '>' + nvertx_4 + '\n' + sequence_fasta_4
		if nvertx_5 :
			multifasta = multifasta + '\n' + '>' + nvertx_5 + '\n' + sequence_fasta_5

		#this chunk calls the MUSCLE alignment
		#create temp file
    	tempin = tempfile.NamedTemporaryFile()
    	#write the fasta to it
    	tempin.write(multifasta)
    	#create a temporary outfile
    	tempout = tempfile.NamedTemporaryFile()
    	tempin.seek(0)
    	#tempin.read()
    	wd = os.path.dirname(os.path.realpath('muscle3.8.31'))
    	os.system(wd + "/srv/www/htdocs/ER/NvER_plotter_django/nemVec_ER/muscle3.8.31 -in " + tempin.name + " -out " + tempout.name +" -html -quiet")
    	#close and remove file
    	tempout.seek(0)
    	out = tempout.read()
    	tempout.close()
    	tempin.close()

	return render(request, 'ER_plotter/nvertxResults.html', locals())

def mfuzz(request):
	gene_search_form = Gene_searchForm(request.POST or None)
	if gene_search_form.is_valid():
		gene_name = gene_search_form.cleaned_data['gene_name']
		gene_search = True
	nvertx_form = NvERTxForm(request.POST or None)
	if nvertx_form.is_valid():
		nvertx_1 = nvertx_form.cleaned_data['nvertx_1']
		if nvertx_1[0] != 'N' :
			nvertx_1 = 'NvERTx.4.' + nvertx_1
		nvertx_2 = nvertx_form.cleaned_data['nvertx_2']
		if nvertx_2 and nvertx_2[0] != 'N' :
			nvertx_2 = 'NvERTx.4.' + nvertx_2
		nvertx_3 = nvertx_form.cleaned_data['nvertx_3']
		if nvertx_3 and nvertx_3[0] != 'N' :
			nvertx_3 = 'NvERTx.4.' + nvertx_3
		nvertx_4 = nvertx_form.cleaned_data['nvertx_4']
		if nvertx_4 and nvertx_4[0] != 'N' :
			nvertx_4 = 'NvERTx.4.' + nvertx_4
		nvertx_5 = nvertx_form.cleaned_data['nvertx_5']
		if nvertx_5 and nvertx_5[0] != 'N' :
			nvertx_5 = 'NvERTx.4.' + nvertx_5
		log2 = nvertx_form.cleaned_data['log2']
		nvertx_search = True
	clusters_list = Mfuzz.objects.all()
	mfuzz_all = Annotation.objects.all()
	mfuzz1_all = mfuzz_all.filter(mfuzz_clust=1)
	return render(request, 'ER_plotter/mfuzz.html', locals())

### Mfuzz home page
def mfuzzHome(request):
	#Side panel forms
	gene_search_form = Gene_searchForm(request.POST or None)
	nvertx_form = NvERTxForm(request.POST or None)
	convert_form = ConvertForm(request.POST or None)

	#list of clusters. Used to create the buttons
	clusters_list = Mfuzz.objects.all()
	#clusters_list_regen = cluster_list.startswith( 'R' )
	clusters_list_regen = clusters_list.filter(mfuzz_cluster_nb__startswith='R')
	clusters_list_embryo = clusters_list.filter(mfuzz_cluster_nb__startswith='E')

	return render(request, 'ER_plotter/mfuzzHome.html', locals())


### Mfuzz page when a cluster has been selected (button clicked)
def mfuzzResults(request,mfuzz_nb):
	#Side panel forms
	gene_search_form = Gene_searchForm(request.POST or None)
	nvertx_form = NvERTxForm(request.POST or None)
	convert_form = ConvertForm(request.POST or None)

	#list of clusters. Used to create the buttons
	if mfuzz_nb[0] == 'E' :
		clusters_list = Mfuzz.objects.all().filter(mfuzz_cluster_nb__startswith='E')
		clusters_annotation = Annotation.objects.filter(Mfuzz_E_Clust=mfuzz_nb).order_by('-Mfuzz_E_Score')
	else :
		clusters_list = Mfuzz.objects.all().filter(mfuzz_cluster_nb__startswith='R')
		clusters_annotation = Annotation.objects.filter(Mfuzz_R_Clust=mfuzz_nb).order_by('-Mfuzz_R_Score')
	
	#Title and plots for the selected cluster
	#mfuzz_cluster_nb = clusters_list.get(mfuzz_cluster_nb=mfuzz_nb).mfuzz_cluster_nb
	mfuzz_graph = clusters_list.get(mfuzz_cluster_nb=mfuzz_nb).cluster_image
	mfuzz_bp_plot = clusters_list.get(mfuzz_cluster_nb=mfuzz_nb).bp_plot_image

	#Details of the selected cluster
	#cluster_annotation = Annotation.objects.filter(mfuzz_clust=mfuzz_nb)
	for elem in clusters_annotation :
		try :
			split = elem.Top_nr_hit_eval.split('|',2)
			elem.ncbi_wo_link_beg = split[0]
			elem.ncbi_link = split[1]
			elem.ncbi_wo_link_end = split[2]
		except :
			pass

	#This is the original pagination i replaced to make the infinite scroll
	
	#pagination for the details
	#page = request.GET.get('page', 1)
	#paginator = diggPaginator.DiggPaginator(clusters_annotation, 50, body=5)
	#try:
	#	cluster_table = paginator.page(page)
	#except PageNotAnInteger:
	#	cluster_table = paginator.page(1)
	#except EmptyPage:
	#	cluster_table = paginator.page(paginator.num_pages)

	#return render(request, 'ER_plotter/mfuzzResults.html', locals())

	page = request.GET.get('page', 1)
	paginator = Paginator(clusters_annotation, 50)
	try:
		cluster_table = paginator.page(page)
	except PageNotAnInteger:
		cluster_table = paginator.page(1)
	except EmptyPage:
		ncluster_table = paginator.page(paginator.num_pages)
	return render(request, 'ER_plotter/mfuzzResults.html', locals())

def home(request):
	gene_search_form = Gene_searchForm(request.GET or None)
	nvertx_form = NvERTxForm(request.POST or None)
	convert_form = ConvertForm(request.POST or None)
	
	form = BlastForm(initial={'sequence_in_form': '', 'evalue_in_form': EVALUE_BLAST_DEFAULT})

	return render(request, 'ER_plotter/home.html', locals())

def searchResults(request):
	gene_search_form = Gene_searchForm(request.GET or None)
	nvertx_form = NvERTxForm(request.POST or None)
	convert_form = ConvertForm(request.POST or None)

	if gene_search_form.is_valid():
		search_query = gene_search_form.cleaned_data['gene_name']

	search_result_all = Annotation.objects.filter(nvertx_id__icontains=search_query) | Annotation.objects.filter(Nemve1_tophit__icontains=search_query) | Annotation.objects.filter(Nemve1_e_val__icontains=search_query) | Annotation.objects.filter(Uniprot_ID__icontains=search_query) | Annotation.objects.filter(Uniprot_Description__icontains=search_query) | Annotation.objects.filter(Top_nr_hit_eval__icontains=search_query) | Annotation.objects.filter(Other_nr_hits__icontains=search_query)

	for elem in search_result_all :
		try :
			split = elem.Top_nr_hit_eval.split('|',2)
			elem.ncbi_wo_link_beg = split[0]
			elem.ncbi_link = split[1]
			elem.ncbi_wo_link_end = split[2]
		except :
			pass
	
	#pagination for the details
	page = request.GET.get('page', 1)
	paginator = Paginator(search_result_all, 50)
	try:
		search_result = paginator.page(page)
	except PageNotAnInteger:
		search_result = paginator.page(1)
	except EmptyPage:
		search_result = paginator.page(paginator.num_pages)

	return render(request, 'ER_plotter/searchResults.html', locals())

	#pagination for the details
	#page = request.GET.get('page', 1)
	#paginator = diggPaginator.DiggPaginator(search_result_all, 50, body=5)
	#try:
	#	search_result = paginator.page(page)
	#except PageNotAnInteger:
	#	search_result = paginator.page(1)
	#except EmptyPage:
	#	search_result = paginator.page(paginator.num_pages)
#
#	return render(request, 'ER_plotter/searchResults.html', locals())

def volcano(request):
	gene_search_form = Gene_searchForm(request.GET or None)
	nvertx_form = NvERTxForm(request.POST or None)
	convert_form = ConvertForm(request.POST or None)
	
	comparison_form = comparisonForm(request.GET or None)

	comparison = False
	
	if comparison_form.is_valid():
		comparison = comparison_form.cleaned_data['comparison']
	
	if comparison:
		try:
			#build up the comparison variables
			pvalue = comparison + '_pvalue'
			fc = comparison + '_logfc'
			fdr = comparison + '_neglogfdr'
			sig = comparison + '_sig'
			
			#this fetches the correct fields and outputs to tuples
			de = de_table.objects.values_list('nvertx_id','Nemve1_tophit','Uniprot_Description','Top_nr_hit_eval',fc,pvalue,fdr,sig)

			#the significance testing was done in R, the script making_DEtable.R
			# the cutoff was ifelse(fischer_7_8$FDR <= 0.05 & abs(fischer_7_8$logFC) >= 2,'Y','N')
			# this splits the data based on the significance test
			sig_de = de.filter(**{ sig+'__iexact': 'Y'})
			not_sig_de = de.filter(**{ sig+'__iexact': 'N'})
			
			#these two blocks extract the individual entries of the tuples and convert them to utf-8
			sig_de_id = [x[0].encode("utf-8") for x in sig_de]
			sig_de_nemve1 = [x[1] for x in sig_de]
			sig_de_uprot = [x[2].encode("utf-8") for x in sig_de]
			sig_de_nr = [x[3].encode("utf-8") for x in sig_de]
			sig_de_fc = [x[4] for x in sig_de]
			sig_de_pvalue = [x[5] for x in sig_de]
			sig_de_fdr = [x[6] for x in sig_de]
			
			not_sig_de_id = [x[0].encode("utf-8") for x in not_sig_de]
			not_sig_de_nemve1 = [x[1] for x in not_sig_de]
			not_sig_de_uprot = [x[2].encode("utf-8") for x in not_sig_de]
			not_sig_de_nr = [x[3].encode("utf-8") for x in not_sig_de]
			not_sig_de_fc = [x[4] for x in not_sig_de]
			not_sig_de_pvalue = [x[5] for x in not_sig_de]
			not_sig_de_fdr = [x[6] for x in not_sig_de]

		except :
			pass
	return render(request, 'ER_plotter/volcano.html', locals())

def about(request):
	gene_search_form = Gene_searchForm(request.GET or None)
	nvertx_form = NvERTxForm(request.POST or None)
	convert_form = ConvertForm(request.POST or None)
	
	return render(request, 'ER_plotter/about.html', locals())
	
def faq(request):
	gene_search_form = Gene_searchForm(request.GET or None)
	nvertx_form = NvERTxForm(request.POST or None)
	convert_form = ConvertForm(request.POST or None)
	
	return render(request, 'ER_plotter/faq.html', locals())
