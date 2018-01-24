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
from django.http import JsonResponse

'''
views.py, part of the nemvecER/ER_plotter module of NvERTx.
Written in python 2.7

by Jacob Warner & Vincent Guerlais

email questions to warner.jacob at gmail.com

This is the main view for the bulk of the application.
Meant to be used in the larger ER_plotter Django app.
Each page request is defined by a function below.
Note that several points use an infinite scroll module in place of Digg pagination.

Most elements exist in try blocks that on exception set an object invalid to True.
These invalidations are used in the template and elsewhere.
'''

class AbsoluteValue(Transform):
    lookup_name = 'abs'
    function = 'ABS'
    
from django.db.models import IntegerField
IntegerField.register_lookup(AbsoluteValue)

class nvertx_results(object):
	'''
	nvertx_results():
	Class to query the database and return a series of dictionaries with the results
	usage: nvertx_1_results = nvertx_results(nvertx, log2)
	
	Arguments: 
	nvertx = the nvertx id, this is used as the key for querying the db
	log2 = boolean to determine whether or not to convert results to log2
	
	dicionary results:
	invalid = series of boolean operators. If a 'Try' block fails, it sets to True. These are used as tests in the template
	sequence_fasta = nucelotide sequence of the transcript. Not in fasta, just raw sequence.
	annotation = holds several fields that contain the annotation objects
	counts = holds several fields for the various datsets
	links = holds the links that are built from annotation	
	'''
	def __init__(self, nvertx, log2):	
		self.invalid = {}
		self.sequence_fasta = {}
		self.annotation = {}
		self.counts = {}
		self.links = {}
			
		self.invalid['regen'] = False
		self.invalid['embryo_warner'] = False
		self.invalid['embryo_fischer'] = False
		self.invalid['embryo_helm'] = False
		self.invalid['embryo_mean'] = False
		self.invalid['links'] =  False
		self.invalid['annotation'] = False
		
		##
		## Sequence Block
		##
		
		try :
			sequence_fasta_query = Fasta.objects.get(nvertx_id=nvertx).fasta_sequence
		except :
			sequence_fasta_query = 'No fasta sequence'
		
		self.sequence_fasta['sequence_fasta'] = sequence_fasta_query
		
		##
		## counts block
		##
				
		try :
			regen = []
			regen_se = []	
							
			#extract regen counts and standard errors and convert to log as needed
			if not log2:
				regen_query = Regen_cpm.objects.all().filter(nvertx_id=nvertx).values_list()
				for timepoint in range(1,17,1):
					regen += [round(math.log(x[timepoint]+1,2),2) for x in regen_query ]
					
				regen_se_query = Regen_log_SE.objects.all().filter(nvertx_id=nvertx).values_list()
				for timepoint in range(1,17,1):
					regen_se += [x[timepoint] for x in regen_se_query ]

			else:
				regen_query = Regen_cpm.objects.all().filter(nvertx_id=nvertx).values_list()
				for timepoint in range(1,17,1):
					regen += [round(x[timepoint],2) for x in regen_query ]
					
				regen_se_query = Regen_SE.objects.all().filter(nvertx_id=nvertx).values_list()
				for timepoint in range(1,17,1):
					regen_se += [x[timepoint] for x in regen_se_query ]

			self.counts['regen'] = regen
			self.counts['regen_se'] = regen_se
						
		except :
			self.invalid['regen'] =  True	
			
		try :
			embryo_warner = []
			embryo_warner_se = []
			embryo_fischer = []
			embryo_fischer_se = []
			embryo_helm = []
			embryo_helm_se = []
			embryo_mean = []
			
			#get all embryo counts
			embryo_query = Embryo_cpm.objects.all().filter(nvertx_id=nvertx).values_list()
			embryo_se_query = Embryo_SE.objects.all().filter(nvertx_id=nvertx).values_list()
			
			#			
			# warner block
			#
			for timepoint in range(1,9,1):
				embryo_warner += [round(x[timepoint],2) for x in embryo_query]
				embryo_warner_se += [x[timepoint] for x in embryo_se_query]

			self.counts['embryo_warner'] = embryo_warner
			self.counts['embryo_warner_se'] = embryo_warner_se
			
			#
			# fischer block
			#
			for timepoint in range(9,29,1):
				embryo_fischer += [round(x[timepoint],2) for x in embryo_query]
				embryo_fischer_se += [x[timepoint] for x in embryo_se_query]
				
			#drop the none that appears at fischer se t=10hpf
			embryo_fischer_se = [0 if x is None else x for x in embryo_fischer_se]
			
			self.counts['embryo_fischer'] = embryo_fischer
			self.counts['embryo_fischer_se'] = embryo_fischer_se
			
			#
			# helm block
			#
			for timepoint in range(29,35,1):
				embryo_helm += [round(x[timepoint],2) for x in embryo_query]
				embryo_helm_se += [x[timepoint] for x in embryo_se_query]

			self.counts['embryo_helm'] = embryo_helm
			self.counts['embryo_helm_se'] = embryo_helm_se
											
		except :
			self.invalid['embryo_warner'] = True
			self.invalid['embryo_fischer'] = True
			self.invalid['embryo_helm'] = True
				
		#
		# averages block
		#
		if not self.invalid['embryo_warner'] and not self.invalid['embryo_fischer'] and not self.invalid['embryo_helm'] :
			try :
				# Don't need to round this since it's not reported in the table
				for timepoint in [9,10,35,12,13,14,15,36,17,18,19,20,37,22,23,24,25,26,27,28,38,2,3,4,39,6,7,8,34]:
					embryo_mean += [round(x[timepoint],2) for x in embryo_query ]
		
				self.counts['embryo_mean'] = embryo_mean

			except :
				self.invalid['embryo_mean'] = False
		else :
			self.invalid['embryo_mean'] = False
		
		##
		## Annotation block
		##
		try :
			self.annotation['nemve1_tophit'] = Annotation.objects.get(nvertx_id=nvertx).Nemve1_tophit
			
			nemve1_e_val = round(Annotation.objects.get(nvertx_id=nvertx).Nemve1_e_val,150)
			if nemve1_e_val == 99 :
				nemve1_e_val = "NA"
			self.annotation['nemve1_e_val'] = nemve1_e_val
			
			self.annotation['mfuzz_r_clust'] = Annotation.objects.get(nvertx_id=nvertx).Mfuzz_R_Clust
			
			mfuzz_r_score = round(Annotation.objects.get(nvertx_id=nvertx).Mfuzz_R_Score,2)
			if mfuzz_r_score == -1 :
				mfuzz_r_score = "NA"
			self.annotation['mfuzz_r_score'] = mfuzz_r_score
			
			self.annotation['mfuzz_e_clust'] = Annotation.objects.get(nvertx_id=nvertx).Mfuzz_E_Clust
			
			mfuzz_e_score = round(Annotation.objects.get(nvertx_id=nvertx).Mfuzz_E_Score,2)
			if mfuzz_e_score == -1 :
				mfuzz_e_score = "NA"
			self.annotation['mfuzz_e_score'] = mfuzz_e_score
			
			self.annotation['uniprot_id'] = Annotation.objects.get(nvertx_id=nvertx).Uniprot_ID
			
			self.annotation['uniprot_description'] = Annotation.objects.get(nvertx_id=nvertx).Uniprot_Description
			
			top_nr_hit_eval = Annotation.objects.get(nvertx_id=nvertx).Top_nr_hit_eval
			self.annotation['top_nr_hit_eval'] = top_nr_hit_eval
			
			if top_nr_hit_eval != "NA" :
				top_nr_hit_eval_split = top_nr_hit_eval.split('|',4)
				self.annotation['nr_beg'] = top_nr_hit_eval_split[0] + "|" + top_nr_hit_eval_split[1] + "|" + top_nr_hit_eval_split[2]
				self.annotation['nr_link'] = top_nr_hit_eval_split[3]
				self.annotation['nr_end'] = top_nr_hit_eval_split[4]
				self.annotation['other_nr_hits'] = Annotation.objects.get(nvertx_id=nvertx).Other_nr_hits
				self.annotation['nr_hit_graph'] = re.search('[\[\- \w]+\]', top_nr_hit_eval).group(0)
			
			else :
				self.invalid['links'] = True				
		except :
			self.invalid['annotation'] = True
		try :
			self.links['ncbi'] = top_nr_hit_eval.split('|')[1]
			self.links['prot'] = top_nr_hit_eval.split('|')[3]
		except :
			self.invalid['links'] = True
	
def results(request):
	gene_search_form = Gene_searchForm(request.POST or None)
	nvertx_form = NvERTxForm(request.POST or None)
	convert_form = ConvertForm(request.POST or None)
	
	nvertx_1 = False
	nvertx_2 = False
	nvertx_3 = False
	nvertx_4 = False
	nvertx_5 = False
	log2 = False

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
	
	if nvertx_1:
		nvertx_1_results = nvertx_results(nvertx_1,log2)
	if nvertx_2:
		nvertx_2_results = nvertx_results(nvertx_2,log2)
	if nvertx_3:
		nvertx_3_results = nvertx_results(nvertx_3,log2)
	if nvertx_4:
		nvertx_4_results = nvertx_results(nvertx_4,log2)
	if nvertx_5:
		nvertx_5_results = nvertx_results(nvertx_5,log2)

	#This chunk cats the multifastas
	multifasta = ''
	if nvertx_1 :
		multifasta = '>' + nvertx_1 + '\n' + nvertx_1_results.sequence_fasta['sequence_fasta']
	if nvertx_2 :
		multifasta += '\n' + '>' + nvertx_2 + '\n' + nvertx_2_results.sequence_fasta['sequence_fasta']
	if nvertx_3 :
		multifasta += '\n' + '>' + nvertx_3 + '\n' + nvertx_3_results.sequence_fasta['sequence_fasta']
	if nvertx_4 :
		multifasta += '\n' + '>' + nvertx_4 + '\n' + nvertx_4_results.sequence_fasta['sequence_fasta']
	if nvertx_5 :
		multifasta += '\n' + '>' + nvertx_5 + '\n' + nvertx_5_results.sequence_fasta['sequence_fasta']

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
	mfuzz_graph = clusters_list.get(mfuzz_cluster_nb=mfuzz_nb).cluster_image
	mfuzz_bp_plot = clusters_list.get(mfuzz_cluster_nb=mfuzz_nb).bp_plot_image

	#Details of the selected cluster
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

			# the significance testing was done in R, the script making_DEtable.R
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
	
def api(request):
	#get the query tag from the url
	query = request.META['QUERY_STRING']
	api_out = {}
	if not query == '':
		if query[0] != 'N' :
			query = 'NvERTx.4.' + query
		if not re.match('NvERTx.4.[0-9]+$', query):
			query = 'NA'
		log2 = True
		
		#query the db
		api_results = nvertx_results(query,log2)
		
		#build up the dictionaries for the response
		#these will be converted to json
		api_out['nvertx_id'] = query
		api_out['regen']= [['counts',api_results.counts['regen']],['error',api_results.counts['regen_se']],['time.hpa',[-1,0,2,4,8,12,16,20,24,36,48,60,72,96,120,144]]]
		api_out['embryo_warner']= [['counts',api_results.counts['embryo_warner']],['error',api_results.counts['embryo_warner_se']],['time.hpf',[24,48,72,96,120,144,168,192]]]
		api_out['embryo_helm']= [['counts',api_results.counts['embryo_helm']],['error',api_results.counts['embryo_helm_se']],['time.hpf',[2,7,12,24,120,240]]]	
		api_out['embryo_fischer']= [['counts',api_results.counts['embryo_fischer']],['error',api_results.counts['embryo_fischer_se']],['time.hpf',[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]]]

		#This block enables COR on the response
		response = JsonResponse(api_out)
		response["Access-Control-Allow-Origin"] = "*"
		response["Access-Control-Allow-Methods"] = "GET, OPTIONS"
		response["Access-Control-Max-Age"] = "1000"
		response["Access-Control-Allow-Headers"] = "X-Requested-With, Content-Type"
		return response
	
	#error if there is no query	
	else :
		api_out['sorry'] = 'no query... add: ?NvERTx.4.100038 or another number to the url query api'
		
		response = JsonResponse(api_out)
		response["Access-Control-Allow-Origin"] = "*"
		response["Access-Control-Allow-Methods"] = "GET, OPTIONS"
		response["Access-Control-Max-Age"] = "1000"
		response["Access-Control-Allow-Headers"] = "X-Requested-With, Content-Type"
		return response
