#DBfill allows to fill the NvERTx Database ;
#Use it on files with columns separated with tabulations (export table from R for example)
#how to use it : 
#	1. run the shell via python (in the /dER_plotter/nemVec_ER dir) : python manage.py shell
#	2. import this script in python : import DBfill
#	3. execute with the name of the file and the name of the table : DBfill.DBFill('fileName','table')
#	the tables are : Regen_cpm ; Fasta ; Embryo_cpm ; Annotation ; Regen_SE ; Embryo_SE ; you can check for their name in the admin panel
#	4. It may take a while depending on the size of your file 'done filling the DB' should print when done
#The tables Annotation and Embryo_cpm can handle missing data. they should be replaced by 'NA'

##################
### Var
##################

import os
from ER_plotter.models import Fasta, Regen_cpm, Embryo_cpm, Annotation, Regen_SE, Embryo_SE, Regen_log_SE, de_table

##################
### Functions
##################

#Split the data on the tabs
def dataSplit(dataLine) :
	return dataLine[:-1].split('\t')

#Create the query to fill the database
def queryCreate(splitLine, DBTableName):
	if DBTableName == 'Regen_cpm' :
		query = "Regen_cpm(nvertx_id = '" + splitLine[0] +"'"
		if splitLine[1] != 'NA' :
			query += ", regen_anc_UC = " + splitLine[1]
		if splitLine[2] != 'NA' :
			query += ", regen_anc_0HPA = " + splitLine[2]
		if splitLine[3] != 'NA' :
			query += ", regen_anc_2HPA = " + splitLine[3]
		if splitLine[4] != 'NA' :
			query += ", regen_anc_4HPA = " + splitLine[4]
		if splitLine[5] != 'NA' :
			query += ", regen_anc_8HPA = " + splitLine[5]
		if splitLine[6] != 'NA' :
			query += ", regen_anc_12HPA = " + splitLine[6]
		if splitLine[7] != 'NA' :
			query += ", regen_anc_16HPA = " + splitLine[7]
		if splitLine[8] != 'NA' :
			query += ", regen_anc_20HPA = " + splitLine[8]
		if splitLine[9] != 'NA' :
			query += ", regen_anc_24HPA = " + splitLine[9]
		if splitLine[10] != 'NA' :
			query += ", regen_anc_36HPA = " + splitLine[10]
		if splitLine[11] != 'NA' :
			query += ", regen_anc_48HPA = " + splitLine[11]
		if splitLine[12] != 'NA' :
			query += ", regen_anc_60HPA = " + splitLine[12]
		if splitLine[13] != 'NA' :
			query += ", regen_anc_72HPA = " + splitLine[13]
		if splitLine[14] != 'NA' :
			query += ", regen_anc_96HPA = " + splitLine[14]
		if splitLine[15] != 'NA' :
			query += ", regen_anc_120HPA = " + splitLine[15]
		if splitLine[16] != 'NA' :
			query += ", regen_anc_144HPA = " + splitLine[16]
		query += ").save()"

	elif DBTableName == 'Fasta' :
		query = "%s(nvertx_id = '%s', fasta_sequence= '%s').save()" % (DBTableName,splitLine[0],splitLine[1])
		#if splitLine[0] == '>' :

	elif DBTableName == 'Embryo_cpm' :
		query = "Embryo_cpm(nvertx_id = '" + splitLine[0] +"'"
		if splitLine[1] != 'NA' :
			query += ", warner_anc_24HPF = " + splitLine[1]
		if splitLine[2] != 'NA' :
			query += ", warner_anc_48HPF = " + splitLine[2]
		if splitLine[3] != 'NA' :
			query += ", warner_anc_72HPF = " + splitLine[3]
		if splitLine[4] != 'NA' :
			query += ", warner_anc_96HPF = " + splitLine[4]
		if splitLine[5] != 'NA' :
			query += ", warner_anc_120HPF = " + splitLine[5]
		if splitLine[6] != 'NA' :
			query += ", warner_anc_144HPF = " + splitLine[6]
		if splitLine[7] != 'NA' :
			query += ", warner_anc_168HPF = " + splitLine[7]
		if splitLine[8] != 'NA' :
			query += ", warner_anc_192HPF = " + splitLine[8]
		if splitLine[9] != 'NA' :
			query += ", fischer_anc_0HPF = " + splitLine[9]
		if splitLine[10] != 'NA' :
			query += ", fischer_anc_1HPF = " + splitLine[10]
		if splitLine[11] != 'NA' :
			query += ", fischer_anc_2HPF = " + splitLine[11]
		if splitLine[12] != 'NA' :
			query += ", fischer_anc_3HPF = " + splitLine[12]
		if splitLine[13] != 'NA' :
			query += ", fischer_anc_4HPF = " + splitLine[13]
		if splitLine[14] != 'NA' :
			query += ", fischer_anc_5HPF = " + splitLine[14]
		if splitLine[15] != 'NA' :
			query += ", fischer_anc_6HPF = " + splitLine[15]
		if splitLine[16] != 'NA' :
			query += ", fischer_anc_7HPF = " + splitLine[16]
		if splitLine[17] != 'NA' :
			query += ", fischer_anc_8HPF = " + splitLine[17]
		if splitLine[18] != 'NA' :
			query += ", fischer_anc_9HPF = " + splitLine[18]
		if splitLine[19] != 'NA' :
			query += ", fischer_anc_10HPF = " + splitLine[19]
		if splitLine[20] != 'NA' :
			query += ", fischer_anc_11HPF = " + splitLine[20]
		if splitLine[21] != 'NA' :
			query += ", fischer_anc_12HPF = " + splitLine[21]
		if splitLine[22] != 'NA' :
			query += ", fischer_anc_13HPF = " + splitLine[22]
		if splitLine[23] != 'NA' :
			query += ", fischer_anc_14HPF = " + splitLine[23]
		if splitLine[24] != 'NA' :
			query += ", fischer_anc_15HPF = " + splitLine[24]
		if splitLine[25] != 'NA' :
			query += ", fischer_anc_16HPF = " + splitLine[25]
		if splitLine[26] != 'NA' :
			query += ", fischer_anc_17HPF = " + splitLine[26]
		if splitLine[27] != 'NA' :
			query += ", fischer_anc_18HPF = " + splitLine[27]
		if splitLine[28] != 'NA' :
			query += ", fischer_anc_19HPF = " + splitLine[28]
		if splitLine[29] != 'NA' :
			query += ", helm_anc_2HPF = " + splitLine[29]
		if splitLine[30] != 'NA' :
			query += ", helm_anc_7HPF = " + splitLine[30]
		if splitLine[31] != 'NA' :
			query += ", helm_anc_12HPF = " + splitLine[31]
		if splitLine[32] != 'NA' :
			query += ", helm_anc_24HPF = " + splitLine[32]
		if splitLine[33] != 'NA' :
			query += ", helm_anc_120HPF = " + splitLine[33]
		if splitLine[34] != 'NA' :
			query += ", helm_anc_240HPF = " + splitLine[34]
		if splitLine[35] != 'NA' :
			query += ", mean_2HPF = " + splitLine[35]
		if splitLine[36] != 'NA' :
			query += ", mean_7HPF = " + splitLine[36]
		if splitLine[37] != 'NA' :
			query += ", mean_12HPF = " + splitLine[37]
		if splitLine[38] != 'NA' :
			query += ", mean_24HPF = " + splitLine[38]
		if splitLine[39] != 'NA' :
			query += ", mean_120HPF = " + splitLine[39]
		query += ").save()"

	elif DBTableName == "Annotation" :
		query = 'Annotation(nvertx_id = "' + splitLine[0] +'"'
		query += ', Nemve1_tophit = "' + splitLine[2] +'"'
		if splitLine[3] != "NA" :
			query += ', Nemve1_e_val = ' + splitLine[3]
		else :
			query += ', Nemve1_e_val = 99'
		query += ', Mfuzz_R_Clust = "' + splitLine[4] +'"'
		if splitLine[5] != "NA" :
			query += ', Mfuzz_R_Score = ' + splitLine[5]
		else :
			query += ', Mfuzz_R_Score = -1'
		query += ', Mfuzz_E_Clust = "' + splitLine[6] +'"'
		if splitLine[7] != "NA" :
			query += ', Mfuzz_E_Score = ' + splitLine[7]
		else :
			query += ', Mfuzz_E_Score = -1'
		query += ', Uniprot_ID = "' + splitLine[8] +'"'
		query += ', Uniprot_Description = "' + splitLine[9] +'"'
		query += ', Top_nr_hit_eval = "' + splitLine[10] +'"'
		query += ', Other_nr_hits = "' + splitLine[11] +'"'
		query += ').save()'
	
	elif DBTableName == 'Regen_SE' :
		query = "Regen_SE(nvertx_id = '" + splitLine[0] +"'"
		if splitLine[1] != 'NA' :
			query += ", regen_se_UC = " + splitLine[1]
		if splitLine[2] != 'NA' :
			query += ", regen_se_0HPA = " + splitLine[2]
		if splitLine[3] != 'NA' :
			query += ", regen_se_2HPA = " + splitLine[3]
		if splitLine[4] != 'NA' :
			query += ", regen_se_4HPA = " + splitLine[4]
		if splitLine[5] != 'NA' :
			query += ", regen_se_8HPA = " + splitLine[5]
		if splitLine[6] != 'NA' :
			query += ", regen_se_12HPA = " + splitLine[6]
		if splitLine[7] != 'NA' :
			query += ", regen_se_16HPA = " + splitLine[7]
		if splitLine[8] != 'NA' :
			query += ", regen_se_20HPA = " + splitLine[8]
		if splitLine[9] != 'NA' :
			query += ", regen_se_24HPA = " + splitLine[9]
		if splitLine[10] != 'NA' :
			query += ", regen_se_36HPA = " + splitLine[10]
		if splitLine[11] != 'NA' :
			query += ", regen_se_48HPA = " + splitLine[11]
		if splitLine[12] != 'NA' :
			query += ", regen_se_60HPA = " + splitLine[12]
		if splitLine[13] != 'NA' :
			query += ", regen_se_72HPA = " + splitLine[13]
		if splitLine[14] != 'NA' :
			query += ", regen_se_96HPA = " + splitLine[14]
		if splitLine[15] != 'NA' :
			query += ", regen_se_120HPA = " + splitLine[15]
		if splitLine[16] != 'NA' :
			query += ", regen_se_144HPA = " + splitLine[16]
		query += ").save()"
	
	elif DBTableName == 'Regen_log_SE' :
		query = "Regen_log_SE(nvertx_id = '" + splitLine[0] +"'"
		if splitLine[1] != 'NA' :
			query += ", regen_log_se_UC = " + splitLine[1]
		if splitLine[2] != 'NA' :
			query += ", regen_log_se_0HPA = " + splitLine[2]
		if splitLine[3] != 'NA' :
			query += ", regen_log_se_2HPA = " + splitLine[3]
		if splitLine[4] != 'NA' :
			query += ", regen_log_se_4HPA = " + splitLine[4]
		if splitLine[5] != 'NA' :
			query += ", regen_log_se_8HPA = " + splitLine[5]
		if splitLine[6] != 'NA' :
			query += ", regen_log_se_12HPA = " + splitLine[6]
		if splitLine[7] != 'NA' :
			query += ", regen_log_se_16HPA = " + splitLine[7]
		if splitLine[8] != 'NA' :
			query += ", regen_log_se_20HPA = " + splitLine[8]
		if splitLine[9] != 'NA' :
			query += ", regen_log_se_24HPA = " + splitLine[9]
		if splitLine[10] != 'NA' :
			query += ", regen_log_se_36HPA = " + splitLine[10]
		if splitLine[11] != 'NA' :
			query += ", regen_log_se_48HPA = " + splitLine[11]
		if splitLine[12] != 'NA' :
			query += ", regen_log_se_60HPA = " + splitLine[12]
		if splitLine[13] != 'NA' :
			query += ", regen_log_se_72HPA = " + splitLine[13]
		if splitLine[14] != 'NA' :
			query += ", regen_log_se_96HPA = " + splitLine[14]
		if splitLine[15] != 'NA' :
			query += ", regen_log_se_120HPA = " + splitLine[15]
		if splitLine[16] != 'NA' :
			query += ", regen_log_se_144HPA = " + splitLine[16]
		query += ").save()"
	
	elif DBTableName == 'Embryo_SE' :
		query = "Embryo_SE(nvertx_id = '" + splitLine[0] +"'"
		if splitLine[1] != 'NA' :
			query += ", warner_se_24HPF = " + splitLine[1]
		if splitLine[2] != 'NA' :
			query += ", warner_se_48HPF = " + splitLine[2]
		if splitLine[3] != 'NA' :
			query += ", warner_se_72HPF = " + splitLine[3]
		if splitLine[4] != 'NA' :
			query += ", warner_se_96HPF = " + splitLine[4]
		if splitLine[5] != 'NA' :
			query += ", warner_se_120HPF = " + splitLine[5]
		if splitLine[6] != 'NA' :
			query += ", warner_se_144HPF = " + splitLine[6]
		if splitLine[7] != 'NA' :
			query += ", warner_se_168HPF = " + splitLine[7]
		if splitLine[8] != 'NA' :
			query += ", warner_se_192HPF = " + splitLine[8]
		if splitLine[9] != 'NA' :
			query += ", fischer_se_0HPF = " + splitLine[9]
		if splitLine[10] != 'NA' :
			query += ", fischer_se_1HPF = " + splitLine[10]
		if splitLine[11] != 'NA' :
			query += ", fischer_se_2HPF = " + splitLine[11]
		if splitLine[12] != 'NA' :
			query += ", fischer_se_3HPF = " + splitLine[12]
		if splitLine[13] != 'NA' :
			query += ", fischer_se_4HPF = " + splitLine[13]
		if splitLine[14] != 'NA' :
			query += ", fischer_se_5HPF = " + splitLine[14]
		if splitLine[15] != 'NA' :
			query += ", fischer_se_6HPF = " + splitLine[15]
		if splitLine[16] != 'NA' :
			query += ", fischer_se_7HPF = " + splitLine[16]
		if splitLine[17] != 'NA' :
			query += ", fischer_se_8HPF = " + splitLine[17]
		if splitLine[18] != 'NA' :
			query += ", fischer_se_9HPF = " + splitLine[18]
		if splitLine[19] != 'NA' :
			query += ", fischer_se_10HPF = " + splitLine[19]
		if splitLine[20] != 'NA' :
			query += ", fischer_se_11HPF = " + splitLine[20]
		if splitLine[21] != 'NA' :
			query += ", fischer_se_12HPF = " + splitLine[21]
		if splitLine[22] != 'NA' :
			query += ", fischer_se_13HPF = " + splitLine[22]
		if splitLine[23] != 'NA' :
			query += ", fischer_se_14HPF = " + splitLine[23]
		if splitLine[24] != 'NA' :
			query += ", fischer_se_15HPF = " + splitLine[24]
		if splitLine[25] != 'NA' :
			query += ", fischer_se_16HPF = " + splitLine[25]
		if splitLine[26] != 'NA' :
			query += ", fischer_se_17HPF = " + splitLine[26]
		if splitLine[27] != 'NA' :
			query += ", fischer_se_18HPF = " + splitLine[27]
		if splitLine[28] != 'NA' :
			query += ", fischer_se_19HPF = " + splitLine[28]
		if splitLine[29] != 'NA' :
			query += ", helm_se_2HPF = " + splitLine[29]
		if splitLine[30] != 'NA' :
			query += ", helm_se_7HPF = " + splitLine[30]
		if splitLine[31] != 'NA' :
			query += ", helm_se_12HPF = " + splitLine[31]
		if splitLine[32] != 'NA' :
			query += ", helm_se_24HPF = " + splitLine[32]
		if splitLine[33] != 'NA' :
			query += ", helm_se_120HPF = " + splitLine[33]
		if splitLine[34] != 'NA' :
			query += ", helm_se_240HPF = " + splitLine[34]
		query += ").save()"

	elif DBTableName == 'de_table' :
		query = "de_table(nvertx_id = '" + splitLine[1] +"'"
		query += ', Nemve1_tophit = "' + splitLine[0] +'"'
		query += ', Uniprot_Description = "' + splitLine[2] +'"'
		query += ', Top_nr_hit_eval = "' + splitLine[3] +'"'
		if splitLine[4] != 'NA' :
			query += ", regen_0_2_logfc = " + splitLine[4]
		if splitLine[5] != 'NA' :
			query += ", regen_0_2_logcpm = " + splitLine[5]
		if splitLine[6] != 'NA' :
			query += ", regen_0_2_pvalue = " + splitLine[6]
		if splitLine[7] != 'NA' :
			query += ", regen_0_2_neglogfdr = " + splitLine[7]
		if splitLine[8] != 'NA' :
			query += ", regen_0_2_sig = '" + splitLine[8]+"'"
		if splitLine[9] != 'NA' :
			query += ", regen_0_4_logfc = " + splitLine[9]
		if splitLine[10] != 'NA' :
			query += ", regen_0_4_logcpm = " + splitLine[10]
		if splitLine[11] != 'NA' :
			query += ", regen_0_4_pvalue = " + splitLine[11]
		if splitLine[12] != 'NA' :
			query += ", regen_0_4_neglogfdr = " + splitLine[12]
		if splitLine[13] != 'NA' :
			query += ", regen_0_4_sig = '" + splitLine[13]+"'"
		if splitLine[14] != 'NA' :
			query += ", regen_0_8_logfc = " + splitLine[14]
		if splitLine[15] != 'NA' :
			query += ", regen_0_8_logcpm = " + splitLine[15]
		if splitLine[16] != 'NA' :
			query += ", regen_0_8_pvalue = " + splitLine[16]
		if splitLine[17] != 'NA' :
			query += ", regen_0_8_neglogfdr = " + splitLine[17]
		if splitLine[18] != 'NA' :
			query += ", regen_0_8_sig = '" + splitLine[18]+"'"
		if splitLine[19] != 'NA' :
			query += ", regen_0_12_logfc = " + splitLine[19]
		if splitLine[20] != 'NA' :
			query += ", regen_0_12_logcpm = " + splitLine[20]
		if splitLine[21] != 'NA' :
			query += ", regen_0_12_pvalue = " + splitLine[21]
		if splitLine[22] != 'NA' :
			query += ", regen_0_12_neglogfdr = " + splitLine[22]
		if splitLine[23] != 'NA' :
			query += ", regen_0_12_sig = '" + splitLine[23]+"'"
		if splitLine[24] != 'NA' :
			query += ", regen_0_16_logfc = " + splitLine[24]
		if splitLine[25] != 'NA' :
			query += ", regen_0_16_logcpm = " + splitLine[25]
		if splitLine[26] != 'NA' :
			query += ", regen_0_16_pvalue = " + splitLine[26]
		if splitLine[27] != 'NA' :
			query += ", regen_0_16_neglogfdr = " + splitLine[27]
		if splitLine[28] != 'NA' :
			query += ", regen_0_16_sig = '" + splitLine[28]+"'"
		if splitLine[29] != 'NA' :
			query += ", regen_0_20_logfc = " + splitLine[29]
		if splitLine[30] != 'NA' :
			query += ", regen_0_20_logcpm = " + splitLine[30]
		if splitLine[31] != 'NA' :
			query += ", regen_0_20_pvalue = " + splitLine[31]
		if splitLine[32] != 'NA' :
			query += ", regen_0_20_neglogfdr = " + splitLine[32]
		if splitLine[33] != 'NA' :
			query += ", regen_0_20_sig = '" + splitLine[33]+"'"
		if splitLine[34] != 'NA' :
			query += ", regen_0_24_logfc = " + splitLine[34]
		if splitLine[35] != 'NA' :
			query += ", regen_0_24_logcpm = " + splitLine[35]
		if splitLine[36] != 'NA' :
			query += ", regen_0_24_pvalue = " + splitLine[36]
		if splitLine[37] != 'NA' :
			query += ", regen_0_24_neglogfdr = " + splitLine[37]
		if splitLine[38] != 'NA' :
			query += ", regen_0_24_sig = '" + splitLine[38]+"'"
		if splitLine[39] != 'NA' :
			query += ", regen_0_36_logfc = " + splitLine[39]
		if splitLine[40] != 'NA' :
			query += ", regen_0_36_logcpm = " + splitLine[40]
		if splitLine[41] != 'NA' :
			query += ", regen_0_36_pvalue = " + splitLine[41]
		if splitLine[42] != 'NA' :
			query += ", regen_0_36_neglogfdr = " + splitLine[42]
		if splitLine[43] != 'NA' :
			query += ", regen_0_36_sig = '" + splitLine[43]+"'"
		if splitLine[44] != 'NA' :
			query += ", regen_0_48_logfc = " + splitLine[44]
		if splitLine[45] != 'NA' :
			query += ", regen_0_48_logcpm = " + splitLine[45]
		if splitLine[46] != 'NA' :
			query += ", regen_0_48_pvalue = " + splitLine[46]
		if splitLine[47] != 'NA' :
			query += ", regen_0_48_neglogfdr = " + splitLine[47]
		if splitLine[48] != 'NA' :
			query += ", regen_0_48_sig = '" + splitLine[48]+"'"
		if splitLine[49] != 'NA' :
			query += ", regen_0_60_logfc = " + splitLine[49]
		if splitLine[50] != 'NA' :
			query += ", regen_0_60_logcpm = " + splitLine[50]
		if splitLine[51] != 'NA' :
			query += ", regen_0_60_pvalue = " + splitLine[51]
		if splitLine[52] != 'NA' :
			query += ", regen_0_60_neglogfdr = " + splitLine[52]
		if splitLine[53] != 'NA' :
			query += ", regen_0_60_sig = '" + splitLine[53]+"'"
		if splitLine[54] != 'NA' :
			query += ", regen_0_72_logfc = " + splitLine[54]
		if splitLine[55] != 'NA' :
			query += ", regen_0_72_logcpm = " + splitLine[55]
		if splitLine[56] != 'NA' :
			query += ", regen_0_72_pvalue = " + splitLine[56]
		if splitLine[57] != 'NA' :
			query += ", regen_0_72_neglogfdr = " + splitLine[57]
		if splitLine[58] != 'NA' :
			query += ", regen_0_72_sig = '" + splitLine[58]+"'"
		if splitLine[59] != 'NA' :
			query += ", regen_0_96_logfc = " + splitLine[59]
		if splitLine[60] != 'NA' :
			query += ", regen_0_96_logcpm = " + splitLine[60]
		if splitLine[61] != 'NA' :
			query += ", regen_0_96_pvalue = " + splitLine[61]
		if splitLine[62] != 'NA' :
			query += ", regen_0_96_neglogfdr = " + splitLine[62]
		if splitLine[63] != 'NA' :
			query += ", regen_0_96_sig = '" + splitLine[63]+"'"
		if splitLine[64] != 'NA' :
			query += ", regen_0_120_logfc = " + splitLine[64]
		if splitLine[65] != 'NA' :
			query += ", regen_0_120_logcpm = " + splitLine[65]
		if splitLine[66] != 'NA' :
			query += ", regen_0_120_pvalue = " + splitLine[66]
		if splitLine[67] != 'NA' :
			query += ", regen_0_120_neglogfdr = " + splitLine[67]
		if splitLine[68] != 'NA' :
			query += ", regen_0_120_sig = '" + splitLine[68]+"'"
		if splitLine[69] != 'NA' :
			query += ", regen_0_144_logfc = " + splitLine[69]
		if splitLine[70] != 'NA' :
			query += ", regen_0_144_logcpm = " + splitLine[70]
		if splitLine[71] != 'NA' :
			query += ", regen_0_144_pvalue = " + splitLine[71]
		if splitLine[72] != 'NA' :
			query += ", regen_0_144_neglogfdr = " + splitLine[72]
		if splitLine[73] != 'NA' :
			query += ", regen_0_144_sig = '" + splitLine[73]+"'"
		if splitLine[74] != 'NA' :
			query += ", fischer_7_8_logfc = " + splitLine[74]
		if splitLine[75] != 'NA' :
			query += ", fischer_7_8_logcpm = " + splitLine[75]
		if splitLine[76] != 'NA' :
			query += ", fischer_7_8_pvalue = " + splitLine[76]
		if splitLine[77] != 'NA' :
			query += ", fischer_7_8_neglogfdr = " + splitLine[77]
		if splitLine[78] != 'NA' :
			query += ", fischer_7_8_sig = '" + splitLine[78]+"'"
		if splitLine[79] != 'NA' :
			query += ", fischer_7_9_logfc = " + splitLine[79]
		if splitLine[80] != 'NA' :
			query += ", fischer_7_9_logcpm = " + splitLine[80]
		if splitLine[81] != 'NA' :
			query += ", fischer_7_9_pvalue = " + splitLine[81]
		if splitLine[82] != 'NA' :
			query += ", fischer_7_9_neglogfdr = " + splitLine[82]
		if splitLine[83] != 'NA' :
			query += ", fischer_7_9_sig = '" + splitLine[83]+"'"
		if splitLine[84] != 'NA' :
			query += ", fischer_7_10_logfc = " + splitLine[84]
		if splitLine[85] != 'NA' :
			query += ", fischer_7_10_logcpm = " + splitLine[85]
		if splitLine[86] != 'NA' :
			query += ", fischer_7_10_pvalue = " + splitLine[86]
		if splitLine[87] != 'NA' :
			query += ", fischer_7_10_neglogfdr = " + splitLine[87]
		if splitLine[88] != 'NA' :
			query += ", fischer_7_10_sig = '" + splitLine[88]+"'"
		if splitLine[89] != 'NA' :
			query += ", fischer_7_11_logfc = " + splitLine[89]
		if splitLine[90] != 'NA' :
			query += ", fischer_7_11_logcpm = " + splitLine[90]
		if splitLine[91] != 'NA' :
			query += ", fischer_7_11_pvalue = " + splitLine[91]
		if splitLine[92] != 'NA' :
			query += ", fischer_7_11_neglogfdr = " + splitLine[92]
		if splitLine[93] != 'NA' :
			query += ", fischer_7_11_sig = '" + splitLine[93]+"'"
		if splitLine[94] != 'NA' :
			query += ", fischer_7_12_logfc = " + splitLine[94]
		if splitLine[95] != 'NA' :
			query += ", fischer_7_12_logcpm = " + splitLine[95]
		if splitLine[96] != 'NA' :
			query += ", fischer_7_12_pvalue = " + splitLine[96]
		if splitLine[97] != 'NA' :
			query += ", fischer_7_12_neglogfdr = " + splitLine[97]
		if splitLine[98] != 'NA' :
			query += ", fischer_7_12_sig = '" + splitLine[98]+"'"
		if splitLine[99] != 'NA' :
			query += ", fischer_7_13_logfc = " + splitLine[99]
		if splitLine[100] != 'NA' :
			query += ", fischer_7_13_logcpm = " + splitLine[100]
		if splitLine[101] != 'NA' :
			query += ", fischer_7_13_pvalue = " + splitLine[101]
		if splitLine[102] != 'NA' :
			query += ", fischer_7_13_neglogfdr = " + splitLine[102]
		if splitLine[103] != 'NA' :
			query += ", fischer_7_13_sig = '" + splitLine[103]+"'"
		if splitLine[104] != 'NA' :
			query += ", fischer_7_14_logfc = " + splitLine[104]
		if splitLine[105] != 'NA' :
			query += ", fischer_7_14_logcpm = " + splitLine[105]
		if splitLine[106] != 'NA' :
			query += ", fischer_7_14_pvalue = " + splitLine[106]
		if splitLine[107] != 'NA' :
			query += ", fischer_7_14_neglogfdr = " + splitLine[107]
		if splitLine[108] != 'NA' :
			query += ", fischer_7_14_sig = '" + splitLine[108]+"'"
		if splitLine[109] != 'NA' :
			query += ", fischer_7_15_logfc = " + splitLine[109]
		if splitLine[110] != 'NA' :
			query += ", fischer_7_15_logcpm = " + splitLine[110]
		if splitLine[111] != 'NA' :
			query += ", fischer_7_15_pvalue = " + splitLine[111]
		if splitLine[112] != 'NA' :
			query += ", fischer_7_15_neglogfdr = " + splitLine[112]
		if splitLine[113] != 'NA' :
			query += ", fischer_7_15_sig = '" + splitLine[113]+"'"
		if splitLine[114] != 'NA' :
			query += ", fischer_7_16_logfc = " + splitLine[114]
		if splitLine[115] != 'NA' :
			query += ", fischer_7_16_logcpm = " + splitLine[115]
		if splitLine[116] != 'NA' :
			query += ", fischer_7_16_pvalue = " + splitLine[116]
		if splitLine[117] != 'NA' :
			query += ", fischer_7_16_neglogfdr = " + splitLine[117]
		if splitLine[118] != 'NA' :
			query += ", fischer_7_16_sig = '" + splitLine[118]+"'"
		if splitLine[119] != 'NA' :
			query += ", fischer_7_17_logfc = " + splitLine[119]
		if splitLine[120] != 'NA' :
			query += ", fischer_7_17_logcpm = " + splitLine[120]
		if splitLine[121] != 'NA' :
			query += ", fischer_7_17_pvalue = " + splitLine[121]
		if splitLine[122] != 'NA' :
			query += ", fischer_7_17_neglogfdr = " + splitLine[122]
		if splitLine[123] != 'NA' :
			query += ", fischer_7_17_sig = '" + splitLine[123]+"'"
		if splitLine[124] != 'NA' :
			query += ", fischer_7_18_logfc = " + splitLine[124]
		if splitLine[125] != 'NA' :
			query += ", fischer_7_18_logcpm = " + splitLine[125]
		if splitLine[126] != 'NA' :
			query += ", fischer_7_18_pvalue = " + splitLine[126]
		if splitLine[127] != 'NA' :
			query += ", fischer_7_18_neglogfdr = " + splitLine[127]
		if splitLine[128] != 'NA' :
			query += ", fischer_7_18_sig = '" + splitLine[128]+"'"
		if splitLine[129] != 'NA' :
			query += ", fischer_7_19_logfc = " + splitLine[129]
		if splitLine[130] != 'NA' :
			query += ", fischer_7_19_logcpm = " + splitLine[130]
		if splitLine[131] != 'NA' :
			query += ", fischer_7_19_pvalue = " + splitLine[131]
		if splitLine[132] != 'NA' :
			query += ", fischer_7_19_neglogfdr = " + splitLine[132]
		if splitLine[133] != 'NA' :
			query += ", fischer_7_19_sig = '" + splitLine[133]+"'"
		if splitLine[134] != 'NA' :
			query += ", helm_7_12_logfc = " + splitLine[134]
		if splitLine[135] != 'NA' :
			query += ", helm_7_12_logcpm = " + splitLine[135]
		if splitLine[136] != 'NA' :
			query += ", helm_7_12_pvalue = " + splitLine[136]
		if splitLine[137] != 'NA' :
			query += ", helm_7_12_neglogfdr = " + splitLine[137]
		if splitLine[138] != 'NA' :
			query += ", helm_7_12_sig = '" + splitLine[138]+"'"
		if splitLine[139] != 'NA' :
			query += ", helm_7_24_logfc = " + splitLine[139]
		if splitLine[140] != 'NA' :
			query += ", helm_7_24_logcpm = " + splitLine[140]
		if splitLine[141] != 'NA' :
			query += ", helm_7_24_pvalue = " + splitLine[141]
		if splitLine[142] != 'NA' :
			query += ", helm_7_24_neglogfdr = " + splitLine[142]
		if splitLine[143] != 'NA' :
			query += ", helm_7_24_sig = '" + splitLine[143]+"'"
		if splitLine[144] != 'NA' :
			query += ", helm_7_120_logfc = " + splitLine[144]
		if splitLine[145] != 'NA' :
			query += ", helm_7_120_logcpm = " + splitLine[145]
		if splitLine[146] != 'NA' :
			query += ", helm_7_120_pvalue = " + splitLine[146]
		if splitLine[147] != 'NA' :
			query += ", helm_7_120_neglogfdr = " + splitLine[147]
		if splitLine[148] != 'NA' :
			query += ", helm_7_120_sig = '" + splitLine[148]+"'"
		if splitLine[149] != 'NA' :
			query += ", helm_7_240_logfc = " + splitLine[149]
		if splitLine[150] != 'NA' :
			query += ", helm_7_240_logcpm = " + splitLine[150]
		if splitLine[151] != 'NA' :
			query += ", helm_7_240_pvalue = " + splitLine[151]
		if splitLine[152] != 'NA' :
			query += ", helm_7_240_neglogfdr = " + splitLine[152]
		if splitLine[153] != 'NA' :
			query += ", helm_7_240_sig = '" + splitLine[153]+"'"
		if splitLine[154] != 'NA' :
			query += ", warner_24_48_logfc = " + splitLine[154]
		if splitLine[155] != 'NA' :
			query += ", warner_24_48_logcpm = " + splitLine[155]
		if splitLine[156] != 'NA' :
			query += ", warner_24_48_pvalue = " + splitLine[156]
		if splitLine[157] != 'NA' :
			query += ", warner_24_48_neglogfdr = " + splitLine[157]
		if splitLine[158] != 'NA' :
			query += ", warner_24_48_sig = '" + splitLine[158]+"'"
		if splitLine[159] != 'NA' :
			query += ", warner_24_72_logfc = " + splitLine[159]
		if splitLine[160] != 'NA' :
			query += ", warner_24_72_logcpm = " + splitLine[160]
		if splitLine[161] != 'NA' :
			query += ", warner_24_72_pvalue = " + splitLine[161]
		if splitLine[162] != 'NA' :
			query += ", warner_24_72_neglogfdr = " + splitLine[162]
		if splitLine[163] != 'NA' :
			query += ", warner_24_72_sig = '" + splitLine[163]+"'"
		if splitLine[164] != 'NA' :
			query += ", warner_24_96_logfc = " + splitLine[164]
		if splitLine[165] != 'NA' :
			query += ", warner_24_96_logcpm = " + splitLine[165]
		if splitLine[166] != 'NA' :
			query += ", warner_24_96_pvalue = " + splitLine[166]
		if splitLine[167] != 'NA' :
			query += ", warner_24_96_neglogfdr = " + splitLine[167]
		if splitLine[168] != 'NA' :
			query += ", warner_24_96_sig = '" + splitLine[168]+"'"
		if splitLine[169] != 'NA' :
			query += ", warner_24_120_logfc = " + splitLine[169]
		if splitLine[170] != 'NA' :
			query += ", warner_24_120_logcpm = " + splitLine[170]
		if splitLine[171] != 'NA' :
			query += ", warner_24_120_pvalue = " + splitLine[171]
		if splitLine[172] != 'NA' :
			query += ", warner_24_120_neglogfdr = " + splitLine[172]
		if splitLine[173] != 'NA' :
			query += ", warner_24_120_sig = '" + splitLine[173]+"'"
		if splitLine[174] != 'NA' :
			query += ", warner_24_144_logfc = " + splitLine[174]
		if splitLine[175] != 'NA' :
			query += ", warner_24_144_logcpm = " + splitLine[175]
		if splitLine[176] != 'NA' :
			query += ", warner_24_144_pvalue = " + splitLine[176]
		if splitLine[177] != 'NA' :
			query += ", warner_24_144_neglogfdr = " + splitLine[177]
		if splitLine[178] != 'NA' :
			query += ", warner_24_144_sig = '" + splitLine[178]+"'"
		if splitLine[179] != 'NA' :
			query += ", warner_24_168_logfc = " + splitLine[179]
		if splitLine[180] != 'NA' :
			query += ", warner_24_168_logcpm = " + splitLine[180]
		if splitLine[181] != 'NA' :
			query += ", warner_24_168_pvalue = " + splitLine[181]
		if splitLine[182] != 'NA' :
			query += ", warner_24_168_neglogfdr = " + splitLine[182]
		if splitLine[183] != 'NA' :
			query += ", warner_24_168_sig = '" + splitLine[183]+"'"
		if splitLine[184] != 'NA' :
			query += ", warner_24_192_logfc = " + splitLine[184]
		if splitLine[185] != 'NA' :
			query += ", warner_24_192_logcpm = " + splitLine[185]
		if splitLine[186] != 'NA' :
			query += ", warner_24_192_pvalue = " + splitLine[186]
		if splitLine[187] != 'NA' :
			query += ", warner_24_192_neglogfdr = " + splitLine[187]
		if splitLine[188] != 'NA' :
			query += ", warner_24_192_sig = '" + splitLine[188]+"'"
		query += ").save()"
		
	else :
		print "error in DBTableName. Valid inputs are : Regen_cpm ; Fasta ; Embryo_cpm ; Annotation ; Regen_SE ; Embryo_SE; de_table"
	return query

#execute the code
def save(DBQueryList) :
	for queries in DBQueryList :
		exec(queries)

##################
### Main
##################

#read the file and process :
def DBFill(fileName,table):
	#file to handle ; temp : to be replaced by command-line args
	dataFile = fileName
	DBQuery = []
	readDataFile = open(dataFile,'r')
	readDataFile.readline() #this only applies if there is a header to skip
	if table == 'Fasta' :
		splitLine = [readDataFile.readline()[1:-1],'']
		for line in readDataFile :
			if line[0] == '>' :
				DBQuery.append(queryCreate(splitLine,table))
				splitLine = [line[1:-1]]
				splitLine.append('')
			else :
				splitLine[1] += line.strip()
		DBQuery.append(queryCreate(splitLine,table))
	else :
		for lines in readDataFile :
			splitLine = dataSplit(lines)
			DBQuery.append(queryCreate(splitLine,table))
	readDataFile.close()
	#save the modifications
	save(DBQuery)
	return 'done filling the DB'