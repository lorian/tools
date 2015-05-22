# Module containing frequently-used functions

import csv
import string
import os
import sys

def get_arg_w_spaces(args):
	"""
	Input sys.args[1:]
	Combine members of a list with spaces between, and strip any initial spaces
	Meant to handle command line arguments with spaces
	"""
	arg = ' '.join(args)
	arg.lstrip(" ")
	return arg

def get_script_path():
	return os.path.dirname(os.path.realpath(sys.argv[0]))

def genome_name_cleanup(raw_names):
	"""
	Takes a list of messy genome names such as found in the Martin_etal
	multifasta and returns cleaned-up, normalized versions.
	"""

#	i100_species = ['acinetobacter_baumannii_sdf', 'alkalilimnicola_ehrlichii_mlhe1', 'alkaliphilus_metalliredigens_qymf', 'anabaena_variabilis_atcc_29413', 'bacillus_anthracis_str_ames', 'bacillus_cereus_atcc_10987', 'bacillus_cereus_atcc_14579', 'bacillus_cereus_e33l', 'bacillus_clausii_ksmk16', 'bacillus_halodurans_c125', 'bacillus_subtilis_subsp_subtilis_str_168', 'bacillus_thuringiensis_str_al_hakam', 'bacteroides_thetaiotaomicron_vpi5482', 'bordetella_bronchiseptica_rb50', 'borrelia_garinii_pbi', 'buchnera_aphidicola_bcc', 'burkholderia_phymatum_stm815', 'burkholderia_pseudomallei_668', 'campylobacter_concisus_13826', 'candidatus_blochmannia_pennsylvanicus_str_bpen', 'chlamydia_trachomatis_l2buch1proctitis', 'chlamydophila_caviae_gpic', 'chlamydophila_pneumoniae_cwl029', 'chlorobium_luteolum_dsm_273', 'chromobacterium_violaceum_atcc_12472', 'clostridium_beijerinckii_ncimb_8052', 'clostridium_botulinum_b1_str_okra', 'clostridium_novyi_nt', 'clostridium_thermocellum_atcc_27405', 'corynebacterium_efficiens_ys314', 'corynebacterium_urealyticum_dsm_7109', 'coxiella_burnetii_rsa_493', 'cyanothece_sp_atcc_51142', 'cytophaga_hutchinsonii_atcc_33406', 'erwinia_tasmaniensis', 'escherichia_coli_apec_o1', 'escherichia_coli_e24377a', 'escherichia_coli_str_k12_substr_mg1655', 'escherichia_coli_str_k12_substr_w3110', 'francisella_tularensis_subsp_tularensis_fsc198', 'geobacter_sulfurreducens_pca', 'haemophilus_influenzae_rd_kw20', 'haemophilus_somnus_129pt', 'haloquadratum_walsbyi_dsm_16790', 'herpetosiphon_aurantiacus_dsm_785', 'hyphomonas_neptunium_atcc_15444', 'idiomarina_loihiensis_l2tr', 'ignicoccus_hospitalis_kin4i', 'lactobacillus_delbrueckii_subsp_bulgaricus_atcc_baa365', 'lactobacillus_salivarius_ucc118', 'lactococcus_lactis_subsp_cremoris_sk11', 'lactococcus_lactis_subsp_lactis_il1403', 'lawsonia_intracellularis_phemn100', 'leptospira_borgpetersenii_serovar_hardjobovis_jb197', 'listeria_welshimeri_serovar_6b_str_slcc5334', 'methanococcus_maripaludis_c7', 'mycobacterium_avium_subsp_paratuberculosis_k10', 'mycobacterium_bovis_af212297', 'mycobacterium_marinum_m', 'mycobacterium_sp_jls', 'mycobacterium_tuberculosis_f11', 'neisseria_gonorrhoeae_fa_1090', 'neisseria_meningitidis_mc58', 'nostoc_punctiforme_pcc_73102', 'ochrobactrum_anthropi_atcc_49188', 'prochlorococcus_marinus_str_mit_9313', 'prochlorococcus_marinus_str_natl1a', 'pseudomonas_putida_w619', 'pseudomonas_stutzeri_a1501', 'psychrobacter_cryohalolentis_k5', 'psychromonas_ingrahamii_37', 'ralstonia_eutropha_h16', 'rhodococcus_jostii_rha1', 'rhodopseudomonas_palustris_bisa53', 'rhodopseudomonas_palustris_cga009', 'rickettsia_prowazekii_str_madrid_e', 'salmonella_enterica_subsp_enterica_serovar_paratyphi_b_str_spb7', 'shewanella_sp_mr7', 'sinorhizobium_meliloti_1021', "sodalis_glossinidius_str_morsitans", 'staphylococcus_aureus_subsp_aureus_jh1', 'staphylococcus_aureus_subsp_aureus_mrsa252', 'staphylococcus_aureus_subsp_aureus_mu3', 'staphylococcus_aureus_subsp_aureus_mw2', 'staphylococcus_aureus_subsp_aureus_nctc_8325', 'staphylococcus_aureus_subsp_aureus_str_newman', 'streptococcus_agalactiae_2603vr', 'streptococcus_pneumoniae_r6', 'streptococcus_pyogenes_mgas315', 'streptococcus_pyogenes_mgas5005', 'synechococcus_elongatus_pcc_6301', 'synechococcus_elongatus_pcc_7942', 'syntrophomonas_wolfei_subsp_wolfei_str_goettingen', 'thermosipho_melanesiensis_bi429', 'thermus_thermophilus_hb8', 'thiobacillus_denitrificans_atcc_25259', 'treponema_denticola_atcc_35405', 'ureaplasma_parvum_serovar_3_str_atcc_700970', 'wigglesworthia_glossinidia_endosymbiont_of_glossina_brevipalpis', 'yersinia_pestis_co92']

	# If BACT_ abbreviations are found, see if they match large genome database
	try:
		min( i for i, sp in enumerate(raw_names) if 'BACT_' in sp )
	except:
		clean_names = raw_names
	else:
		with open(os.path.join(get_script_path(),'M_biggenomenames.txt'), 'r') \
					as biggenomefile:
			big_genome = dict(csv.reader(biggenomefile))
		clean_names = [big_genome[s.partition('|')[0]]+"_"+s.partition('|')[2]
						if s.partition('|')[0] in big_genome.keys()
						else s for s in raw_names]

	# Force clean_names to be lowercase, underscored, and without punctuation
	bad_punct = "!\"#$%&'()*+,-./:;<=>@[\]^`{}~" # string.punctuation w/no _ | and ? (latter used to separate synonyms)
	clean_names = [s.lower().translate(string.maketrans("",""),bad_punct)
					.replace(" ","_") for s in clean_names]

	'''
	# Look for matches in the i100 genome list, to remove fragment numbers
	for i100 in i100_species:
		for i,s in enumerate(clean_names):
			if s.startswith(i100):
				clean_names[i] = i100
	'''
	return clean_names
