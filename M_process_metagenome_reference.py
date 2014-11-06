'''
Process large metagenomic reference dataset:
	Make names fully readable by replacing spaces with _ (done)
	Discard virus genomes (done)
	Chop genome records into equal-size chunks (done)
	Assemble plasmids?? (no)
	Split into files of less than 3.6 billion characters for bowtie2 index (done)
		note that the plasmid file never gets too big, so I don't bother checking if I should split it
	Remove duplicate genomes (done)

	Write unit tests
	check that all lines begin as expected: with > or with nucleotides
'''

import sys
import os
import math
import lanthpy
import pprint

# Is the strain already present in the i100 fasta?
def is_duplicate(strain):
	i100_strains = ['acinetobacter_baumannii_sdf', 'alkalilimnicola_ehrlichii_mlhe1', 'alkaliphilus_metalliredigens_qymf', 'anabaena_variabilis_atcc_29413', 'bacillus_anthracis_str_ames', 'bacillus_cereus_atcc_10987', 'bacillus_cereus_atcc_14579', 'bacillus_cereus_e33l', 'bacillus_clausii_ksmk16', 'bacillus_halodurans_c125', 'bacillus_subtilis_subsp_subtilis_str_168', 'bacillus_thuringiensis_str_al_hakam', 'bacteroides_thetaiotaomicron_vpi5482', 'bordetella_bronchiseptica_rb50', 'borrelia_garinii_pbi', 'buchnera_aphidicola_bcc', 'burkholderia_phymatum_stm815', 'burkholderia_pseudomallei_668', 'campylobacter_concisus_13826', 'candidatus_blochmannia_pennsylvanicus_str_bpen', 'chlamydia_trachomatis_l2buch1proctitis', 'chlamydophila_caviae_gpic', 'chlamydophila_pneumoniae_cwl029', 'chlorobium_luteolum_dsm_273', 'chromobacterium_violaceum_atcc_12472', 'clostridium_beijerinckii_ncimb_8052', 'clostridium_botulinum_b1_str_okra', 'clostridium_novyi_nt', 'clostridium_thermocellum_atcc_27405', 'corynebacterium_efficiens_ys314', 'corynebacterium_urealyticum_dsm_7109', 'coxiella_burnetii_rsa_493', 'cyanothece_sp_atcc_51142', 'cytophaga_hutchinsonii_atcc_33406', 'erwinia_tasmaniensis', 'escherichia_coli_apec_o1', 'escherichia_coli_e24377a', 'escherichia_coli_str_k12_substr_mg1655', 'escherichia_coli_str_k12_substr_w3110', 'francisella_tularensis_subsp_tularensis_fsc198', 'geobacter_sulfurreducens_pca', 'haemophilus_influenzae_rd_kw20', 'haemophilus_somnus_129pt', 'haloquadratum_walsbyi_dsm_16790', 'herpetosiphon_aurantiacus_dsm_785', 'hyphomonas_neptunium_atcc_15444', 'idiomarina_loihiensis_l2tr', 'ignicoccus_hospitalis_kin4i', 'lactobacillus_delbrueckii_subsp_bulgaricus_atcc_baa365', 'lactobacillus_salivarius_ucc118', 'lactococcus_lactis_subsp_cremoris_sk11', 'lactococcus_lactis_subsp_lactis_il1403', 'lawsonia_intracellularis_phemn100', 'leptospira_borgpetersenii_serovar_hardjobovis_jb197', 'listeria_welshimeri_serovar_6b_str_slcc5334', 'methanococcus_maripaludis_c7', 'mycobacterium_avium_subsp_paratuberculosis_k10', 'mycobacterium_bovis_af212297', 'mycobacterium_marinum_m', 'mycobacterium_sp_jls', 'mycobacterium_tuberculosis_f11', 'neisseria_gonorrhoeae_fa_1090', 'neisseria_meningitidis_mc58', 'nostoc_punctiforme_pcc_73102', 'ochrobactrum_anthropi_atcc_49188', 'prochlorococcus_marinus_str_mit_9313', 'prochlorococcus_marinus_str_natl1a', 'pseudomonas_putida_w619', 'pseudomonas_stutzeri_a1501', 'psychrobacter_cryohalolentis_k5', 'psychromonas_ingrahamii_37', 'ralstonia_eutropha_h16', 'rhodococcus_jostii_rha1', 'rhodopseudomonas_palustris_bisa53', 'rhodopseudomonas_palustris_cga009', 'rickettsia_prowazekii_str_madrid_e', 'salmonella_enterica_subsp_enterica_serovar_paratyphi_b_str_spb7', 'shewanella_sp_mr7', 'sinorhizobium_meliloti_1021', "sodalis_glossinidius_str_morsitans", 'staphylococcus_aureus_subsp_aureus_jh1', 'staphylococcus_aureus_subsp_aureus_mrsa252', 'staphylococcus_aureus_subsp_aureus_mu3', 'staphylococcus_aureus_subsp_aureus_mw2', 'staphylococcus_aureus_subsp_aureus_nctc_8325', 'staphylococcus_aureus_subsp_aureus_str_newman', 'streptococcus_agalactiae_2603vr', 'streptococcus_pneumoniae_r6', 'streptococcus_pyogenes_mgas315', 'streptococcus_pyogenes_mgas5005', 'synechococcus_elongatus_pcc_6301', 'synechococcus_elongatus_pcc_7942', 'syntrophomonas_wolfei_subsp_wolfei_str_goettingen', 'thermosipho_melanesiensis_bi429', 'thermus_thermophilus_hb8', 'thiobacillus_denitrificans_atcc_25259', 'treponema_denticola_atcc_35405', 'ureaplasma_parvum_serovar_3_str_atcc_700970', 'wigglesworthia_glossinidia_endosymbiont_of_glossina_brevipalpis', 'yersinia_pestis_co92']

	if any(dupe in lanthpy.genome_name_cleanup([strain])[0] for dupe in i100_strains):
		return True
	return False

# Calculate most balanced length of genome fragments, to avoid having particularly short fragments
def find_fragment_size(genome_length,ideal_size):
	num_of_fragments = max(round(float(genome_length)/ideal_size),1) # rounded to nearest int, minimum 1
	length_of_fragments = int(math.ceil(genome_length/float(num_of_fragments))) # rounded up to nearest bp
	return length_of_fragments

# Cut fasta entry into pieces of specified size
def chop_genome(genome_sequence,species,IDEAL_FRAGMENT_SIZE):

	fragment_size = find_fragment_size(len(genome_sequence), IDEAL_FRAGMENT_SIZE)

	chopped_lines = ""
	for i in range(0,len(genome_sequence),fragment_size):
		chopped_lines += species[:-1] +'_'+ str(i) +'\n'+ genome_sequence[i:i+fragment_size] +'\n'

	chopped_lines = chopped_lines.replace('\n\n','\n') # remove extraneous double newlines
	return chopped_lines

def process_file(filename,SPLIT_SIZE,IDEAL_FRAGMENT_SIZE,skip_dupes=False):
	mfa = open(filename,'r')

	file_count_g = 0
	file_count_p = 0
	line_count = 0

	# Open files
	genomes_string = filename.rsplit('.', 1)[0] + "_fragmented_genomes_{0}.mfa" #mfa file format is for compatibility with mass-indexing script
	plasmids_string = filename.rsplit('.', 1)[0] + "_fragmented_plasmids_{0}.mfa"

	fa_genomes = open(genomes_string.format(file_count_g), 'wt')
	fa_genomes.seek(0) #overwrite if it exists
	fa_plasmids = open(plasmids_string.format(file_count_p), 'wt')
	fa_plasmids.seek(0) #overwrite if it exists

	prev_species = "."
	last_type = ""
	prev_genome = ""
	skip_current = False

	list_of_strains = []
	for line in mfa:
		if line[0] == '>': # fasta name
			list_of_strains.append(line[:-1])
			skip_current = False
			if prev_genome: # if we have a genome in buffer, we need to chop it and write it before handling the new entry

				fa_genomes.write(chop_genome(prev_genome,prev_species,IDEAL_FRAGMENT_SIZE))
				prev_genome = "" # wipe last genome once it's written

			line = line.replace (" ", "_")
			if line.find('VIRL') != -1: #ignore viruses
				last_type = 'virus'
			elif line.find('plasmid') != -1: #just copy plasmids
				last_type = 'plasmid'
				fa_plasmids.write(line)
			else: # prepare to chop genomes up into pieces
				last_type = 'genome'
				if skip_dupes:
					skip_current = is_duplicate(line)

				if not skip_current:
					prev_species = line # for purposes of labeling the genome fragments

					# Split file past size limit, but only between species entries
					fa_genomes.flush()
					filesize = os.fstat(fa_genomes.fileno())
					if filesize.st_size > SPLIT_SIZE:
						fa_genomes.close()
						file_count_g += 1
						fa_genomes = open(genomes_string.format(file_count_g), 'wt')

		else: #actual fasta content
			if last_type == 'genome' and not skip_current: # first store entire genome in variable
				prev_genome += line

			elif last_type == 'plasmid':
				fa_plasmids.write(line)

	# Catch the last genome if it's the last item in the file
	if prev_genome: # if we have a genome in buffer, we need to chop it and write it before handling the new entry
		fa_genomes.write(chop_genome(prev_genome,prev_species,IDEAL_FRAGMENT_SIZE))
		prev_genome = "" # wipe last genome once it's written

	# close files
	fa_genomes.truncate()
	fa_genomes.close()
	fa_plasmids.truncate()
	fa_plasmids.close()

def main():
	SPLIT_SIZE = 3500000000 # 3.5gb
	IDEAL_FRAGMENT_SIZE = 500000 # bp

	#Note: File needs to be sorted!

	filename = ""
	iterarg = iter(sys.argv)
	next(iterarg) #skip name of function
	for arg in iterarg:
		if filename == "":
			filename = arg #avoid extra space at beginning
		else:
			filename = filename + " " + arg

	skip_dupes = filename.startswith('Martin')

	process_file(filename,SPLIT_SIZE,IDEAL_FRAGMENT_SIZE,skip_dupes)


if __name__ == '__main__':
	main()
