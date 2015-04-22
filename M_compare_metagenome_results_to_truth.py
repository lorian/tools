# Compares output of metagenomic analysis tools with ground truth of dataset

import sys
import csv
import numpy
import string
import matplotlib
import math
import pickle
import lanthpy
import pprint

numpy.set_printoptions(precision=4)

class Dataset():

	def __init__(self, species=[], abundance=[], counts=[], size=[]):
		self.species = list(species)
		self.abundance = list(abundance)
		self.counts = list(counts)
		self.size = list(size) # can be empty
		self.clean_names()

	def add_record(self, new_species, new_abundance, new_counts, new_size=None):
		self.species.extend(lanthpy.genome_name_cleanup([new_species]))
		self.abundance.append(new_abundance)
		self.counts.append(new_counts)
		self.size.append(new_size)

	def set_by_array(self, array):
		# Warning: this overwrites all existing data
		try:
			self.species, self.abundance, self.counts, self.size = zip(*array)
		except:
			self.species, self.abundance, self.counts = zip(*array)
		self.species = list(self.species)
		self.abundance = list(self.abundance)
		self.counts = list(self.counts)
		self.size = list(self.size)

	def print_record(self, species):
		print "{0}:\n\t{1} ({2} reads)".format(species,self.lookup_abundance(species),self.lookup_count(species))

	def get_array(self):
		return zip(self.species,self.abundance,self.counts)

	def lookup_abundance(self, species):
		try:
			return self.abundance[self.species.index(species)]
		except:
			return 0

	def lookup_count(self, species):
		try:
			return self.counts[self.species.index(species)]
		except:
			return 0

	def lookup_length(self):
		assert len(self.species) == len(self.abundance) == len(self.counts)
		return len(self.species)

	def clean_names(self):
		self.species = lanthpy.genome_name_cleanup(self.species)

	def sort_by_name(self):
		# Alphabetical by species
		self.set_by_array(sorted(self.get_array(),key=lambda x:x[0]))

	def set_threshold(self, threshold=0.001):
		# Set very low estimated abundances to 0
		self.set_by_array([r if r[1] > threshold else (r[0],0,r[2]) for r in self.get_array()])
		print "Number of non-zero abundances: {0}".format(self.lookup_length() - self.abundance.count(0))

	def convert_to_percentage(self):
		# Normalize TPM/FPKM abundances
		total_ab = math.fsum(self.abundance)
		self.abundance = [100*v/total_ab for v in self.abundance]
		assert math.fsum(self.abundance) == 100

	def remove_matches(self, target):
		self.set_by_array([r for r in self.get_array() if (r[0].find(target) == -1)])
		print "Number of entries after removing {0}: {1}".format(target,len(self.species))


def collapse_strains(strains,abundance): #NOT WORKING
	""" Group strains together by species and genus """

	just_genus = [s.split("_")[0] for s in strains] # Note: ignores synonyms (probably safe)
	just_species = [s.split("_")[0]+s.split("_")[1] if (s.count("_") > 1) else s for s in strains]
	#just_species = [s.partition(',')[0].split("_")[0]+"_"+s.partition(',')[0].split("_")[1]+","+s.partition(',')[2].split("_")[0]+"_"+s.partition(',')[2].split("_")[1] if s.find(',') > -1 else s.split("_")[0]+"_"+s.split("_")[1] for s in strains]

	species_combo = zip(just_species,abundance)
	species_dict = {}
	for s,a in species_combo:
		species_dict.setdefault(s,0)
		species_dict[s] = species_dict[s] + a

	genus_combo = zip(just_genus,abundance)
	genus_dict = {}
	for s,a in genus_combo:
		genus_dict.setdefault(s,0)
		genus_dict[s] = genus_dict[s] + a

	temp1,temp2 = zip(*sorted(species_dict.items()))
	temp3,temp4 = zip(*sorted(genus_dict.items()))
	return list(temp1),numpy.array(temp2),list(temp3),numpy.array(temp4)

def collapse_duplicates(raw_data):

	# Create dictionary of lists of duplicates
	dup_data = raw_data.get_array()
	set_sp = {}
	set_ab = {}
	set_co = {}
	for sp,ab,co in dup_data:
		name = sp.partition('_gi')[0] #the prepended strain name
		set_sp.setdefault(name,[]).append(sp)
		set_ab.setdefault(name,[]).append(ab)
		set_co.setdefault(name,[]).append(co)
	assert(set_ab.keys() == set_co.keys() == set_sp.keys())

	# New, clean dataset for data without duplicates
	undupe = Dataset()

	for k,v in set_sp.items():
		if len(v) == 1: # just add record directly if it has no duplicates
			undupe.add_record(k,set_ab[k][0],set_co[k][0])
		else:
			if all('chromosome' in dupe for dupe in v):
				# Average abundances of separate chromosomes
				undupe.add_record(k,math.fsum(set_ab[k])/len(v),math.fsum(set_co[k]))
			elif all('complete' in dupe for dupe in v):
				# Sum abundances of duplicate genomes
				undupe.add_record(k,math.fsum(set_ab[k]),math.fsum(set_co[k]))
			else:
				# Should insert more complicated handling here; just summing for now
				print "Possible error in {0}: neither chromosome nor genome?".format(k)
				undupe.add_record(k,math.fsum(set_ab[k]),math.fsum(set_co[k]))

	print "Number of entries after combining duplicates: {0}".format(undupe.lookup_length())

	return undupe

def process_input(filename,size,fragmented=False):
	""" Pull species names, abundance, and counts out of input gasic or express file """

	suffix = filename.rpartition('.')[2] #xprs or txt file

	input_file = open(filename,'r')
	input_csv = csv.reader(input_file, 'excel-tab')
	input_data = [r for r in input_csv]
	input_data = input_data[1:] #remove header row

	raw_est = Dataset()
	if suffix == 'xprs': #express file
		raw_est.species = zip(*input_data)[1] # express species names are in second column
		raw_est.counts = [float(i) for i in zip(*input_data)[6]] # used when dataset abundances are given in raw counts
		raw_est.abundance = [float(i) for i in zip(*input_data)[10]] # used when dataset abundances are given as percentages

	elif filename.endswith('expression.txt'): # kallisto file
		raw_est.species = zip(*input_data)[0] # kallisto species names are in first column
		raw_est.counts = [float(i) for i in zip(*input_data)[4]]
		raw_est.abundance = [float(i) for i in zip(*input_data)[3]]

	elif suffix == 'txt': #gasic file - NOT WORKING
		raw_est.species = zip(*input_data)[0] # gasic species names are in first column
		raw_est.counts = [float(i) for i in zip(*input_data)[2]]
		raw_est.abundance = raw_est.counts/(numpy.sum(raw_est.counts)*numpy.array(size)) # used when dataset abundances are given as percentages

	else:
		print "File is not supported input type"

	print "Number of raw entries: {0}".format(raw_est.lookup_length())
	raw_est.clean_names()
	raw_est.remove_matches('plasmid')
	# Also want to remove ribosomal_RNA_gene, mitochondrion, and chloroplast

	est = collapse_duplicates(raw_est)
	est.convert_to_percentage()
	est.sort_by_name()
	est.set_threshold()

	with open('non-zero_abundances.txt','w') as f:
		for i,species in enumerate(est.species):
			if est.abundance[i] != 0:
				f.write('{0}\t{1}\n'.format(species,est.abundance[i]))

	est_species_alone, est_species_ab, est_genus_alone, est_genus_ab = collapse_strains(est.species,est.abundance)
	return est,est_species_alone,est_species_ab,est_genus_alone,est_genus_ab

def dataset_truth(dataset):
	""" truth is in format |species|abundance|counts|genome size| where undefined counts is 0
	 and species may include alternate species names separated by a comma """

	# ? used to separate synonyms
	truth = Dataset()
	if dataset == 'i100':
		truth.set_by_array([("acinetobacter_baumannii_sdf",0.9918569894,0,3477999),("alkalilimnicola_ehrlichii_mlhe-1",0.9534103189,0,3275944),("alkaliphilus_metalliredigens_qymf",1.0337468099,0,4929566),("anabaena_variabilis_atcc_29413",0.8701233754,0,7068604),("bacillus_anthracis_str._ames",1.0629349921,0,5227293),("bacillus_cereus_atcc_10987",0.9077576661,0,5432653),("bacillus_cereus_atcc_14579",1.4619387586,0,5427084),("bacillus_cereus_e33l",0.9147432601,0,5843240),("bacillus_clausii_ksm-k16",0.9028911617,0,4303871),("bacillus_halodurans_c-125",0.8847367596,0,4202352),("bacillus_subtilis_subsp._subtilis_str._168",0.8765090628,0,4214547),("bacillus_thuringiensis_str._al_hakam",0.9505918962,0,5313031),("bacteroides_thetaiotaomicron_vpi-5482",0.8661600805,0,6293400),("bordetella_bronchiseptica_rb50",1.1977086557,0,5339179),("borrelia_garinii_pbi",1.1270363889,0,986916),("buchnera_aphidicola_bcc",0.8652027125,0,422435),("burkholderia_phymatum_stm815",1.0474017422,0,8676565),("burkholderia_pseudomallei_668",0.8953944405,0,7040404),("campylobacter_concisus_13826",1.0403648471,0,2099415),("candidatus_blochmannia_pennsylvanicus_str._bpen",0.9129394805,0,791654),("chlamydia_trachomatis_l2b/uch-1/proctitis",0.9593251688,0,1038863),("chlamydophila_caviae_gpic",0.8711493582,0,1181357),("chlamydophila_pneumoniae_cwl029",0.8810983686,0,1230230),("chlorobium_luteolum_dsm_273",0.9264853736,0,2364842),("chromobacterium_violaceum_atcc_12472",0.9624326898,0,4751080),("clostridium_beijerinckii_ncimb_8052",0.9656484821,0,6000632),("clostridium_botulinum_b1_str._okra",0.9478594593,0,4107014),("clostridium_novyi_nt",1.1580488951,0,2547720),("clostridium_thermocellum_atcc_27405",0.9330381139,0,3843301),("corynebacterium_efficiens_ys-314",1.0275058892,0,3219507),("corynebacterium_urealyticum_dsm_7109",1.0008453082,0,2369219),("coxiella_burnetii_rsa_493",0.9184749477,0,2032675),("cyanothece_sp._atcc_51142",0.9401343144,0,5460382),("cytophaga_hutchinsonii_atcc_33406",0.9563196367,0,4433218),("erwinia_tasmaniensis",1.1018606017,0,4067869),("escherichia_coli_apec_o1",0.9061021819,0,5497655),("escherichia_coli_e24377a",0.9244094549,0,5249294),("escherichia_coli_str._k-12_substr._mg1655",0.9836152016,0,4639637),("escherichia_coli_str._k-12_substr._w3110",1.2861538913,0,4646332),("francisella_tularensis_subsp._tularensis_fsc198",0.8799240389,0,1892616),("geobacter_sulfurreducens_pca",0.8822915619,0,3814128),("haemophilus_influenzae_rd_kw20",0.8925873316,0,1830138),("haemophilus_somnus_129pt",0.9353392713,0,2012879),("haloquadratum_walsbyi_dsm_16790",0.9724329129,0,3179362),("herpetosiphon_aurantiacus_dsm_785",0.8671305022,0,6785432),("hyphomonas_neptunium_atcc_15444",0.9797413592,0,3705021),("idiomarina_loihiensis_l2tr",0.9111749381,0,2839318),("ignicoccus_hospitalis_kin4/i",0.9286138232,0,1297538),("lactobacillus_delbrueckii_subsp._bulgaricus_atcc_baa-365",1.1765523477,0,1856951),("lactobacillus_salivarius_ucc118",0.9013332759,0,2133980),("lactococcus_lactis_subsp._cremoris_sk11",0.9377036228,0,2598353),("lactococcus_lactis_subsp._lactis_il1403",1.3855635178,0,2365589),("lawsonia_intracellularis_phe/mn1-00",0.9204063033,0,1719017),("leptospira_borgpetersenii_serovar_hardjo-bovis_jb197",0.8633259013,0,3876236),("listeria_welshimeri_serovar_6b_str._slcc5334",0.8859899334,0,2814130),("methanococcus_maripaludis_c7",0.9452084692,0,1772694),("mycobacterium_avium_subsp._paratuberculosis_k-10",0.8872643124,0,4829781),("mycobacterium_bovis_af2122/97",0.8998056707,0,4345492),("mycobacterium_marinum_m",0.8776298824,0,6660145),("mycobacterium_sp._jls",1.7694735192,0,6048425),("mycobacterium_tuberculosis_f11",0.9044804156,0,4424435),("neisseria_gonorrhoeae_fa_1090",0.9876497049,0,2153922),("neisseria_meningitidis_mc58",0.9307972132,0,2272360),("nostoc_punctiforme_pcc_73102",0.8983073115,0,9059196),("ochrobactrum_anthropi_atcc_49188",1.0160158194,0,5205782),("prochlorococcus_marinus_str._mit_9313",1.0107078939,0,2410873),("prochlorococcus_marinus_str._natl1a",0.9689793951,0,1864731),("pseudomonas_putida_w619",1.2512067245,0,5774330),("pseudomonas_stutzeri_a1501",0.8912213302,0,4567418),("psychrobacter_cryohalolentis_k5",2.2276838861,0,3101098),("psychromonas_ingrahamii_37",0.8691118092,0,4559598),("ralstonia_eutropha_h16",0.8885605469,0,3364647),("rhodococcus_jostii_rha1",0.8835041681,0,9702740),("rhodopseudomonas_palustris_bisa53",1.0715558086,0,5505494),("rhodopseudomonas_palustris_cga009",1.0056581462,0,5476069),("rickettsia_prowazekii_str._madrid_e",0.8939780971,0,1111523),("salmonella_enterica_subsp._enterica_serovar_paratyphi_b_str._spb7",1.113841943,0,4858887),("shewanella_sp._mr-7",0.896837214,0,4799110),("sinorhizobium_meliloti_1021",0.8743176015,0,6691696),("sodalis_glossinidius_str._'morsitans'",0.9426347213,0,4292505),("staphylococcus_aureus_subsp._aureus_jh1",1.1416733976,0,2936937),("staphylococcus_aureus_subsp._aureus_mrsa252",0.8754051191,0,2902619),("staphylococcus_aureus_subsp._aureus_mu3",0.8642580865,0,2880168),("staphylococcus_aureus_subsp._aureus_mw2",1.0909113027,0,2820462),("staphylococcus_aureus_subsp._aureus_nctc_8325",0.8787680456,0,2821361),("staphylococcus_aureus_subsp._aureus_str._newman",0.9760172336,0,2878897),("streptococcus_agalactiae_2603v/r",1.054906113,0,2160267),("streptococcus_pneumoniae_r6",1.2222455073,0,2038615),("streptococcus_pyogenes_mgas315",0.8732460774,0,1900521),("streptococcus_pyogenes_mgas5005",1.0216060126,0,1838554),("synechococcus_elongatus_pcc_6301",0.8681143003,0,2696255),("synechococcus_elongatus_pcc_7942",1.5752103822,0,2796473),("syntrophomonas_wolfei_subsp._wolfei_str._goettingen",0.9094481382,0,2936195),("thermosipho_melanesiensis_bi429",0.9223838063,0,1915238),("thermus_thermophilus_hb8",0.9962505203,0,2197216),("thiobacillus_denitrificans_atcc_25259",1.0808489742,0,2909809),("treponema_denticola_atcc_35405",0.8898793162,0,2843201),("ureaplasma_parvum_serovar_3_str._atcc_700970",0.87219013,0,751719),("wigglesworthia_glossinidia_endosymbiont_of_glossina_brevipalpis",1.3295506979,0,703005),("yersinia_pestis_co92",0.9165878579,0,4829858)])
	elif dataset == 'i400':
		truth.set_by_array([("acaryochloris_marina_mbic11017",0.6580017152,0,5300000),("acholeplasma_laidlawii_pg-8a",0.5226579129,0,3370000),("acidiphilium_cryptum_jf-5",0.4652774748,0,6100000),("acidothermus_cellulolyticus_11b",0.4092605674,0,860000),("acidovorax_sp._js42",0.3927157911,0,1110000),("acinetobacter_baumannii_atcc_17978",0.3798974674,0,6180000),("acinetobacter_sp._adp1",0.3695749545,0,3100000),("actinobacillus_pleuropneumoniae_l20",0.3610205403,0,750000),("actinobacillus_succinogenes_130z",0.3537729723,0,4600000),("aeromonas_hydrophila_subsp._hydrophila_atcc_7966",0.3475239317,0,6220000),("aeromonas_salmonicida_subsp._salmonicida_a449",0.3420584778,0,7210000),("alcanivorax_borkumensis_sk2",0.3372215683,0,7100000),("Alkalilimnicola ehrlichii MLHE-1?Alkalilimnicola ehrlichei MLHE-1",0.3328981646,0,2000000),("alkaliphilus_oremlandii_ohilas",0.3290008576,0,2230000),("anaeromyxobacter_dehalogenans_2cp-c",0.3254618711,0,7260000),("anaeromyxobacter_sp._fw109-5",0.3222277239,0,1200000),("anaplasma_marginale_str._st._maries",0.3192555655,0,5310000),("anaplasma_phagocytophilum_hz",0.3165105985,0,4590000),("aquifex_aeolicus_vf5",0.3139642264,0,7310000),("archaeoglobus_fulgidus_dsm_4304",0.3115926977,0,1690000),("arcobacter_butzleri_rm4018",0.3093760956,0,3900000),("arthrobacter_aurescens_tc1",0.307297574,0,4600000),("arthrobacter_sp._fb24",0.3053427725,0,1600000),("aster_yellows_witches'-broom_phytoplasma_aywb",0.3034993617,0,8530000),("Azoarcus aromaticum EbN1?Azoarcus sp. EbN1",0.3017566867,0,2100000),("azorhizobium_caulinodans_ors_571",0.3001054844,0,1800000),("Bacillus amyloliquefaciens subsp. plantarum str. FZB42?Bacillus amyloliquefaciens FZB42",0.2985376569,0,3340000),("bacillus_anthracis_str._'ames_ancestor'",0.2970460886,0,2040000),("bacillus_cereus_e33l",0.2956244975,0,3200000),("bacillus_halodurans_c-125",0.2942673129,0,2940000),("bacillus_licheniformis_atcc_14580",0.2929695745,0,6800000),("bacillus_pumilus_safr-032",0.2917268486,0,1010000),("bacillus_thuringiensis_str._al_hakam",0.290535158,0,4750000),("bacillus_weihenstephanensis_kbab4",0.2893909225,0,5230000),("bacteroides_fragilis_ych46",0.2882909095,0,6620000),("bacteroides_vulgatus_atcc_8482",0.2872321915,0,5850000),("bartonella_bacilliformis_kc583",0.2862121093,0,2300000),("Bartonella henselae strain Houston-1?Bartonella henselae str. Houston-1",0.2852282415,0,3300000),("bartonella_quintana_str._toulouse",0.2842783774,0,780000),("bartonella_tribocorum_cip_105476",0.283360494,0,1180000),("bifidobacterium_adolescentis_atcc_15703",0.2824727355,0,5800000),("bifidobacterium_longum_djo10a",0.2816133963,0,4300000),("bordetella_avium_197n",0.2807809052,0,6460000),("bordetella_parapertussis_strain_12822",0.2799738122,0,1510000),("Bordetella parapertussis strain 12822?Bordetella parapertussis 12822",0.2791907765,0,5430000),("Bordetella petrii strain DSM 12804?Bordetella petrii DSM 12804",0.2784305562,0,5190000),("borrelia_afzelii_pko",0.277691999,0,3400000),("borrelia_burgdorferi_b31",0.2769740339,0,1590000),("borrelia_garinii_pbi",0.276275664,0,6090000),("bradyrhizobium_japonicum_usda_110",0.27559596,0,8360000),("bradyrhizobium_sp._btai1",0.2749340544,0,2640000),("brucella_canis_atcc_23365",0.274289136,0,5900000),("Brucella melitensis bv. 1 str. 16M?Brucella melitensis 16M",0.2736604456,0,1300000),("brucella_ovis_atcc_25840",0.2730472715,0,5700000),("brucella_suis_1330",0.2724489459,0,3600000),("buchnera_aphidicola_str._aps_(acyrthosiphon_pisum)",0.2718648413,0,2170000),("burkholderia_cenocepacia_mc0-3",0.2712943675,0,7100000),("burkholderia_mallei_nctc_10247",0.2707369687,0,5500000),("burkholderia_multivorans_atcc_17616",0.270192121,0,1240000),("burkholderia_phymatum_stm815",0.2696593299,0,5030000),("burkholderia_pseudomallei_1710b",0.2691381286,0,2400000),("burkholderia_thailandensis_e264",0.2686280753,0,3630000),("burkholderia_vietnamiensis_g4",0.2681287524,0,5900000),("burkholderia_xenovorans_lb400",0.2676397641,0,2700000),("caldicellulosiruptor_saccharolyticus_dsm_8903",0.2671607352,0,4940000),("caldivirga_maquilingensis_ic-167",0.2666913096,0,5250000),("campylobacter_curvus_525.92",0.2662311494,0,1470000),("campylobacter_hominis_atcc_baa-381",0.2657799333,0,1900000),("campylobacter_jejuni_rm1221",0.2653373557,0,3270000),("candidatus_blochmannia_pennsylvanicus_str._bpen",0.2649031259,0,1600000),("candidatus_carsonella_ruddii_pv",0.264476967,0,6910000),("candidatus_desulforudis_audaxviator_mp104c",0.2640586148,0,5100000),("candidatus_korarchaeum_cryptofilum_opf8",0.2636478177,0,2100000),("Candidatus Koribacter versatilis Ellin345?Acidobacteria bacterium Ellin345",0.43181989,0,800000),("candidatus_pelagibacter_ubique_htcc1062",0.2632443353,0,5270000),("candidatus_ruthia_magnifica_str._cm_(calyptogena_magnifica)",0.2624584069,0,2200000),("candidatus_vesicomyosocius_okutanii_ha",0.2620755319,0,5910000),("carboxydothermus_hydrogenoformans_z-2901",0.2616991124,0,2700000),("caulobacter_crescentus_cb15",0.2613289565,0,1900000),("chlamydia_muridarum_nigg",0.26096488,0,730000),("chlamydia_trachomatis_duw-3cx",0.2606067066,0,3850000),("chlamydophila_caviae_gpic",0.2602542674,0,5800000),("chlamydophila_pneumoniae_ar39",0.2599074,0,660000),("Chlorobium luteolum DSM 273?Pelodictyon luteolum DSM 273",0.2322340449,0,3800000),("chlorobium_tepidum_tls",0.2595659487,0,6500000),("chloroflexus_aurantiacus_j-10-fl",0.2592297639,0,2000000),("chromobacterium_violaceum_atcc_12472",0.258898702,0,2400000),("chromohalobacter_salexigens_dsm_3043",0.2585726249,0,3710000),("citrobacter_koseri_atcc_baa-895",0.2582513996,0,2560000),("clostridium_acetobutylicum_atcc_824",0.2579348983,0,2750000),("clostridium_beijerinckii_ncimb_8052",0.2576229981,0,4550000),("clostridium_botulinum_a3_str._loch_maree",0.2573155803,0,2100000),("clostridium_difficile_630",0.2570125308,0,3900000),("clostridium_kluyveri_dsm_555",0.2567137396,0,4500000),("clostridium_novyi_nt",0.2564191006,0,3400000),("clostridium_perfringens_atcc_13124",0.2561285114,0,1700000),("clostridium_phytofermentans_isdg",0.2558418733,0,2310000),("clostridium_tetani_e88",0.2555590909,0,6720000),("clostridium_thermocellum_atcc_27405",0.2552800722,0,5960000),("colwellia_psychrerythraea_34h",0.2550047282,0,3800000),("corynebacterium_diphtheriae_nctc_13129",0.2547329729,0,5100000),("corynebacterium_efficiens_ys-314",0.2544647232,0,4060000),("corynebacterium_jeikeium_k411",0.2541998989,0,1300000),("corynebacterium_urealyticum_dsm_7109",0.2539384221,0,3810000),("coxiella_burnetii_dugway_5j108-111",0.2536802176,0,1930000),("cupriavidus_taiwanensis",0.2534252127,0,3150000),("cyanothece_sp._atcc_51142",0.2531733367,0,5300000),("cytophaga_hutchinsonii_atcc_33406",0.2529245213,0,2700000),("dechloromonas_aromatica_rcb",0.2526787005,0,4300000),("dehalococcoides_ethenogenes_195",0.25243581,0,3970000),("dehalococcoides_sp._cbdb1",0.2521957877,0,1900000),("deinococcus_radiodurans_r1",0.2519585732,0,3700000),("delftia_acidovorans_sph-1",0.251724108,0,1110000),("desulfococcus_oleovorans_hxd3",0.2514923354,0,2100000),("desulfotalea_psychrophila_lsv54",0.2512632003,0,9100000),("desulfotomaculum_reducens_mi-1",0.2510366492,0,3600000),("Desulfovibrio alaskensis G20?Desulfovibrio desulfuricans subsp. desulfuricans str. G20",0.2508126301,0,9010000),("desulfovibrio_vulgaris_subsp._vulgaris_str._hildenborough",0.2505910926,0,5200000),("dichelobacter_nodosus_vcs1703a",0.2503719877,0,1900000),("dinoroseobacter_shibae_dfl_12",0.2501552677,0,1600000),("ehrlichia_chaffeensis_str._arkansas",0.2499408863,0,4540000),("ehrlichia_ruminantium_str._welgevonden",0.2497287985,0,1910000),("enterobacter_sp._638",0.2495189605,0,2100000),("erwinia_tasmaniensis",0.2493113297,0,4700000),("erythrobacter_litoralis_htcc2594",0.2491058647,0,5900000),("exiguobacterium_sibiricum_255-15",0.248902525,0,5010000),("finegoldia_magna_atcc_29328",0.2487012715,0,2300000),("flavobacterium_johnsoniae_uw101",0.2485020659,0,7900000),("flavobacterium_psychrophilum_jip0286",0.248304871,0,4150000),("francisella_philomiragia_subsp._philomiragia_atcc_25017",0.2481096505,0,4730000),("Frankia alni str. ACN14a?Frankia alni ACN14a",0.2479163691,0,4600000),("fusobacterium_nucleatum_subsp._nucleatum_atcc_25586",0.2477249923,0,2650000),("geobacillus_kaustophilus_hta426",0.2475354867,0,1770000),("geobacillus_thermodenitrificans_ng80-2",0.2473478194,0,2600000),("geobacter_lovleyi_sz",0.2471619587,0,2300000),("geobacter_sulfurreducens_pca",0.2469778733,0,2620000),("geobacter_uraniireducens_rf4",0.246795533,0,2800000),("gloeobacter_violaceus_pcc_7421",0.2466149081,0,5200000),("gluconacetobacter_diazotrophicus_pal_5",0.2464359699,0,2000000),("gluconobacter_oxydans_621h",0.24625869,0,4310000),("gramella_forsetii_kt0803",0.246083041,0,1900000),("granulibacter_bethesdensis_cgdnih1",0.245908996,0,980000),("Haemophilus ducreyi strain 35000HP?Haemophilus ducreyi 35000HP",0.2457365288,0,3590000),("haemophilus_influenzae_pittgg",0.2455656137,0,2660000),("haemophilus_somnus_2336",0.2453962259,0,1900000),("haloarcula_marismortui_atcc_43049",0.2452283407,0,2940000),("halobacterium_salinarum_r1",0.2450619343,0,5800000),("halobacterium_sp._nrc-1",0.2448969835,0,2590000),("halorhodospira_halophila_sl1",0.2447334653,0,2130000),("helicobacter_acinonychis_str._sheeba",0.2445713574,0,2500000),("helicobacter_hepaticus_atcc_51449",0.2444106382,0,4690000),("helicobacter_pylori_26695",0.2442512861,0,2400000),("heliobacterium_modesticaldum_ice1",0.2440932805,0,4850000),("herminiimonas_arsenicoxydans",0.2439366008,0,2810000),("hyperthermus_butylicus_dsm_5456",0.2437812271,0,2580000),("hyphomonas_neptunium_atcc_15444",0.2436271399,0,880000),("idiomarina_loihiensis_l2tr",0.2434743201,0,1550000),("ignicoccus_hospitalis_kin4i",0.243322749,0,2150000),("jannaschia_sp._ccs1",0.2431724083,0,7800000),("janthinobacterium_sp._marseille",0.24302328,0,1760000),("kineococcus_radiotolerans_srs30216",0.2428753465,0,5000000),("klebsiella_pneumoniae_subsp._pneumoniae_mgh_78578",0.2427285908,0,2370000),("lactobacillus_acidophilus_ncfm",0.242582996,0,6600000),("lactobacillus_brevis_atcc_367",0.2424385455,0,2570000),("lactobacillus_casei_atcc_334",0.2422952233,0,490000),("lactobacillus_delbrueckii_subsp._bulgaricus_atcc_baa-365",0.2421530134,0,1520000),("lactobacillus_fermentum_ifo_3956",0.2420119005,0,3980000),("lactobacillus_gasseri_atcc_33323",0.2418718692,0,6380000),("lactobacillus_helveticus_dpc_4571",0.2417329048,0,2200000),("lactobacillus_plantarum_wcfs1",0.2415949925,0,1800000),("Lactobacillus reuteri DSM 20016?Lactobacillus reuteri F275",0.2414581181,0,2690000),("Lactobacillus sakei strain 23K?Lactobacillus sakei subsp. sakei 23K",0.2413222676,0,4080000),("lactococcus_lactis_subsp._cremoris_sk11",0.2411874271,0,1900000),("lawsonia_intracellularis_phemn1-00",0.2410535832,0,1800000),("legionella_pneumophila_str._paris",0.2409207226,0,1230000),("leifsonia_xyli_subsp._xyli_str._ctcb07",0.2407888324,0,3300000),("leptospira_biflexa_serovar_patoc_strain_'patoc_1_(paris)'",0.2406578997,0,4300000),("leptospira_borgpetersenii_serovar_hardjo-bovis_l550",0.2405279122,0,2300000),("leptospira_interrogans_serovar_lai_str._56601",0.2403988574,0,2920000),("leptothrix_cholodnii_sp-6",0.2402707234,0,2000000),("leuconostoc_citreum_km20",0.2401434984,0,6000000),("leuconostoc_mesenteroides_subsp._mesenteroides_atcc_8293",0.2400171708,0,2500000),("listeria_innocua_clip11262",0.2398917291,0,2350000),("listeria_welshimeri_serovar_6b_str._slcc5334",0.2397671622,0,2400000),("magnetospirillum_magneticum_amb-1",0.2396434591,0,7500000),("mannheimia_succiniciproducens_mbel55e",0.2395206091,0,3800000),("maricaulis_maris_mcs10",0.2393986014,0,3600000),("marinobacter_aquaeolei_vt8",0.2392774258,0,6000000),("mesoplasma_florum_l1",0.2391570721,0,2500000),("mesorhizobium_loti_maff303099",0.2390375301,0,2900000),("Chelativorans sp. BNC1?Mesorhizobium sp. BNC1",0.23891879,0,3300000),("methanobrevibacter_smithii_atcc_35061",0.2388008422,0,1900000),("methanocaldococcus_jannaschii_dsm_2661",0.2386836771,0,1200000),("Methanocella arvoryzae MRE50?uncultured methanogenic archaeon RC-I",0.224524087,0,4090000),("methanococcoides_burtonii_dsm_6242",0.2385672853,0,4000000),("methanococcus_aeolicus_nankai-3",0.2384516577,0,2840000),("methanococcus_maripaludis_c5",0.2383367852,0,2560000),("methanococcus_vannielii_sb",0.238222659,0,5350000),("methanocorpusculum_labreanum_z",0.2381092702,0,4600000),("methanoculleus_marisnigri_jr1",0.2379966104,0,9100000),("methanopyrus_kandleri_av19",0.237884671,0,1800000),("methanosaeta_thermophila_pt",0.2377734438,0,5200000),("methanosarcina_barkeri_str._fusaro",0.2376629206,0,4020000),("methanosarcina_mazei_go1",0.2375530934,0,5800000),("methanosphaera_stadtmanae_dsm_3091",0.2374439542,0,1500000),("methanospirillum_hungatei_jf-1",0.2373354953,0,3100000),("methanothermobacter_thermautotrophicus_str._delta_h",0.237227709,0,1770000),("methylobacillus_flagellatus_kt",0.2371205878,0,2870000),("methylobacterium_extorquens_pa1",0.2370141244,0,5300000),("methylobacterium_radiotolerans_jcm_2831",0.2369083114,0,1600000),("methylobacterium_sp._4-46",0.2368031416,0,1270000),("methylococcus_capsulatus_str._bath",0.2366986081,0,3700000),("microcystis_aeruginosa_nies-843",0.2365947038,0,2300000),("moorella_thermoacetica_atcc_39073",0.236491422,0,2150000),("Mycobacterium bovis BCG Pasteur 1173P2?Mycobacterium bovis BCG str. Pasteur 1173P2",0.2363887558,0,790000),("mycobacterium_gilvum_pyr-gck",0.2362866987,0,1040000),("mycobacterium_leprae_tn",0.2361852442,0,960000),("mycobacterium_marinum_m",0.2360843858,0,3700000),("mycobacterium_smegmatis_str._mc2_155",0.2359841173,0,4430000),("mycobacterium_sp._kms",0.2358844323,0,5500000),("mycobacterium_ulcerans_agy99",0.2357853247,0,3950000),("mycobacterium_vanbaalenii_pyr-1",0.2356867886,0,4900000),("mycoplasma_agalactiae_pg2",0.2355888179,0,5220000),("mycoplasma_capricolum_subsp._capricolum_atcc_27343",0.2354914069,0,2360000),("Mycoplasma gallisepticum str. R(low)?Mycoplasma gallisepticum R",0.2353945496,0,1180000),("mycoplasma_genitalium_g37",0.2352982405,0,8400000),("mycoplasma_mobile_163k",0.235202474,0,2490000),("mycoplasma_mycoides_subsp._mycoides_sc_str._pg1",0.2351072445,0,3640000),("mycoplasma_penetrans_hf-2",0.2350125465,0,1600000),("mycoplasma_pulmonis_uab_ctip",0.2349183748,0,9670000),("mycoplasma_synoviae_53",0.234824724,0,5690000),("myxococcus_xanthus_dk_1622",0.2347315889,0,4390000),("nanoarchaeum_equitans_kin4-m",0.2346389645,0,2300000),("natronomonas_pharaonis_dsm_2160",0.2345468456,0,5080000),("neisseria_gonorrhoeae_fa_1090",0.2344552272,0,3010000),("neisseria_meningitidis_mc58",0.2343641044,0,1400000),("Neorickettsia sennetsu strain Miyayama?Neorickettsia sennetsu str. Miyayama",0.2342734725,0,3240000),("nitrobacter_hamburgensis_x14",0.2341833266,0,4400000),("nitrobacter_winogradskyi_nb-255",0.2340936619,0,3600000),("nitrosococcus_oceani_atcc_19707",0.234004474,0,4280000),("nitrosomonas_europaea_atcc_19718",0.233915758,0,1200000),("nitrosomonas_eutropha_c91",0.2338275097,0,5500000),("nitrosopumilus_maritimus_scm1",0.2337397244,0,1400000),("nitrosospira_multiformis_atcc_25196",0.2336523979,0,1520000),("nocardia_farcinica_ifm_10152",0.2335655257,0,5500000),("nocardioides_sp._js614",0.2334791036,0,1800000),("nostoc_punctiforme_pcc_73102",0.2333931274,0,2180000),("nostoc_sp._pcc_7120",0.2333075929,0,4660000),("novosphingobium_aromaticivorans_dsm_12444",0.2332224961,0,2200000),("oceanobacillus_iheyensis_hte831",0.2331378327,0,1470000),("ochrobactrum_anthropi_atcc_49188",0.233053599,0,3300000),("oenococcus_oeni_psu-1",0.2329697908,0,1700000),("onion_yellows_phytoplasma_oy-m",0.2328864044,0,3730000),("opitutus_terrae_pb90-1",0.2328034359,0,1580000),("orientia_tsutsugamushi_str._boryong",0.2327208815,0,1600000),("parabacteroides_distasonis_atcc_8503",0.2326387374,0,5360000),("Parachlamydia-related symbiont UWE25?Candidatus Protochlamydia amoebophila UWE25",0.2628479382,0,2250000),("paracoccus_denitrificans_pd1222",0.232557,0,5260000),("parvibaculum_lavamentivorans_ds-1",0.2324756656,0,5630000),("pasteurella_multocida_subsp._multocida_str._pm70",0.2323947307,0,5200000),("pelobacter_carbinolicus_dsm_2380",0.2323141916,0,2030000),("petrotoga_mobilis_sj95",0.2321542871,0,3540000),("photobacterium_profundum_ss9",0.2320749147,0,2100000),("photorhabdus_luminescens_subsp._laumondii_tto1",0.2319959245,0,5310000),("picrophilus_torridus_dsm_9790",0.2319173131,0,1800000),("polaromonas_naphthalenivorans_cj2",0.2318390772,0,3100000),("polaromonas_sp._js666",0.2317612135,0,7780000),("Polynucleobacter necessarius subsp. asymbioticus QLW-P1DMWA-1?Polynucleobacter sp. QLW-P1DMWA-1",0.2316065901,0,8200000),("Polynucleobacter necessarius subsp. necessarius STIR1?Polynucleobacter necessarius STIR1",0.2316837189,0,5450000),("porphyromonas_gingivalis_atcc_33277",0.2315298241,0,2200000),("propionibacterium_acnes_kpa171202",0.2314534178,0,3200000),("pseudoalteromonas_atlantica_t6c",0.231377368,0,2060000),("Pseudoalteromonas haloplanktis str. TAC125?Pseudoalteromonas haloplanktis TAC125",0.2313016718,0,5000000),("pseudomonas_aeruginosa_pa7",0.2312263263,0,1410000),("Pseudomonas protegens Pf-5?Pseudomonas fluorescens Pf-5",0.2311513284,0,7790000),("pseudomonas_putida_kt2440",0.2310766752,0,4910000),("pseudomonas_stutzeri_a1501",0.2310023639,0,3090000),("pseudomonas_syringae_pv._tomato_str._dc3000",0.2309283917,0,3700000),("psychrobacter_arcticus_273-4",0.2308547558,0,4130000),("psychrobacter_cryohalolentis_k5",0.2307814533,0,4840000),("psychrobacter_sp._prwf-1",0.2307084815,0,160000),("psychromonas_ingrahamii_37",0.2306358378,0,4200000),("pyrobaculum_aerophilum_str._im2",0.2305635194,0,9090000),("pyrobaculum_arsenaticum_dsm_13514",0.2304915238,0,2730000),("pyrobaculum_calidifontis_jcm_11548",0.2304198483,0,4710000),("pyrobaculum_islandicum_dsm_4184",0.2303484903,0,2560000),("pyrococcus_abyssi_ge5",0.2302774472,0,1990000),("pyrococcus_furiosus_dsm_3638",0.2302067167,0,5510000),("pyrococcus_horikoshii_ot3",0.2301362961,0,5930000),("ralstonia_eutropha_jmp134",0.230066183,0,1580000),("Cupriavidus metallidurans CH34?Ralstonia metallidurans CH34",0.2299963749,0,2100000),("ralstonia_solanacearum_gmi1000",0.2299268695,0,3230000),("renibacterium_salmoninarum_atcc_33209",0.2298576644,0,3000000),("rhizobium_etli_cfn_42",0.2297887573,0,6530000),("rhizobium leguminosarum bv. viciae chromosome complete genome, strain 3841?Rhizobium leguminosarum bv. viciae 3841",0.2297201457,0,5600000),("rhodobacter_sphaeroides_2.4.1",0.2296518274,0,2200000),("Rhodococcus jostii RHA1?Rhodococcus sp. RHA1",0.2295838002,0,1700000),("rhodoferax_ferrireducens_t118",0.2295160617,0,2820000),("rhodopirellula_baltica_sh_1",0.2294486098,0,4100000),("rhodopseudomonas_palustris_bisa53",0.2293814422,0,2910000),("rhodospirillum_rubrum_atcc_11170",0.2293145568,0,2200000),("rickettsia_akari_str._hartford",0.2292479514,0,4410000),("rickettsia_bellii_rml369-c",0.229181624,0,1800000),("rickettsia_canadensis_str._mckiel",0.2291155723,0,4860000),("rickettsia_felis_urrwxcal2",0.2290497944,0,4610000),("rickettsia_massiliae_mtu5",0.2289842881,0,5690000),("rickettsia_prowazekii_str._madrid_e",0.2289190514,0,580000),("rickettsia_rickettsii_str._'sheila_smith'",0.2288540823,0,4800000),("rickettsia_typhi_str._wilmington",0.2287893789,0,1000000),("roseiflexus_castenholzii_dsm_13941",0.228724939,0,3540000),("roseiflexus_sp._rs-1",0.2286607609,0,6820000),("roseobacter_denitrificans_och_114",0.2285968426,0,3320000),("rubrobacter_xylanophilus_dsm_9941",0.2285331821,0,3400000),(	"Ruegeria pomeroyi DSS-3?Silicibacter pomeroyi DSS-3",0.2271942665,0,6990000),("saccharophagus_degradans_2-40",0.2284697775,0,3920000),("Saccharopolyspora erythraea NRRL2338?Saccharopolyspora erythraea NRRL 2338",0.2284066271,0,1200000),("salinibacter_ruber_dsm_13855",0.2283437289,0,1810000),("salinispora_arenicola_cns-205",0.2282810811,0,3290000),("salinispora_tropica_cnb-440",0.2282186819,0,4700000),("Salmonella enterica subsp. enterica serovar typhimurium str. LT2?Salmonella typhimurium LT2",0.2281565296,0,3000000),("serratia_proteamaculans_568",0.2280946223,0,1830000),("shewanella_amazonensis_sb2b",0.2280329583,0,7000000),("shewanella_baltica_os195",0.2279715358,0,4990000),("shewanella_denitrificans_os217",0.2279103532,0,4200000),("shewanella_frigidimarina_ncimb_400",0.2278494087,0,3300000),("shewanella_halifaxensis_haw-eb4",0.2277887006,0,4100000),("shewanella_loihica_pv-4",0.2277282274,0,4770000),("shewanella_pealeana_atcc_700345",0.2276679872,0,2200000),("shewanella_putrefaciens_cn-32",0.2276079786,0,2930000),("shewanella_sediminis_haw-eb3",0.2275481999,0,3100000),("shewanella_sp._ana-3",0.2274886494,0,3900000),("shewanella_woodyi_atcc_51908",0.2274293257,0,6740000),("shigella_boydii_cdc_3083-94",0.227370227,0,3590000),("shigella_dysenteriae_sd197",0.227311352,0,1900000),("shigella_flexneri_2a_str._301",0.227252699,0,3230000),("silicibacter_sp._tm1040",0.227136053,0,5370000),("sinorhizobium_medicae_wsm419",0.2270780571,0,4900000),("sinorhizobium_meliloti_1021",0.2270202772,0,4990000),("sphingomonas_wittichii_rw1",0.2269627119,0,860000),("sphingopyxis_alaskensis_rb2256",0.2269053597,0,4670000),("staphylococcus_aureus_subsp._aureus_usa300",0.2268482192,0,2260000),("staphylococcus_epidermidis_rp62a",0.226791289,0,4800000),("staphylococcus_haemolyticus_jcsc1435",0.2267345677,0,7600000),("staphylococcus_saprophyticus_subsp._saprophyticus_atcc_15305",0.2266780539,0,9770000),("staphylothermus_marinus_f1",0.2266217461,0,5700000),("streptococcus_agalactiae_nem316",0.2265656431,0,2000000),("streptococcus_gordonii_str._challis_substr._ch1",0.2265097435,0,1140000),("streptococcus_mutans_ua159",0.2264540459,0,4270000),("streptococcus_pneumoniae_cgsp14",0.226398549,0,6540000),("streptococcus_pyogenes_mgas10750",0.2263432515,0,1700000),("streptococcus_sanguinis_sk36",0.2262881521,0,1800000),("streptococcus_suis_05zyh33",0.2262332495,0,1800000),("streptococcus_thermophilus_lmd-9",0.2261785424,0,2400000),("streptomyces_avermitilis_ma-4680",0.2261240296,0,1670000),("streptomyces_coelicolor_a3(2)",0.2260697098,0,930000),("sulfolobus_acidocaldarius_dsm_639",0.2260155818,0,790000),("sulfolobus_solfataricus_p2",0.2259616443,0,5280000),("sulfurihydrogenibium_sp._yo3aop1",0.2259078961,0,5620000),("sulfurimonas_denitrificans_dsm_1251",0.225854336,0,8700000),("sulfurovum_sp._nbc37-1",0.2258009628,0,1300000),("symbiobacterium_thermophilum_iam_14863",0.2257477753,0,1000000),("synechococcus_elongatus_pcc_7942",0.2256947724,0,2410000),("synechococcus_sp._pcc_7002",0.2256419528,0,4400000),("syntrophomonas_wolfei_subsp._wolfei_str._goettingen",0.2255893155,0,3050000),("thermoanaerobacter_pseudethanolicus_atcc_33223",0.2255368592,0,1360000),("thermoanaerobacter_sp._x514",0.2254845829,0,2400000),("thermoanaerobacter_tengcongensis_mb4",0.2254324854,0,2750000),("thermobifida_fusca_yx",0.2253805656,0,4100000),("thermofilum_pendens_hrk_5",0.2253288225,0,1900000),("thermoplasma_acidophilum_dsm_1728",0.2252772548,0,4250000),("thermoplasma_volcanium_gss1",0.2252258616,0,6920000),("thermoproteus_neutrophilus_v24sta",0.2251746417,0,4970000),("thermosipho_melanesiensis_bi429",0.2251235941,0,5170000),("thermosynechococcus_elongatus_bp-1",0.2250727177,0,1080000),("thermotoga_lettingae_tmo",0.2250220116,0,4820000),("thermotoga_petrophila_rku-1",0.2249714746,0,4230000),("thermotoga_sp._rq2",0.2249211057,0,6290000),("thermus_thermophilus_hb27",0.2248709039,0,4660000),("thiobacillus_denitrificans_atcc_25259",0.2248208683,0,9090000),("thiomicrospira_crunogena_xcl-2",0.2247709977,0,5400000),("treponema_denticola_atcc_35405",0.2247212913,0,4800000),("treponema_pallidum_subsp._pallidum_str._nichols",0.2246717479,0,1600000),("trichodesmium_erythraeum_ims101",0.2246223667,0,3370000),("tropheryma_whipplei_str._twist",0.2245731467,0,1760000),("ureaplasma_parvum_serovar_3_str._atcc_27815",0.2244751865,0,1940000),("verminephrobacter_eiseniae_ef01-2",0.2244264443,0,2480000),("vibrio_cholerae_o395",0.2243778595,0,3020000),("vibrio_fischeri_es114",0.2243294313,0,2100000),("vibrio_harveyi_atcc_baa-1116",0.2242811585,0,3640000),("vibrio_parahaemolyticus_rimd_2210633",0.2242330405,0,2410000),("vibrio_vulnificus_yj016",0.2241850761,0,3000000),("wolbachia_endosymbiont_of_drosophila_melanogaster",0.2241372646,0,2580000),("wolinella_succinogenes_dsm_1740",0.2240896051,0,3600000),("xanthobacter_autotrophicus_py2",0.2240420967,0,2800000),("xanthomonas_axonopodis_pv._citri_str._306",0.2239947385,0,1700000),("xanthomonas_oryzae_pv._oryzae_maff_311018",0.2239475296,0,1200000),("xylella_fastidiosa_9a5c",0.2239004692,0,2100000),("yersinia_enterocolitica_subsp._enterocolitica_8081",0.2238535564,0,1590000),("yersinia_pseudotuberculosis_ip_31758",0.2238067905,0,2600000),("zymomonas_mobilis_subsp._mobilis_zm4",0.2237601705,0,1400000)])
	elif dataset == 'simLC':
		truth.set_by_array([("Actinobacillus_succinogenes_130Z",0.49,252,2319663),("Alkalilimnicola_ehrlichii_MLHE-1?Alkalilimnicola_ehrlichei_MLHE-1",0.52,373,3275944),("Alkaliphillus_metalliredigenes_UNDEF",0.45,489,4929566),("Anabaena_variabilis_ATCC_29413",0.61,855,6365727),("Anaeromyxobacter_dehalogenans_2CP-C",0.53,584,5013479),("Arthrobacter_sp._FB24",0.55,570,4698945),("Azotobacter_vinelandii_AvOP",0.55,650,5365318),("Bacillus_cereus_NVH391-98",0.58,520,4087024),("Bifidobacterium_longum_DJO10A",0.55,288,2375792),("Bradyrhizobium_BTAi1,Bradyrhizobium_sp._BTAi1",5.08,9277,8264687),("Brevibacterium_linens_BL2",0.56,542,4367044),("Burkholderia_ambifaria_AMMD",0.58,955,7484988),("Burkholderia_cenocepacia_AU_1054",0.55,879,7279118),("Burkholderia_cenocepacia_HI2424",0.57,956,7537985),("Burkholderia_sp._sp.strain_383?Burkholderia_sp._383",1.35,1074,3587082),("Burkholderia_vietnamiensis_G4",0.61,992,7305582),("Burkholderia_xenovorans_LB40?Burkholderia_xenovorans_LB400",0.53,1149,9731140),("Caldicellulosiruptor_accharolyticus_UNDEF",0.56,367,2970275),("Chlorobium_limicola_DSMZ_245T",0.62,381,2763181),("Chlorobium_phaeobacteroides_DSM_266",0.52,359,3133902),("Chlorobium_vvibrioforme_f._thiosulfatophilum_DSMZ_265T",0.46,199,1966858),("Chloroflexus_aurantiacus_J-10-fl",0.58,679,5258541),("Chromohalobacter_salexigens_DSM3043",0.56,454,3696649),("Clostridium_beijerincki_NCIMB_8052",0.56,737,6000632),("Clostridium_thermocellum_ATCC_27405",0.54,461,3843301),("Crocosphaera_watsonii_WH_8501",0.59,812,6238478),("Cytophaga_hutchinsonii_ATCC_33406",5.27,5168,4433218),("Dechloromonas_aromatica_RCB",0.54,537,4501104),("Deinococcus_geothermalis_DSM_11300?Deinococcus_geothermalis_DSM11300",0.76,415,2467205),("Desulfitobacterium_hafniense_DCB-2",0.66,769,5279134),("Desulfovibrio_desulfuricans_G20",0.59,484,3730232),("Ehrlichia_canis_Jake",0.67,196,1315030),("Ehrlichia_chaffeensis_sapulpa",0.67,150,1005936),("Enterococcus_faecium_DO",0.6,359,2698137),("Exiguobacterium_UNDEF_255-15",0.56,377,3034136),("Ferroplasma_acidarmanus_fer1",0.55,238,1947953),("Frankia_sp._CcI3",0.54,645,5433628),("Frankia_sp._EAN1pec",0.56,1109,8982042),("Geobacter_metallireducens_GS-15",0.58,515,3997420),("Haemophilus_somnus_129PT",0.52,232,2007700),("Jannaschia_sp._CCS1",0.57,543,4317977),("Kineococcus_radiotolerans_SRS30216",0.54,566,4761183),("Lactobacillus_brevis_ATCC_367",0.35,177,2291220),("Lactobacillus_casei_ATCC_334",0.57,362,2895264),("Lactobacillus_delbrueckii_bulgaricus_ATCC_BAA-365",0.48,195,1856951),("Lactobacillus_gasseri_ATCC_33323",0.58,244,1894360),("Lactococcus_lactis_cremoris_SK11",0.56,301,2438589),("Leuconostoc_mesenteroides_mesenteroides_ATCC_8293",0.52,235,2038396),("Magnetococcus_sp._MC-1",0.48,504,4719581),("Marinobacter_aquaeolei_VT8",0.57,547,4326849),("Mesorhizobium_sp._BNC1",0.58,567,4412446),("Methanococcoides_burtonii_DSM6242",0.47,268,2575032),("Methanosarcina_barkeri_Fusaro",0.51,545,4837408),("Methanospirillum_hungatei_JF-1",0.55,429,3544738),("Methylobacillus_flagellatus_strain_KT",0.56,365,2971517),("Moorella_thermoacetica_ATCC_39073",1.16,674,2628784),("Nitrobacter_hamburgensis_UNDEF",0.65,630,4406967),("Nitrobacter_winogradskyi_Nb-255",0.57,427,3402093),("Nitrosococcus_oceani_UNDEF",0.53,409,3481691),("Nitrosomonas_eutropha_C71",0.53,314,2661057),("Nitrosospira_multiformis_ATCC_25196",0.54,378,3184243),("Nocardioides_sp._JS614",0.58,636,4985871),("Novosphingobium_aromaticivorans_DSM_12444_F199",0.66,520,3561584),("Oenococcus_oeni_PSU-1",0.46,182,1780517),("Paracoccus_denitrificans_PD1222",0.58,585,4582380),("Pediococcus_pentosaceus_ATCC_25745",0.54,217,1832387),("Pelobacter_carbinolicus_DSM_2380",0.6,489,3665893),("Pelobacter_propionicus_DSM_2379",0.57,508,4008000),("Pelodictyon_luteolum_UNDEF",0.48,250,2364842),("Pelodictyon_phaeoclathratiforme_BU-1_DSMZ_5477T",0.6,402,3018238),("Polaromonas_sp._JS666",0.64,733,5200264),("Prochlorococcus_marinus_str._MIT_9312",0.48,183,1709204),("Prochlorococcus_sp._NATL2A",0.62,253,1842899),("Prosthecochloris_aestuarii_SK413/DSMZ_271t",0.51,282,2512923),("Prosthecochloris_sp._BS1",0.8,483,2736403),("Pseudoalteromonas_atlantica_T6c",0.51,588,5187005),("Pseudomonas_fluorescens_PfO-1",0.51,730,6438405),("Pseudomonas_putida_F1",0.51,675,5959964),("Pseudomonas_syringae_B728a",0.55,746,6093698),("Psychrobacter_arcticum_273-4",0.56,327,2650701),("Psychrobacter_cryopegella_UNDEF",0.62,422,3059876),("Rhodobacter_sphaeroides_2.4.1",0.56,514,4131626),("Rhodoferax_ferrireducens_UNDEF",0.58,599,4712337),("Rhodopseudomonas_palustris_BisA53",0.52,636,5505494),("Rhodopseudomonas_palustris_BisB18",0.57,699,5513844),("Rhodopseudomonas_palustris_BisB5",0.53,575,4892717),("Rhodopseudomonas_palustris_HaA2",24.49,28861,5331656),("Rhodospirillum_rubrum_ATCC_11170",0.58,559,4352825),("Rubrobacter_xylanophilus_DSM_9941",0.57,409,3225748),("Saccharophagus_degradans_2-40",0.52,582,5057531),("Shewanella_amazonensis_SB2B",0.56,536,4306142),("Shewanella_baltica_OS155",0.55,621,5127376),("Shewanella_frigidimarina_NCMB400",0.51,551,4845257),("Shewanella_putefaciens_UNDEF",0.55,565,4659220),("Shewanella_sp._ANA-3",0.6,664,4972204),("Shewanella_sp._MR-7",0.54,568,4792610),("Shewanella_sp._PV-4",0.52,524,4602594),("Shewanella_sp._W3-18-1",0.51,533,4708380),("Silicibacter_sp._TM1040",0.66,469,3200938),("Sphingopyxis_alaskensis_RB2256",0.59,438,3345170),("Streptococcus_suis_89/1591",0.56,263,2143334),("Streptococcus_thermophilus_LMD-9",0.43,178,1856368),("Synechococcus_sp._PCC_7942_elongatus",0.53,316,2695903),("Syntrophobacter_fumaroxidans_MPOB",0.55,606,4990251),("Syntrophomonas_wolfei_Goettingen",0.48,314,2936195),("Thermoanaerobacter_ethanolicus_39E",0.6,315,2362816),("Thermobifida_fusca_YX",0.54,434,3642249),("Thiobacillus_denitrificans_ATCC_25259",0.61,395,2909809),("Thiomicrospira_crunogena_XCL-2",0.51,274,2427734),("Thiomicrospira_denitrificans_ATCC_338890?Thiomicrospira_denitrificans_ATCC_33889",0.57,277,2201561),("Trichodesmium_erythraeum_IMS101",0.57,977,7750108),("Xylella_fastidiosa_Ann-1",2.51,2836,5115560),("Xylella_fastidiosa_Dixon",1.04,601,2622359)])
	else:
		print "Dataset not recognized. Input format is <filename> <dataset>."
		return 0

	truth.clean_names()
	#true_species = [s.lower().translate(string.maketrans("",""),bad_punct).replace(" ","_") for s in truth[0]]

	true_species_alone,true_species_ab,true_genus_alone,true_genus_ab = collapse_strains(truth.species,truth.abundance)

	return truth,true_species_alone,true_species_ab,true_genus_alone,true_genus_ab

def calc_error(truth,est):
	""" Calculate errors between true and estimated abundances """

	adjusted_abundance = numpy.zeros((truth.lookup_length()))

	if truth.species != est.species:
		if truth.lookup_length() == est.lookup_length():
			adjusted_abundance = est.abundance
			for ind,sp in enumerate(est.species):
				if not (est.species[ind] == truth.species[ind] or est.species[ind] == truth.species[ind].partition(',')[0] or est.species[ind] == truth.species[ind].partition(',')[2]): # does not match species or a listed synonym
					print "Error: species match not found between estimated {0} and actual {1}. Pretending they are the same for the sake of further calculations.".format(est.species[ind],truth.species[ind])
		else:
			print "Mapped samples don't match actual samples"
			for ind,sp in enumerate(truth.species):
				try:
					sp_match = est.species.index(sp)
				except ValueError: # most likely due to synonym
					try:
						sp_match = est.species.index(truth.species[ind].partition(',')[0])
					except ValueError:
						try:
							sp_match = est.species.index(truth.species[ind].partition(',')[2])
						except ValueError: # now we've really failed
							print "Error: species match for {0} not found in estimated data. Setting abundance and counts to 0 for {0}.".format(sp)
							adjusted_abundance[ind] = 0
						else:
							adjusted_abundance[ind] = est.abundance[sp_match]
					else:
						adjusted_abundance[ind] = est.abundance[sp_match]
				else:
					adjusted_abundance[ind] = est.abundance[sp_match]
	else:
		adjusted_abundance = est.abundance

	print "Total number of samples: {0}".format(len(adjusted_abundance))
	diff = 100*(numpy.array(adjusted_abundance) - numpy.array(truth.abundance))/numpy.array(truth.abundance)
	print "Average relative error: {0}".format(numpy.mean(abs(diff)))
	diff_sq = [d*d/10000 for d in diff]
	print "Relative root mean squared error: {0}".format(numpy.mean(diff_sq) ** (0.5) * 100)

	return diff,adjusted_abundance

def graph_error(true_species, true_abundance, est_species, est_abundance,
				adjusted_abundance, diff, expname, tier, showgraphs=True):
	# These imports are here so script can run on server w/out graphics
	import lanthplot

	true_sp = [x.replace('_',' ') for x in true_species]
	all_species = list(set(est_species + true_species))

	filtered_ab = [r for r in est_abundance if r != 0] # exclude 0's from est_abundance to get a more sensible mean
	mean_ab = sum(filtered_ab)/len(filtered_ab)

	# pickle raw diffs for all species
	all_est = [est_abundance[est_species.index(sp)] if sp in est_species else 0 for sp in all_species]
	all_true = [true_abundance[true_species.index(sp)] if sp in true_species else 0 for sp in all_species]
	all_diff = []
	for i,a in enumerate(all_est):
		try:
			all_diff.append(100*(a - all_true[i]) / max(a,all_true[i]))
		except ZeroDivisionError:
			all_diff.append(0) # if both the estimate and the actual are 0, we're good here
	pickle.dump(zip(all_species,all_diff),open(expname+"_diffs.pickle",'w'))

	# graph true abundances
	if len(all_species) == len(true_species):
		xmax = len(true_species)
		x = numpy.array(range(0,xmax))

		ab_filter = zip(true_sp,true_abundance,adjusted_abundance)
		ab_filter.sort( key=lambda x: x[1],reverse=True )
		ab_species,ab_true,ab_adjusted = zip(*ab_filter)

		lanthplot.plot_setup_pre("Graph of true {1}-level abundances in {0}"
			.format(expname,tier), xlabels = ab_species, xticks = range(0,xmax), xrotation = -90)
		lanthplot.plot(x, ab_true, color='blue', label='True')
		lanthplot.plot(x,ab_adjusted, color='red', label='Estimated')
		lanthplot.plot_setup_post(save_file = expname +'_'+ tier +'_abundances.png')

	# graph abundances for all species, not just the true ones, if above mean
	else:
		print "Filtering out any estimated results under the mean abundance of {0}".format(mean_ab)
		present_true = []
		present_est = []
		present_species = []
		for i,sp in enumerate(all_species):
			if (sp in true_species) or (sp in est_species and est_abundance[est_species.index(sp)]>mean_ab): #only include species if it has a non-zero estimate or non-zero actual abundance
				present_species.append(sp)
				try:
					present_est.append(est_abundance[est_species.index(sp)])
				except:
					present_est.append(0)
				try:
					present_true.append(true_abundance[true_species.index(sp)])
				except:
					present_true.append(0)

		xmax = len(present_species)
		present_sp = [x.replace('_',' ') for x in present_species]
		x = numpy.array(range(0,xmax))

		# sort based on true abundance
		all_filter = zip(present_sp,present_true,present_est)
		all_filter.sort( key=lambda x: x[1],reverse=True )
		fil_sp,fil_true,fil_est = zip(*all_filter)

		lanthplot.plot_setup_pre(
			"Graph of all above-average estimated abundances in {0} at {1}-level"
			.format(expname,tier), xlabels = fil_sp, xticks = range(0,xmax),
			xrotation = -90)
		lanthplot.plot(x, fil_true, color='blue', label="True")
		lanthplot.plot(x, fil_est, color='red', label="Estimated")
		lanthplot.plot_setup_post(save_file = expname +'_'+ tier +'_sigabundances.png')

		# sort based on name
		all_filter = zip(present_sp,present_true,present_est)
		all_filter.sort( key=lambda x: x[0],reverse=True )
		fil_sp,fil_true,fil_est = zip(*all_filter)

		lanthplot.plot_setup_pre(
			"Graph of all above-average estimated abundances in {0} at {1}-level"
			.format(expname,tier), xlabels = fil_sp, xticks = range(0,xmax),
			xrotation = -90)
		lanthplot.plot(x, fil_true, color='blue', label="True")
		lanthplot.plot(x, fil_est, color='red', label="Estimated")
		lanthplot.plot_setup_post(save_file = expname +'_'+ tier +'_sigabundances2.png')

	# graph diffs
	if len(all_species) == len(true_species):
		diff_combo = zip(true_sp, diff)
		diff_filter = diff_combo
		diff_filter.sort( key=lambda x: x[1], reverse=True )
		try:
			diff_species,diff_ab = zip(*diff_filter)
		except:
			pass
		else:
			xmax = len(diff_species)
			if xmax < 200:
				x = numpy.array(range(0, xmax))
				lanthplot.plot_setup_pre(
					"Graph of % error in present species in {0}-level {1}"
					.format(tier, expname), xlabels = diff_species, xticks = range(0,xmax),
					xrotation = -90)
				lanthplot.plot(x, diff_ab, plot_type='bar', width=0.8, color='purple')
				lanthplot.plot_setup_post(save_file = expname +'_'+ tier +'_errors.png',legend=False)

	else:
		total_diff = []
		for i,a in enumerate(present_est):
			total_diff.append(100*(a-present_true[i])/max(a,present_true[i]))

		diff_combo = zip(present_sp,total_diff)
		diff_filter = diff_combo
		diff_filter.sort( key=lambda x: x[1],reverse=True )
		try:
			diff_species,diff_ab = zip(*diff_filter)
		except:
			pass
		else:
			xmax = len(diff_species)
			if xmax < 200:
				x = numpy.array(range(0,xmax))

				lanthplot.plot_setup_pre(
					"Graph of % error in all true or estimated species in {0}-level {1}"
					.format(tier, expname), xlabels = diff_species, xticks = range(0,xmax),
					xrotation = -90)
				lanthplot.plot(x, diff_ab, plot_type='bar', width=0.8, color='purple')
				lanthplot.plot_setup_post(save_file = expname +'_'+ tier +'_allerrors.png',legend=False)

	return

def main(argv=sys.argv):
	"""
	Command line usage: python M_compare_metagenomic_results_to_truth.py
						[filename] [dataset] <fragmented data? (defaults to false)
						<show graphs? (defaults to true)>
	"""

	filename = argv[1]
	exp_name = filename.rpartition('_')[0] # will be used for graph-naming purposes
	dataset = argv[2] #i100 or simLC

	if len(argv) > 3 and (argv[3].lower() == 'true' or argv[3].lower() == 't'):
		fragmented = True
	else:
		fragmented = False

	if len(argv) > 4 and (argv[4].lower() == 'false' or argv[4].lower() == 'f'):
		show_graphs = False
	else:
		show_graphs = True

	truth,true_species_alone,true_species_ab,true_genus_alone,true_genus_ab = dataset_truth(dataset)
	est,est_species_alone,est_species_ab,est_genus_alone,est_genus_ab = process_input(filename,truth.size)

	print "Strain-level error:"
	diff, adjusted_abundance = calc_error(truth,est)
	if show_graphs:
		graph_error(true_species, true_abundance, est_species, est_abundance,
					adjusted_abundance, diff, exp_name, 'strain', fragmented)
	'''
	print "Species-level error:"
	diff, adjusted_abundance = calc_error(true_species_alone, true_species_ab,
											est_species_alone, est_species_ab)
	if show_graphs:
		graph_error(true_species_alone, true_species_ab, est_species_alone,
					est_species_ab, adjusted_abundance, diff, exp_name, 'species', fragmented)

	print "Genus-level error:"
	diff, adjusted_abundance = calc_error(true_genus_alone, true_genus_ab,
											est_genus_alone, est_genus_ab)
	if show_graphs:
		graph_error(true_genus_alone, true_genus_ab, est_genus_alone,
					est_genus_ab, adjusted_abundance, diff, exp_name,'genus', fragmented)
	'''

if __name__ == "__main__":
    main()
