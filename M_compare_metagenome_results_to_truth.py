import sys
import csv
import numpy as np

np.set_printoptions(precision=4)

def process_input(filename,size):
	suffix = filename.rsplit('.')[1] #xprs or txt file

	input_file = open(filename,'r')
	input_csv = csv.reader(input_file, 'excel-tab')
	input_data = [r for r in input_csv]
	input_data = input_data[1:] #remove header row

	if suffix == 'xprs': #express file
		input_data.sort(key=lambda x:x[1]) # sort alphabetically by species
		species = zip(*input_data)[1] # express species names are in second column
		est_counts = np.array([float(i) for i in zip(*input_data)[6]]) # used when dataset abundances are given in raw counts
		fpkm = np.array([float(i) for i in zip(*input_data)[10]]) # used when dataset abundances are given as percentages
		fpkm_total = np.sum(fpkm)
		est_abundance = fpkm*100/fpkm_total
	elif suffix == 'txt': #gasic file
		input_data.sort(key=lambda x:x[0]) # sort alphabetically by species
		species = zip(*input_data)[0] # gasic species names are in first column
		est_counts = np.array([float(i) for i in zip(*input_data)[2]]) # used when dataset abundances are given in raw counts

		est_abundance = est_counts/(np.sum(est_counts)*np.array(size)) # used when dataset abundances are given as percentages
		est_abundance = 100*est_abundance/np.sum(est_abundance) # normalize

	# Force species to be lowercase and underscored
	species = [s.lower().replace(" ","_") for s in species]

	return species,est_abundance,est_counts

def dataset_truth(dataset):
	# truth is in format |species|abundance|counts|genome size| where undefined counts is 0
	# and species may include alternate species names separated by a comma

	if dataset == 'i100':
		truth = np.array(zip(*[("acinetobacter_baumannii_sdf",0.9918569894,0,3477999),("alkalilimnicola_ehrlichii_mlhe-1",0.9534103189,0,3275944),("alkaliphilus_metalliredigens_qymf",1.0337468099,0,4929566),("anabaena_variabilis_atcc_29413",0.8701233754,0,7068604),("bacillus_anthracis_str._ames",1.0629349921,0,5227293),("bacillus_cereus_atcc_10987",0.9077576661,0,5432653),("bacillus_cereus_atcc_14579",1.4619387586,0,5427084),("bacillus_cereus_e33l",0.9147432601,0,5843240),("bacillus_clausii_ksm-k16",0.9028911617,0,4303871),("bacillus_halodurans_c-125",0.8847367596,0,4202352),("bacillus_subtilis_subsp._subtilis_str._168",0.8765090628,0,4214547),("bacillus_thuringiensis_str._al_hakam",0.9505918962,0,5313031),("bacteroides_thetaiotaomicron_vpi-5482",0.8661600805,0,6293400),("bordetella_bronchiseptica_rb50",1.1977086557,0,5339179),("borrelia_garinii_pbi",1.1270363889,0,986916),("buchnera_aphidicola_bcc",0.8652027125,0,422435),("burkholderia_phymatum_stm815",1.0474017422,0,8676565),("burkholderia_pseudomallei_668",0.8953944405,0,7040404),("campylobacter_concisus_13826",1.0403648471,0,2099415),("candidatus_blochmannia_pennsylvanicus_str._bpen",0.9129394805,0,791654),("chlamydia_trachomatis_l2b/uch-1/proctitis",0.9593251688,0,1038863),("chlamydophila_caviae_gpic",0.8711493582,0,1181357),("chlamydophila_pneumoniae_cwl029",0.8810983686,0,1230230),("chlorobium_luteolum_dsm_273",0.9264853736,0,2364842),("chromobacterium_violaceum_atcc_12472",0.9624326898,0,4751080),("clostridium_beijerinckii_ncimb_8052",0.9656484821,0,6000632),("clostridium_botulinum_b1_str._okra",0.9478594593,0,4107014),("clostridium_novyi_nt",1.1580488951,0,2547720),("clostridium_thermocellum_atcc_27405",0.9330381139,0,3843301),("corynebacterium_efficiens_ys-314",1.0275058892,0,3219507),("corynebacterium_urealyticum_dsm_7109",1.0008453082,0,2369219),("coxiella_burnetii_rsa_493",0.9184749477,0,2032675),("cyanothece_sp._atcc_51142",0.9401343144,0,5460382),("cytophaga_hutchinsonii_atcc_33406",0.9563196367,0,4433218),("erwinia_tasmaniensis",1.1018606017,0,4067869),("escherichia_coli_apec_o1",0.9061021819,0,5497655),("escherichia_coli_e24377a",0.9244094549,0,5249294),("escherichia_coli_str._k-12_substr._mg1655",0.9836152016,0,4639637),("escherichia_coli_str._k-12_substr._w3110",1.2861538913,0,4646332),("francisella_tularensis_subsp._tularensis_fsc198",0.8799240389,0,1892616),("geobacter_sulfurreducens_pca",0.8822915619,0,3814128),("haemophilus_influenzae_rd_kw20",0.8925873316,0,1830138),("haemophilus_somnus_129pt",0.9353392713,0,2012879),("haloquadratum_walsbyi_dsm_16790",0.9724329129,0,3179362),("herpetosiphon_aurantiacus_dsm_785",0.8671305022,0,6785432),("hyphomonas_neptunium_atcc_15444",0.9797413592,0,3705021),("idiomarina_loihiensis_l2tr",0.9111749381,0,2839318),("ignicoccus_hospitalis_kin4/i",0.9286138232,0,1297538),("lactobacillus_delbrueckii_subsp._bulgaricus_atcc_baa-365",1.1765523477,0,1856951),("lactobacillus_salivarius_ucc118",0.9013332759,0,2133980),("lactococcus_lactis_subsp._cremoris_sk11",0.9377036228,0,2598353),("lactococcus_lactis_subsp._lactis_il1403",1.3855635178,0,2365589),("lawsonia_intracellularis_phe/mn1-00",0.9204063033,0,1719017),("leptospira_borgpetersenii_serovar_hardjo-bovis_jb197",0.8633259013,0,3876236),("listeria_welshimeri_serovar_6b_str._slcc5334",0.8859899334,0,2814130),("methanococcus_maripaludis_c7",0.9452084692,0,1772694),("mycobacterium_avium_subsp._paratuberculosis_k-10",0.8872643124,0,4829781),("mycobacterium_bovis_af2122/97",0.8998056707,0,4345492),("mycobacterium_marinum_m",0.8776298824,0,6660145),("mycobacterium_sp._jls",1.7694735192,0,6048425),("mycobacterium_tuberculosis_f11",0.9044804156,0,4424435),("neisseria_gonorrhoeae_fa_1090",0.9876497049,0,2153922),("neisseria_meningitidis_mc58",0.9307972132,0,2272360),("nostoc_punctiforme_pcc_73102",0.8983073115,0,9059196),("ochrobactrum_anthropi_atcc_49188",1.0160158194,0,5205782),("prochlorococcus_marinus_str._mit_9313",1.0107078939,0,2410873),("prochlorococcus_marinus_str._natl1a",0.9689793951,0,1864731),("pseudomonas_putida_w619",1.2512067245,0,5774330),("pseudomonas_stutzeri_a1501",0.8912213302,0,4567418),("psychrobacter_cryohalolentis_k5",2.2276838861,0,3101098),("psychromonas_ingrahamii_37",0.8691118092,0,4559598),("ralstonia_eutropha_h16",0.8885605469,0,3364647),("rhodococcus_jostii_rha1",0.8835041681,0,9702740),("rhodopseudomonas_palustris_bisa53",1.0715558086,0,5505494),("rhodopseudomonas_palustris_cga009",1.0056581462,0,5476069),("rickettsia_prowazekii_str._madrid_e",0.8939780971,0,1111523),("salmonella_enterica_subsp._enterica_serovar_paratyphi_b_str._spb7",1.113841943,0,4858887),("shewanella_sp._mr-7",0.896837214,0,4799110),("sinorhizobium_meliloti_1021",0.8743176015,0,6691696),("sodalis_glossinidius_str._'morsitans'",0.9426347213,0,4292505),("staphylococcus_aureus_subsp._aureus_jh1",1.1416733976,0,2936937),("staphylococcus_aureus_subsp._aureus_mrsa252",0.8754051191,0,2902619),("staphylococcus_aureus_subsp._aureus_mu3",0.8642580865,0,2880168),("staphylococcus_aureus_subsp._aureus_mw2",1.0909113027,0,2820462),("staphylococcus_aureus_subsp._aureus_nctc_8325",0.8787680456,0,2821361),("staphylococcus_aureus_subsp._aureus_str._newman",0.9760172336,0,2878897),("streptococcus_agalactiae_2603v/r",1.054906113,0,2160267),("streptococcus_pneumoniae_r6",1.2222455073,0,2038615),("streptococcus_pyogenes_mgas315",0.8732460774,0,1900521),("streptococcus_pyogenes_mgas5005",1.0216060126,0,1838554),("synechococcus_elongatus_pcc_6301",0.8681143003,0,2696255),("synechococcus_elongatus_pcc_7942",1.5752103822,0,2796473),("syntrophomonas_wolfei_subsp._wolfei_str._goettingen",0.9094481382,0,2936195),("thermosipho_melanesiensis_bi429",0.9223838063,0,1915238),("thermus_thermophilus_hb8",0.9962505203,0,2197216),("thiobacillus_denitrificans_atcc_25259",1.0808489742,0,2909809),("treponema_denticola_atcc_35405",0.8898793162,0,2843201),("ureaplasma_parvum_serovar_3_str._atcc_700970",0.87219013,0,751719),("wigglesworthia_glossinidia_endosymbiont_of_glossina_brevipalpis",1.3295506979,0,703005),("yersinia_pestis_co92",0.9165878579,0,4829858)]))
	elif dataset == 'simLC':
		truth = np.array(zip(*[("Actinobacillus_succinogenes_130Z",0.49,252,2319663),("Alkalilimnicola_ehrlichii_MLHE-1,Alkalilimnicola_ehrlichei_MLHE-1",0.52,373,3275944),("Alkaliphillus_metalliredigenes_UNDEF",0.45,489,4929566),("Anabaena_variabilis_ATCC_29413",0.61,855,6365727),("Anaeromyxobacter_dehalogenans_2CP-C",0.53,584,5013479),("Arthrobacter_sp._FB24",0.55,570,4698945),("Azotobacter_vinelandii_AvOP",0.55,650,5365318),("Bacillus_cereus_NVH391-98",0.58,520,4087024),("Bifidobacterium_longum_DJO10A",0.55,288,2375792),("Bradyrhizobium_BTAi1,Bradyrhizobium_sp._BTAi1",5.08,9277,8264687),("Brevibacterium_linens_BL2",0.56,542,4367044),("Burkholderia_ambifaria_AMMD",0.58,955,7484988),("Burkholderia_cenocepacia_AU_1054",0.55,879,7279118),("Burkholderia_cenocepacia_HI2424",0.57,956,7537985),("Burkholderia_sp._sp.strain_383,Burkholderia_sp._383",1.35,1074,3587082),("Burkholderia_vietnamiensis_G4",0.61,992,7305582),("Burkholderia_xenovorans_LB40,Burkholderia_xenovorans_LB400",0.53,1149,9731140),("Caldicellulosiruptor_accharolyticus_UNDEF",0.56,367,2970275),("Chlorobium_limicola_DSMZ_245T",0.62,381,2763181),("Chlorobium_phaeobacteroides_DSM_266",0.52,359,3133902),("Chlorobium_vvibrioforme_f._thiosulfatophilum_DSMZ_265T",0.46,199,1966858),("Chloroflexus_aurantiacus_J-10-fl",0.58,679,5258541),("Chromohalobacter_salexigens_DSM3043",0.56,454,3696649),("Clostridium_beijerincki_NCIMB_8052",0.56,737,6000632),("Clostridium_thermocellum_ATCC_27405",0.54,461,3843301),("Crocosphaera_watsonii_WH_8501",0.59,812,6238478),("Cytophaga_hutchinsonii_ATCC_33406",5.27,5168,4433218),("Dechloromonas_aromatica_RCB",0.54,537,4501104),("Deinococcus_geothermalis_DSM_11300,Deinococcus_geothermalis_DSM11300",0.76,415,2467205),("Desulfitobacterium_hafniense_DCB-2",0.66,769,5279134),("Desulfovibrio_desulfuricans_G20",0.59,484,3730232),("Ehrlichia_canis_Jake",0.67,196,1315030),("Ehrlichia_chaffeensis_sapulpa",0.67,150,1005936),("Enterococcus_faecium_DO",0.6,359,2698137),("Exiguobacterium_UNDEF_255-15",0.56,377,3034136),("Ferroplasma_acidarmanus_fer1",0.55,238,1947953),("Frankia_sp._CcI3",0.54,645,5433628),("Frankia_sp._EAN1pec",0.56,1109,8982042),("Geobacter_metallireducens_GS-15",0.58,515,3997420),("Haemophilus_somnus_129PT",0.52,232,2007700),("Jannaschia_sp._CCS1",0.57,543,4317977),("Kineococcus_radiotolerans_SRS30216",0.54,566,4761183),("Lactobacillus_brevis_ATCC_367",0.35,177,2291220),("Lactobacillus_casei_ATCC_334",0.57,362,2895264),("Lactobacillus_delbrueckii_bulgaricus_ATCC_BAA-365",0.48,195,1856951),("Lactobacillus_gasseri_ATCC_33323",0.58,244,1894360),("Lactococcus_lactis_cremoris_SK11",0.56,301,2438589),("Leuconostoc_mesenteroides_mesenteroides_ATCC_8293",0.52,235,2038396),("Magnetococcus_sp._MC-1",0.48,504,4719581),("Marinobacter_aquaeolei_VT8",0.57,547,4326849),("Mesorhizobium_sp._BNC1",0.58,567,4412446),("Methanococcoides_burtonii_DSM6242",0.47,268,2575032),("Methanosarcina_barkeri_Fusaro",0.51,545,4837408),("Methanospirillum_hungatei_JF-1",0.55,429,3544738),("Methylobacillus_flagellatus_strain_KT",0.56,365,2971517),("Moorella_thermoacetica_ATCC_39073",1.16,674,2628784),("Nitrobacter_hamburgensis_UNDEF",0.65,630,4406967),("Nitrobacter_winogradskyi_Nb-255",0.57,427,3402093),("Nitrosococcus_oceani_UNDEF",0.53,409,3481691),("Nitrosomonas_eutropha_C71",0.53,314,2661057),("Nitrosospira_multiformis_ATCC_25196",0.54,378,3184243),("Nocardioides_sp._JS614",0.58,636,4985871),("Novosphingobium_aromaticivorans_DSM_12444_F199",0.66,520,3561584),("Oenococcus_oeni_PSU-1",0.46,182,1780517),("Paracoccus_denitrificans_PD1222",0.58,585,4582380),("Pediococcus_pentosaceus_ATCC_25745",0.54,217,1832387),("Pelobacter_carbinolicus_DSM_2380",0.6,489,3665893),("Pelobacter_propionicus_DSM_2379",0.57,508,4008000),("Pelodictyon_luteolum_UNDEF",0.48,250,2364842),("Pelodictyon_phaeoclathratiforme_BU-1_DSMZ_5477T",0.6,402,3018238),("Polaromonas_sp._JS666",0.64,733,5200264),("Prochlorococcus_marinus_str._MIT_9312",0.48,183,1709204),("Prochlorococcus_sp._NATL2A",0.62,253,1842899),("Prosthecochloris_aestuarii_SK413/DSMZ_271t",0.51,282,2512923),("Prosthecochloris_sp._BS1",0.8,483,2736403),("Pseudoalteromonas_atlantica_T6c",0.51,588,5187005),("Pseudomonas_fluorescens_PfO-1",0.51,730,6438405),("Pseudomonas_putida_F1",0.51,675,5959964),("Pseudomonas_syringae_B728a",0.55,746,6093698),("Psychrobacter_arcticum_273-4",0.56,327,2650701),("Psychrobacter_cryopegella_UNDEF",0.62,422,3059876),("Rhodobacter_sphaeroides_2.4.1",0.56,514,4131626),("Rhodoferax_ferrireducens_UNDEF",0.58,599,4712337),("Rhodopseudomonas_palustris_BisA53",0.52,636,5505494),("Rhodopseudomonas_palustris_BisB18",0.57,699,5513844),("Rhodopseudomonas_palustris_BisB5",0.53,575,4892717),("Rhodopseudomonas_palustris_HaA2",24.49,28861,5331656),("Rhodospirillum_rubrum_ATCC_11170",0.58,559,4352825),("Rubrobacter_xylanophilus_DSM_9941",0.57,409,3225748),("Saccharophagus_degradans_2-40",0.52,582,5057531),("Shewanella_amazonensis_SB2B",0.56,536,4306142),("Shewanella_baltica_OS155",0.55,621,5127376),("Shewanella_frigidimarina_NCMB400",0.51,551,4845257),("Shewanella_putefaciens_UNDEF",0.55,565,4659220),("Shewanella_sp._ANA-3",0.6,664,4972204),("Shewanella_sp._MR-7",0.54,568,4792610),("Shewanella_sp._PV-4",0.52,524,4602594),("Shewanella_sp._W3-18-1",0.51,533,4708380),("Silicibacter_sp._TM1040",0.66,469,3200938),("Sphingopyxis_alaskensis_RB2256",0.59,438,3345170),("Streptococcus_suis_89/1591",0.56,263,2143334),("Streptococcus_thermophilus_LMD-9",0.43,178,1856368),("Synechococcus_sp._PCC_7942_elongatus",0.53,316,2695903),("Syntrophobacter_fumaroxidans_MPOB",0.55,606,4990251),("Syntrophomonas_wolfei_Goettingen",0.48,314,2936195),("Thermoanaerobacter_ethanolicus_39E",0.6,315,2362816),("Thermobifida_fusca_YX",0.54,434,3642249),("Thiobacillus_denitrificans_ATCC_25259",0.61,395,2909809),("Thiomicrospira_crunogena_XCL-2",0.51,274,2427734),("Thiomicrospira_denitrificans_ATCC_338890,Thiomicrospira_denitrificans_ATCC_33889",0.57,277,2201561),("Trichodesmium_erythraeum_IMS101",0.57,977,7750108),("Xylella_fastidiosa_Ann-1",2.51,2836,5115560),("Xylella_fastidiosa_Dixon",1.04,601,2622359)]))

	true_species = [s.lower().replace(" ","_") for s in truth[0]]
	true_abundance = np.array([float(i) for i in truth[1]])
	true_counts = np.array([float(i) for i in truth[2]])
	true_size = np.array([float(i) for i in truth[3]])

	return true_species,true_abundance,true_counts,true_size

def calc_error(true_species,true_abundance,est_species,est_abundance):
	adjusted_abundance = np.zeros((len(true_abundance)))
	if true_species != est_species:
		print "Possible species mismatch."
		if len(true_species) == len(est_species): # name synonyms
			print "Looking for species synonyms..."
			adjusted_abundance = est_abundance
			for ind,sp in enumerate(est_species):
				if not (est_species[ind] == true_species[ind] or est_species[ind] == true_species[ind].partition(',')[0] or est_species[ind] == true_species[ind].partition(',')[2]): # does not match species or a listed synonym
					print "Error: species match not found between {0} and {1}. Pretending they are the same for the sake of further calculations.".format(est_species[ind],true_species[ind])
		else: # mapped against more species than actually in sample
			print "Mapped against larger sample of species than in dataset."
			for ind,sp in enumerate(true_species):
				try:
					sp_match = est_species.index(sp)
				except ValueError: # most likely due to synonym
					try:
						sp_match = est_species.index(true_species[ind].partition(',')[0])
					except ValueError:
						try:
							sp_match = est_species.index(true_species[ind].partition(',')[2])
						except ValueError: # now we've really failed
							print "Error: species match for {0} not found in data. Setting abundance and counts to 0 for {0}.".format(sp)
							adjusted_abundance[ind] = 0
						else:
							adjusted_abundance[ind] = est_abundance[sp_match]
					else:
						adjusted_abundance[ind] = est_abundance[sp_match]
				else:
					adjusted_abundance[ind] = est_abundance[sp_match]
	else:
		adjusted_abundance = est_abundance

	diff = 100*abs(true_abundance - adjusted_abundance)/true_abundance

	# graph diff
	print "Error: {0}".format(list(diff))

	print "Average relative error: {0}".format(np.mean(diff))
	diff_sq = [d*d/10000 for d in diff]
	print "Relative root mean squared error: {0}".format(np.mean(diff_sq) ** (0.5) * 100)

	return diff,np.mean(diff),np.mean(diff_sq) ** (0.5) * 100

def main(argv=sys.argv):
	filename = argv[1]
	dataset = argv[2] #i100 or simLC

	true_species,true_abundance,true_counts,true_size = dataset_truth(dataset)
	est_species,est_abundance,est_counts = process_input(filename,true_size)

	calc_error(true_species,true_abundance,est_species,est_abundance)

if __name__ == "__main__":
    main()
