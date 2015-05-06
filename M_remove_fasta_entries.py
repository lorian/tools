# Remove genomes from a fasta found in a separate list

import sys
import os
import math
import lanthpy
import pprint

# Is the strain already present in the i100 fasta?
def is_duplicate(strain,remove):
	i100_strains = ['acinetobacter_baumannii_sdf', 'alkalilimnicola_ehrlichii_mlhe1', 'alkaliphilus_metalliredigens_qymf', 'anabaena_variabilis_atcc_29413', 'bacillus_anthracis_str_ames', 'bacillus_cereus_atcc_10987', 'bacillus_cereus_atcc_14579', 'bacillus_cereus_e33l', 'bacillus_clausii_ksmk16', 'bacillus_halodurans_c125', 'bacillus_subtilis_subsp_subtilis_str_168', 'bacillus_thuringiensis_str_al_hakam', 'bacteroides_thetaiotaomicron_vpi5482', 'bordetella_bronchiseptica_rb50', 'borrelia_garinii_pbi', 'buchnera_aphidicola_bcc', 'burkholderia_phymatum_stm815', 'burkholderia_pseudomallei_668', 'campylobacter_concisus_13826', 'candidatus_blochmannia_pennsylvanicus_str_bpen', 'chlamydia_trachomatis_l2buch1proctitis', 'chlamydophila_caviae_gpic', 'chlamydophila_pneumoniae_cwl029', 'chlorobium_luteolum_dsm_273', 'chromobacterium_violaceum_atcc_12472', 'clostridium_beijerinckii_ncimb_8052', 'clostridium_botulinum_b1_str_okra', 'clostridium_novyi_nt', 'clostridium_thermocellum_atcc_27405', 'corynebacterium_efficiens_ys314', 'corynebacterium_urealyticum_dsm_7109', 'coxiella_burnetii_rsa_493', 'cyanothece_sp_atcc_51142', 'cytophaga_hutchinsonii_atcc_33406', 'erwinia_tasmaniensis', 'escherichia_coli_apec_o1', 'escherichia_coli_e24377a', 'escherichia_coli_str_k12_substr_mg1655', 'escherichia_coli_str_k12_substr_w3110', 'francisella_tularensis_subsp_tularensis_fsc198', 'geobacter_sulfurreducens_pca', 'haemophilus_influenzae_rd_kw20', 'haemophilus_somnus_129pt', 'haloquadratum_walsbyi_dsm_16790', 'herpetosiphon_aurantiacus_dsm_785', 'hyphomonas_neptunium_atcc_15444', 'idiomarina_loihiensis_l2tr', 'ignicoccus_hospitalis_kin4i', 'lactobacillus_delbrueckii_subsp_bulgaricus_atcc_baa365', 'lactobacillus_salivarius_ucc118', 'lactococcus_lactis_subsp_cremoris_sk11', 'lactococcus_lactis_subsp_lactis_il1403', 'lawsonia_intracellularis_phemn100', 'leptospira_borgpetersenii_serovar_hardjobovis_jb197', 'listeria_welshimeri_serovar_6b_str_slcc5334', 'methanococcus_maripaludis_c7', 'mycobacterium_avium_subsp_paratuberculosis_k10', 'mycobacterium_bovis_af212297', 'mycobacterium_marinum_m', 'mycobacterium_sp_jls', 'mycobacterium_tuberculosis_f11', 'neisseria_gonorrhoeae_fa_1090', 'neisseria_meningitidis_mc58', 'nostoc_punctiforme_pcc_73102', 'ochrobactrum_anthropi_atcc_49188', 'prochlorococcus_marinus_str_mit_9313', 'prochlorococcus_marinus_str_natl1a', 'pseudomonas_putida_w619', 'pseudomonas_stutzeri_a1501', 'psychrobacter_cryohalolentis_k5', 'psychromonas_ingrahamii_37', 'ralstonia_eutropha_h16', 'rhodococcus_jostii_rha1', 'rhodopseudomonas_palustris_bisa53', 'rhodopseudomonas_palustris_cga009', 'rickettsia_prowazekii_str_madrid_e', 'salmonella_enterica_subsp_enterica_serovar_paratyphi_b_str_spb7', 'shewanella_sp_mr7', 'sinorhizobium_meliloti_1021', "sodalis_glossinidius_str_morsitans", 'staphylococcus_aureus_subsp_aureus_jh1', 'staphylococcus_aureus_subsp_aureus_mrsa252', 'staphylococcus_aureus_subsp_aureus_mu3', 'staphylococcus_aureus_subsp_aureus_mw2', 'staphylococcus_aureus_subsp_aureus_nctc_8325', 'staphylococcus_aureus_subsp_aureus_str_newman', 'streptococcus_agalactiae_2603vr', 'streptococcus_pneumoniae_r6', 'streptococcus_pyogenes_mgas315', 'streptococcus_pyogenes_mgas5005', 'synechococcus_elongatus_pcc_6301', 'synechococcus_elongatus_pcc_7942', 'syntrophomonas_wolfei_subsp_wolfei_str_goettingen', 'thermosipho_melanesiensis_bi429', 'thermus_thermophilus_hb8', 'thiobacillus_denitrificans_atcc_25259', 'treponema_denticola_atcc_35405', 'ureaplasma_parvum_serovar_3_str_atcc_700970', 'wigglesworthia_glossinidia_endosymbiont_of_glossina_brevipalpis', 'yersinia_pestis_co92']
	i400_strains = ["acaryochloris_marina_mbic11017","acholeplasma_laidlawii_pg-8a","acidiphilium_cryptum_jf-5","acidothermus_cellulolyticus_11b","acidovorax_sp._js42","acinetobacter_baumannii_atcc_17978","acinetobacter_sp._adp1","actinobacillus_pleuropneumoniae_l20","actinobacillus_succinogenes_130z","aeromonas_hydrophila_subsp._hydrophila_atcc_7966","aeromonas_salmonicida_subsp._salmonicida_a449","alcanivorax_borkumensis_sk2","Alkalilimnicola ehrlichii MLHE-1","Alkalilimnicola ehrlichei MLHE-1","alkaliphilus_oremlandii_ohilas","anaeromyxobacter_dehalogenans_2cp-c","anaeromyxobacter_sp._fw109-5","anaplasma_marginale_str._st._maries","anaplasma_phagocytophilum_hz","aquifex_aeolicus_vf5","archaeoglobus_fulgidus_dsm_4304","arcobacter_butzleri_rm4018","arthrobacter_aurescens_tc1","arthrobacter_sp._fb24","aster_yellows_witches'-broom_phytoplasma_aywb","Azoarcus aromaticum EbN1","Azoarcus sp. EbN1","azorhizobium_caulinodans_ors_571","Bacillus amyloliquefaciens subsp. plantarum str. FZB42","Bacillus amyloliquefaciens FZB42","bacillus_anthracis_str._'ames_ancestor'","bacillus_cereus_e33l","bacillus_halodurans_c-125","bacillus_licheniformis_atcc_14580","bacillus_pumilus_safr-032","bacillus_thuringiensis_str._al_hakam","bacillus_weihenstephanensis_kbab4","bacteroides_fragilis_ych46","bacteroides_vulgatus_atcc_8482","bartonella_bacilliformis_kc583","Bartonella henselae strain Houston-1","Bartonella henselae str. Houston-1","bartonella_quintana_str._toulouse","bartonella_tribocorum_cip_105476","bifidobacterium_adolescentis_atcc_15703","bifidobacterium_longum_djo10a","bordetella_avium_197n","bordetella_parapertussis_strain_12822","bordetella_pertussis_tohama_i","Bordetella parapertussis 12822","Bordetella petrii strain DSM 12804","Bordetella petrii DSM 12804","borrelia_afzelii_pko","borrelia_burgdorferi_b31","borrelia_garinii_pbi","bradyrhizobium_japonicum_usda_110","bradyrhizobium_sp._btai1","brucella_canis_atcc_23365","Brucella melitensis bv. 1 str. 16M","Brucella melitensis 16M","brucella_ovis_atcc_25840","brucella_suis_1330","buchnera_aphidicola_str._aps_(acyrthosiphon_pisum)","burkholderia_cenocepacia_mc0-3","burkholderia_mallei_nctc_10247","burkholderia_multivorans_atcc_17616","burkholderia_phymatum_stm815","burkholderia_pseudomallei_1710b","burkholderia_thailandensis_e264","burkholderia_vietnamiensis_g4","burkholderia_xenovorans_lb400","caldicellulosiruptor_saccharolyticus_dsm_8903","caldivirga_maquilingensis_ic-167","campylobacter_curvus_525.92","campylobacter_hominis_atcc_baa-381","campylobacter_jejuni_rm1221","candidatus_blochmannia_pennsylvanicus_str._bpen","candidatus_carsonella_ruddii_pv","candidatus_desulforudis_audaxviator_mp104c","candidatus_korarchaeum_cryptofilum_opf8","Candidatus Koribacter versatilis Ellin345","Acidobacteria bacterium Ellin345","candidatus_pelagibacter_ubique_htcc1062","candidatus_ruthia_magnifica_str._cm_(calyptogena_magnifica)","candidatus_vesicomyosocius_okutanii_ha","carboxydothermus_hydrogenoformans_z-2901","caulobacter_crescentus_cb15","chlamydia_muridarum_nigg","chlamydia_trachomatis_duw-3cx","chlamydophila_caviae_gpic","chlamydophila_pneumoniae_ar39","Chlorobium luteolum DSM 273","Pelodictyon luteolum DSM 273","chlorobium_tepidum_tls","chloroflexus_aurantiacus_j-10-fl","chromobacterium_violaceum_atcc_12472","chromohalobacter_salexigens_dsm_3043","citrobacter_koseri_atcc_baa-895","clostridium_acetobutylicum_atcc_824","clostridium_beijerinckii_ncimb_8052","clostridium_botulinum_a3_str._loch_maree","clostridium_difficile_630","clostridium_kluyveri_dsm_555","clostridium_novyi_nt","clostridium_perfringens_atcc_13124","clostridium_phytofermentans_isdg","clostridium_tetani_e88","clostridium_thermocellum_atcc_27405","colwellia_psychrerythraea_34h","corynebacterium_diphtheriae_nctc_13129","corynebacterium_efficiens_ys-314","corynebacterium_jeikeium_k411","corynebacterium_urealyticum_dsm_7109","coxiella_burnetii_dugway_5j108-111","cupriavidus_taiwanensis","cyanothece_sp._atcc_51142","cytophaga_hutchinsonii_atcc_33406","dechloromonas_aromatica_rcb","dehalococcoides_ethenogenes_195","dehalococcoides_sp._cbdb1","deinococcus_radiodurans_r1","delftia_acidovorans_sph-1","desulfococcus_oleovorans_hxd3","desulfotalea_psychrophila_lsv54","desulfotomaculum_reducens_mi-1","Desulfovibrio alaskensis G20","Desulfovibrio desulfuricans subsp. desulfuricans str. G20","desulfovibrio_vulgaris_subsp._vulgaris_str._hildenborough","dichelobacter_nodosus_vcs1703a","dinoroseobacter_shibae_dfl_12","ehrlichia_chaffeensis_str._arkansas","ehrlichia_ruminantium_str._welgevonden","enterobacter_sp._638","erwinia_tasmaniensis","erythrobacter_litoralis_htcc2594","exiguobacterium_sibiricum_255-15","finegoldia_magna_atcc_29328","flavobacterium_johnsoniae_uw101","flavobacterium_psychrophilum_jip0286","francisella_philomiragia_subsp._philomiragia_atcc_25017","Frankia alni str. ACN14a","Frankia alni ACN14a","fusobacterium_nucleatum_subsp._nucleatum_atcc_25586","geobacillus_kaustophilus_hta426","geobacillus_thermodenitrificans_ng80-2","geobacter_lovleyi_sz","geobacter_sulfurreducens_pca","geobacter_uraniireducens_rf4","gloeobacter_violaceus_pcc_7421","gluconacetobacter_diazotrophicus_pal_5","gluconobacter_oxydans_621h","gramella_forsetii_kt0803","granulibacter_bethesdensis_cgdnih1","Haemophilus ducreyi strain 35000HP","Haemophilus ducreyi 35000HP","haemophilus_influenzae_pittgg","haemophilus_somnus_2336","haloarcula_marismortui_atcc_43049","halobacterium_salinarum_r1","halobacterium_sp._nrc-1","halorhodospira_halophila_sl1","helicobacter_acinonychis_str._sheeba","helicobacter_hepaticus_atcc_51449","helicobacter_pylori_26695","heliobacterium_modesticaldum_ice1","herminiimonas_arsenicoxydans","hyperthermus_butylicus_dsm_5456","hyphomonas_neptunium_atcc_15444","idiomarina_loihiensis_l2tr","ignicoccus_hospitalis_kin4i","jannaschia_sp._ccs1","janthinobacterium_sp._marseille","kineococcus_radiotolerans_srs30216","klebsiella_pneumoniae_subsp._pneumoniae_mgh_78578","lactobacillus_acidophilus_ncfm","lactobacillus_brevis_atcc_367","lactobacillus_casei_atcc_334","lactobacillus_delbrueckii_subsp._bulgaricus_atcc_baa-365","lactobacillus_fermentum_ifo_3956","lactobacillus_gasseri_atcc_33323","lactobacillus_helveticus_dpc_4571","lactobacillus_plantarum_wcfs1","Lactobacillus reuteri DSM 20016","Lactobacillus reuteri F275","Lactobacillus sakei strain 23K","Lactobacillus sakei subsp. sakei 23K","lactococcus_lactis_subsp._cremoris_sk11","lawsonia_intracellularis_phemn1-00","legionella_pneumophila_str._paris","leifsonia_xyli_subsp._xyli_str._ctcb07","leptospira_biflexa_serovar_patoc_strain_'patoc_1_(paris)'","leptospira_borgpetersenii_serovar_hardjo-bovis_l550","leptospira_interrogans_serovar_lai_str._56601","leptothrix_cholodnii_sp-6","leuconostoc_citreum_km20","leuconostoc_mesenteroides_subsp._mesenteroides_atcc_8293","listeria_innocua_clip11262","listeria_welshimeri_serovar_6b_str._slcc5334","magnetospirillum_magneticum_amb-1","mannheimia_succiniciproducens_mbel55e","maricaulis_maris_mcs10","marinobacter_aquaeolei_vt8","mesoplasma_florum_l1","mesorhizobium_loti_maff303099","Chelativorans sp. BNC1","Mesorhizobium sp. BNC1","methanobrevibacter_smithii_atcc_35061","methanocaldococcus_jannaschii_dsm_2661","Methanocella arvoryzae MRE50","uncultured methanogenic archaeon RC-I","methanococcoides_burtonii_dsm_6242","methanococcus_aeolicus_nankai-3","methanococcus_maripaludis_c5","methanococcus_vannielii_sb","methanocorpusculum_labreanum_z","methanoculleus_marisnigri_jr1","methanopyrus_kandleri_av19","methanosaeta_thermophila_pt","methanosarcina_barkeri_str._fusaro","methanosarcina_mazei_go1","methanosphaera_stadtmanae_dsm_3091","methanospirillum_hungatei_jf-1","methanothermobacter_thermautotrophicus_str._delta_h","methylobacillus_flagellatus_kt","methylobacterium_extorquens_pa1","methylobacterium_radiotolerans_jcm_2831","methylobacterium_sp._4-46","methylococcus_capsulatus_str._bath","microcystis_aeruginosa_nies-843","moorella_thermoacetica_atcc_39073","Mycobacterium bovis BCG Pasteur 1173P2","Mycobacterium bovis BCG str. Pasteur 1173P2","mycobacterium_gilvum_pyr-gck","mycobacterium_leprae_tn","mycobacterium_marinum_m","mycobacterium_smegmatis_str._mc2_155","mycobacterium_sp._kms","mycobacterium_ulcerans_agy99","mycobacterium_vanbaalenii_pyr-1","mycoplasma_agalactiae_pg2","mycoplasma_capricolum_subsp._capricolum_atcc_27343","Mycoplasma gallisepticum str. R(low)","Mycoplasma gallisepticum R","mycoplasma_genitalium_g37","mycoplasma_mobile_163k","mycoplasma_mycoides_subsp._mycoides_sc_str._pg1","mycoplasma_penetrans_hf-2","mycoplasma_pulmonis_uab_ctip","mycoplasma_synoviae_53","myxococcus_xanthus_dk_1622","nanoarchaeum_equitans_kin4-m","natronomonas_pharaonis_dsm_2160","neisseria_gonorrhoeae_fa_1090","neisseria_meningitidis_mc58","Neorickettsia sennetsu strain Miyayama","Neorickettsia sennetsu str. Miyayama","nitrobacter_hamburgensis_x14","nitrobacter_winogradskyi_nb-255","nitrosococcus_oceani_atcc_19707","nitrosomonas_europaea_atcc_19718","nitrosomonas_eutropha_c91","nitrosopumilus_maritimus_scm1","nitrosospira_multiformis_atcc_25196","nocardia_farcinica_ifm_10152","nocardioides_sp._js614","nostoc_punctiforme_pcc_73102","nostoc_sp._pcc_7120","novosphingobium_aromaticivorans_dsm_12444","oceanobacillus_iheyensis_hte831","ochrobactrum_anthropi_atcc_49188","oenococcus_oeni_psu-1","onion_yellows_phytoplasma_oy-m","opitutus_terrae_pb90-1","orientia_tsutsugamushi_str._boryong","parabacteroides_distasonis_atcc_8503","Parachlamydia-related symbiont UWE25","Candidatus Protochlamydia amoebophila UWE25","paracoccus_denitrificans_pd1222","parvibaculum_lavamentivorans_ds-1","pasteurella_multocida_subsp._multocida_str._pm70","pelobacter_carbinolicus_dsm_2380","petrotoga_mobilis_sj95","photobacterium_profundum_ss9","photorhabdus_luminescens_subsp._laumondii_tto1","picrophilus_torridus_dsm_9790","polaromonas_naphthalenivorans_cj2","polaromonas_sp._js666","Polynucleobacter necessarius subsp. asymbioticus QLW-P1DMWA-1","Polynucleobacter sp. QLW-P1DMWA-1","Polynucleobacter necessarius subsp. necessarius STIR1","Polynucleobacter necessarius STIR1","porphyromonas_gingivalis_atcc_33277","propionibacterium_acnes_kpa171202","pseudoalteromonas_atlantica_t6c","Pseudoalteromonas haloplanktis str. TAC125","Pseudoalteromonas haloplanktis TAC125","pseudomonas_aeruginosa_pa7","Pseudomonas protegens Pf-5","Pseudomonas fluorescens Pf-5","pseudomonas_putida_kt2440","pseudomonas_stutzeri_a1501","pseudomonas_syringae_pv._tomato_str._dc3000","psychrobacter_arcticus_273-4","psychrobacter_cryohalolentis_k5","psychrobacter_sp._prwf-1","psychromonas_ingrahamii_37","pyrobaculum_aerophilum_str._im2","pyrobaculum_arsenaticum_dsm_13514","pyrobaculum_calidifontis_jcm_11548","pyrobaculum_islandicum_dsm_4184","pyrococcus_abyssi_ge5","pyrococcus_furiosus_dsm_3638","pyrococcus_horikoshii_ot3","ralstonia_eutropha_jmp134","Cupriavidus metallidurans CH34","Ralstonia metallidurans CH34","ralstonia_solanacearum_gmi1000","renibacterium_salmoninarum_atcc_33209","rhizobium_etli_cfn_42","rhizobium leguminosarum bv. viciae chromosome complete genome","rhodobacter_sphaeroides_2.4.1","Rhodococcus jostii RHA1","Rhodococcus sp. RHA1","rhodoferax_ferrireducens_t118","rhodopirellula_baltica_sh_1","rhodopseudomonas_palustris_bisa53","rhodospirillum_rubrum_atcc_11170","rickettsia_akari_str._hartford","rickettsia_bellii_rml369-c","rickettsia_canadensis_str._mckiel","rickettsia_felis_urrwxcal2","rickettsia_massiliae_mtu5","rickettsia_prowazekii_str._madrid_e","rickettsia_rickettsii_str._'sheila_smith'","rickettsia_typhi_str._wilmington","roseiflexus_castenholzii_dsm_13941","roseiflexus_sp._rs-1","roseobacter_denitrificans_och_114","rubrobacter_xylanophilus_dsm_9941","Ruegeria pomeroyi DSS-3","Silicibacter pomeroyi DSS-3""","saccharophagus_degradans_2-40","Saccharopolyspora erythraea NRRL2338","Saccharopolyspora erythraea NRRL 2338","salinibacter_ruber_dsm_13855","salinispora_arenicola_cns-205","salinispora_tropica_cnb-440","Salmonella enterica subsp. enterica serovar typhimurium str. LT2","Salmonella typhimurium LT2","serratia_proteamaculans_568","shewanella_amazonensis_sb2b","shewanella_baltica_os195","shewanella_denitrificans_os217","shewanella_frigidimarina_ncimb_400","shewanella_halifaxensis_haw-eb4","shewanella_loihica_pv-4","shewanella_pealeana_atcc_700345","shewanella_putrefaciens_cn-32","shewanella_sediminis_haw-eb3","shewanella_sp._ana-3","shewanella_woodyi_atcc_51908","shigella_boydii_cdc_3083-94","shigella_dysenteriae_sd197","shigella_flexneri_2a_str._301","silicibacter_sp._tm1040","sinorhizobium_medicae_wsm419","sinorhizobium_meliloti_1021","sphingomonas_wittichii_rw1","sphingopyxis_alaskensis_rb2256","staphylococcus_aureus_subsp._aureus_usa300","staphylococcus_epidermidis_rp62a","staphylococcus_haemolyticus_jcsc1435","staphylococcus_saprophyticus_subsp._saprophyticus_atcc_15305","staphylothermus_marinus_f1","streptococcus_agalactiae_nem316","streptococcus_gordonii_str._challis_substr._ch1","streptococcus_mutans_ua159","streptococcus_pneumoniae_cgsp14","streptococcus_pyogenes_mgas10750","streptococcus_sanguinis_sk36","streptococcus_suis_05zyh33","streptococcus_thermophilus_lmd-9","streptomyces_avermitilis_ma-4680","streptomyces_coelicolor_a3(2)","sulfolobus_acidocaldarius_dsm_639","sulfolobus_solfataricus_p2","sulfurihydrogenibium_sp._yo3aop1","sulfurimonas_denitrificans_dsm_1251","sulfurovum_sp._nbc37-1","symbiobacterium_thermophilum_iam_14863","synechococcus_elongatus_pcc_7942","synechococcus_sp._pcc_7002","syntrophomonas_wolfei_subsp._wolfei_str._goettingen","thermoanaerobacter_pseudethanolicus_atcc_33223","thermoanaerobacter_sp._x514","thermoanaerobacter_tengcongensis_mb4","thermobifida_fusca_yx","thermofilum_pendens_hrk_5","thermoplasma_acidophilum_dsm_1728","thermoplasma_volcanium_gss1","thermoproteus_neutrophilus_v24sta","thermosipho_melanesiensis_bi429","thermosynechococcus_elongatus_bp-1","thermotoga_lettingae_tmo","thermotoga_petrophila_rku-1","thermotoga_sp._rq2","thermus_thermophilus_hb27","thiobacillus_denitrificans_atcc_25259","thiomicrospira_crunogena_xcl-2","treponema_denticola_atcc_35405","treponema_pallidum_subsp._pallidum_str._nichols","trichodesmium_erythraeum_ims101","tropheryma_whipplei_str._twist","ureaplasma_parvum_serovar_3_str._atcc_27815","verminephrobacter_eiseniae_ef01-2","vibrio_cholerae_o395","vibrio_fischeri_es114","vibrio_harveyi_atcc_baa-1116","vibrio_parahaemolyticus_rimd_2210633","vibrio_vulnificus_yj016","wolbachia_endosymbiont_of_drosophila_melanogaster","wolinella_succinogenes_dsm_1740","xanthobacter_autotrophicus_py2","xanthomonas_axonopodis_pv._citri_str._306","xanthomonas_oryzae_pv._oryzae_maff_311018","xylella_fastidiosa_9a5c","yersinia_enterocolitica_subsp._enterocolitica_8081","yersinia_pseudotuberculosis_ip_31758","zymomonas_mobilis_subsp._mobilis_zm4"]

	if remove == 'i100':
		remove_strains = lanthpy.genome_name_cleanup(i100_strains)
	elif remove == 'i400':
		remove_strains = lanthpy.genome_name_cleanup(i400_strains)
	else:
		sys.exit("Don't have list of {0} strains".format(remove))

	if any(dupe in lanthpy.genome_name_cleanup([strain])[0] for dupe in lanthpy.genome_name_cleanup(remove_strains)):
		return True
	print "{0} is not a duplicate".format(strain)
	return False

def process_file(filename,remove):
	mfa = open(filename,'r')

	file_count_g = 0
	file_count_p = 0
	line_count = 0

	# Open output file
	output_filename = filename.rsplit('.', 1)[0] + "_dedupe.fasta"
	output = open(output_filename, 'wt')
	output.seek(0) #overwrite if it exists

	skip_current = False
	list_of_strains = []
	for line in mfa:
		if line[0] == '>': # fasta name
			list_of_strains.append(line[:-1])
			skip_current = False
			line = line.replace (" ", "_")

			if 'VIRL' in line: #ignore viruses
				skip_current = True
			else: # check if the genome makes the bad list
				skip_current = is_duplicate(line,remove)

				if not skip_current:
					output.write(line)

		else: #actual fasta content
			if not skip_current: # first store entire genome in variable
				output.write(line)

	# close files
	output.truncate()
	output.close()

def main():
	''' python M_remove_fasta_entries.py <filename> <i100 or i400> '''

	filename = sys.argv[1]
	remove = sys.argv[2]
	process_file(filename,remove)

if __name__ == '__main__':
	main()