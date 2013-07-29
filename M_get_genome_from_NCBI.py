# Download genomes of listed species, avoiding duplicates

import urllib2
import string
import os.path

def parse_line(line,species,entries):

	if line.startswith('>'): # new fasta entry
		name = line[line.find(" ")+1:].lower()[:-1]
		if name.find(species) != -1: # matches species
			if name in entries or (name.find('complete genome') != -1 and entries != []): # duplicate name or duplicate 'complete genome'
				print "\tSkipping duplicate {0}".format(name)
				return False
			else: # valid new fasta that matches species
				entries.append(name)
				print "Added {0}".format(name)
				return '>' + species.replace (" ", "_") + " " + line[1:] #add name to beginning of ID line
		else: # does not match species
			print "\tSkipping {0}".format(name)
			return False
	elif line.startswith('<'): # xml instead of fasta
		print "ERROR: EMPTY PAGE"
		return False

	return line # normal line of fasta

def get_genome(species_orig):
	species = species_orig.lower()
	filename = species.replace (" ", "_").replace('/',"") + '.mfa'

	if os.path.isfile(filename) and os.path.getsize(filename) > 100: # don't re-run if file exists and is sensible size
		return
	else:
		print "SPECIES: {0}".format(species_orig)

	mfa = open(filename) # new fasta file for every species
	mfa.seek(0)

	# get genome ID
	if not 'ID_list' in globals(): #if I don't already have a list
		page_genome = urllib2.urlopen('http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=genome&term=' + urllib2.quote(species))
		genome_id = ""
		for line in page_genome:
			if line.startswith("<Id>"):
				if genome_id != "":
					print "Duplicate IDs for {0}; skipping.".format(species)
					return
				genome_id = line[4:-6]
	else:
		genome_id = ID_list[species_list.index(species_orig)]

	# get query key and web env for chromosome
	page_nucleotide = urllib2.urlopen('http://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=genome&db=nuccore&id={0}&term=gene+in+chromosome[prop]&cmd=neighbor_history'.format(genome_id))
	for line in page_nucleotide:
		line = string.replace(line,"\t","")
		if line.startswith("<QueryKey>"):
			query_key = line[10:-12]
		if line.startswith("<WebEnv>"):
			web_env = line[8:-10]

	# get fasta for chromosome
	try:
		page_fasta = urllib2.urlopen('http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&query_key=' + query_key + '&WebEnv=' + web_env + '&rettype=fasta&retmode=text')
	except:
		print "FAILED TO GET PAGE FOR {0}".format(species)
		return 0

	# copy chromosome fastas, trying to skip duplicates
	fasta = ""
	skip = False
	entries = []
	for line in page_fasta:
		if not skip or line.startswith('>'): # if skip, avoid entire loop until next new record
			skip = False
			nl = parse_line(line,species,entries)
			if nl:
				fasta = fasta + nl
			else:
				skip = True

	# get query key and web env for plasmid
	page_nucleotide = urllib2.urlopen('http://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=genome&db=nuccore&id=' + genome_id + '&term=gene+in+plasmid[prop]&cmd=neighbor_history')
	for line in page_nucleotide:
		line = string.replace(line,"\t","")
		if line.startswith("<QueryKey>"):
			query_key = line[10:-12]
		if line.startswith("<WebEnv>"):
			web_env = line[8:-10]

	# get fasta for plasmid
	page_fasta = urllib2.urlopen('http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&query_key=' + query_key + '&WebEnv=' + web_env + '&rettype=fasta&retmode=text')

	# copy plasmid fastas, trying to skip duplicates
	fasta = fasta + "\n"
	skip = False
	entries = []
	for line in page_fasta:
		if not skip or line.startswith('>'): # if skip, avoid entire loop until next new record
			skip = False
			nl = parse_line(line,species,entries)
			if nl:
				fasta = fasta + nl
			else:
				skip = True

	fasta = fasta + "\n"
	mfa.write(fasta) # write fastas to file
	mfa.truncate()
	mfa.close()

i100_species_list = [
'Psychrobacter cryohalolentis K5',
'Mycobacterium sp. JLS',
'Synechococcus elongatus PCC 7942',
'Bacillus cereus ATCC 14579',
'Lactococcus lactis subsp. lactis Il1403',
'Wigglesworthia glossinidia endosymbiont of Glossina brevipalpis',
'Pseudomonas putida W619',
'Streptococcus pneumoniae R6',
'Bordetella bronchiseptica RB50',
'Lactobacillus delbrueckii subsp. bulgaricus ATCC BAA-365',
'Clostridium novyi NT',
'Staphylococcus aureus subsp. aureus JH1',
'Borrelia garinii PBi',
'Salmonella enterica subsp. enterica serovar Paratyphi B str. SPB7',
'Erwinia tasmaniensis',
'Staphylococcus aureus subsp. aureus MW2',
'Thiobacillus denitrificans ATCC 25259',
'Rhodopseudomonas palustris BisA53',
'Bacillus anthracis str. Ames',
'Streptococcus agalactiae 2603V/R',
'Burkholderia phymatum STM815',
'Campylobacter concisus 13826',
'Alkaliphilus metalliredigens QYMF',
'Corynebacterium efficiens YS-314',
'Streptococcus pyogenes MGAS5005',
'Ochrobactrum anthropi ATCC 49188',
'Prochlorococcus marinus str. MIT 9313',
'Rhodopseudomonas palustris CGA009',
'Corynebacterium urealyticum DSM 7109',
'Thermus thermophilus HB8',
'Acinetobacter baumannii SDF',
'Neisseria gonorrhoeae FA 1090',
'Hyphomonas neptunium ATCC 15444',
'Staphylococcus aureus subsp. aureus str. Newman',
'Haloquadratum walsbyi DSM 16790',
'Prochlorococcus marinus str. NATL1A',
'Clostridium beijerinckii NCIMB 8052',
'Chromobacterium violaceum ATCC 12472',
'Chlamydia trachomatis L2b/UCH-1/proctitis',
'Cytophaga hutchinsonii ATCC 33406',
'Alkalilimnicola ehrlichii MLHE-1', #fixed spelling
'Bacillus thuringiensis str. Al Hakam',
'Clostridium botulinum B1 str. Okra',
'Methanococcus maripaludis C7',
'Sodalis glossinidius str. \'morsitans\'',
'Cyanothece sp. ATCC 51142',
'Lactococcus lactis subsp. cremoris SK11',
'Haemophilus somnus 129PT',
'Clostridium thermocellum ATCC 27405',
'Neisseria meningitidis MC58',
'Ignicoccus hospitalis KIN4/I',
'Escherichia coli E24377A',
'Thermosipho melanesiensis BI429',
'Lawsonia intracellularis PHE/MN1-00',
'Coxiella burnetii RSA 493',
'Yersinia pestis CO92',
'Bacillus cereus E33L',
'Candidatus Blochmannia pennsylvanicus str. BPEN',
'Idiomarina loihiensis L2TR',
'Syntrophomonas wolfei subsp. wolfei str. Goettingen',
'Bacillus cereus ATCC 10987',
'Escherichia coli APEC O1',
'Mycobacterium tuberculosis F11',
'Bacillus clausii KSM-K16',
'Lactobacillus salivarius UCC118',
'Mycobacterium bovis AF2122/97',
'Nostoc punctiforme PCC 73102',
'Shewanella sp. MR-7',
'Burkholderia pseudomallei 668',
'Rickettsia prowazekii str. Madrid E',
'Haemophilus influenzae Rd KW20',
'Pseudomonas stutzeri A1501',
'Treponema denticola ATCC 35405',
'Ralstonia eutropha H16',
'Mycobacterium avium subsp. paratuberculosis K-10',
'Listeria welshimeri serovar 6b str. SLCC5334',
'Bacillus halodurans C-125',
'Geobacter sulfurreducens PCA',
'Chlamydophila pneumoniae CWL029',
'Francisella tularensis subsp. tularensis FSC198',
'Staphylococcus aureus subsp. aureus NCTC 8325',
'Mycobacterium marinum M',
'Bacillus subtilis subsp. subtilis str. 168',
'Staphylococcus aureus subsp. aureus MRSA252',
'Sinorhizobium meliloti 1021',
'Streptococcus pyogenes MGAS315',
'Ureaplasma parvum serovar 3 str. ATCC 700970',
'Chlamydophila caviae GPIC',
'Anabaena variabilis ATCC 29413',
'Psychromonas ingrahamii 37',
'Synechococcus elongatus PCC 6301',
'Bacteroides thetaiotaomicron VPI-5482',
'Staphylococcus aureus subsp. aureus Mu3',
'Leptospira borgpetersenii serovar Hardjo-bovis JB197',
'Herpetosiphon aurantiacus DSM 785', #new strain name
'Buchnera aphidicola BCc', #new strain name
'Rhodococcus jostii RHA1', #new species name
'Escherichia coli str. K-12 substr. W3110', #added -
'Escherichia coli str. K-12 substr. MG1655', #added -
'Chlorobium luteolum DSM 273', #new genus
'pseudomonas fluorescens pf0-1',
'Rhodopseudomonas palustris HaA2',
'Burkholderia ambifaria AMMD']

ignore_ID_list = [ #NCBI taxonomy id
316058,
288000,
269798,
155920,
264732,
155919,
269484,
266117,
292415,
315749,
269483,
351627,
324602,
290402,
165597,
332415,
59920,
335284,
316056,
94122,
60481,
292414,
340099,
339671,
339670,
290315,
319795,
290400,
266940,
265072,
323098,
279238,
324925,
296591,
342610,
316057,
317655,
326298,
203124,
187272,
240292,
290397,
290399,
322710,
331271,
331272,
269482,
203119,
272564,
207559,
262543,
106370,
298653,
269799,
321967,
324831,
351348,
323259,
323097,
323261,
335283,
323848,
196162,
338963,
338966,
351746,
205918,
338969,
316055,
269796,
326297,
325240,
318167,
319224,
323850,
286604,
335543,
293826,
205913,
321955,
290317,
159087,
333146,
205914,
272622,
203120,
156889,
318586,
278197,
319225,
74546,
290512,
331678,
203122,
351745,
335541,
317025,
266265,
290318,
290398,
333849,
321956,
266779,
259564,
269797,
205922,
259536,
1140,
269800,
203123,
272943,
387344,
322159]

species_list = [
"Rhodopseudomonas palustris HaA2",
"Bradyrhizobium sp. BTAi1",
"Cytophaga hutchinsonii ATCC 33406",
"Xylella fastidiosa Ann-1",
"Moorella thermoacetica ATCC 39073",
"Xylella fastidiosa Dixon",
"Ehrlichia canis Jake",
"Rubrobacter xylanophilus DSM 9941",
"Thiobacillus denitrificans ATCC 25259",
"Bacillus cereus NVH391-98",
"Burkholderia sp. sp.strain 383",
"Caldicellulosiruptor accharolyticus UNDEF",
"Chloroflexus aurantiacus J-10-fl",
"Clostridium beijerincki NCIMB 8052",
"Crocosphaera watsonii WH 8501",
"Ehrlichia chaffeensis sapulpa",
"Prochlorococcus sp. NATL2A",
"Psychrobacter cryopegella UNDEF",
"Rhodopseudomonas palustris BisB18",
"Shewanella sp. ANA-3",
"Shewanella sp. MR-7",
"Silicibacter sp. TM1040",
"Thermoanaerobacter ethanolicus 39E",
"Actinobacillus succinogenes 130Z",
"Burkholderia ambifaria AMMD",
"Chlorobium limicola DSMZ 245(T)",
"Deinococcus geothermalis DSM11300",
"Jannaschia sp. CCS1",
"Kineococcus radiotolerans SRS30216",
"Methylobacillus flagellatus strain KT",
"Nitrobacter winogradskyi Nb-255",
"Novosphingobium aromaticivorans DSM 12444 (F199)",
"Pelodictyon phaeoclathratiforme BU-1 (DSMZ 5477(T))",
"Polaromonas sp. JS666",
"Pseudoalteromonas atlantica T6c",
"Rhodopseudomonas palustris BisB5",
"Sphingopyxis alaskensis RB2256",
"Thiomicrospira denitrificans ATCC 33889",
"Trichodesmium erythraeum IMS101",
"Alkalilimnicola ehrlichei MLHE-1",
"Anabaena variabilis ATCC 29413",
"Anaeromyxobacter dehalogenans 2CP-C",
"Arthrobacter sp. FB24",
"Azotobacter vinelandii AvOP",
"Burkholderia cenocepacia AU 1054",
"Burkholderia cenocepacia HI2424",
"Burkholderia vietnamiensis G4",
"Clostridium thermocellum ATCC 27405",
"Desulfitobacterium hafniense DCB-2",
"Desulfovibrio desulfuricans G20",
"Exiguobacterium UNDEF 255-15",
"Frankia sp. CcI3",
"Frankia sp. EAN1pec",
"Geobacter metallireducens GS-15",
"Lactobacillus casei ATCC 334",
"Lactobacillus gasseri ATCC 33323",
"Marinobacter aquaeolei VT8",
"Methanospirillum hungatei JF-1",
"Nitrobacter hamburgensis UNDEF",
"Nitrosococcus oceani UNDEF",
"Nitrosomonas eutropha C71",
"Nitrosospira multiformis ATCC 25196",
"Nocardioides sp. JS614",
"Pelobacter carbinolicus DSM 2380",
"Pelobacter propionicus DSM 2379",
"Pseudomonas putida F1",
"Pseudomonas syringae B728a",
"Rhodoferax ferrireducens UNDEF",
"Rhodopseudomonas palustris BisA53",
"Rhodospirillum rubrum ATCC 11170",
"Shewanella amazonensis SB2B",
"Shewanella baltica OS155",
"Shewanella frigidimarina NCMB400",
"Shewanella putefaciens UNDEF",
"Shewanella sp. PV-4",
"Streptococcus suis 89/1591",
"Syntrophobacter fumaroxidans MPOB",
"Alkaliphillus metalliredigenes UNDEF",
"Bifidobacterium longum DJO10A",
"Brevibacterium linens BL2",
"Chlorobium phaeobacteroides DSM 266",
"Dechloromonas aromatica RCB",
"Ferroplasma acidarmanus fer1",
"Haemophilus somnus 129PT",
"Lactococcus lactis cremoris SK11",
"Leuconostoc mesenteroides mesenteroides ATCC 8293",
"Magnetococcus sp. MC-1",
"Paracoccus denitrificans PD1222",
"Pediococcus pentosaceus ATCC 25745",
"Pelodictyon luteolum UNDEF",
"Prochlorococcus marinus str. MIT 9312",
"Prosthecochloris aestuarii SK413/DSMZ 271(t)",
"Prosthecochloris sp. BS1",
"Saccharophagus degradans 2-40",
"Shewanella sp. W3-18-1",
"Syntrophomonas wolfei Goettingen",
"Thiomicrospira crunogena XCL-2",
"Burkholderia xenovorans LB400",
"Chlorobium vvibrioforme f. thiosulfatophilum DSMZ 265(T)",
"Chromohalobacter salexigens DSM3043",
"Enterococcus faecium DO",
"Lactobacillus delbrueckii bulgaricus ATCC BAA-365",
"Mesorhizobium sp. BNC1",
"Methanococcoides burtonii DSM6242",
"Methanosarcina barkeri Fusaro",
"Pseudomonas fluorescens PfO-1",
"Psychrobacter arcticum 273-4",
"Synechococcus sp. PCC 7942 (elongatus)",
"Thermobifida fusca YX",
"Oenococcus oeni PSU-1",
"Rhodobacter sphaeroides 2.4.1",
"Lactobacillus brevis ATCC 367",
"Streptococcus thermophilus LMD-9"]

for s in species_list:
	get_genome(s)

