# Download genomes of listed species, avoiding duplicates

import urllib2
import string
import os.path
import sys

def replace_spaces(line):
	return line.replace(" ", "_")

def has_genome(entries): #is there a complete genome in the list?
	for e in entries:
		if e.find('plasmid') == -1 and e.find('complete genome') != -1:
			return True
	return False

def check_duplicate(entries, name): #do we already have this entry in our fasta?
	# first standardize chromosome names:
	if (name in entries or (has_genome([name]) and has_genome(entries))):
		#print "\tSkipping {0}".format(name)
		return True #duplicate
	# filter based on whole genome project
	if name.find('whole genome') != -1 and name.find('project') != -1:
		return True
	# filter out specific genes
	if name.find('gene ') != -1 or name.find('rna ') != -1: # space at end is to try to catch only words, not interior sequences
		print "\tSkipping gene {0}".format(name)
		return True
	print "Adding {0}".format(name)
	return False #not a dupe

def find_best(entries):
	print entries
	if len(entries) == 1:
		return entries
	best = [r for r in entries if (r.find('plasmid') == -1)] # start off by adding all plasmids

	return best

def parse_line(line,species,entries):

	if line.startswith('>'): # new fasta entry
		name = line[line.find(" ")+1:].lower()[:-1].replace('chromosome ii','chromosome 2').replace('chromosome i','chromosome 1') #standardize name
		if name.find(species) != -1: # matches species
			if check_duplicate(entries,name):
				return False
			else: # valid new fasta that matches species
				entries.append(name)
				#print "Added {0}".format(name)
				return '>' + replace_spaces(species) + '_' + replace_spaces(line[1:]) #add name to beginning of ID line
		else: # does not match species
			#print "\tSkipping {0}".format(name)
			return False
	elif line.startswith('<'): # xml instead of fasta
		print "ERROR: FASTA EMPTY PAGE"
		return False

	return line # normal line of fasta


def write_fasta(filename,fasta):
	mfa = open(filename, 'w') # overwrite fasta file
	mfa.seek(0)
	mfa.write(fasta) # write fastas to file
	mfa.truncate()
	mfa.close()


def parse_fasta(page_fasta,species):
	# copy fastas, trying to skip duplicates
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

	return entries, fasta +'\n'


def get_genome(species_orig):
	species = species_orig.lower()
	filename = species.replace(" ", "_").replace('/',"") + '.mfa'
	print "\nLooking up {0}".format(species_orig)

	if os.path.isfile(filename) and os.path.getsize(filename) > 200: # don't re-run if file exists and is sensible size
		return

	if 'tax_ID_list' in globals():
		# species ID
		tax_id = tax_ID_list[species_list.index(species_orig)]
		print "\tID: {0}".format(tax_id)

		# get species name
		try:
			page_tax = urllib2.urlopen('http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id={0}'.format(tax_id))
		except:
			#print "FAILED TO GET TAXONOMY PAGE FOR {0}".format(species)
			return 0

		for line in page_tax:
			line = string.replace(line,"\t","") # unfortunately this xml page uses spaces
			if line.startswith('    <ScientificName>'):
				species = line[20:-18].lower()
				#print 'Actual species name: {0}'.format(species)
				break

	# get genome ID
	if 'refseq_ID_list' in globals():
		chr_id = refseq_ID_list[species_list.index(species_orig)]
		#print "\tREFSEQ: {0}".format(chr_id)
		page_genome = urllib2.urlopen('http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=genome&term={0}&retmax=100000'.format(chr_id))
	else:
		page_genome = urllib2.urlopen('http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=genome&term={0}&retmax=100000'.format(urllib2.quote(species)))

	genome_id = ""
	for line in page_genome:
		line = string.replace(line,"\t","")
		if line.startswith("<Id>"):
			if genome_id != "":
				#print "Duplicate IDs for {0}; skipping.".format(species)
				return
			genome_id = line[4:-6]

	# get query key and web env
	page_nucleotide = urllib2.urlopen('http://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=genome&db=nuccore&id={0}&term=gene+in+chromosome[prop]+OR+gene+in+plasmid[prop]+OR+wgs[prop]&cmd=neighbor_history'.format(genome_id))
	for line in page_nucleotide:
		line = string.replace(line,"\t","")
		if line.startswith("<QueryKey>"):
			query_key = line[10:-12]
		if line.startswith("<WebEnv>"):
			web_env = line[8:-10]

	# get fastas
	try:
		page_fasta = urllib2.urlopen('http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&query_key={0}&WebEnv={1}&rettype=fasta&retmode=text&retmax=100000'.format(query_key,web_env))
		print "Got fasta page"
	except:
		print "Failed to get fasta page"
	else:
		entries, fasta = parse_fasta(page_fasta,species)

		mfa = open(filename, 'w') # new fasta file for every species
		mfa.seek(0)
		mfa.write(fasta) # write fastas to file
		mfa.truncate()
		mfa.close()

	# If above method doesn't work, get chr directly
	if not os.path.isfile(filename) or os.path.getsize(filename) < 200:
		if 'refseq_ID_list' in globals():
			chr_id = refseq_ID_list[species_list.index(species_orig)]

			page_chr = urllib2.urlopen('http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore&term={0}'.format(chr_id))
			genome_id = ""
			for line in page_chr:
				line = string.replace(line,"\t","")
				if line.startswith("<Id>"):
					if genome_id != "":
						#print "Duplicate IDs for {0}; skipping.".format(species)
						return
					genome_id = line[4:-6]

			page_fasta = urllib2.urlopen('http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={0}&rettype=fasta&retmode=text&retmax=100000'.format(genome_id))
			entries, fasta = parse_fasta(page_fasta,species)

			mfa = open(filename, 'w') # overwrite fasta file
			mfa.seek(0)
			mfa.write(fasta) # write fastas to file
			mfa.truncate()
			mfa.close()
		else: # try looking up nucleotide directly
			page_list = urllib2.urlopen('http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi/?db=nuccore&term={0}&retmax=100000'.format(urllib2.quote(species)))
			id_list = []
			for line in page_list:
				line = string.replace(line,"\t","")
				if line.startswith("<Id>"):
					id_list.append(line[4:-6])

			all_ids = ','.join(id_list)
			print all_ids
			fasta_page = urllib2.urlopen('http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={0}&rettype=fasta&retmode=text&retmax=100000'.format(all_ids))

			entries, fasta = parse_fasta(fasta_page,species)

			write_fasta(filename,fasta)

	if not os.path.isfile(filename) or os.path.getsize(filename) < 200:
		# complete failure
		try:
			print "UTTERLY FAILED: {0} {1}".format(refseq_ID_list[species_list.index(species_orig)],species)
		except:
			print "UTTERLY FAILED: {0}".format(species)


species_list = [
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
	'Prochlorococcus marinus str. MIT9313', # removed space after MIT
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
	'Mycobacterium avium subsp. paratuberculosis K10', # removed -
	'Listeria welshimeri serovar 6b str. SLCC5334',
	'Bacillus halodurans C-125',
	'Geobacter sulfurreducens PCA',
	'Chlamydophila pneumoniae CWL029',
	'Francisella tularensis subsp. tularensis FSC 198', # added space after FSC
	'Staphylococcus aureus subsp. aureus NCTC 8325',
	'Mycobacterium marinum M',
	'Bacillus subtilis subsp. subtilis str. 168',
	'Staphylococcus aureus subsp. aureus MRSA 252', # added space after MRSA
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
	'Escherichia coli str. K12 substr. W3110',
	'Escherichia coli str. K-12 substr. MG1655', #added -
	'Chlorobium luteolum DSM 273', #new genus
	'pseudomonas fluorescens pf0-1',
	'Rhodopseudomonas palustris HaA2',
	'Burkholderia ambifaria AMMD']

'''
tax_id_list = [ #NCBI taxonomy id
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
'''
#i400

i400_refseq_ID_list = [
	"NC_009925",
	"NC_010163",
	"NC_009484",
	"NC_008009",
	"NC_008578",
	"NC_008782",
	"NC_009085",
	"NC_005966",
	"NC_009053",
	"NC_009655",
	"NC_008570",
	"NC_009348",
	"NC_008260",
	"NC_008340",
	"NC_009922",
	"NC_007760",
	"NC_009675",
	"NC_004842",
	"NC_007797",
	"NC_000918",
	"NC_000917",
	"NC_009850",
	"NC_008711",
	"NC_008541",
	"NC_007716",
	"NC_006513",
	"NC_009937",
	"NC_009725",
	"NC_007530",
	"NC_006274",
	"NC_002570",
	"NC_006270",
	"NC_009848",
	"NC_008600",
	"NC_010184",
	"NC_006347",
	"NC_009614",
	"NC_008783",
	"NC_005956",
	"NC_005955",
	"NC_010161",
	"NC_008618",
	"NC_010816",
	"NC_010645",
	"NC_002928",
	"NC_002929",
	"NC_010170",
	"NC_008277",
	"NC_001318",
	"NC_006156",
	"NC_004463",
	"NC_009485",
	"NC_010103",
	"NC_003317",
	"NC_009505",
	"NC_004310",
	"NC_002528",
	"NC_010508",
	"NC_009080",
	"NC_010084",
	"NC_010622",
	"NC_007434",
	"NC_007650",
	"NC_009256",
	"NC_007951",
	"NC_009437",
	"NC_009954",
	"NC_009715",
	"NC_009714",
	"NC_003912",
	"NC_007292",
	"NC_008512",
	"NC_010424",
	"NC_010482",
	"NC_007205",
	"NC_005861",
	"NC_008610",
	"NC_009465",
	"NC_007503",
	"NC_002696",
	"NC_002620",
	"NC_000117",
	"NC_003361",
	"NC_002179",
	"NC_002932",
	"NC_010175",
	"NC_005085",
	"NC_007963",
	"NC_009792",
	"NC_003030",
	"NC_009617",
	"NC_010520",
	"NC_009089",
	"NC_009706",
	"NC_008593",
	"NC_008261",
	"NC_010001",
	"NC_004557",
	"NC_009012",
	"NC_003910",
	"NC_002935",
	"NC_004369",
	"NC_007164",
	"NC_010545",
	"NC_009727",
	"NC_010528",
	"NC_010546",
	"NC_008255",
	"NC_007298",
	"NC_002936",
	"NC_007356",
	"NC_001263",
	"NC_010002",
	"NC_009943",
	"NC_006138",
	"NC_009253",
	"NC_007519",
	"NC_002937",
	"NC_009446",
	"NC_009952",
	"NC_007799",
	"NC_006832",
	"NC_009436",
	"NC_010694",
	"NC_007722",
	"NC_010556",
	"NC_010376",
	"NC_009441",
	"NC_009613",
	"NC_010336",
	"NC_008278",
	"NC_003454",
	"NC_006510",
	"NC_009328",
	"NC_010814",
	"NC_002939",
	"NC_009483",
	"NC_005125",
	"NC_010125",
	"NC_006677",
	"NC_008571",
	"NC_008343",
	"NC_002940",
	"NC_009567",
	"NC_010519",
	"NC_006396",
	"NC_010364",
	"NC_002607",
	"NC_008789",
	"NC_008229",
	"NC_004917",
	"NC_000915",
	"NC_010337",
	"NC_009138",
	"NC_008818",
	"NC_008358",
	"NC_006512",
	"NC_009776",
	"NC_007802",
	"NC_009659",
	"NC_009664",
	"NC_009648",
	"NC_006814",
	"NC_008497",
	"NC_008526",
	"NC_008529",
	"NC_010610",
	"NC_008530",
	"NC_010080",
	"NC_004567",
	"NC_009513",
	"NC_007576",
	"NC_008527",
	"NC_008011",
	"NC_006368",
	"NC_006087",
	"NC_010602",
	"NC_008508",
	"NC_004342",
	"NC_010524",
	"NC_010471",
	"NC_008531",
	"NC_003212",
	"NC_008555",
	"NC_007626",
	"NC_006300",
	"NC_008347",
	"NC_008740",
	"NC_006055",
	"NC_002678",
	"NC_008254",
	"NC_009515",
	"NC_000909",
	"NC_007955",
	"NC_009635",
	"NC_009135",
	"NC_009634",
	"NC_008942",
	"NC_009051",
	"NC_003551",
	"NC_008553",
	"NC_007355",
	"NC_003901",
	"NC_007681",
	"NC_007796",
	"NC_000916",
	"NC_007947",
	"NC_010172",
	"NC_010505",
	"NC_010511",
	"NC_002977",
	"NC_010296",
	"NC_007644",
	"NC_008769",
	"NC_009338",
	"NC_002677",
	"NC_010612",
	"NC_008596",
	"NC_008705",
	"NC_008611",
	"NC_008726",
	"NC_009497",
	"NC_007633",
	"NC_004829",
	"NC_000908",
	"NC_006908",
	"NC_005364",
	"NC_004432",
	"NC_002771",
	"NC_007294",
	"NC_008095",
	"NC_005213",
	"NC_007426",
	"NC_002946",
	"NC_003112",
	"NC_007798",
	"NC_007964",
	"NC_007406",
	"NC_007484",
	"NC_004757",
	"NC_008344",
	"NC_010085",
	"NC_007614",
	"NC_006361",
	"NC_008699",
	"NC_010628",
	"NC_003272",
	"NC_007794",
	"NC_004193",
	"NC_009667",
	"NC_008528",
	"NC_005303",
	"NC_010571",
	"NC_009488",
	"NC_009615",
	"NC_008686",
	"NC_009719",
	"NC_002663",
	"NC_007498",
	"NC_007512",
	"NC_010003",
	"NC_006370",
	"NC_005126",
	"NC_005877",
	"NC_008781",
	"NC_007948",
	"NC_010531",
	"NC_009379",
	"NC_010729",
	"NC_006085",
	"NC_008228",
	"NC_007481",
	"NC_009656",
	"NC_004129",
	"NC_002947",
	"NC_009434",
	"NC_004578",
	"NC_007204",
	"NC_007969",
	"NC_009524",
	"NC_008709",
	"NC_003364",
	"NC_009376",
	"NC_009073",
	"NC_008701",
	"NC_000868",
	"NC_003413",
	"NC_000961",
	"NC_007347",
	"NC_007973",
	"NC_003295",
	"NC_010168",
	"NC_007761",
	"NC_008380",
	"NC_007493",
	"NC_008268",
	"NC_007908",
	"NC_005027",
	"NC_008435",
	"NC_007643",
	"NC_009881",
	"NC_007940",
	"NC_009879",
	"NC_007109",
	"NC_009900",
	"NC_000963",
	"NC_009882",
	"NC_006142",
	"NC_009767",
	"NC_009523",
	"NC_008209",
	"NC_008148",
	"NC_007912",
	"NC_009142",
	"NC_007677",
	"NC_009953",
	"NC_009380",
	"NC_003197",
	"NC_009832",
	"NC_008700",
	"NC_009997",
	"NC_007954",
	"NC_008345",
	"NC_010334",
	"NC_009092",
	"NC_009901",
	"NC_009438",
	"NC_009831",
	"NC_008577",
	"NC_010506",
	"NC_010658",
	"NC_007606",
	"NC_004337",
	"NC_003911",
	"NC_008044",
	"NC_009636",
	"NC_003047",
	"NC_009511",
	"NC_008048",
	"NC_007793",
	"NC_002976",
	"NC_007168",
	"NC_007350",
	"NC_009033",
	"NC_004368",
	"NC_009785",
	"NC_004350",
	"NC_010582",
	"NC_008024",
	"NC_009009",
	"NC_009442",
	"NC_008532",
	"NC_003155",
	"NC_003888",
	"NC_007181",
	"NC_002754",
	"NC_010730",
	"NC_007575",
	"NC_009663",
	"NC_006177",
	"NC_007604",
	"NC_010475",
	"NC_008346",
	"NC_010321",
	"NC_010320",
	"NC_003869",
	"NC_007333",
	"NC_008698",
	"NC_002578",
	"NC_002689",
	"NC_010525",
	"NC_009616",
	"NC_004113",
	"NC_009828",
	"NC_009486",
	"NC_010483",
	"NC_005835",
	"NC_007404",
	"NC_007520",
	"NC_002967",
	"NC_000919",
	"NC_008312",
	"NC_004572",
	"NC_009464",
	"NC_010503",
	"NC_008786",
	"NC_009456",
	"NC_006840",
	"NC_009783",
	"NC_004603",
	"NC_005139",
	"NC_002978",
	"NC_005090",
	"NC_009720",
	"NC_003919",
	"NC_007705",
	"NC_002488",
	"NC_008800",
	"NC_009708",
	"NC_006526"
]

i400_species_list = [
	"Acaryochloris marina MBIC11017",
	"Acholeplasma laidlawii PG-8A",
	"Acidiphilium cryptum JF-5",
	"Candidatus Koribacter versatilis Ellin345", #"Acidobacteria bacterium Ellin345",
	"Acidothermus cellulolyticus 11B",
	"Acidovorax sp. JS42",
	"Acinetobacter baumannii ATCC 17978",
	"Acinetobacter sp. ADP1",
	"Actinobacillus pleuropneumoniae L20",
	"Actinobacillus succinogenes 130Z",
	"Aeromonas hydrophila subsp. hydrophila ATCC 7966",
	"Aeromonas salmonicida subsp. salmonicida A449",
	"Alcanivorax borkumensis SK2",
	"Alkalilimnicola ehrlichii MLHE-1", #"Alkalilimnicola ehrlichei MLHE-1",
	"Alkaliphilus oremlandii OhILAs",
	"Anaeromyxobacter dehalogenans 2CP-C",
	"Anaeromyxobacter sp. Fw109-5",
	"Anaplasma marginale str. St. Maries",
	"Anaplasma phagocytophilum HZ",
	"Aquifex aeolicus VF5",
	"Archaeoglobus fulgidus DSM 4304",
	"Arcobacter butzleri RM4018",
	"Arthrobacter aurescens TC1",
	"Arthrobacter sp. FB24",
	"Aster yellows witches'-broom phytoplasma AYWB",
	"Azoarcus aromaticum EbN1", #"Azoarcus sp. EbN1",
	"Azorhizobium caulinodans ORS 571",
	"Bacillus amyloliquefaciens subsp. plantarum str. FZB42", #"Bacillus amyloliquefaciens FZB42",
	"Bacillus anthracis str. 'Ames Ancestor'",
	"Bacillus cereus E33L",
	"Bacillus halodurans C-125",
	"Bacillus licheniformis ATCC 14580",
	"Bacillus pumilus SAFR-032",
	"Bacillus thuringiensis str. Al Hakam",
	"Bacillus weihenstephanensis KBAB4",
	"Bacteroides fragilis YCH46",
	"Bacteroides vulgatus ATCC 8482",
	"Bartonella bacilliformis KC583",
	"Bartonella henselae strain Houston-1", #"Bartonella henselae str. Houston-1",
	"Bartonella quintana str. Toulouse",
	"Bartonella tribocorum CIP 105476",
	"Bifidobacterium adolescentis ATCC 15703",
	"Bifidobacterium longum DJO10A",
	"Bordetella avium 197N",
	"Bordetella parapertussis strain 12822", #"Bordetella parapertussis 12822",
	"Bordetella pertussis Tohama I",
	"Bordetella petrii strain DSM 12804", #"Bordetella petrii DSM 12804",
	"Borrelia afzelii PKo",
	"Borrelia burgdorferi B31",
	"Borrelia garinii PBi",
	"Bradyrhizobium japonicum USDA 110",
	"Bradyrhizobium sp. BTAi1",
	"Brucella canis ATCC 23365",
	"Brucella melitensis bv. 1 str. 16M", #"Brucella melitensis 16M",
	"Brucella ovis ATCC 25840",
	"Brucella suis 1330",
	"Buchnera aphidicola str. APS (Acyrthosiphon pisum)",
	"Burkholderia cenocepacia MC0-3",
	"Burkholderia mallei NCTC 10247",
	"Burkholderia multivorans ATCC 17616",
	"Burkholderia phymatum STM815",
	"Burkholderia pseudomallei 1710b",
	"Burkholderia thailandensis E264",
	"Burkholderia vietnamiensis G4",
	"Burkholderia xenovorans LB400",
	"Caldicellulosiruptor saccharolyticus DSM 8903",
	"Caldivirga maquilingensis IC-167",
	"Campylobacter curvus 525.92",
	"Campylobacter hominis ATCC BAA-381",
	"Campylobacter jejuni RM1221",
	"Candidatus Blochmannia pennsylvanicus str. BPEN",
	"Candidatus Carsonella ruddii PV",
	"Candidatus Desulforudis audaxviator MP104C",
	"Candidatus Korarchaeum cryptofilum OPF8",
	"Candidatus Pelagibacter ubique HTCC1062",
	"Parachlamydia-related symbiont UWE25", #"Candidatus Protochlamydia amoebophila UWE25",
	"Candidatus Ruthia magnifica str. Cm (Calyptogena magnifica)",
	"Candidatus Vesicomyosocius okutanii HA",
	"Carboxydothermus hydrogenoformans Z-2901",
	"Caulobacter crescentus CB15",
	"Chlamydia muridarum Nigg",
	"Chlamydia trachomatis D/UW-3/CX",
	"Chlamydophila caviae GPIC",
	"Chlamydophila pneumoniae AR39",
	"Chlorobium tepidum TLS",
	"Chloroflexus aurantiacus J-10-fl",
	"Chromobacterium violaceum ATCC 12472",
	"Chromohalobacter salexigens DSM 3043",
	"Citrobacter koseri ATCC BAA-895",
	"Clostridium acetobutylicum ATCC 824",
	"Clostridium beijerinckii NCIMB 8052",
	"Clostridium botulinum A3 str. Loch Maree",
	"Clostridium difficile 630",
	"Clostridium kluyveri DSM 555",
	"Clostridium novyi NT",
	"Clostridium perfringens ATCC 13124",
	"Clostridium phytofermentans ISDg",
	"Clostridium tetani E88",
	"Clostridium thermocellum ATCC 27405",
	"Colwellia psychrerythraea 34H",
	"Corynebacterium diphtheriae NCTC 13129",
	"Corynebacterium efficiens YS-314",
	"Corynebacterium jeikeium K411",
	"Corynebacterium urealyticum DSM 7109",
	"Coxiella burnetii Dugway 5J108-111",
	"Cupriavidus taiwanensis",
	"Cyanothece sp. ATCC 51142",
	"Cytophaga hutchinsonii ATCC 33406",
	"Dechloromonas aromatica RCB",
	"Dehalococcoides ethenogenes 195",
	"Dehalococcoides sp. CBDB1",
	"Deinococcus radiodurans R1",
	"Delftia acidovorans SPH-1",
	"Desulfococcus oleovorans Hxd3",
	"Desulfotalea psychrophila LSv54",
	"Desulfotomaculum reducens MI-1",
	"Desulfovibrio alaskensis G20", #"Desulfovibrio desulfuricans subsp. desulfuricans str. G20",
	"Desulfovibrio vulgaris subsp. vulgaris str. Hildenborough",
	"Dichelobacter nodosus VCS1703A",
	"Dinoroseobacter shibae DFL 12",
	"Ehrlichia chaffeensis str. Arkansas",
	"Ehrlichia ruminantium str. Welgevonden",
	"Enterobacter sp. 638",
	"Erwinia tasmaniensis",
	"Erythrobacter litoralis HTCC2594",
	"Exiguobacterium sibiricum 255-15",
	"Finegoldia magna ATCC 29328",
	"Flavobacterium johnsoniae UW101",
	"Flavobacterium psychrophilum JIP02/86",
	"Francisella philomiragia subsp. philomiragia ATCC 25017",
	"Frankia alni str. ACN14a", #"Frankia alni ACN14a",
	"Fusobacterium nucleatum subsp. nucleatum ATCC 25586",
	"Geobacillus kaustophilus HTA426",
	"Geobacillus thermodenitrificans NG80-2",
	"Geobacter lovleyi SZ",
	"Geobacter sulfurreducens PCA",
	"Geobacter uraniireducens Rf4",
	"Gloeobacter violaceus PCC 7421",
	"Gluconacetobacter diazotrophicus PAl 5",
	"Gluconobacter oxydans 621H",
	"Gramella forsetii KT0803",
	"Granulibacter bethesdensis CGDNIH1",
	"Haemophilus ducreyi strain 35000HP", #"Haemophilus ducreyi 35000HP",
	"Haemophilus influenzae PittGG",
	"Haemophilus somnus 2336",
	"Haloarcula marismortui ATCC 43049",
	"Halobacterium salinarum R1",
	"Halobacterium sp. NRC-1",
	"Halorhodospira halophila SL1",
	"Helicobacter acinonychis str. Sheeba",
	"Helicobacter hepaticus ATCC 51449",
	"Helicobacter pylori 26695",
	"Heliobacterium modesticaldum Ice1",
	"Herminiimonas arsenicoxydans",
	"Hyperthermus butylicus DSM 5456",
	"Hyphomonas neptunium ATCC 15444",
	"Idiomarina loihiensis L2TR",
	"Ignicoccus hospitalis KIN4/I",
	"Jannaschia sp. CCS1",
	"Janthinobacterium sp. Marseille",
	"Kineococcus radiotolerans SRS30216",
	"Klebsiella pneumoniae subsp. pneumoniae MGH 78578",
	"Lactobacillus acidophilus NCFM",
	"Lactobacillus brevis ATCC 367",
	"Lactobacillus casei ATCC 334",
	"Lactobacillus delbrueckii subsp. bulgaricus ATCC BAA-365",
	"Lactobacillus fermentum IFO 3956",
	"Lactobacillus gasseri ATCC 33323",
	"Lactobacillus helveticus DPC 4571",
	"Lactobacillus plantarum WCFS1",
	"Lactobacillus reuteri DSM 20016", #"Lactobacillus reuteri F275",
	"Lactobacillus sakei strain 23K", #"Lactobacillus sakei subsp. sakei 23K",
	"Lactococcus lactis subsp. cremoris SK11",
	"Lawsonia intracellularis PHE/MN1-00",
	"Legionella pneumophila str. Paris",
	"Leifsonia xyli subsp. xyli str. CTCB07",
	"Leptospira biflexa serovar Patoc strain 'Patoc 1 (Paris)'",
	"Leptospira borgpetersenii serovar Hardjo-bovis L550",
	"Leptospira interrogans serovar Lai str. 56601",
	"Leptothrix cholodnii SP-6",
	"Leuconostoc citreum KM20",
	"Leuconostoc mesenteroides subsp. mesenteroides ATCC 8293",
	"Listeria innocua Clip11262",
	"Listeria welshimeri serovar 6b str. SLCC5334",
	"Magnetospirillum magneticum AMB-1",
	"Mannheimia succiniciproducens MBEL55E",
	"Maricaulis maris MCS10",
	"Marinobacter aquaeolei VT8",
	"Mesoplasma florum L1",
	"Mesorhizobium loti MAFF303099",
	"Chelativorans sp. BNC1", #"Mesorhizobium sp. BNC1",
	"Methanobrevibacter smithii ATCC 35061",
	"Methanocaldococcus jannaschii DSM 2661",
	"Methanococcoides burtonii DSM 6242",
	"Methanococcus aeolicus Nankai-3",
	"Methanococcus maripaludis C5",
	"Methanococcus vannielii SB",
	"Methanocorpusculum labreanum Z",
	"Methanoculleus marisnigri JR1",
	"Methanopyrus kandleri AV19",
	"Methanosaeta thermophila PT",
	"Methanosarcina barkeri str. Fusaro",
	"Methanosarcina mazei Go1",
	"Methanosphaera stadtmanae DSM 3091",
	"Methanospirillum hungatei JF-1",
	"Methanothermobacter thermautotrophicus str. Delta H",
	"Methylobacillus flagellatus KT",
	"Methylobacterium extorquens PA1",
	"Methylobacterium radiotolerans JCM 2831",
	"Methylobacterium sp. 4-46",
	"Methylococcus capsulatus str. Bath",
	"Microcystis aeruginosa NIES-843",
	"Moorella thermoacetica ATCC 39073",
	"Mycobacterium bovis BCG Pasteur 1173P2", #"Mycobacterium bovis BCG str. Pasteur 1173P2",
	"Mycobacterium gilvum PYR-GCK",
	"Mycobacterium leprae TN",
	"Mycobacterium marinum M",
	"Mycobacterium smegmatis str. MC2 155",
	"Mycobacterium sp. KMS",
	"Mycobacterium ulcerans Agy99",
	"Mycobacterium vanbaalenii PYR-1",
	"Mycoplasma agalactiae PG2",
	"Mycoplasma capricolum subsp. capricolum ATCC 27343",
	"Mycoplasma gallisepticum str. R(low)", #"Mycoplasma gallisepticum R",
	"Mycoplasma genitalium G37",
	"Mycoplasma mobile 163K",
	"Mycoplasma mycoides subsp. mycoides SC str. PG1",
	"Mycoplasma penetrans HF-2",
	"Mycoplasma pulmonis UAB CTIP",
	"Mycoplasma synoviae 53",
	"Myxococcus xanthus DK 1622",
	"Nanoarchaeum equitans Kin4-M",
	"Natronomonas pharaonis DSM 2160",
	"Neisseria gonorrhoeae FA 1090",
	"Neisseria meningitidis MC58",
	"Neorickettsia sennetsu strain Miyayama", #"Neorickettsia sennetsu str. Miyayama",
	"Nitrobacter hamburgensis X14",
	"Nitrobacter winogradskyi Nb-255",
	"Nitrosococcus oceani ATCC 19707",
	"Nitrosomonas europaea ATCC 19718",
	"Nitrosomonas eutropha C91",
	"Nitrosopumilus maritimus SCM1",
	"Nitrosospira multiformis ATCC 25196",
	"Nocardia farcinica IFM 10152",
	"Nocardioides sp. JS614",
	"Nostoc punctiforme PCC 73102",
	"Nostoc sp. PCC 7120",
	"Novosphingobium aromaticivorans DSM 12444",
	"Oceanobacillus iheyensis HTE831",
	"Ochrobactrum anthropi ATCC 49188",
	"Oenococcus oeni PSU-1",
	"Onion yellows phytoplasma OY-M",
	"Opitutus terrae PB90-1",
	"Orientia tsutsugamushi str. Boryong",
	"Parabacteroides distasonis ATCC 8503",
	"Paracoccus denitrificans PD1222",
	"Parvibaculum lavamentivorans DS-1",
	"Pasteurella multocida subsp. multocida str. Pm70",
	"Pelobacter carbinolicus DSM 2380",
	"Chlorobium luteolum DSM 273", #"Pelodictyon luteolum DSM 273",
	"Petrotoga mobilis SJ95",
	"Photobacterium profundum SS9",
	"Photorhabdus luminescens subsp. laumondii TTO1",
	"Picrophilus torridus DSM 9790",
	"Polaromonas naphthalenivorans CJ2",
	"Polaromonas sp. JS666",
	"Polynucleobacter necessarius subsp. necessarius STIR1", #"Polynucleobacter necessarius STIR1",
	"Polynucleobacter necessarius subsp. asymbioticus QLW-P1DMWA-1", #"Polynucleobacter sp. QLW-P1DMWA-1",
	"Porphyromonas gingivalis ATCC 33277",
	"Propionibacterium acnes KPA171202",
	"Pseudoalteromonas atlantica T6c",
	"Pseudoalteromonas haloplanktis str. TAC125", #"Pseudoalteromonas haloplanktis TAC125",
	"Pseudomonas aeruginosa PA7",
	"Pseudomonas protegens Pf-5", #"Pseudomonas fluorescens Pf-5",
	"Pseudomonas putida KT2440",
	"Pseudomonas stutzeri A1501",
	"Pseudomonas syringae pv. tomato str. DC3000",
	"Psychrobacter arcticus 273-4",
	"Psychrobacter cryohalolentis K5",
	"Psychrobacter sp. PRwf-1",
	"Psychromonas ingrahamii 37",
	"Pyrobaculum aerophilum str. IM2",
	"Pyrobaculum arsenaticum DSM 13514",
	"Pyrobaculum calidifontis JCM 11548",
	"Pyrobaculum islandicum DSM 4184",
	"Pyrococcus abyssi GE5",
	"Pyrococcus furiosus DSM 3638",
	"Pyrococcus horikoshii OT3",
	"Ralstonia eutropha JMP134",
	"Cupriavidus metallidurans CH34", #"Ralstonia metallidurans CH34",
	"Ralstonia solanacearum GMI1000",
	"Renibacterium salmoninarum ATCC 33209",
	"Rhizobium etli CFN 42",
	"rhizobium leguminosarum bv. viciae chromosome complete genome, strain 3841", #"Rhizobium leguminosarum bv. viciae 3841",
	"Rhodobacter sphaeroides 2.4.1",
	"Rhodococcus jostii RHA1", #"Rhodococcus sp. RHA1",
	"Rhodoferax ferrireducens T118",
	"Rhodopirellula baltica SH 1",
	"Rhodopseudomonas palustris BisA53",
	"Rhodospirillum rubrum ATCC 11170",
	"Rickettsia akari str. Hartford",
	"Rickettsia bellii RML369-C",
	"Rickettsia canadensis str. McKiel",
	"Rickettsia felis URRWXCal2",
	"Rickettsia massiliae MTU5",
	"Rickettsia prowazekii str. Madrid E",
	"Rickettsia rickettsii str. 'Sheila Smith'",
	"Rickettsia typhi str. Wilmington",
	"Roseiflexus castenholzii DSM 13941",
	"Roseiflexus sp. RS-1",
	"Roseobacter denitrificans OCh 114",
	"Rubrobacter xylanophilus DSM 9941",
	"Saccharophagus degradans 2-40",
	"Saccharopolyspora erythraea NRRL2338", #"Saccharopolyspora erythraea NRRL 2338",
	"Salinibacter ruber DSM 13855",
	"Salinispora arenicola CNS-205",
	"Salinispora tropica CNB-440",
	"Salmonella enterica subsp. enterica serovar typhimurium str. LT2", #"Salmonella typhimurium LT2",
	"Serratia proteamaculans 568",
	"Shewanella amazonensis SB2B",
	"Shewanella baltica OS195",
	"Shewanella denitrificans OS217",
	"Shewanella frigidimarina NCIMB 400",
	"Shewanella halifaxensis HAW-EB4",
	"Shewanella loihica PV-4",
	"Shewanella pealeana ATCC 700345",
	"Shewanella putrefaciens CN-32",
	"Shewanella sediminis HAW-EB3",
	"Shewanella sp. ANA-3",
	"Shewanella woodyi ATCC 51908",
	"Shigella boydii CDC 3083-94",
	"Shigella dysenteriae Sd197",
	"Shigella flexneri 2a str. 301",
	"Ruegeria pomeroyi DSS-3", #"Silicibacter pomeroyi DSS-3",
	"Silicibacter sp. TM1040",
	"Sinorhizobium medicae WSM419",
	"Sinorhizobium meliloti 1021",
	"Sphingomonas wittichii RW1",
	"Sphingopyxis alaskensis RB2256",
	"Staphylococcus aureus subsp. aureus USA300",
	"Staphylococcus epidermidis RP62A",
	"Staphylococcus haemolyticus JCSC1435",
	"Staphylococcus saprophyticus subsp. saprophyticus ATCC 15305",
	"Staphylothermus marinus F1",
	"Streptococcus agalactiae NEM316",
	"Streptococcus gordonii str. Challis substr. CH1",
	"Streptococcus mutans UA159",
	"Streptococcus pneumoniae CGSP14",
	"Streptococcus pyogenes MGAS10750",
	"Streptococcus sanguinis SK36",
	"Streptococcus suis 05ZYH33",
	"Streptococcus thermophilus LMD-9",
	"Streptomyces avermitilis MA-4680",
	"Streptomyces coelicolor A3(2)",
	"Sulfolobus acidocaldarius DSM 639",
	"Sulfolobus solfataricus P2",
	"Sulfurihydrogenibium sp. YO3AOP1",
	"Sulfurimonas denitrificans DSM 1251",
	"Sulfurovum sp. NBC37-1",
	"Symbiobacterium thermophilum IAM 14863",
	"Synechococcus elongatus PCC 7942",
	"Synechococcus sp. PCC 7002",
	"Syntrophomonas wolfei subsp. wolfei str. Goettingen",
	"Thermoanaerobacter pseudethanolicus ATCC 33223",
	"Thermoanaerobacter sp. X514",
	"Thermoanaerobacter tengcongensis MB4",
	"Thermobifida fusca YX",
	"Thermofilum pendens Hrk 5",
	"Thermoplasma acidophilum DSM 1728",
	"Thermoplasma volcanium GSS1",
	"Thermoproteus neutrophilus V24Sta",
	"Thermosipho melanesiensis BI429",
	"Thermosynechococcus elongatus BP-1",
	"Thermotoga lettingae TMO",
	"Thermotoga petrophila RKU-1",
	"Thermotoga sp. RQ2",
	"Thermus thermophilus HB27",
	"Thiobacillus denitrificans ATCC 25259",
	"Thiomicrospira crunogena XCL-2",
	"Treponema denticola ATCC 35405",
	"Treponema pallidum subsp. pallidum str. Nichols",
	"Trichodesmium erythraeum IMS101",
	"Tropheryma whipplei str. Twist",
	"Methanocella arvoryzae MRE50", #"uncultured methanogenic archaeon RC-I",
	"Ureaplasma parvum serovar 3 str. ATCC 27815",
	"Verminephrobacter eiseniae EF01-2",
	"Vibrio cholerae O395",
	"Vibrio fischeri ES114",
	"Vibrio harveyi ATCC BAA-1116",
	"Vibrio parahaemolyticus RIMD 2210633",
	"Vibrio vulnificus YJ016",
	"Wolbachia endosymbiont of Drosophila melanogaster",
	"Wolinella succinogenes DSM 1740",
	"Xanthobacter autotrophicus Py2",
	"Xanthomonas axonopodis pv. citri str. 306",
	"Xanthomonas oryzae pv. oryzae MAFF 311018",
	"Xylella fastidiosa 9a5c",
	"Yersinia enterocolitica subsp. enterocolitica 8081",
	"Yersinia pseudotuberculosis IP 31758",
	"Zymomonas mobilis subsp. mobilis ZM4"
	]

def main(argv=sys.argv):
	"""
	Command line usage: python M_get_genome_from_NCBI.py
						[species (optional)]
	"""
	if len(argv) > 1:
		incoming_arg = argv[1] # can enter single species for lookup on command line

		if os.path.exists(incoming_arg):
			print "Taking species to look up from {0}".format(incoming_arg)
			species_source = []
			with open(incoming_arg,'r') as sp_file:
				for l in sp_file:
					species_source.append(l.rstrip())
		else:
			species_source = [incoming_arg]
			print "Species to look up set as {0}".format(species_source)
	else:
		species_source = species_list

	for s in species_source:
		get_genome(s)

if __name__ == "__main__":
    main()
