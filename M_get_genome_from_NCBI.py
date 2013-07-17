# Download genomes of listed species, avoiding duplicates

import urllib2
import string

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
		print "ERROR: Empty page"
		return False

	return line # normal line of fasta

def get_genome(species):
	print "SPECIES: {0}".format(species)
	species = species.lower()

	mfa = open(species.replace (" ", "_").replace('/',"") + '.mfa','w') # new fasta file for every species
	mfa.seek(0)

	# get genome ID
	page_genome = urllib2.urlopen('http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=genome&term=' + urllib2.quote(species))
	genome_id = ""
	for line in page_genome:
		if line.startswith("<Id>"):
			if genome_id != "":
				print "Duplicate IDs for {0}; skipping.".format(species)
				return
			genome_id = line[4:-6]

	# get query key and web env for chromosome
	page_nucleotide = urllib2.urlopen('http://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=genome&db=nuccore&id=' + genome_id + '&term=gene+in+chromosome[prop]&cmd=neighbor_history')
	for line in page_nucleotide:
		line = string.replace(line,"\t","")
		if line.startswith("<QueryKey>"):
			query_key = line[10:-12]
		if line.startswith("<WebEnv>"):
			web_env = line[8:-10]

	# get fasta for chromosome
	page_fasta = urllib2.urlopen('http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&query_key=' + query_key + '&WebEnv=' + web_env + '&rettype=fasta&retmode=text')

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

species= [
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
'Leptospira borgpetersenii serovar Hardjo-bovis JB197']

more_species = [
'Herpetosiphon aurantiacus DSM 785', #new strain name
'Buchnera aphidicola BCc', #new strain name
'Rhodococcus jostii RHA1', #new species name
'Escherichia coli str. K-12 substr. W3110', #added -
'Escherichia coli str. K-12 substr. MG1655', #added -
'Chlorobium luteolum DSM 273'] #new genus

species_list = [
'pseudomonas fluorescens pf0-1',
'Rhodopseudomonas palustris HaA2',
'Burkholderia ambifaria AMMD']

for s in species_list:
	get_genome(s)

