# Python3 script which takes in an annotation file(gtf/gff3) and a transcriptomic fasta file
# and produces an sqlite file which can be uploaded to Trips-Viz
# All co-ordinates produced are 1 based
# All start codon positions (including cds_start) should be at the first nucleotide of the codon
# All stop codon positions (including cds_stop) should be at the last nucleotide of the codon
import sys
import re
import sqlite3
from intervaltree import Interval, IntervalTree
import itertools
from sqlitedict import SqliteDict
import os

# python create_annotation_sqlite.py /data2/Jack/homo_sapiens/Homo_sapiens.GRCh38.108_ensembl.gtf \
#  /data2/Jack/homo_sapiens/transcript_ensembl_108.fa \
# 0 ENST00000673477 DDX11L17 \
#  -o /data2/Jack/homo_sapiens/homo_sapiens.ensembl108.l108.sqlite

organism = sys.argv[1]
#This should be a GTF or GFF3 file
annotation_file = open(sys.argv[2],"r")
#This needs to be the transcriptomic fasta file
fasta_file = open(sys.argv[3],"r")
#This value will be added used to create UTRs of this length, useful when looking at transcriptomes without annotated UTRs
pseudo_utr_len = int(sys.argv[4])
#An example of a transcript_id from the annotation file, e.g ENST000000123456
user_transcript_id = sys.argv[5]
#An example of a gene name from the annotation file
user_gene_name = sys.argv[6]
# Set to true if transcript version is included in transcript_id, e.g: ENST000000123456.1
TRAN_VERSION = False


if os.path.isfile("{}.sqlite".format(organism)):
	print ("{}.sqlite already exists".format(organism))
	sys.exit()


# old_exons = SqliteDict("/home/DATA/www/tripsviz/tripsviz/trips_annotations/mus_musculus/transcriptomic_to_genomic.sqlite")



delimiters = {"transcripts":{"before":[],"after":['.'],"annot_types": ["cds","utr"]},
			  "genes":{"before":[],"after":['"'],"annot_types": ["lnc_rna"]}}

punctuation = [";"," ","-",":","-",".","=","\t"]
#Find delimiters in the annotation and fasta files using the user_transcript_id
#and user_gene_name examples given by user.
for line in annotation_file:
	if user_transcript_id in line:
		tabsplitline = line.split("\t")
		annot_type = tabsplitline[2]
		if annot_type not in delimiters["transcripts"]["annot_types"]:
			delimiters["transcripts"]["annot_types"].append(annot_type.lower())
		splitline = line.split(user_transcript_id)
		before_delimiter = splitline[0]
		for item in punctuation:
			if item in before_delimiter:
				if len(before_delimiter.split(item)[-1]) >= 5:
					before_delimiter = before_delimiter.split(item)[-1]
		after_delimiter = splitline[1][:2]
		if before_delimiter not in delimiters["transcripts"]["before"] and len(before_delimiter) >= 5:
			delimiters["transcripts"]["before"].append(before_delimiter)
		if after_delimiter not in delimiters["transcripts"]["after"]:
			delimiters["transcripts"]["after"].append(after_delimiter)
	if user_gene_name in line:
		tabsplitline = line.split("\t")
		annot_type = tabsplitline[2]
		print ("ANNOT TYPE",annot_type)
		if annot_type not in delimiters["genes"]["annot_types"]:
			delimiters["genes"]["annot_types"].append(annot_type.lower())
		splitline = line.split(user_gene_name)
		before_delimiter = splitline[0]
		for item in punctuation:
			if item in before_delimiter:
				if len(before_delimiter.split(item)[-1]) >= 5:
					before_delimiter = before_delimiter.split(item)[-1]
		after_delimiter = splitline[1][0]
		if before_delimiter not in delimiters["genes"]["before"] and len(before_delimiter) >= 5:
			delimiters["genes"]["before"].append(before_delimiter)
		if after_delimiter not in delimiters["genes"]["after"]:
			if after_delimiter in punctuation:
				delimiters["genes"]["after"].append(after_delimiter)

print ("delimeters[genes]",delimiters["transcripts"]["annot_types"])

for line in fasta_file:
	if user_transcript_id in line:
		splitline = line.split(user_transcript_id)
		before_delimiter = splitline[0]
		for item in punctuation:
			if item in before_delimiter:
				if len(before_delimiter.split(item)[1]) >= 5:
					before_delimiter = before_delimiter.split(item)[1]
		after_delimiter = splitline[1][0]
		if before_delimiter not in delimiters["transcripts"]["before"] and len(before_delimiter) >= 5:
			delimiters["transcripts"]["before"].append(before_delimiter)
		if after_delimiter not in delimiters["transcripts"]["after"]:
			delimiters["transcripts"]["after"].append(after_delimiter)
fasta_file.close()
annotation_file.close()





if delimiters["transcripts"]["before"] == []:
	print ("ERROR: No transcript_id with the name {} could be found in the annotation file".format(user_transcript_id))
	sys.exit()
if delimiters["genes"]["before"] == []:
	print ("ERROR: No gene with the name {} could be found in the annotation file".format(user_gene_name))
	sys.exit()

master_dict = {}
coding_dict = {}
notinfasta = open("notinfasta.csv","w")

#Given a nucleotide sequence returns the positions of all start and stop codons.
def get_start_stops(transcript_sequence):
	transcript_sequence = transcript_sequence.upper()
	start_codons = ['ATG']
	stop_codons = ['TAA', 'TAG', 'TGA']
	seq_frames = {'starts': [], 'stops': []}
	for codons, positions in ((start_codons, 'starts'),(stop_codons, 'stops')):
		if len(codons) > 1:
			pat = re.compile('|'.join(codons))
		else:
			pat = re.compile(codons[0])
		for m in re.finditer(pat, transcript_sequence):
			# Increment position by 1, Frame 1 starts at position 1 not 0,
			# if it's a stop codon add another 2 so it points to the last nuc of the codon
			if positions == "starts":
				start = m.start() + 1
			else:
				start = m.start() + 3
			seq_frames[positions].append(start)
	return seq_frames

#parse fasta to get the nucleotide sequence of transcripts and the positions of start/stop codons.
fasta_file = open(sys.argv[3],"r")
read_fasta = fasta_file.read()
split_fasta = read_fasta.split(">")
for entry in split_fasta[1:]:
	newline_split = entry.split("\n")
	tran = newline_split[0]
	for item in delimiters["transcripts"]["after"]:
		if item in tran:
			tran = tran.split(item)[0]
	tran = tran.replace("-","_").replace("(","").replace(")","")
	seq = ("".join(newline_split[1:]))
	if "_PAR_Y" in tran:
		tran += "_chrY"
	elif "_PAR_X" in tran:
		tran += "_chrX"
	tran = tran.upper()
	starts_stops = get_start_stops(seq)
	if tran not in master_dict:
		master_dict[tran] = {"utr":[], "cds":[], "exon":[],"start_codon":[],"stop_codon":[],"start_list":starts_stops["starts"],
							"stop_list":starts_stops["stops"],"transcript":[], "strand":"" ,"gene_name":"","chrom":"","seq":seq,"cds_start":"NULL","cds_stop":"NULL",
							"length":len(seq),"principal":0,"version":"NULL"}



def to_ranges(iterable):
	tup_list = []
	iterable = sorted(set(iterable))
	for key, group in itertools.groupby(enumerate(iterable),lambda t: t[1] - t[0]):
		group = list(group)
		tup_list.append((group[0][1], group[-1][1]))
	return tup_list

#parse annotation file to get chromsome, exon location and CDS info for each transcript
def parse_gtf_file(annotation_file):
	for line in annotation_file:
		if line == "\n":
			continue
		if line[0] != '#':
			splitline = (line.replace("\n","")).split("\t")
			chrom = splitline[0]
			try:
				annot_type = splitline[2].lower()
			except:
				print ("ERROR tried to index to second item in splitline: ",splitline, line)
				sys.exit()
			#if annot_type not in ["cds", "utr", "exon", "transcript","five_prime_utr", "three_prime_utr","stop_codon","start_codon"]:
			#	continue
			if annot_type not in delimiters["transcripts"]["annot_types"] and annot_type not in delimiters["genes"]["annot_types"]:
				continue
			if annot_type == "five_prime_utr" or annot_type == "three_prime_utr":
				annot_type = "utr"
			strand = splitline[6]
			if strand == "+":
				start = int(splitline[3])
				end = int(splitline[4])
			else:
				start = int(splitline[3])+1
				end = int(splitline[4])+1
			desc = splitline[8]
			tran = desc
			gene = desc
			for item in delimiters["transcripts"]["before"]:
				if item in tran:
					tran = tran.split(item)[1]
			for item in delimiters["transcripts"]["after"]:
				if item in tran:
					tran = tran.split(item)[0]
			if "." in tran and TRAN_VERSION == True:
				#print ("raw tran",tran)
				tran = tran.split(".")
				version = int(tran[-1].split("_")[0])
				tran = tran[0]
			else:
				version = "NULL"
			tran = tran.replace("-","_").replace(".","_")
			tran = tran.replace("(","").replace(")","")
			tran = tran.replace(" ","").replace("\t","")
			tran = tran.upper()
			tran = tran.replace("GENE_","").replace("ID_","")
			if "_PAR_Y" in desc:
				#print ("adding _PAR_Y to tran")
				tran = tran+"_PAR_Y"
				#print ("New tran ", tran)
			#if "PAR_Y" in line:
			#	print (line)
			#	#sys.exit()
			#print ("tran",tran,version)
			#if tran == "ENST00000316448":
			#	print ("annot type",annot_type)
			#	print ("appending exon to tran", start, end,line)

			gene_found = False

			if annot_type in delimiters["genes"]["annot_types"]:
				for item in delimiters["genes"]["before"]:
					if item in gene:
						gene_found = True
						gene = gene.split(item)[1]
				for item in delimiters["genes"]["after"]:
					if item in gene:
						gene = gene.split(item)[0]
				gene = gene.replace("'","''")
				gene = gene.replace("GENE_","")
				gene = gene.replace("ID_","")
				gene = gene.upper()
			if tran in master_dict:
				master_dict[tran]["strand"] = strand
				if strand == "+":
					if annot_type in master_dict[tran]:
        	                                master_dict[tran][annot_type].append((start, end))
				else:
					if annot_type in master_dict[tran]:
						master_dict[tran][annot_type].append((start, end))
				master_dict[tran]["chrom"] = chrom
				master_dict[tran]["version"] = version
				if gene_found == True:
					master_dict[tran]["gene_name"] = gene
			else:
				notinfasta.write("{}\n".format(tran))

annotation_file = open(sys.argv[2],"r")
parse_gtf_file(annotation_file)


#remove transcripts that were in fasta file but not in annotation_file
notinannotation = []
for tran in master_dict:
	if master_dict[tran]["chrom"] == "":
		#print ("tran {} has no chrom :(".format(tran))
		notinannotation.append(tran)
for tran in notinannotation:
	del master_dict[tran]
#Dictionary to store the coding status of a gene, if any transcript of this gene is coding, the value will be True
coding_genes_dict = {}
#parse master_dict to calculate length, cds_start/stop and exon junction positions
for tran in master_dict:
	try:
		transeq = master_dict[tran]["seq"]
	except Exception as e:
		notinfasta.write("{}\n".format(tran))
		continue
	exon_junctions = []
	total_length = len(transeq)
	three_len = 1
	five_len = 1
	strand = master_dict[tran]["strand"]
	if master_dict[tran]["gene_name"] == "":
		master_dict[tran]["gene_name"] = tran
	gene = master_dict[tran]["gene_name"]
	if gene not in coding_genes_dict:
		coding_genes_dict[gene] = False

	if master_dict[tran]["cds"] == []:
		tran_type = "noncoding"
		cds_start = 'NULL'
		cds_stop = 'NULL'
	else:
		#get utr lengths from annotation
		tran_type = "coding"
		coding_genes_dict[gene] = True
		sorted_exons = sorted(master_dict[tran]["exon"])
		sorted_cds = sorted(master_dict[tran]["cds"])
	
		min_cds = sorted_cds[0][0]
		#Some annotation files do not have utr annotation types, so fix that here if thats the case
		if master_dict[tran]["utr"] == []:
			for exon_tup in master_dict[tran]["exon"]:
				if exon_tup not in master_dict[tran]["cds"]:
					# Now check if this overlaps with any of the CDS exons
					overlap = False
					for cds_tup in master_dict[tran]["cds"]:
						if exon_tup[0] == cds_tup[0] and exon_tup[1] != cds_tup[1]:
							master_dict[tran]["utr"].append((cds_tup[1],exon_tup[1]))
							overlap = True
						if exon_tup[0] != cds_tup[0] and exon_tup[1] == cds_tup[1]:
							master_dict[tran]["utr"].append((exon_tup[0],cds_tup[0]))
							overlap = True
					if overlap == False:
						master_dict[tran]["utr"].append(exon_tup)

		for tup in sorted(master_dict[tran]["utr"]):
			if tup[0] < min_cds:
				five_len += (tup[1]-tup[0])+1
			elif tup[0] > min_cds:
				three_len += (tup[1] - tup[0])+1
			else:
				pass
		if strand == "+":
			if len(sorted_exons) > 1:
				sorted_exons[0] = (sorted_exons[0][0]-pseudo_utr_len, sorted_exons[0][1])
				sorted_exons[-1] = (sorted_exons[-1][0], sorted_exons[-1][1]+pseudo_utr_len)
			else:
				sorted_exons[0] = (sorted_exons[0][0]-pseudo_utr_len, sorted_exons[0][1]+pseudo_utr_len)
			master_dict[tran]["exon"] = sorted_exons
			cds_start = (five_len+pseudo_utr_len)
			cds_stop = ((total_length - three_len)-pseudo_utr_len)+4
		elif strand == "-":
			if len(sorted_exons) > 1:
				sorted_exons[0] = ((sorted_exons[0][0]-pseudo_utr_len), sorted_exons[0][1])
				sorted_exons[-1] = (sorted_exons[-1][0], (sorted_exons[-1][1]+pseudo_utr_len))
			else:
				sorted_exons[0] = ((sorted_exons[0][0]-pseudo_utr_len), (sorted_exons[0][1]+pseudo_utr_len))
			master_dict[tran]["exon"] = sorted_exons
			cds_start = (three_len+pseudo_utr_len)
			cds_stop  = ((total_length - (five_len))-pseudo_utr_len)+4
			#if tran == "ENST00000381401":
			#	print ("cds start, cds stop, five_len, three_len",cds_start,cds_stop,five_len,three_len)
			#	#sys.exit()
		else:
			print ("strand is unknown: {}".format(strand))
			sys.exit()
	#get exon junctions, cds is easy just get end of each tuple except last, same for utr except for if same as cds start/stop
	total_intronic = 0
	try:
		if strand == "+":
			tx_start = min(sorted(master_dict[tran]["exon"]))[0]
			prev_end = tx_start
			for tup in sorted(master_dict[tran]["exon"])[:-1]:
				total_intronic += tup[0]-prev_end
				exon_junctions.append(((tup[1])-tx_start)-total_intronic)
				prev_end = tup[1]
		elif strand == "-":
			tx_start = max(sorted(master_dict[tran]["exon"]))[-1]
			prev_end = tx_start
			for tup in (sorted(master_dict[tran]["exon"])[1:])[::-1]:
				total_intronic += (tup[0]+1)-prev_end
				exon_junctions.append(((tup[1])-tx_start)-total_intronic)
				prev_end = tup[1]
	except:
		if strand == "+":
			tx_start = min(sorted(master_dict[tran]["cds"]))[0]
			prev_end = tx_start
			for tup in sorted(master_dict[tran]["cds"])[:-1]:
				total_intronic += tup[0]-prev_end
				exon_junctions.append(((tup[1])-tx_start)-total_intronic)
				prev_end = tup[1]
		elif strand == "-":
			tx_start = max(sorted(master_dict[tran]["cds"]))[-1]
			prev_end = tx_start
			for tup in (sorted(master_dict[tran]["cds"])[1:])[::-1]:
				total_intronic += (tup[0]+1)-prev_end
				exon_junctions.append(((tup[1])-tx_start)-total_intronic)
				prev_end = tup[1]
	#This can happen when a coding transcript doesn't have a properly annotated 3' trailer
	if cds_stop != "NULL":
		if cds_stop > total_length:
			cds_stop = total_length
	if strand == "+" and cds_start != "NULL":
		master_dict[tran]["cds_start"] = cds_start
		master_dict[tran]["cds_stop"] = cds_stop
	elif strand == "-" and cds_start != "NULL":
		master_dict[tran]["cds_start"] = cds_start
		master_dict[tran]["cds_stop"] = cds_stop
	
	master_dict[tran]["strand"] = strand
	master_dict[tran]["tran_type"] = tran_type
	master_dict[tran]["exon_junctions"] = exon_junctions

longest_tran_dict = {}
for tran in master_dict:
	try:
		gene = master_dict[tran]["gene_name"]
	except:
		continue
	if coding_genes_dict[gene] == True:
		if "cds_start" in master_dict[tran]:
			if master_dict[tran]["cds_stop"] != "NULL" and master_dict[tran]["cds_start"] != "NULL":
				cds_length = master_dict[tran]["cds_stop"]- master_dict[tran]["cds_start"]
				if gene not in longest_tran_dict:
					longest_tran_dict[gene] = {"tran":tran,"length":cds_length}
				else:
					if cds_length > longest_tran_dict[gene]["length"]:
						longest_tran_dict[gene] = {"tran":tran,"length":cds_length}
					if cds_length == longest_tran_dict[gene]["length"]:
						if master_dict[tran]["length"] > master_dict[longest_tran_dict[gene]["tran"]]["length"]:
							longest_tran_dict[gene] = {"tran":tran,"length":cds_length}
	else:
		length = master_dict[tran]["length"]
		if gene not in longest_tran_dict:
			longest_tran_dict[gene] = {"tran":tran,"length":length}
		elif length > longest_tran_dict[gene]["length"]:
			longest_tran_dict[gene] = {"tran":tran,"length":length}





for gene in longest_tran_dict:
	longest_tran = longest_tran_dict[gene]["tran"]
	master_dict[longest_tran]["principal"] = 1

gene_sample = []
for key in list(master_dict)[:10]:
	try:
		gene_sample.append(master_dict[key]["gene_name"])
	except:
		pass
print ("Here is a sample of the transcript ids: {}".format(list(master_dict)[:10]))
print ("Here is a sample of the gene names: {}".format(gene_sample))


# Takes a transcript, transcriptomic position and a master_dict (see ribopipe scripts) and returns the genomic position, positions should be passed 1 at a time.
def tran_to_genome(tran, start_pos, end_pos, master_dict):
	pos_list = []
	for i in range(start_pos,end_pos+1):
		pos_list.append(i)
	genomic_pos_list = []
	if tran in master_dict:
		transcript_info = master_dict[tran]
	else:
		return ("Null", [])

	chrom = transcript_info["chrom"]
	strand = transcript_info["strand"]
	exons = sorted(transcript_info["exon"])
	#print ("chrom,strand,exons",chrom,strand,exons)
	for pos in pos_list:
		#print ("pos",pos)
		if strand == "+":
			exon_start = 0
			for tup in exons:
				#print ("tup",tup)
				exon_start = tup[0]
				exonlen = tup[1] - tup[0]
				if pos > exonlen:
					pos = (pos - exonlen)-1
				else:
					break
			#print ("appending exon_start-pos", exon_start, pos, exon_start+pos)
			genomic_pos_list.append((exon_start+pos)-1)
		elif strand == "-":
			exon_start = 0
			for tup in exons[::-1]:
				#print ("tup",tup)
				exon_start = tup[1]
				exonlen = tup[1] - tup[0]
				#print ("exonlen",exonlen)
				if pos > exonlen:
					#print ("pos is greater")
					pos = (pos - exonlen)-1
					#print ("new pos",pos)
				else:
					break
			#print ("appending exon_start-pos", exon_start, pos, exon_start-pos)
			genomic_pos_list.append((exon_start-pos)+1)
	return (chrom, genomic_pos_list)




orf_dict = {"uorf":{},
	"ouorf":{},
	"cds":{},
	"nested":{},
	"odorf":{},
	"dorf":{},
	"minusone":{},
	"readthrough":{},
	"plusone":{},
	"noncoding":{},
	"extension":{},
	"inframe_stop":{}
		}

start_codons = ['ATG','CTG','GTG','TTG','ATC','ATA','ATT','ACG','AAG','AGG']

stop_codons = ["TAG","TAA","TGA"]


# Keep track of the longest transcript for each noncoding gene, append this to transcript list later
longest_noncoding = {}


tran_count = 0
# This section is used to gather all cds regions, convert them to genomic regions and store them in a dictionary to check against later (all transcript contribute to this not just those
# in the transcript list)
genomic_cds_dict = {}
tree_dict = {}
for transcript in master_dict:
	#print (transcript, master_dict[transcript]["tran_type"])
	tran_count += 1
	if "seq" not in master_dict[transcript]:
		continue
	chrom = master_dict[transcript]["chrom"]
	if chrom not in genomic_cds_dict:
		genomic_cds_dict[chrom] = []
	if "cds_start" in master_dict[transcript]:
		cds_start = master_dict[transcript]["cds_start"]
		cds_stop = master_dict[transcript]["cds_stop"]
		if cds_start != "NULL":
			cds_pos = []
			for i in range(cds_start, cds_stop+1):
				cds_pos.append(i)

			for tup in master_dict[transcript]["cds"]:
				if tup[0] != tup[1]:
					if tup not in genomic_cds_dict[chrom]:
						genomic_cds_dict[chrom].append(tup)

print ("genomic cds dict built")
for chrom in genomic_cds_dict:
	tree_dict[chrom] = IntervalTree.from_tuples(genomic_cds_dict[chrom])

#print (list(tree_dict))


connection = sqlite3.connect("{}.sqlite".format(organism))
cursor = connection.cursor()
cursor.execute("CREATE TABLE IF NOT EXISTS transcripts (transcript VARCHAR(50), gene VARCHAR(50), length INT(6), cds_start INT(6), cds_stop INT(6), sequence VARCHAR(50000), strand CHAR(1), stop_list VARCHAR(10000), start_list VARCHAR(10000), exon_junctions VARCHAR(1000), tran_type INT(1), gene_type INT(1), principal INT(1), version INT(2),gc INT(3),five_gc INT(3), cds_gc INT(3), three_gc INT(3), chrom VARCHAR(20));")
cursor.execute("CREATE TABLE IF NOT EXISTS coding_regions (transcript VARCHAR(50), coding_start INT(6), coding_stop INT(6));")
cursor.execute("CREATE TABLE IF NOT EXISTS exons (transcript VARCHAR(50), exon_start INT(6), exon_stop INT(6));")
cursor.execute("CREATE TABLE IF NOT EXISTS uorf (transcript VARCHAR(300), start_codon VARCHAR(10), length INT(6), start INT(6), stop INT(6), cds_coverage FLOAT(20));")
cursor.execute("CREATE TABLE IF NOT EXISTS ouorf (transcript VARCHAR(300), start_codon VARCHAR(10), length INT(6), start INT(6), stop INT(6), cds_coverage FLOAT(20));")
cursor.execute("CREATE TABLE IF NOT EXISTS cds (transcript VARCHAR(300), start_codon VARCHAR(10), length INT(6), start INT(6), stop INT(6), cds_coverage FLOAT(20));")
cursor.execute("CREATE TABLE IF NOT EXISTS nested (transcript VARCHAR(300), start_codon VARCHAR(10), length INT(6), start INT(6), stop INT(6), cds_coverage FLOAT(20));")
cursor.execute("CREATE TABLE IF NOT EXISTS odorf (transcript VARCHAR(300), start_codon VARCHAR(10), length INT(6), start INT(6), stop INT(6), cds_coverage FLOAT(20));")
cursor.execute("CREATE TABLE IF NOT EXISTS dorf (transcript VARCHAR(300), start_codon VARCHAR(10), length INT(6), start INT(6), stop INT(6), cds_coverage FLOAT(20));")
cursor.execute("CREATE TABLE IF NOT EXISTS minusone(transcript VARCHAR(300), start_codon VARCHAR(10), length INT(6), start INT(6), stop INT(6), cds_coverage FLOAT(20));")
cursor.execute("CREATE TABLE IF NOT EXISTS readthrough (transcript VARCHAR(300), start_codon VARCHAR(10), length INT(6), start INT(6), stop INT(6), cds_coverage FLOAT(20));")
cursor.execute("CREATE TABLE IF NOT EXISTS plusone (transcript VARCHAR(300), start_codon VARCHAR(10), length INT(6), start INT(6), stop INT(6), cds_coverage FLOAT(20));")
cursor.execute("CREATE TABLE IF NOT EXISTS noncoding (transcript VARCHAR(300), start_codon VARCHAR(10), length INT(6), start INT(6), stop INT(6), cds_coverage FLOAT(20));")
cursor.execute("CREATE TABLE IF NOT EXISTS extension (transcript VARCHAR(300), start_codon VARCHAR(10), length INT(6), start INT(6), stop INT(6), cds_coverage FLOAT(20));")
cursor.execute("CREATE TABLE IF NOT EXISTS inframe_stop (transcript VARCHAR(300), start_codon VARCHAR(10), length INT(6), start INT(6), stop INT(6), cds_coverage FLOAT(20));")
connection.commit();


print ("Finding ORFs")
transcript_count = 0
total_transcripts = len(list(master_dict))
for transcript in master_dict:
	#print ("transcript",transcript)
	#if transcript != "ENST00000316448":
	#	continue
	transcript_count += 1
	if transcript_count%100 == 0:
		print ("Transcripts complete: {}/{}".format(transcript_count,total_transcripts))
	if "seq" not in master_dict[transcript]:
		print ("transcript {} has no sequence".format(transcript))
		continue
	seq = master_dict[transcript]["seq"]
	cds_start = "NULL"
	cds_stop = "NULL"
	transcript_len = len(seq)
	if "cds_start" in master_dict[transcript]:
		cds_start = master_dict[transcript]["cds_start"]
		cds_stop = master_dict[transcript]["cds_stop"]

	#Find out what regions of this transcript overlap with any other coding regions
	coding_positions = []
	if cds_start != "NULL":
		#If this is a coding transcript don't bother checking the CDS
		for i in range(cds_start,cds_stop):
			coding_positions.append(i)
		#check 5' leader
		chrom, pos_list = tran_to_genome(transcript, 0, cds_start, master_dict)
		for i in range(0,cds_start):
			genomic_pos = pos_list[i]
			overlap = tree_dict[chrom][genomic_pos]
			if len(overlap) != 0:
				coding_positions.append(i)
		#check 3' trailer
		chrom, pos_list = tran_to_genome(transcript, cds_stop, transcript_len, master_dict)
		for i in range(cds_stop,transcript_len+1):
			#print ("i",i)
			genomic_pos = pos_list[i-cds_stop]
			#print ("genomic position",genomic_pos)
			overlap = tree_dict[chrom][genomic_pos]
			if len(overlap) != 0:
				#print ("overlap not empty appending i",overlap)
				coding_positions.append(i)
	else:
		#check entire transcript
		chrom, pos_list = tran_to_genome(transcript, 0, transcript_len, master_dict)
		for i in range(0,transcript_len):
			genomic_pos = pos_list[i]
			overlap = tree_dict[chrom][genomic_pos]
			if len(overlap) != 0:
				coding_positions.append(i)
	coding_positions_tuple = to_ranges(coding_positions)
	coding_dict[transcript] = coding_positions_tuple
	coding_positions = set(coding_positions)
	#if this is a coding transcript find the minusone, readhtrough, plusone co-ordinates
	if cds_start != "NULL":
		#print ("transcript", transcript)
		#if pseudo_utr_len != 0:
		#	cds_stop -= 3 # take 3 from stop so we can match it with orf_stop, do it here rather than above in case cds_stop is null
		recoding_dict = {2:"minusone",0:"readthrough",1:"plusone"}
		for addition in recoding_dict:
			orftype = recoding_dict[addition]
			for i in range(cds_stop+addition,transcript_len,3):
				if seq[i:i+3] in stop_codons:
					#orf_seq = seq[cds_stop:i+3]
					orf_start = cds_stop
					if orftype == "readthrough":
						orf_start -= 2
					if orftype == "plusone":
						orf_start -= 1
					orf_stop = i+3 # +2 so it refers to the end of the stop codon
					start_codon = None
					length = (i+3)-cds_stop
					orf_pos_list = []
					#determine how many nucleotides in this orf overlap with an annotated coding region
					cds_cov_count = 0.0
					for position in range(orf_start,orf_stop):
						orf_pos_list.append(position)
					for pos in range(orf_start, orf_stop+1):
						if pos in coding_positions:
							cds_cov_count += 1
					cds_cov = cds_cov_count/length
					#print ("orftype, start, stop", orftype, orf_start, orf_stop)
					cursor.execute("INSERT INTO {} VALUES('{}','{}',{},{},{},{});".format(orftype, transcript, start_codon, length,orf_start,orf_stop,cds_cov))
					break
		#sys.exit()
	for frame in [0,1,2]:
		for i in range(frame,transcript_len,3):
			if seq[i:i+3] in start_codons:
				for x in range(i, transcript_len,3):
					if seq[x:x+3] in stop_codons:
						#orf_seq = seq[i:x+3]
						orf_start = i+1
						orf_stop = x+3 # +2 so it refers to the end of the stop codon
						start_codon = seq[i:i+3]
						length = (x+3)-i
						orf_pos_list = []
						#determine how many nucleotides in this orf overlap with an annotated coding region
						cds_cov_count = 0.0
						for pos in range(orf_start, orf_stop+1):
							if pos in coding_positions:
								cds_cov_count += 1
						cds_cov = float(cds_cov_count)/float(length)
						# Now determine orf type
						if cds_start == "NULL":
							orftype = "noncoding"
						else:
							#print ("cds start is not null :{}:{}".format(cds_start,cds_stop))
							#print ("orf start, orf stop", orf_start, orf_stop)
							if orf_start == cds_start and orf_stop == cds_stop:
								orftype = "cds"
								#print ("orf type is cds")
							elif orf_start < cds_start and orf_stop == cds_stop:
								orftype = "extension"
								#special case for extensions, we only take from the orf_start to the cds_start, and re-calculate cds coverage
								orf_stop = cds_start
								cds_cov_count = 0.0
								for pos in range(orf_start, orf_stop+1):
									if pos in coding_positions:
										cds_cov_count += 1
								cds_cov = float(cds_cov_count)/float(length)
								#orf_seq = seq[orf_start:cds_start]
								length = cds_start-orf_start
								#print ("orf type is extension")
							elif orf_start < cds_start and orf_stop <= cds_start:
								orftype = ("uorf")
								#print ("orf type is uorf")
							elif orf_start < cds_start and orf_stop > cds_start:
								orftype = "ouorf"
								#print ("orf type is ouorf")
								#sys.exit()
							elif orf_start >= cds_start and orf_start <= cds_stop and orf_stop <= cds_stop:
								if orf_stop == cds_stop:
									break
								#print ("Tran, cds_start, cds_stop, orf_start, orf_stop, tranlen",tran, cds_start, cds_stop, orf_start, orf_stop,transcript_len)
								if (orf_stop < transcript_len and orf_stop%3 == cds_stop%3) or (cds_start != 1 and orf_stop%3 == (cds_start+2)%3):
									#print ("Transcript {} has an inframe stop codon".format(transcript))
									break
								orftype = "nested"
								#print ("orf type is nested")
							elif orf_start >= cds_start and orf_start <= cds_stop and orf_stop > cds_stop:
								orftype = "odorf"
								#print ("orftype is odorf")
							elif orf_start > cds_stop and orf_stop > cds_stop:
								orftype = "dorf"
								#print ("orftype is dorf")
							#if orf_stop > cds_start and orf_stop < cds_stop:
							#	if (orf_stop+1)%3 == cds_start%3:
							#		orftype = "inframe_stop"
							#		print ("inframe stop, transcript, orf_stop", transcript, orf_stop)
							#		sys.exit()
							#		if transcript not in orf_dict:
							#			orf_dict[orftype][transcript] = []
							#	#print ("some weird stop or something")
						cursor.execute("INSERT INTO {} VALUES('{}','{}',{},{},{},{});".format(orftype, transcript, start_codon, length,orf_start,orf_stop,cds_cov))
						break
# Used to keep track of the codons at cds_start and cds_stop positions,
# If there is an issue with the cds co-ordinates the starts and stops counts will
# be much lower than the other count, start with 1 to prevent division by 0
nuc_dict = {"stops":{"stops":1,"other":0}, "starts":{"starts":1,"other":0}}

def calcgc(seq):
	if seq == "":
		return "NULL"
	g_count = 0
	c_count = 0
	a_count = 0
	t_count = 0
	for char in seq:
		if char == "A":
			a_count += 1
		if char == "T":
			t_count += 1
		if char == "G":
			g_count += 1
		if char == "C":
			c_count += 1
		gc = ((g_count+c_count)/float(len(seq)))*100
	return round(gc,2)




for transcript in master_dict:
	#print ("transcripts", transcript)
	length = master_dict[transcript]["length"]
	cds_start = master_dict[transcript]["cds_start"]
	cds_stop = master_dict[transcript]["cds_stop"]
	seq = master_dict[transcript]["seq"]
	strand = master_dict[transcript]["strand"]
	chrom = master_dict[transcript]["chrom"]
	gene = master_dict[transcript]["gene_name"]
	gc = calcgc(seq)
	five_gc = "NULL"
	cds_gc = "NULL"
	three_gc = "NULL"
	if cds_start != "NULL":
		five_gc = calcgc(seq[:cds_start])
		cds_gc = calcgc(seq[cds_start:cds_stop])
		three_gc = calcgc(seq[cds_stop:])
		# check that the nucleotide cds_start points to is the first of the start codon
		# take one becase cds_start is 1 based but python indexing is 0 based
		start_nuc = seq[cds_start-1:cds_start+2]
		#print ("start nuc",start_nuc)
		if start_nuc == "ATG":
			nuc_dict["starts"]["starts"] += 1
		else:
			nuc_dict["starts"]["other"] += 1
		stop_nuc = seq[cds_stop-3:cds_stop]
		#print ("stop_nuc",stop_nuc)
		if stop_nuc in ["TAG","TAA","TGA"]:
			nuc_dict["stops"]["stops"] += 1
		else:
			nuc_dict["stops"]["other"] += 1
	tran_type = master_dict[transcript]["tran_type"]
	if coding_genes_dict[gene] == True:
		gene_type = 1
	else:
		gene_type = 0
	#print ("tran type before",tran_type)
	if tran_type == "coding":
		tran_type = 1
	else:
		tran_type = 0
	#print ("tran type after",tran_type)
	start_list = str(master_dict[transcript]["start_list"]).replace(" ","").strip("[]")
	stop_list = str(master_dict[transcript]["stop_list"]).replace(" ","").strip("[]")
	exon_junctions = str(master_dict[transcript]["exon_junctions"]).replace(" ","").strip("[]")
	principal = master_dict[transcript]["principal"]
	version = master_dict[transcript]["version"]
	#print (master_dict[transcript])
	#print (tran_type)
	#print (gene_type)
	#print (principal)
	#print (version)
	#print ("INSERT INTO transcripts VALUES('{}','{}',{},{},{},'{}','{}','{}','{}','{}',{},{},{},{});".format(transcript, gene, length, cds_start, cds_stop, seq, strand,stop_list, start_list, exon_junctions, tran_type,gene_type,principal,version))
	cursor.execute("INSERT INTO transcripts VALUES('{}','{}',{},{},{},'{}','{}','{}','{}','{}',{},{},{},{},{},{},{},{},'{}');".format(transcript, gene, length, cds_start, cds_stop, seq, strand,stop_list, start_list, exon_junctions, tran_type,gene_type,principal,version,gc,five_gc,cds_gc,three_gc,chrom))

	for tup in master_dict[transcript]["exon"]:
		cursor.execute("INSERT INTO exons VALUES('{}',{},{});".format(transcript,tup[0],tup[1]))
	if transcript in coding_dict:
		for tup in coding_dict[transcript]:
			cursor.execute("INSERT INTO coding_regions VALUES('{}',{},{});".format(transcript,tup[0],tup[1]))

connection.commit()
connection.close()

print("delim", delimiters)
if (nuc_dict["starts"]["other"]/nuc_dict["starts"]["starts"]) > 0.05:
	print ("Warning: {} transcripts do not have a an AUG at the CDS start position".format(nuc_dict["starts"]["other"]))
if (nuc_dict["stops"]["other"]/nuc_dict["stops"]["stops"]) > 0.05:
	print ("Warning: {} transcripts do not have a an stop codon at the CDS stop position".format(nuc_dict["stops"]["other"]))
if len(notinannotation) >0:
	print ("Warning: {} transcripts were in the fasta file, but not the annotation file, these will be discarded".format(len(notinannotation)))
