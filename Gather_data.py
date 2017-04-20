from Bio import SeqUtils
from Bio import Entrez
from urllib2 import HTTPError
import time
import sqlite3    

Entrez.email ="ambat.jacob@gmail.com"#login to Entrez#

sqlite_file = '/my_db.sqlite'#generate db file in the current directory#

#database global variables#

sqlite_file = 'fastadb.sqlite'   
table_name2 = 'Gather_data' 
new_field = 'No'
field_type = 'INTEGER'
id_column = 'gi_accession'
column_type2 = 'TEXT' 
description_column = 'description'
column_type3 = 'TEXT'
seq_column = 'sequence'
column_type4 = 'TEXT'
comp_seq_column = 'Complementary_sequence'
column_type5 = 'TEXT'
PAM_column1 = 'PAMsites_exons'
column_type6 = 'VARCHAR'
PAM_column2 = 'PAMsites_complementary'
column_type7 = 'VARCHAR'

#creating gene information in sqlite db#
#function used to create db available in the current directory#
def createdb():
	gis = [100753385, 100689306, 100751648]	
	accession = []
	description = []
	sequence = []
	
	request = Entrez.epost("nucleotide",id=",".join(map(str,gis)))
	result = Entrez.read(request)
	webEnv = result["WebEnv"]
	queryKey = result["QueryKey"]
	handle = Entrez.efetch(db="nucleotide",retmode="xml", webenv=webEnv, query_key=queryKey)
	for r in Entrez.parse(handle):
		# Grab the GI# 
		try:
			gi=int([x for x in r['GBSeq_other-seqids'] if "gi" in x][0].split("|")[1])
		except ValueError:
			gi=None
		fastaseq = ">GI ",gi," "+r["GBSeq_primary-accession"]+" "+r["GBSeq_definition"]+"\n"+r["GBSeq_sequence"][0:20]
		accession.append(''.join(fastaseq[0].strip() + str(fastaseq[1])))
		description.append(' '.join(fastaseq[2].split()[0:3]))
		sequence.append(fastaseq[2].split()[-1].upper())
	
	alt_map = {'ins':'0'}
	complement = {'A':'T','G':'C','T':'A','C':'G'}
	
	# getting the complementary sequence#
	def reverse_complement(seq):    
	    for k,v in alt_map.iteritems():
	        seq = seq.replace(k,v)
	    bases = list(seq) 
	    bases = reversed([complement.get(base,base) for base in bases])
	    bases = ''.join(bases)
	    for k,v in alt_map.iteritems():
	        bases = bases.replace(v,k)
	    return bases
	
	complementary_sequence = [reverse_complement(seq) for seq in sequence]
	
	
	#print sequence,complementary_sequence#
	
	#fetching the positions of 'GG' from the sequence
	exon = []
	comp_exon = []
	pattern = 'GG'
	for exons in sequence:
		
		exon_search = str(SeqUtils.nt_search(exons, pattern))
		exon.append(exon_search)
		
	for comp in complementary_sequence:
		
		comp_exon_search = str(SeqUtils.nt_search(comp, pattern))
		comp_exon.append(comp_exon_search)
	
	#print exon
	#print comp_exon
	
	conn = sqlite3.connect(sqlite_file)
	c = conn.cursor()
	
	c.execute('CREATE TABLE {tn} ({nf} {ft} PRIMARY KEY)'\
			.format (tn=table_name2, nf=new_field, ft=field_type))
	c.execute("ALTER TABLE {tn} ADD COLUMN '{cn}' {ct}"\
	        .format(tn=table_name2, cn=id_column, ct=column_type2))
	c.execute("ALTER TABLE {tn} ADD COLUMN '{cn}' {ct}"\
	        .format(tn=table_name2, cn=description_column, ct=column_type3))
	c.execute("ALTER TABLE {tn} ADD COLUMN '{cn}' {ct}"\
	        .format(tn=table_name2, cn=seq_column, ct=column_type4))
	c.execute("ALTER TABLE {tn} ADD COLUMN '{cn}' {ct}"\
	        .format(tn=table_name2, cn=comp_seq_column, ct=column_type5))       
	c.execute("ALTER TABLE {tn} ADD COLUMN '{cn}' {ct}"\
	        .format(tn=table_name2, cn=PAM_column1, ct=column_type6))
	c.execute("ALTER TABLE {tn} ADD COLUMN '{cn}' {ct}"\
	        .format(tn=table_name2, cn=PAM_column2, ct=column_type7))       
	        
	c.execute('''INSERT INTO Gather_data(No, gi_accession, description, sequence,Complementary_sequence,PAMsites_exons,PAMsites_complementary) VALUES(?,?,?,?,?,?,?)''', (1, accession[0], description[0], sequence[0],complementary_sequence[0],exon[0],comp_exon[0]))
	c.execute('''INSERT INTO Gather_data(No, gi_accession, description, sequence,Complementary_sequence,PAMsites_exons,PAMsites_complementary) VALUES(?,?,?,?,?,?,?)''', (2, accession[1], description[1], sequence[1],complementary_sequence[0],exon[1],comp_exon[1]))
	c.execute('''INSERT INTO Gather_data(No, gi_accession, description, sequence,Complementary_sequence,PAMsites_exons,PAMsites_complementary) VALUES(?,?,?,?,?,?,?)''', (3, accession[2], description[1], sequence[2],complementary_sequence[0],exon[2],comp_exon[2]))
	conn.commit()
	conn.close()

########################################################################
#this function retrives values of PAM sites forthe given Gene name,
#@param gene_name : gene name sent from user
#@return : values of PAM positions
########################################################################
def get_names(gene_name):
	conn = sqlite3.connect(sqlite_file)
	c = conn.cursor()
	c.execute('SELECT PAMsites_exons,PAMsites_complementary FROM Gather_data WHERE gi_accession = ?'\
			.format(),(gene_name,))

	row = c.fetchone()
	
	val3 = ''
	if row != None:
		
		val1 = row[0]	
		val2 = row[1]
		val3 = 'PAMsites_exons = ' + val1+' PAMsites_complementary = '+val2 	
		print val1,val2
	else:
		val3 = 'no data found'
	conn.close()
	return val3



