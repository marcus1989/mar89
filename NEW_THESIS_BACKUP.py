########################################################################
						     MODULE1
########################################################################
#######################--blast_parser_tool--############################
import csv
import time
import re
import os
import sys
from collections import Counter
import operator
from fractions import *
import glob
import ntpath
from collections import defaultdict



path = open('config.txt').read().splitlines()[0].split('=')[-1]

rootDir = '.'		
blast_files = []
curdir = os.getcwd()
curdir_up = '/'.join(curdir.split('/')[:-1])


for dirName, subdirList, fileList in os.walk(rootDir, topdown = False):
	
	for fname in fileList:
		
		if fname.startswith("S.A"):
			
			fname = os.path.join(dirName, fname)
			blast_files.append(fname)

print 'Module1'
print '		step 1.1 : Parsing the input Blastp files'	


for blastfiles in blast_files[:]:

	if 'Prot' not in blastfiles:
		
		qids=[]
		query_lengths = []
		counter = 0
		seqsas = []
		file1 = open(blastfiles,'r').read()
		queries = file1.split('Query=')
		datas = queries[1:]	
		
		for item in datas[:]:
			
				lines = item.split('\n')
				qid = item.split()[0]
				qids.append(qid)
				
				for line in lines[:]:
					
					if line.startswith('Length='):
						
						query_lengths.append(int(line.split('=')[-1]))
						break
					
		for i,data in enumerate(datas[:]):
			
			lines = data.split('\n')
			record = False
			
			for line in lines[:]:
				
				if line.startswith(">") :
					
					tmp = line.split(">")
					tmp_name = tmp[1]
					tmp_name1 = tmp_name.split("[")	
					tmp_hit = ''.join(tmp_name[0:-1])
					
					if 'Staphylococcus' in line:
						
						record  = True
						
					else:
						
						record = False
						
				if line.startswith(" Score") and record:
					
					tmp = line.strip().split()
					tmp_score_s = tmp[2]
					tmp_score = float(tmp_score_s)
					tmp_evalue = float(tmp[7].replace(",",""))
					seqsas.append([qids[i],tmp_hit,tmp_score,tmp_evalue])
		
				if line.startswith(" Identities")and counter <len(seqsas) and record:
					
					tmp = line.strip().split() 
					tmp_id = tmp[3]
					tmp_ids = tmp_id.replace('(','').replace(')','').replace('%','').replace(',','')
					ids = int(tmp_ids)
					tmp_po = tmp[7]
					tmp_pos = tmp_po.replace('(','').replace(')','').replace('%','').replace(',','')
					pos = int(tmp_pos)
					tmp_gap = tmp[11]
					tmp_gaps = tmp_gap.replace('(','').replace(')','').replace('%','').replace(',','')
					gaps_percent = int(tmp_gaps)
					gap_number = int(tmp[10].split('/')[0])
					alignment_length = int(tmp[10].split('/')[-1])
					coverage_percent = round(float((alignment_length - gap_number))/query_lengths[i] * 100, 2)
					seqsas[counter].append(ids)
					seqsas[counter].append(pos)
					seqsas[counter].append(gaps_percent)
					seqsas[counter].append(gap_number)
					seqsas[counter].append(alignment_length)
					seqsas[counter].append(coverage_percent)
					counter+=1
					
		path1 = '%s/RESULT/MODULE1/P1' % curdir_up

		if not os.path.exists(path1):
			
			os.makedirs(path1)
		file_name = ntpath.basename('blast_out1%s' % blastfiles) + '.txt'
		
		with open(os.path.join(path1,file_name),'w') as out1:
			
			for item in seqsas[:]:
				
				item = '\t'.join([str(x) for x in item])
				out1.write('%s\n' %item)
				
			out1.close()
					
	else:
		
		strsas = []
		qids=[]
		query_lengths = []
		counter = 0
		file2 = open(blastfiles,'r').read()
		queries = file2.split('Query=')
		datas = queries[1:]	
		
		for item in datas[:]:
			
				lines = item.split('\n')
				qid = item.split()[0]
				qids.append(qid)
				
				for line in lines[:]:
					
					if line.startswith('Length='):
						
						query_lengths.append(int(line.split('=')[-1]))
						break
						
		for i,data in enumerate(datas[:]):
			
			lines = data.split('\n')
			record = False
			
			for line in lines[:]:
				
				if line.startswith(">") :
					
					tmp = line.split(">")
					tmp_name = tmp[1]
					tmp_hit = tmp_name.split("|")[0]
			
					
				if line.startswith(" Score") :
					
					tmp = line.strip().split()
					tmp_score_s = tmp[2]
					tmp_score = float(tmp_score_s)
					tmp_evalue = float(tmp[7].replace(",",""))
					strsas.append([qids[i],tmp_hit,tmp_score,tmp_evalue])
		
				if line.startswith(" Identities") and counter < len(strsas):
					
					tmp = line.strip().split()
					tmp_id = tmp[3]
					tmp_ids = tmp_id.replace('(','').replace(')','').replace('%','').replace(',','')
					ids = int(tmp_ids)
					tmp_po = tmp[7]
					tmp_pos = tmp_po.replace('(','').replace(')','').replace('%','').replace(',','')
					pos = int(tmp_pos)
					tmp_gap = tmp[11]
					tmp_gaps = tmp_gap.replace('(','').replace(')','').replace('%','').replace(',','')
					gaps_percent = int(tmp_gaps)
					gap_number_1 = Fraction(tmp[10])
					gap_number = int(tmp[10].split('/')[0])
					alignment_length = int(tmp[10].split('/')[-1])
					coverage_percent = round(float((alignment_length - gap_number))/query_lengths[i] * 100, 2)
					strsas[counter].append(ids)
					strsas[counter].append(pos)
					strsas[counter].append(gaps_percent)
					strsas[counter].append(gap_number)
					strsas[counter].append(alignment_length)
					strsas[counter].append(coverage_percent)
					counter +=1
					
		path1 = '%s/RESULT/MODULE1/P1' %curdir_up
		
		if not os.path.exists(path1):
			
			os.makedirs(path1)
		prot_file_name = ntpath.basename('prot_blast_out1%s' % blastfiles) + '.txt'
		
		with open(os.path.join(path1,prot_file_name),'w') as out2:
			
			for item in strsas[:]:
				
				item = '\t'.join([str(x) for x in item])
				out2.write('%s\n' %item)
				
			out2.close()
		
def parser2():
	
		os.chdir('%s/RESULT/MODULE1/P1' %curdir_up)
		
		for file1 in glob.glob('*.txt'):
			file_s = open(file1).readlines()
			prepsas = []
			
			for item in file_s[:]:
				
				item = item.strip().split('\t')
				hit = item[1]
				e = float(item[3])
				ids = int(item[4])
				cov = float(item[9])
				if e <=1e-10 and ids >= 35 and cov >= 75:
					
					prepsas.append(item)
						
				if len(item) < 10:
					
					print 'not match'
			
			prot_file_name_s = str(file1) 
			
			path2 = '%s/RESULT/MODULE1/P2' %curdir_up
			
			if not os.path.exists(path2):
				
				os.makedirs(path2)
		
			with open(os.path.join(path2,prot_file_name_s),'w') as prepsas1:
				
				for hits in prepsas[:]:
					
					hits = '\t'.join([str(x) for x in hits])
					prepsas1.write('%s\n' %hits)
					
				prepsas1.close()
		
				
def parser3():
	
	os.chdir('%s/RESULT/MODULE1/P2' %curdir_up)
	
	for file2 in glob.glob('*.txt'):
		
		file3 =open(file2).readlines()
		d = {}
		
		for filters in file3[:]:
			
			key, value = filters.strip("\n").split("\t")[0],filters.strip("\n").split("\t")[1:]
			key = key.strip('\t')
			value = [str(x)[0:]for x in value]
			
			if key not in d:
				
				d[key] = [value]
				
			elif key in d and len(d[key]) <= 250:
				
				d[key].append(value)
				
		prot_file_name_s = str(file2) 
		
		path2 = '%s/RESULT/MODULE1/P3' %curdir_up
		
		if not os.path.exists(path2):
			
			os.makedirs(path2)	
				
		with open(os.path.join(path2,prot_file_name_s),'w') as fp:
			
			for item in d.keys()[:]:
				
				line = item
				hits = d[item]

				for hit in hits:
					
					hit2 = ','.join(hit)
					line += '\t%s' % hit2
					
				fp.write("%s\n" % line)

	
parser2()
parser3()

#########################--ftp_downloads.py--###########################	
from Bio import Entrez
from urllib2 import HTTPError
import time
import os
import glob
import ntpath

print "		step 1.2: FTP Download of the NCBI sequences"

curdir = os.getcwd()
curdir_up = '/'.join(curdir.split('/')[:-1])

if os.path.isdir("%s/RESULT/MODULE1/FTP" %curdir_up):
	
	print "		FTP already exists"
	 
else:	
	
	Entrez.email ="ambat.jacob@gmail.com"
	
	path = open('config.txt').read().splitlines()[1].split('=')[-1].strip()

	for root, dirs, files in os.walk(path):

		filepaths = []
		
		for file in files:
			
			filepath = os.path.join(root, file)
			filepaths.append(filepath)
		
		for filepath in filepaths[:]:
			if 'Prot' not in filepath:
				
				p3files = open(filepath,'r').readlines()
				d = {}
				
				for hits in p3files[:]:
					
					key,value = hits.strip("\n").split("\t")[0],hits.strip("\n").split("\t")[1:]
					key = key.strip('\t')
					value = [x.split('|')[0:] for x in value]
					value2 =[x[1] for x in value]
					value3 = [ x for x in value2 if len(x) > 4]

					for items_ids in value3:
						if key not in d:
							d[key] = [items_ids]
						else:
							d[key].append(items_ids)	

			new = ['_ids.txt']
			if 'Prot' not in filepath: 
				
				ftps = ntpath.basename("%s" %filepath)
				newftps = ftps.split('.txt')[:-1]
				filtered_ids_names = ''.join(newftps + new)
				
				with open(os.path.join(path,filtered_ids_names),'w') as output_ids:	
					
					for item in d.keys()[:]:
						
						line = item
						hits = d[item]
						line += '\t%s' %hits
						output_ids.write("%s\n" % line)	
						
					output_ids.close()

				for items in d.keys()[:]:
					
					ftp_path = '%s/RESULT/MODULE1/FTP' %curdir_up
					os.system('mkdir %s' % ftp_path)
					item_path = ftp_path + '/%s' %items
					
					if item_path not in os.listdir(ftp_path):
						
						os.system('mkdir %s' % item_path)
						
					os.chdir(item_path)
					print os.getcwd()
					counter = 0
					values = d[items]
	
					for value in values:
						
						out = open('%s' % value, 'w')
						fetch_success = False
						
						while(fetch_success == False):
							
							try:
								handle = Entrez.efetch(db="protein", id=value, retmode="xml")
								records = Entrez.read(handle)
								fastaseq = value.rstrip()+" "+records[0]["GBSeq_primary-accession"]+" "+records[0]["GBSeq_definition"]+"\n"+records[0]["GBSeq_sequence"]
								out.write(fastaseq)
								out.close()
								fetch_success = True
								
							except:
								
								pass
							time.sleep(1) # to make sure not many requests go per second to ncbi	
							
					counter += 1

#############################--Merge.py--###############################
import os
import fnmatch
from collections import Counter
import glob

print "		step 1.3:  Merge the NCBI sequences to one sequence.fasta file"

path = open('config.txt').read().splitlines()[2].split('=')[-1].strip()

dirpaths = []

for root, dirs, files in os.walk(path):  
	
	for dir in dirs:
		
		path = os.path.join(root, dir)
		dirpaths.append(path)
		

for dirpath in dirpaths[:]:
	
	folders = dirpath.split("/")[-1]
	merged_content = ''
	files =  glob.glob('%s/*' % dirpath)
	newcontent = {}
	
	for file in files[:]:
		
		content = open(file).read()	
		key = ">" + content.split()[0]
		value = content.split()[-1]
		
		if key not in newcontent:
			newcontent[key] = value
		elif key in content:
			newcontent[key].append(value)

	new_ncbi_file = dirpath + "/" + "sequence.fasta"

	with open(new_ncbi_file,"w") as fasta_file:
		
		for merge in newcontent.keys()[:]:
			line = merge
			sequences = newcontent[merge]
			line += '\n'
			for merged in sequences:
				
				merged2 = ' '.join(merged)
				line += '%s' % merged2 
			
			fasta_file.write("%s\n" % line)
	fasta_file.close()

######################--MSA_parser.py--#################################
import os
import re
import fnmatch
from collections import OrderedDict

print "		step 1.4:  Preparing the ClustalO input files"

curdir = os.getcwd()
curdir_up = '/'.join(curdir.split('/')[:-1])

path = open('config.txt').read().splitlines()[2].split('=')[-1].strip()

colfiles = path
lists = os.listdir(colfiles)


seq_input = ["%s/S.AureusCOL" %curdir_up,"%s/S.AureusN315" %curdir_up]

filepaths = []
for root, dirs, files in os.walk(path): 
	  
	for file in files:
		
		if '.fasta' in file:
			
			paths = os.path.join(root, file)
			filepaths.append(paths)



for strains in seq_input[:]:
	
	def seqfilter():
			
		files = open(strains).read().split(">")
		q_ids = []
		qrys = []
		for i,ids in enumerate(files[:]):
			
			query_ids = ids.strip("\n")
			q_ids.append(query_ids)
			query = ids.strip(" ")[0:9]
			querys = query.split(" ")[0]
			
			for filepath in filepaths[:]:
				
				if querys in filepath:
					
					try:
						
						content = open(filepath).read()
					
					
						if ids not in content:
							
							new_content =  ">" + ids  + "\n" + content
							newfile = open(path+'/%s.fasta' % querys, 'w')
							newfile.write(new_content)	
					except:
						
						pass
						
				else:
					
					pass
			
	#########################################################################
	
	def blancofilter():		
				
		blanco  = open('%s/blancoDBproteins' %curdir_up).read().split('>')
		blancos = []
		
		for i,pdbs in enumerate(blanco[:]):
			
			pdb_seqs = pdbs.strip("\n")
			blancos.append(pdb_seqs)
			pdbids = pdbs.strip(" ")[0:6]
			
			path3 = '%s/pdbseq' %curdir_up
			
			if not os.path.exists(path3):
				
				os.makedirs(path3)	
				
			with open(os.path.join(path3,"%s.fasta" % pdbids),'w') as outfile:
				outfile.write(pdbs)
	
	#########################################################################
	
	pdbfiles = '%s/pdbseq' %curdir_up
	
	
	def msafileprep():
		
		dirlistings = os.listdir(pdbfiles)
		editfiles =[]
		structure_names_ext = []
		
		for item in dirlistings:
			
			if ".fasta" in item:
				
				editfiles.append(item)

		new_file = ["%s/RESULT/MODULE1/P3/S.AureusCOL_Prot_Blast.txt" %curdir_up,"%s/RESULT/MODULE1/P3/S.AureusN315_Prot_Blast.txt" %curdir_up]
			
		input_fasta_files = []
		
		for input_file in lists:
			
			if input_file.endswith("fasta"):
				
				input_fasta_files.append(input_file)

		fasta_names = []
		
		for fasta_name_ext in input_fasta_files:

			fasta_name = fasta_name_ext.split(".")[0]
			fasta_names.append(fasta_name)

		prot_names = []
		new_fil_pdbs = []
		str_dict = {}
		
		
		for p3out in new_file[:]:
			
			new_files = open(p3out).readlines()
			
			for i,prots in enumerate(new_files[:]):
				
				values = prots.split('\t')
				key = values[0]
				prot_values = [item.split(',')[0].strip() + '.fasta' for item in values[1:]]
				
				if key not in str_dict:
					
					str_dict[key] = prot_values
					
				elif key in str_dict:
					
					str_dict[key].append(prot_values)
		
		set1 = set(editfiles)
		newlist = [key for key, value in str_dict.iteritems() if str_dict.values() and set1]	
		structure_name_ext = [value for key, value in str_dict.iteritems() if str_dict.values() and set1]
		
		
		for strain_name in newlist[:]:	
			
			if strain_name in fasta_names[:]:
				
				fasta_file = open(colfiles + "/" + strain_name+'.fasta').read()
				pdbs = str_dict[strain_name]
				content = fasta_file
				
				for newseqs in pdbs:
					
					sequence = open(pdbfiles + "/" + newseqs).read()
					sequence_str = '>' + sequence
					content += sequence_str
					fasta_files = open(colfiles + "/" + strain_name + '.fasta','w')
					fasta_files.write(content)
					fasta_files.close()
			
			
		for fastas in fasta_names[:]:
			
			filename = colfiles + "/" + fastas + '.fasta'
			
			if fastas not in newlist[:]:
				
				try:
					os.remove(filename)	
					
				except:
					
					pass
				
	seqfilter()	
	blancofilter()
	msafileprep()

######################--clustalo.py--###################################
import subprocess
import sys
import os
import glob
from Bio.Align.Applications import ClustalOmegaCommandline

print "		step 1.5:  Multiple Sequence Alignment using ClustalO"

curdir = os.getcwd()
curdir_up = '/'.join(curdir.split('/')[:-1])

os.chdir("%s/RESULT/MODULE1/FTP" %curdir_up)

files = glob.glob("*.fasta")

def clustalo():
	
	for fasta_files in files[:]:
		
		fasta_file_s = fasta_files.split('.')[0]
		print fasta_file_s
		cline = ClustalOmegaCommandline('clustalo',infile = fasta_files,outfile = fasta_file_s + '.aln' ,verbose= False, auto=True)
		cline()
		
clustalo()

print "End of Module 1"
########################################################################
							MODULE2
########################################################################
###########################--htmlparse.py--#############################
from bs4 import BeautifulSoup
import urllib2
import csv
import os
import sys
import glob
import shutil

print "Module 2"

print "		step 2.1 : Creating the transporter folders and alignment files"

curdir = os.getcwd()
curdir_up = '/'.join(curdir.split('/')[:-1])

os.chdir("%s/RESULT/MODULE1/FTP" %curdir_up)
aln_files = glob.glob("*.aln")

strains = {}

for i in range(3):
	
	url=('http://www.membranetransport.org/all_type_btab.php?oOID=saur%s' % i)
	header = {'User-Agent': 'Morzilla/5.0'}
	req=urllib2.Request(url, headers = header)
	page = urllib2.urlopen(req)
	soup = BeautifulSoup(page)
	
	
	ORF = ""
	Family_ID = ""
	Family_TC = ""
	Family_Name = ""
	Transporter_Type = ""
	Substrate = ""
	Note = ""
	TC = ""
	
	table = soup.find('body', {"bgcolor" : "#FFFFFF"})
	

	for row in table.findAll("tr"):
		
		cells = row.findAll("td")
		
		if len(cells) == 8:
			
			ORF = cells[0].find(text = True)
			ORF = str(ORF)
			
			if 'ORF' not in ORF:
				
				Family_ID = cells[1].find(text = True)
				Family_ID = str(Family_ID)
				Family_TC = cells[2].find(text = True)
				Family_Name = cells[3].find(text = True)
				Transporter_Type = cells[4].find(text = True)
				Substrate = cells[5].find(text = True)
				Note = cells[6].find(text = True)
				TC = cells[7].find(text = True)
				
				if Family_ID not in strains:
					
					strains[Family_ID] = [ORF]
					
				elif Family_ID in strains:
					
					strains[Family_ID].append(ORF)
				 

file_names = []

for alns in aln_files[:]:
	names  =  alns.split('.aln')[0]
	file_names.append(names)
	

f = open('%s/RESULT/transport_families.txt' %curdir_up, 'wb')

if os.path.isdir("%s/RESULT/MODULE2" %curdir_up):
	
	print "	folder already exists"
	 
else:	
	
	for item in strains.keys()[:]:
		
		line = item
		folders = os.makedirs(curdir_up + "/RESULT/MODULE2/%s" %line)
		hits = strains[item]
		
		for news in file_names:
			
			if news in hits:
				
				news_folders = os.makedirs(curdir_up + "/RESULT/MODULE2/%s/%s" %(line,news))
				src = curdir_up + "/RESULT/MODULE1/FTP/%s.aln" %news
				dest = curdir_up + "/RESULT/MODULE2/%s/%s" %(line,news)
				shutil.copy2(src, dest)
				
		for hit in hits:
			
			line += '\n\t%s\n' % hit
			
		f.write("%s\n" % line)	 


root_dir = "%s/RESULT/MODULE2" %curdir_up
newaln = glob.glob('*.aln')

empty_count = 0

for curdir, subdirs, files in os.walk(root_dir):
	
	if len(subdirs) == 0 and len(files) == 0: 
		
		empty_count += 1 
		os.rmdir(curdir)
#########################--rate4site.py--###############################
import os
import glob

print "		step 2.2 : Rate4Site program"

curdir = os.getcwd()
curdir_up = '/'.join(curdir.split('/')[:-1])
os.chdir(curdir_up + "/RESULT/MODULE2/")

newdir = os.getcwd()

def rate4site():
	
	for curdir, subdirs, files in os.walk(newdir):
	
		for file in files:
			
			rate_file = file.split('.')[0] + '_rate'
			file_path = os.path.join(curdir, file)
			outpath = os.path.join(curdir, rate_file)
			os.system('rate4site -s %s -o %s' %(file_path,outpath))
			
rate4site()	
#######################-rate_check.py--#################################

import os
import re


print "		step 2.3 : Check if the rate files are empty"
   

curdir = os.getcwd()
curdir_up = '/'.join(curdir.split('/')[:-1])
rate_path = "%s/RESULT/MODULE2" %curdir_up

for curdir, subdirs, files in os.walk(rate_path): 
	
	for file in files[:]:
		
		newpath = os.path.join(curdir, file)
		if os.stat(newpath).st_size == 0:
			os.remove(newpath)
			
#######################--rate4siteparser.py--###########################
import os
import re
import scipy.stats as stats
import matplotlib.pyplot as plt
import pylab as plt2
from pylab import *
import numpy as np
from decimal import *
from scipy.stats import norm
from sklearn import preprocessing
from math import floor, ceil 

print "		step 2.3 : Rate4Site Parser and Plotting the Bins"

def rate4site():
	curdir = os.getcwd()
	curdir_up = '/'.join(curdir.split('/')[:-1])
	rate_path = "%s/RESULT/MODULE2" %curdir_up
	
	for curdir, subdirs, files in os.walk(rate_path): 
		
		for file in files[:]:
			
			if '_rate' in file: 
				
				query = file.split('_')[0]
				outpath = os.path.join(curdir, file)
				rate_file = open(outpath)
				confidence = []
				data_reqs = []
				fil_data = []
				scores = [] 
				residues = []
				range1 = []
				range2=[]
				y_count = []
				
				for i,interval in enumerate(rate_file):
					
					if i == 0:
						
						True
						
					else:
						
						if not interval.startswith('#') and interval != '\n':
							
							interval = interval.split(',')
							parts1 = interval[0].split()
							seq = parts1[0]
							res = parts1[1]
							residues.append(res)
							y_count.append(seq)
							score = parts1[2]				
							scores.append(float(score))					
							range1.append(parts1[3][1:])
							range1 = [item for item in range1 if item.strip()]					
							parts2 = interval[-1].split()
							range2.append(parts2[0][:-1])
							range2 = [item for item in range2 if item.strip()]
							range_diff = range2 + range1
							confidence.append(interval)
							data_reqs.append([score,range_diff])
					
					
				fig, ax = plt.subplots(1, 1)
				mean, var, skew, kurt = norm.stats(moments='mvsk')
				rang = [min(scores), max(scores)]
				Long = len(scores)
				Maxim = max(scores) #MaxValue
				Minim = min(scores) #MinValue
				av = np.mean(scores) #Average
				StDev = np.std(scores) #Standard Dev.
				
				
				x = np.linspace(Minim, Maxim, Long)
				ax.plot(x, norm.pdf(x, av, StDev),'r-', lw=3, alpha=0.9, label='RATE SCORES')
				
				weights = np.ones_like(scores)/len(scores)
				normalized = [(s-min(scores))/(max(scores)-min(scores)) for s in scores]
				newpath = os.path.join(curdir)
			
			
				ax.hist(normalized, weights = weights, normed=True, histtype='stepfilled', alpha=0.2,label='NORMALIZED RATE SCORES')
				plt.title('%s' %query +'_Normalized_Rate4Site_Scores')
				plt.xlabel('Rate4Site_Scores', fontsize=14)
				plt.ylabel('Sequence_Count', fontsize=14)
				plt.legend(loc='upper right')
				fig.savefig(newpath +'_%s_normalized.png'%query)
				plt.close("all")

				y_count = map(int, y_count)
				
				color_path = os.path.join(curdir)
				
				bins = np.arange(floor(min(normalized[:])), ceil(max(normalized[:])), 0.10)
				
				colors = ('blue', 'red', 'green','cyan','purple','pink','violet','lime','aqua')
				
				# get the max count for a particular bin for all classes combined
				max_bin = max(np.histogram(normalized[:], bins=bins)[0])
				plt.figure()
				n, bins, patches = plt.hist(normalized[:], bins , alpha=0.3)
				
				for c, p in zip(colors, patches):
					
					plt.setp(p, 'facecolor', c)  
					
				plt.ylim([0, max_bin*1.3])
				plt.title('%s'%query + '_Normalized_Scores_In_9_Bins')
				plt.xlabel('Color_Bins', fontsize=14)
				plt.ylabel('Sequence_Count', fontsize=14)
				plt.legend(loc='upper right')
				plt.savefig(newpath + '%s_colorbins.png' %query)
				plt.close("all")
				
				pairs = [(x,y,z) for x,y,z in zip(normalized, residues, y_count)]
				group1 = []
				group2 = []
				group3 = []
				
				for item in pairs[:]:
					
					if item[0]<= 0.3:
						
						group1.append(item)
						
					elif 0.4 <= item[0] <=0.7:
						
						group2.append(item)
						
					else:
						
						group3.append(item)	
				
					
				conserved_residues = open(newpath + '/%s_conserved_residues.txt' %query, 'w')

				respo = []
				
				for items in group1[:]:
					
					res_po = tuple((items[2],items[1]))
					respo.append(res_po)
					
				strs =" ".join(str(x) for x in respo)
				conserved_residues.write(strs+"\n")
				conserved_residues.close()
			
rate4site()
###################--pdb_folders.py--###################################

import os
import re
import glob
import shutil

print "		step 2.5 : Create the pdb folders with opm files"
   

curdir = os.getcwd()
curdir_up = '/'.join(curdir.split('/')[:-1])
rate_path = "%s/RESULT/MODULE2" %curdir_up


for curdir, subdirs, files in os.walk(rate_path): 

	for file in files[:]:
		
		if '.aln' in file:
				
			newfile = open(os.path.join(curdir,file)).read().split('>')
			
			for ite in newfile[:]:
				
				if "PDBID" in ite:
					
					fold = ite.split('|')[0]
					new_fold = fold.lower()
					new_folds = new_fold.split(":")[0]
					
					if os.path.isdir(os.path.join(curdir,new_folds)):
						
						pass
					else:
						
						os.makedirs(os.path.join(curdir,new_folds))
					
					 
src = os.chdir(curdir_up + "/alpha")
alpha_files = glob.glob("*")

for curdir, subdirs, files in os.walk(rate_path): 

	for subs in subdirs[:]:
		
		pdbs = subs.split(":")[0] + ".pdb"
		
		if pdbs in alpha_files:
			
			src = curdir_up + "/alpha/%s" %pdbs
			dest = curdir + "/%s" %subs
			shutil.copy2(src, dest)

print "End of Module 2"	

########################################################################
							MODULE3
########################################################################
#########################--filter.py--##################################
import os

print "		step 3.1 : Create the coordinate file of the opm files"

curdir = os.getcwd()
curdir_up = '/'.join(curdir.split('/')[:-1])
rate_path = "%s/RESULT/MODULE2" %curdir_up

for curdir, subdirs, files in os.walk(rate_path): 

	for file in files[:]:
		
		if '.pdb' in file:
			
			atom_data = open(os.path.join(curdir, file),'r').readlines()
			outfile = open(os.path.join(curdir, file.split(".")[0]) + "ATOM", 'w')
			for atoms in atom_data:
				if atoms.startswith('ATOM'):
					outfile.write(atoms)
######################--filter_larger.py--##############################
import os

print "		step 3.1 : Create the coordinate file of the opm files"

curdir = os.getcwd()
curdir_up = '/'.join(curdir.split('/')[:-1])
rate_path = "%s/RESULT/MODULE2-LARGER" %curdir_up

for curdir, subdirs, files in os.walk(rate_path): 

	for file in files[:]:
		
		if '.pdb' in file:
			
			atom_data = open(os.path.join(curdir, file),'r').readlines()
			outfile = open(os.path.join(curdir, file.split(".")[0]) + "ATOM", 'w')
			for atoms in atom_data:
				if atoms.startswith('ATOM'):
					outfile.write(atoms)
			outfile.close()
			
for curdir, subdirs, files in os.walk(rate_path): 

	for file in files[:]:
		
		if 'ATOM' in file:
			file = os.path.join(curdir, file)
			os.system('grep CB %s > %sCB' %(file,file))
			print file
print 'done'
					
#########################--box.py--#####################################
#!/usr/bin/python
def box(com_2): # com is list, xyz is int
	
	x = 1
	y = 1
	z = com_2[2]
	xs = []
	ys = []
	zs = []
	points = []
	p1 = [com_2[0]+x, com_2[1], com_2[2]+15]
	points.append(p1)
	p2 = [com_2[0]-x, com_2[1], com_2[2]+15]
	points.append(p2)
	p3 = [com_2[0], com_2[1]+y, com_2[2]+15]
	points.append(p3)
	p4 = [com_2[0], com_2[1]-y, com_2[2]+15]
	points.append(p4)
	p5 = [com_2[0]+x, com_2[1], com_2[2]-15]
	points.append(p5)
	p6 = [com_2[0]-x, com_2[1], com_2[2]-15]
	points.append(p6)
	p7 = [com_2[0], com_2[1]+y, com_2[2]-15]
	points.append(p7)
	p8 = [com_2[0], com_2[1]-y, com_2[2]-15]
	points.append(p8)
	
	for point in points:
		
		xs.append(point[0])
		ys.append(point[1])
		zs.append(point[2])
		
	return xs,ys,zs
#########################--com_2.py--#####################################

#!/usr/bin/python
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
from itertools import product, combinations
import os
import re
import math
from pymol.cgo import *
from pymol import cmd
from pymol.vfont import plain
import box
import glob
from os.path import basename
#import invokepores

print "		step 3.2 : Calculating the center of Mass of the pore and creating the BOX"

curdir = os.getcwd()
curdir_up = '/'.join(curdir.split('/')[:-1])
rate_path = "%s/RESULT/MODULE2" %curdir_up

for curdir, subdirs, files in os.walk(rate_path):

	for file in files[:]:

		if '.pdb' in file:
		
			new = file.split(".")[0]
			pdb = open(os.path.join(curdir , file),'r').readlines()

			charged_res = []
			x =[]
			z =[]
			y =[]
			com =[]
			z_cutoff = []
			
			for item in pdb:
				
				if item.startswith('ATOM'):

					xcoor = item[31:38]
					ycoor = item[39:46]
					zcoor = item[47:54]
					x.append(float(xcoor))
					y.append(float(ycoor))
					z.append(float(zcoor))
					
			for item in z[:]:
				
				if -15.0<= item and item <=15.0:
					z_cutoff.append(item)
					
					
			x_cord = sum(x)/len(x)
			
			y_cord = sum(y)/len(y)
			
			try:
				
				z_cord = sum(z_cutoff)/len(z_cutoff)
				
			except ZeroDivisionError:
				
				continue
				
			com.append([x_cord,y_cord,z_cord])
			out1 = open(os.path.join(curdir,new + 'out' + '_com ' + '.txt'),'w')
			out1.write(str(com))
			out1.close()
			
			one_comp = box.box([x_cord,y_cord,z_cord])
			out2 = open(os.path.join(curdir,new + 'out' + '_box ' + '.txt'),'w')
			out2.write(str(one_comp))
			out2.close()
		
			x_max = max(one_comp[0])
			x_min = min(one_comp[0])
			y_max = max(one_comp[1])
			y_min = min(one_comp[1])
			z_max = max(one_comp[2])
			z_min = min(one_comp[2])
			
			ranges = [x_max,x_min,y_max,y_min,z_max,z_min]
			out3 = open(os.path.join(curdir,new + 'out' + '_range ' + '.txt'),'w')
			out3.write(str(ranges))
			out3.close() 
			
			###########################3d-plot###################################

			fig = plt.figure()
			ax = fig.gca(projection='3d')
			ax.set_aspect("equal")
			centre = ax.scatter(com[0][0],com[0][1],com[0][2],color="g",s=100)
			points = box.box((com[0][0],com[0][1],com[0][2]))
			scatterfile = ax.scatter(points[0], points[1], points[2],color = "r")
			os.chdir(os.path.join(curdir))
			
			plt.savefig('%s_box.png'%new)
			plt.close("all")
			
invokepores.invoke()

#####################--invokepores.py--#################################

import os
import subprocess

print "		step 3.3 : Create the coordinate file of the opm files"

def invoke():
	
	bpath = '/home/jacob/Downloads/PROPORES'
	os.environ['PERLLIB'] = bpath + '/module_lib' 
	os.environ['ROTALIB'] =  bpath + '/ROTA_lib_30'
	subprocess.check_call(['perl', bpath + '/Pore_ID.pl' ,'-f', bpath + '/3rfuATOM', '-r' ,'1.0','-s','1.2', '-c' ,'1.4', '-n' ,'3rfu'], env = os.environ)

###################--get-plrs.py--######################################

import os
import glob
import re 
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import box


fig = plt.figure()
ax = fig.gca(projection='3d')
ax.set_aspect("equal")

print "		step 3.4 : Filtering out the PLRs"

def sorted_nicely(l):
	
    convert = lambda text: int(text) if text.isdigit() else text 
    alphanum_key = lambda key:[convert(c) for c in re.split('([0-9]+)',key)] 
    return sorted(l, key = alphanum_key)

curdir = os.getcwd()
curdir_up = '/'.join(curdir.split('/')[:-1])
rate_path = "%s/RESULT/MODULE2" %curdir_up
ftp_path = "%s/RESULT/MODULE1" %curdir_up
range_paths = []
ftp_fastas = []

for curdir, subdirs, files in os.walk(ftp_path): 
	
		for file in files[:]:	
			
			if '.aln' in file:
				
				name = file.split('.')[0]
				fasta_path = os.path.join(curdir, name + ".fasta")
				ftp_fastas.append(fasta_path)	
			
for curdir, subdirs, files in os.walk(rate_path): 

	for file in files[:]:

		if 'range' in file:
			
			filepath = os.path.join(curdir, file)
			range_paths.append(filepath)			

source_files = []

for paths in range_paths:
	
	srce = paths.split('/')[:-1]
	srce_path = '/'.join(srce)
	source_files.append(srce_path)	

for i,srce in enumerate(source_files[:]):
	
	source_file = os.chdir(srce)	
	files = sorted_nicely(glob.glob('*.PTin'))
	zcors = []
	xcors = []
	ycors=[]
	coins = []
	
	for file in sorted_nicely(files[:]):
		
		cord = open(file,'r').read().splitlines()[2:-1]
		record = False
		for line in cord:
			
			f = line.split()
			x = float(f[0:3][0])
			y = float(f[0:3][1])
			z = float(f[0:3][2])
			plrs = [x,y,z]
			range_file = range_paths[i]

			with open(range_file) as box_range:
				
				cords = box_range.read().strip('[').strip(',').strip(']').split(',')
				
				if z >= float(cords[5]) and z<= float(cords[4]) and y>= float(cords[3]) and y <= float(cords[2]) and x>= float(cords[1]) and x<= float(cords[0]):
					
					record = True
					xcors.append(x)
					ycors.append(y)
					zcors.append(z)
						
		if record:
			
			coins.append(file)

	listcontent = ''
	plr = open('plr.txt', 'w')
	
	for lis in coins:
		
		name = lis[:-4] + 'list'
		content =  open(name).readlines()[:-1]
		content = ''.join(content)
		listcontent += content
		
	plr.write(listcontent)
	plr.close()
	
	plr_list = open('plr_list.txt','w')
	unipid = srce.split('/')[-2]
	
	for fas_path in ftp_fastas[:]:
		
		if unipid in fas_path:

			pdb_chain = open(fas_path).read().split(">")
			
			for pdbids in pdb_chain[:]:
				
				pds_identifier = pdbids.find('PDBID')
				
				if pds_identifier >=1:
					
					chain = pdbids.split("|")[0]
					chain_info = chain.split(":")[-1]

			merge = open('%s/plr.txt'% srce).readlines()
			final = []
			
			for info in merge[:]:
				
				final_lists = list(info.split()[1:])
				final.append(final_lists)
				
			for resi in final[:]:
				
				resi = '\t'.join([str(x) for x in resi])
				plr_list.write("%s\n" %resi)
						
			plr_list.close()	

			ax.scatter([x for x in xcors],[y for y in ycors],[z for z in zcors], alpha = 0.2,color = "r")
	
			comname = range_paths[i].split('/')[-1].split('_')[0]+"_com .txt"
			compath = srce +'/' + comname
			comfiles = open(compath).read().strip('[').strip(']').split(',')
		
			points = box.box([float(comfiles[0]),float(comfiles[1]),float(comfiles[2])])
			plt.savefig('filter.png')
			plt.close("all")

############################--plr_single_lettercode.py--################

import os
import glob

print "		step 3.5 : Filtering out the PLRs"

curdir = os.getcwd()
curdir_up = '/'.join(curdir.split('/')[:-1])
rate_path = "%s/RESULT/MODULE2" %curdir_up

for curdir, subdirs, files in os.walk(rate_path): 
	
		for file in files[:]:	
			
			if 'plr_list.txt' in file:
				
				amino_name = ["ALA","CYS","ASP","GLU","PHE","GLY","HIS","ILE","LYS","LEU","MET","ASN","PRO","GLN","ARG","SER","THR","VAL","TRP","TYR"]
				code = ["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"] 
				
				amino_dict = {}
				
				destination = os.path.join(curdir,file)
					
				for i in range(len(amino_name)):
					amino_dict[amino_name[i]]=code[i]
				
				os.chdir('/'.join(destination.split('/')[:-1]))
				
				plr_files = set(glob.glob('plr_list.txt'))
				
				dictl = {}
				positions = []
				for line in plr_files:
					hold = True
					contents = open(line).readlines()
				
					for  content in contents[:]:
						con = content.split()
						key = con[1]
				 		pos = content.split()[0]
						positions.append(pos)
					
						if key not in dictl:
							dictl[key] = [pos]
						elif key in dictl:
							dictl[key].append(pos)
				
				newdict = {}
				
				for key, val in dictl.iteritems():
				
					new_key = amino_dict[key]
					new_val = val
					
					if new_key not in newdict:
						newdict[new_key] = [new_val]
					elif new_key in newdict:
						newdict[new_key].append(new_val)
						
					 	 
				plrs = open('plrs.txt','w')
				for items in newdict.keys():
					lines = items
					hit = newdict[items]
					for hits in hit:
						lines +='\t%s' % hits
					plrs.write("%s\n" %lines)
	
print "End of Module 3"

########################--plr_position.py--#############################

import os
import re
from collections import namedtuple
import numpy as np
import matplotlib.pyplot as plt 
from pylab import *
from weblogolib import *
import pylab as plt
from matplotlib.legend_handler import HandlerLine2D

print "		step 4 : Mapping the PLRs and Conserved residues on the primary sequence of protein"

curdir = os.getcwd()
curdir_up = '/'.join(curdir.split('/')[:-1])
rate_path = "%s/RESULT/MODULE2" %curdir_up

plr_destinations = []
clustal_files = []
pdb_files = []
fasta_files = []
conserved = []

for curdir, subdirs, files in os.walk(rate_path): 
	
		for file in files[:]:	
			
			if 'plrs.txt' in  file:
				
				plr_destination = os.path.join(curdir,file)
				plr_destinations.append(plr_destination)

			if '.aln' in file:

				clustal_file = os.path.join(curdir,file)
				clustal_files.append(clustal_file)
					
			if 'ATOM' in file:
				
				pdb_file = os.path.join(curdir,file)
				pdb_files.append(pdb_file)
				
			if '.fasta' in file:
				
				fasta_file = os.path.join(curdir,file)
				fasta_files.append(fasta_file)
			
			if '_conserved_residues.txt' in file:
			
				conserved_file = os.path.join(curdir, file)
				conserved.append(conserved_file)
		
for ipdb,pdb_file in enumerate(pdb_files):
	
	curdir2 = '/'.join(pdb_file.split('/')[:-1])
	
	fasta_file = ''
	clustal_file = ''
	conserved_file = ''
	
	plr_destination = plr_destinations[ipdb]
	protein_name = str.upper(pdb_file.split('/')[-1].strip('ATOM'))
	protein_name_lower = pdb_file.split('/')[-1].strip('ATOM')

	strain_name = pdb_file.split('/')[-3]
	
	for fastafile in fasta_files:
		
		if strain_name in fastafile:
			
			fasta_file = fastafile
			
	
	for clustalfile in clustal_files:
		
		if strain_name in clustalfile:
			
			clustal_file = clustalfile
				
	for cons_set in conserved:
		
		if strain_name in cons_set:
			
			conserved_file = cons_set
	
							
	clustals = open(clustal_file).read().split('>')
	pdbs = open(pdb_file).readlines()
	fastas = open(fasta_file).read().split('>')	
	data_conserved_file = open(conserved_file).read()
	
	clustal_dict1 = {}
	
	for qids in clustals[:]:

		if protein_name in qids :
			
			strs = qids.split()
			str_ids = strs[0]
			sequences = ''.join(strs[1:]).split()

			if str_ids not in clustal_dict1:
				
				clustal_dict1[str_ids] = sequences
				
			elif str_ids in clustal_dict1:
				
				clustal_dict1[str_ids].append(sequences)
		else:
			pass
		
	#print clustal_dict1 
	#break
	
	clustal_dict2 = {}
	
	for seqs in clustal_dict1:
		
		sequences = clustal_dict1[seqs]
		for residues in sequences:
			
			ind_resis = list(residues)
			
			pairs = []
			for iinres, ind_resi in enumerate(ind_resis):
				pairs.append([iinres+1, ind_resi])
			clustal_dict2[seqs] = pairs
	
	#print clustal_dict2, clustal_dict1
	#break
	
	pdb_dict1 = {}
	
	for maps in pdbs:
		
		residue_name = maps[17:20]
		residue_sequence_number = int(maps[23:26])
	
		if residue_sequence_number not in pdb_dict1:
			
			 pdb_dict1[residue_sequence_number] = [residue_name]
			 
		elif residue_sequence_number in pdb_dict1:
			
			pdb_dict1[residue_sequence_number].append(residue_name)
	
	#print pdb_dict1,ipdb
	#break
		
	amino_name = ["ALA","CYS","ASP","GLU","PHE","GLY","HIS","ILE","LYS","LEU","MET","ASN","PRO","GLN","ARG","SER","THR","VAL","TRP","TYR"]
	code = ["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"] 
	
	amino_dict = {}
	
	for i in range(len(amino_name)):
		
		amino_dict[amino_name[i]]=code[i]
	pdb_dict2 = {}
	
	for pos in pdb_dict1:
		
		res = pdb_dict1[pos]
		
		for resi in res:
			
			single = amino_dict[resi]
			
			if pos not in pdb_dict2:
				
				pdb_dict2[pos] = [single]
				
			elif pos in pdb_dict2:
				
				pdb_dict2[pos].append(single)	
	
	#print len(pdb_dict2),len(pdb_dict1)
	#break
	
	pdb_dict3 = {}	
	
	for single_items in pdb_dict2:
		
		pdb_dict3[single_items] = pdb_dict2[single_items][0]
	
	#print len(pdb_dict3)
	#break
	
	start_pos = min(pdb_dict3.keys())
	
	start_res = pdb_dict3[start_pos]
	
	##print start_pos, start_res
	##break
	
	new_clustal_dict1 ={}	
			
	for keys in clustal_dict2:
		
		values = clustal_dict2[keys]
		new_clustal_dict1[keys] = []
		start = False
		start_pos2 = start_pos
		
		for items in values:
			
			if items[1] == '-':
				
				pass
				
			else:	
				
				res =  items[1]
		
				if res == start_res and not start:
					
					new_pair  = start_pos2, res
					new_clustal_dict1[keys].append(new_pair)	
					start = True
					start_pos2 += 1
					
				elif start:
					
					new_pair  = start_pos2, res
					new_clustal_dict1[keys].append(new_pair)
					start_pos2 += 1
	
	
	for clus in new_clustal_dict1.keys()[:]:
		
		line = clus+ '\n'
		hits = new_clustal_dict1[clus]
		clustal_dict_file = open(curdir2 + '/%sclustal_file' % clus, 'w')
		content = ''
		
		for hit in hits:
			
			newcontents = str(hit[0]) + ',' + str(hit[1]) + '\n'
			content += newcontents
		clustal_dict_file.write(content)
		
	clustal_dict_file.close()					
			
	#print new_clustal_dict1,clustal_dict2
	#break
	
	new_clustal_dict2={}	
	plr_file = open(plr_destination).readlines()
	plrs_tuple = []
	
	for i,plr in enumerate(plr_file[:]):
		
		residue =  plr.split()[0]
		plr_s =  re.findall(r"'\s*([^']*?)\s*'",plr[1:-1])
		newplrs = [int(x) for x in plr_s]
		
		for i2,items in enumerate(newplrs[:]):
			
			plr = tuple((newplrs[i2],residue))
			plrs_tuple.append(plr)	
			
	vs = new_clustal_dict1.values()
	ks = new_clustal_dict1.keys()
	
	for j,k in enumerate(ks):
		
		new_clustal_dict2[k] = []
		for n in vs[j]:
			
			if n in plrs_tuple:
				
				n2 = (n[0],'p')
				new_clustal_dict2[k].append(n2)

			else:
				
				n = (n[0],'*')
				new_clustal_dict2[k].append(n)
			
	#print new_clustal_dict2, new_clustal_dict1
	#break
	
	fastafile1 = {}
	
	for s,seq in enumerate(fastas[:]):
		
		if protein_name in seq:
			
			fasta_strs = seq.split()
			fasta_ids = fasta_strs[0]
			fasta_seqs = ''.join(fasta_strs[1:]).split()
	
			if fasta_ids not in fastafile1:
				
				fastafile1[fasta_ids] = fasta_seqs
				
			elif fasta_ids in fastafile1:
				
				fastafile1[fasta_ids].append(fasta_seqs)
			
	#print fastafile1
	#break
	
	fastafile2 = {}
	
	for seqs in fastafile1:
		
		sequences = fastafile1[seqs]
		
		for residues in sequences:
			
			ind_resis2 = list(residues)
			pairs = []
	
			for ires, ind_residues in enumerate(ind_resis2):
				
				pairs.append([ires+1, ind_residues])
			fastafile2[seqs] = pairs
			
	#print fastafile1, fastafile2
	#break
	
	fastafile3 = {}
	pre = ''
	post = ''
	for items in new_clustal_dict2.keys():
		
		line = items + '\n'
	
		min_pdb_len = min(new_clustal_dict2[items])
		pos1,res1 = tuple(min_pdb_len)
		max_pdb_len = max(new_clustal_dict2[items])
		pos2,res2 = tuple(max_pdb_len)
		
		#print min_pdb_len,max_pdb_len	
	
		for item in fastafile2.keys():
			
			newline = item + '\n'
			min_newseqdat_len = min(fastafile2[item])
			pos3,res3 = tuple(min_newseqdat_len)
			max_newseqdat_len = max(fastafile2[item])
			pos4,res4 = tuple(max_newseqdat_len)
			
			#print min_newseqdat_len, max_newseqdat_len	
			
			min_len = pos1 - pos3
			max_len = pos4 - pos2
			
			#print min_len, max_len
			
			for plr in new_clustal_dict2.keys():
					
				pre = range(1, min_len + 1)
				post = range(1, max_len +1)
			
		for num in pre:
			
			positions = min_pdb_len[0] - num
			new_pos = (positions, '*')
			new_clustal_dict2[items].insert(0,new_pos)
		
		for num1 in post:
			
			positions1 = max_pdb_len[0] + num1
			new_pos1 = (positions1, '*')
			new_clustal_dict2[items].append(new_pos1)  
		fastafile3[items] = new_clustal_dict2[items]	
	
	#print fastafile3, fastafile2
	#print new_clustal_dict2

	data1 = data_conserved_file.split(')')
	data2 = [ item.replace('(','').strip().split(',') for item in data1[:-1]]
	data3 = [(int(item[0]), item[1].strip().strip("'")) for item in data2]
	
	res = new_clustal_dict1.values()
	cons = new_clustal_dict1.keys()
	conserved_set = {}

	for c, con in enumerate(cons):
		
		conserved_set[con] = []
		
		for s in new_clustal_dict1[con]:

			if s in data3:
				
				s2 = (s[0],'c')
				conserved_set[con].append(s2)
				
			else:
			
				s = (s[0],'*')
				conserved_set[con].append(s)
	
	#print conserved_set
	#print new_clustal_dict1
	#print fastafile3

	new_conserved_set= {}
	pre1 = ''
	post1 = ''
	
	for items in conserved_set.keys():
		
		line = items + '\n'
		min_pdb_len = min(conserved_set[items])
		pos1,res1 = tuple(min_pdb_len)
		max_pdb_len = max(conserved_set[items])
		pos2,res2 = tuple(max_pdb_len)
		
		#print min_pdb_len,max_pdb_len	
	
		for item in fastafile2.keys():
			
			newline = item + '\n'
			min_newseqdat_len = min(fastafile2[item])
			pos3,res3 = tuple(min_newseqdat_len)
			max_newseqdat_len = max(fastafile2[item])
			pos4,res4 = tuple(max_newseqdat_len)
			
			#print min_newseqdat_len, max_newseqdat_len	
			
			min_len = pos1 - pos3
			max_len = pos4 - pos2
				
			for plr in conserved_set.keys():
					
				pre1 = range(1, min_len + 1)
				post1 = range(1, max_len +1)
				
			for num in pre1:
				
				positions = min_pdb_len[0] - num
				new_pos = (positions, '*')
				conserved_set[items].insert(0,new_pos)
			
			for num1 in post1:
				
				positions1 = max_pdb_len[0] + num1
				new_pos1 = (positions1, '*')
				conserved_set[items].append(new_pos1) 
				 
			new_conserved_set[items] = conserved_set[items]
	
	
	#print fastafile3
	#print fastafile3.keys()
	#print new_conserved_set.keys()
	#print pdb_file
	#print new_conserved_set
	
	for fast in fastafile3.keys()[:]:

		fastahits = fastafile3[fast]
		fasta_dict_file = open(curdir2 + '/%sfastafile3' % fast, 'w')
		fasta_content = ''
		
		for hit in fastahits:
			
			fastcontents = str(hit[0]) + ',' + str(hit[1]) + '\n'
			fasta_content += fastcontents
		fasta_dict_file.write(fasta_content)
		
	fasta_dict_file.close()				
	
	for consd in new_conserved_set.keys()[:]:
		
		conservedhits = new_conserved_set[consd]
		conserved_dict_file = open(curdir2 + '/%sconserved_set_file' % consd, 'w')
		conserved_content = ''
		
		for hit in conservedhits:
			
			conscontents = str(hit[0]) + ',' + str(hit[1]) + '\n'
			conserved_content += conscontents
		conserved_dict_file.write(conserved_content)
		
	conserved_dict_file.close()	
	
	plr_pos = open(curdir2 + '/plr_pos', 'w')
	
	#print fastafile3, new_conserved_set
	#print pdb_file
	for item in fastafile3.keys()[:]:
		
		line = item + '\n'
		hits = fastafile3[item]	
		residues = []
		residues1 = []
		residues2 = []
		cosposs = []
		poss = []
		for cos in new_conserved_set.keys()[:]:
			
			for costuples in conserved_set[cos]:
				
				cospos,cosres = costuples
				
				if cos in item:
					
					residues.append(cosres) 
					cosposs.append(cospos)
		##print residues	
		
		for tuples in hits:
			
			pos,residue_name = tuples
			residues1.append(residue_name)
			poss.append(pos)
				
		for res in fastafile1.keys()[:]:	
			
			if res in item:
				
				residues2.append(''.join(fastafile1[res]))
				
		#print len(residues) #,cosposs
		#print 'yo'
		#print len(residues1) #, poss
		#print 'yoho'
		#print len(residues2[0]) #, len(residues2[0])
		#print 'yoyoyo' ,pdb_file
	
		
		for hit in residues:		
			
			line += '%s' %hit
		
		line += '\n'
	
		for hits in residues1:		
			
			line += '%s' %hits
		line += '\n'
		
		for hitss in residues2:		
			
			line += '%s' %hitss
			
		plr_pos.write("%s\n" % line)	
	
	plr_pos.close()
	
print "End of Module 4"	
	
#######################--plr_pos_plot.py--##############################

import os
import re
from collections import namedtuple
import numpy as np
import matplotlib.pyplot as plt 
from pylab import *
from weblogolib import *
import pylab as plt
from matplotlib.legend_handler import HandlerLine2D
fig = plt.figure()
ax = axes()
hold = True

with open('/home/jacob/Desktop/RESULT/MODULE2/DAACS/SA2172/1xfh/plr_pos') as f:
	
	data = f.read()		
	
data = data.split()[3]
labels = [x for x in data]

kytefile = open('/home/jacob/Desktop/RESULT/MODULE2/DAACS/SA2172/1xfh/kyte').readlines()
kyte_scores = {}
k_scores = []

for kscore in kytefile:
	
	if kscore.startswith("Position"):
		
		kytescore = float(kscore.strip(" ").split()[3])
		kyteindex = int(kscore.strip(" ").split()[1])
		k_scores.append(kytescore)
		
		if kyteindex not in kyte_scores:
			
			kyte_scores[kyteindex] = kytescore
			
		elif kyteindex in kyte_scores:
			
			kyte_scores[kyteindex].append(kytescore)
			

	
rate_file = open('/home/jacob/Desktop/RESULT/MODULE2/DAACS/SA2172/SA2172_rate','r')
rate_dict = {}
rate_scores_dict ={} 
r_scores = [] 

query_len = 0
for i,interval in enumerate(rate_file):
	
	if i == 0:
		
		True
		
	else:
		
		if not interval.startswith('#') and interval != '\n':
			
			interval = interval.split(',')
			parts1 = interval[0].split()
			rate_key = 'SA2172'
			rate_pos = int(parts1[0])
			rate_res = parts1[1]
			score = float(parts1[2])				
			r_scores.append(float(score))	
			rate_tuple = tuple((rate_pos,rate_res))
			
			if rate_pos not in rate_dict:
				
				rate_dict[rate_pos] = rate_tuple
				
			elif rate_pos in rate_dict:
				
				rate_dict[rate_pos].append(rate_tuple)
			
			if rate_pos not in rate_scores_dict:
				
				rate_scores_dict[rate_pos] = score
				
			elif rate_pos in rate_scores_dict:
				
				rate_scores_dict[rate_pos].append(score) 	
				
#print k_scores
#print r_scores

clustal_files = open('/home/jacob/Desktop/RESULT/MODULE2/DAACS/SA2172/1xfh/1XFH:A|PDBID|CHAIN|SEQUENCEclustal_file').readlines()
pdb_tuples = []
for datas in clustal_files:
	
	data1 = datas.split(',')
	data2 = int(data1[0]),data1[1].strip('\n')
	pdb_tuples.append(data2)

pdb_tuple  = sorted(pdb_tuples)

rate_tuples  = sorted(rate_dict.values())

#print len(pdb_tuple),len(rate_tuples)

mapped_rate = {}
pdb_po_start = 0
rate_po_start = 0
rate_po_index = 0
stop = False

for tup in pdb_tuple:
	
	if stop == False:
		
		pdbpos,pdbres = tup
		for i,ratetup in enumerate(rate_tuples):
			ratepos, rateres = ratetup
			
			if pdbres == rateres:

				pdb_po_start = pdbpos
				rate_po_start = ratepos
				rate_po_index = i
				stop = True
				break		
				
#print pdb_tuple	
#print rate_tuples	

for tup1 in pdb_tuple:

	try: 
		rate_tuple = rate_tuples[rate_po_index]
		key = tup1[0]
		value = (tup1[0],tup1[1],rate_tuple[0],rate_tuple[1])
		mapped_rate[key] = value
		rate_po_index += 1
		
	except IndexError:
		
		break
				
#print kyte_scores
#print rate_scores_dict

graph_dict = {}

for maps in mapped_rate:
	
	vals = mapped_rate[maps]
	rate_key = vals[2]
	pdb_pos = vals[0]
	pdb_res = vals[1]
	graph_keys = pdb_pos
	
	if pdb_pos in kyte_scores:
		
		r,k= rate_scores_dict[rate_key], kyte_scores[pdb_pos]
		graph_vals = (pdb_pos,pdb_res,r,k)
		graph_dict[graph_keys] = graph_vals

#print graph_dict

fasta_filess = open('/home/jacob/Desktop/RESULT/MODULE2/DAACS/SA2172/1xfh/1XFH:B|PDBID|CHAIN|SEQUENCEfastafile3').readlines()
fasta_tuples = []

for fastadatas in fasta_filess:
	
	fasta_data1 = fastadatas.split(',')
	fasta_data2 = int(fasta_data1[0]),fasta_data1[1].strip('\n')
	fasta_tuples.append(fasta_data2)

conserved_files = open('/home/jacob/Desktop/RESULT/MODULE2/DAACS/SA2172/1xfh/1XFH:B|PDBID|CHAIN|SEQUENCEconserved_set_file').readlines()
conserved_tuples = []

for conserveddatas in conserved_files:
	
	conserveddata1 = conserveddatas.split(',')
	conserveddata2 = int(conserveddata1[0]),conserveddata1[1].strip('\n')
	conserved_tuples.append(conserveddata2)

plr_paired = []
not_plr_paired = []

for tup  in fasta_tuples:

	p,r = tup
	
	if '*' not in r:
		
		plr_pair = (p,r)
		plr_paired.append(plr_pair)
		
	else:
		
		not_plr_pair = (p,r)
		not_plr_paired.append(not_plr_pair)
				
#print len(plr_paired),len(not_plr_paired)

conserved = []

for tup1 in conserved_tuples:
	
	p2,r2 = tup1

	if '*' not in r2:
		
		cons_res = (p2,r2)
		conserved.append(cons_res)
		
#print conserved

conserved_plrset = []
not_conserved_plrset = []

for pls in plr_paired[:]:
	
	p_1,r_1 = pls
	
	for conss in conserved[:]:
		
		cp_1,cr_1 = conss
		
		if cp_1 == p_1:
			
			new_con_set = (cp_1,cr_1)
			conserved_plrset.append(new_con_set)
			
		else:
			
			pass	

conserved_nonplr = list(set(conserved)-set(conserved_plrset))

print  len(conserved), len(conserved_plrset), len(conserved_nonplr)

#print conserved_plrset
	
conserved_plr = open("/home/jacob/Desktop/RESULT/MODULE2/DAACS/SA2172/conserved_plrset_1xfh.txt","w")

for point in graph_dict.keys():
	
	p1,r1,r_s,k_s = graph_dict[point]
	ps_list = []
	
	for pa in conserved_plrset:
		
		ps,rs = pa
		ps_list = []
		
		if ps == p1:
			
			ps_list.append((p1,r1,r_s,k_s))
			plt.scatter(p1,r_s,color='blue')
			plt.scatter(p1,k_s,color='green')
			for i, txt in enumerate(ps_list):
				txt2 = (list(txt[:2]) + [txt[3]])
		
				ax.annotate(txt[1],xy =(txt[0],txt[2]),xytext=(0.0,0.1),textcoords='offset points', ha='center', va='bottom',
			bbox=dict(boxstyle='round,pad=0.1', fc='yellow', alpha=0.1),size = 10)
							
				ax.annotate(txt[1],xy = (txt[0],txt[3]),xytext=(0.0,0.1),textcoords='offset points', ha='center', va='bottom',
			bbox=dict(boxstyle='round,pad=0.1', fc='yellow', alpha=0.1),size = 10)

			for dat in ps_list[:]:
				
				dat1 =', '.join(map(str, dat[:]))
				#print dat1
				conserved_plr.write("%s\n" %dat1)
			
conserved_plr.close()
residue_pos = range(1,len(data)+1)
yrange = ylim(-4.5,4.5)

x = [residue_pos]
y = [yrange]		

plt.title('DAACS_1XFH CONSERVED PLRs')
plt.xlabel('SEQUENCE POSITION', fontsize=14)
plt.ylabel('KYTE-RATE SCORES', fontsize=14)

hB, = plot([1],marker='o', label = 'KYTE-SCORE', color = 'green')
hR, = plot([1],marker='o', label = 'RATE-SCORE', color = 'blue')

plt.legend(handler_map={hB: HandlerLine2D(numpoints=2)})
plt.legend(handler_map={hR: HandlerLine2D(numpoints=2)})

hB.set_visible(False)
hR.set_visible(False)

ax.axhline(0, color='red', lw=0.5)

plt.show()

#########################--kyte_score.py--##############################

import os
import glob
import re
import scipy.stats as stats
from matplotlib.font_manager import FontProperties
import matplotlib.pyplot as plt
from pylab import *
import numpy as np
from decimal import *
from scipy.stats import norm
from sklearn import preprocessing
from math import floor, ceil 

folder_path = '/home/jacob/Desktop/RESULT/COMPARE'


groups = []
for root, subfolders, files in os.walk(folder_path):
	
	for folders in subfolders:
		groups.append(folders)


group1 = folder_path+ '/'+ groups[0]
group2 = folder_path + '/' + groups[1]
group3 = folder_path + '/' + groups[2]

points_path = folder_path + '/Points'
points = []
for root, subfolders, files in os.walk(points_path):
	
	for folders in subfolders:
		points.append(folders)
		
points1 = points_path + '/' + points[0] 
points2 = points_path + '/' + points[1]
points3 = points_path + '/' + points[2]
scatter = [points1,points2,points3]


conserved_kyte_scores = []
conserved_rate_scores = []
for items in scatter[:]:
	
	for root, folders, fil in os.walk(items):
		
		for files in fil:
			
			#print files
			counter_1 = 0
			if 'PLR' in files:
				
				plrfile = open(items + '/' + files).readlines()
				con_kyte = []
				con_rate = []
				for con in plrfile:
					
					consk = con.split()[2]
					consk = consk.strip(',')
					conskyte = float(consk)
					con_kyte.append(conskyte)
					consrate = con.split()[3]
					consr = consrate.strip(',')
					consrate = float(consr)
					con_rate.append(consrate)
			
				conserved_kyte_scores.append(con_kyte)	
				conserved_rate_scores.append(con_rate)	
				
#print conserved_kyte_scores[0]

box_items = [group1, group2, group3]

kyte_scaled_scores = []
rate_scaled_scores=[]
kyte_mean_list = []
kyte_sd_list = []
rate_mean_list = []
rate_sd_list = []
for items in box_items[:]:
	
	for root, folders, fil in os.walk(items):
		
		for files in fil:
			
			#print files
			counter_1 = 0
			if 'kyte' in files:
				
				kytefile = open(items + '/' + files).readlines()
				description = []
				for items in kytefile:
					
					if items.startswith("Position"):
						
						descriptions = items.strip(" ").split()
						description.append(descriptions)
				
				
				kytescores = []
				for kscores in description[:]:
					
					kdscores = float(kscores[3])
					kytescores.append(kdscores)
				
				#print min(kytescores), max(kytescores)
				
				
				newscores = [(x - min(kytescores))/(max(kytescores) - min(kytescores)) for x in kytescores]	
				kyte_scaled_scores.append(kytescores)
				
				av_kyte = np.mean(kytescores)
				kyte_mean_list.append(av_kyte)
				StDev_kyte = np.std(kytescores)
				kyte_sd_list.append(StDev_kyte)
				
				print av_kyte,StDev_kyte

			else:
				pass

for items in box_items[:]:
	
	for root, folders, fil in os.walk(items):
		
		for files in fil:
			
			#print files
			counter_1 = 0
			if 'RATE' in files:	
							
				rate_file = open(items+ '/' + files).readlines()
				#print rate_file
				confidence = []
				data_reqs = []
				fil_data = []
				scores = []
				residues = []
				range1 = []
				range2=[]
				y_count = []
				for i,interval in enumerate(rate_file):
					
					if i == 0:
						
						True
						
					else:
						
						if not interval.startswith('#') and interval != '\n':
							
							interval = interval.split(',')
							parts1 = interval[0].split()
							seq = parts1[0]
							res = parts1[1]
							residues.append(res)
							#print res
							y_count.append(seq)
							score = parts1[2]		
							scores.append(float(score))					
				
				#print min(scores),max(scores)
				normalized = [(s-min(scores))/(max(scores)-min(scores)) for s in scores]
				rate_scaled_scores.append(scores)



				av_rate = np.mean(scores) 
				rate_mean_list.append(av_rate)
				StDev_rate = np.std(scores) #Standard Dev.
				rate_sd_list.append(StDev_rate)
				
				print av_rate,StDev_rate

def setBoxColors(bp):
	
    setp(bp['boxes'][0], color='blue')
    setp(bp['caps'][0], color='blue')
    setp(bp['caps'][1], color='blue')
    setp(bp['whiskers'][0], color='blue')
    setp(bp['whiskers'][1], color='blue')
    setp(bp['fliers'][0], color='blue')
    setp(bp['fliers'][1], color='blue')
    setp(bp['medians'][0], color='blue')

    setp(bp['boxes'][1], color='red')
    setp(bp['caps'][2], color='red')
    setp(bp['caps'][3], color='red')
    setp(bp['whiskers'][2], color='red')
    setp(bp['whiskers'][3], color='red')
    setp(bp['fliers'][2], color='red')
    setp(bp['fliers'][3], color='red')
    setp(bp['medians'][1], color='red')

daacs = [kyte_scaled_scores[0],rate_scaled_scores[0]]
mfs = [kyte_scaled_scores[1],rate_scaled_scores[1]]
pot = [kyte_scaled_scores[2],rate_scaled_scores[2]] 


data_to_plot = [daacs,mfs,pot]

fig = plt.figure()
ax = axes()
hold = True

bp = boxplot(daacs,positions = [1.5, 2.5])
plt.scatter([1.5], [kyte_mean_list[0]],color = 'blue')
plt.scatter([1.5],[kyte_mean_list[0]-kyte_sd_list[0]],s=80,color = 'blue', marker=(5, 2))
plt.scatter([1.5],[kyte_mean_list[0]+kyte_sd_list[0]],s=80,color = 'blue', marker=(5, 2))



plt.scatter([2.5], [rate_mean_list[0]],color = 'red')
plt.scatter([2.5],[rate_mean_list[0]-rate_sd_list[0]],s=80,color = 'red', marker=(5, 2))
plt.scatter([2.5],[rate_mean_list[0]+rate_sd_list[0]],s=80,color = 'red', marker=(5, 2))

setBoxColors(bp)

bp = boxplot(mfs, positions = [4.5, 5.5])
plt.scatter([4.5], [kyte_mean_list[0]],color = 'blue')
plt.scatter([4.5],[kyte_mean_list[0]-kyte_sd_list[0]],s=80,color = 'blue', marker=(5, 2))
plt.scatter([4.5],[kyte_mean_list[0]+kyte_sd_list[0]],s=80,color = 'blue', marker=(5, 2))

plt.scatter([5.5], [rate_mean_list[0]],color = 'red')
plt.scatter([5.5],[rate_mean_list[0]-rate_sd_list[0]],s=80,color = 'red', marker=(5, 2))
plt.scatter([5.5],[rate_mean_list[0]+rate_sd_list[0]],s=80,color = 'red', marker=(5, 2))

setBoxColors(bp)

bp = boxplot(pot, positions = [7.5, 8.5])
plt.scatter([7.5], [kyte_mean_list[0]],color = 'blue')
plt.scatter([7.5],[kyte_mean_list[0]-kyte_sd_list[0]],s=80,color = 'blue', marker=(5, 2))
plt.scatter([7.5],[kyte_mean_list[0]+kyte_sd_list[0]],s=80,color = 'blue', marker=(5, 2))

plt.scatter([8.5], [rate_mean_list[0]],color = 'red')
plt.scatter([8.5],[rate_mean_list[0]-rate_sd_list[0]],s=80,color = 'red', marker=(5, 2))
plt.scatter([8.5],[rate_mean_list[0]+rate_sd_list[0]],s=80,color = 'red', marker=(5, 2))

setBoxColors(bp)



for i in [1.5,2.5,4.5,5.5,7.5,8.5]:
	
	if i == 1.5:
		
		y = conserved_kyte_scores[:][0]
		x = np.random.normal(i, 0.04, size=len(y))
		colors = matplotlib.cm.rainbow(np.linspace(0, 1, len(conserved_kyte_scores[:][0])))
		cs = [colors[i//len(x)] for i in range(len(conserved_kyte_scores[:][0])*len(x))]
		matplotlib.pyplot.scatter(x,y,color=cs)
		
	elif i == 2.5:
		
		y = conserved_rate_scores[:][0]
		x = np.random.normal(i, 0.04, size=len(y))
		colors = matplotlib.cm.rainbow(np.linspace(0, 1, len(conserved_rate_scores[:][0])))
		cs = [colors[i//len(x)] for i in range(len(conserved_rate_scores[:][0])*len(x))]
		matplotlib.pyplot.scatter(x,y,color=cs)
		
	elif i == 4.5:
		
		y = conserved_kyte_scores[:][1]
		x = np.random.normal(i, 0.04, size=len(y))
		colors = matplotlib.cm.rainbow(np.linspace(0, 1, len(conserved_kyte_scores[:][1])))
		cs = [colors[i//len(x)] for i in range(len(conserved_kyte_scores[:][1])*len(x))]
		matplotlib.pyplot.scatter(x,y,color=cs)
		
	elif i == 5.5:
		
		y = conserved_rate_scores[:][1]
		x = np.random.normal(i, 0.04, size=len(y))
		colors = matplotlib.cm.rainbow(np.linspace(0, 1, len(conserved_rate_scores[:][1])))
		cs = [colors[i//len(x)] for i in range(len(conserved_rate_scores[:][1])*len(x))]
		matplotlib.pyplot.scatter(x,y,color=cs)
		
	elif i == 7.5:
		
		y = conserved_kyte_scores[:][2]
		x = np.random.normal(i, 0.04, size=len(y))
		colors = matplotlib.cm.rainbow(np.linspace(0, 1, len(conserved_kyte_scores[:][2])))
		cs = [colors[i//len(x)] for i in range(len(conserved_kyte_scores[:][2])*len(x))]
		matplotlib.pyplot.scatter(x,y,color=cs)
		
	elif i ==8.5:
		
		y = conserved_rate_scores[:][2]
		x = np.random.normal(i, 0.04, size=len(y))
		colors = matplotlib.cm.rainbow(np.linspace(0, 1, len(conserved_rate_scores[:][2])))
		cs = [colors[i//len(x)] for i in range(len(conserved_rate_scores[:][2])*len(x))]
		matplotlib.pyplot.scatter(x,y,color=cs)
						
	else:
		
		pass

ylim(-5.0,5.0)
xlim(0,9)

plt.ylabel('KYTE/RATE SCORES', fontsize=14)
ax.set_xticklabels(['DAACS', 'MFS', 'POT'])
ax.set_xticks([2.0, 5.0, 8.0])

# draw temporary red and blue lines and use them to create a legend
hB, = plot([1, 1],'b-')
hR, = plot([1, 1],'r-')

fontP = FontProperties()
fontP.set_size('large')
legend((hB, hR),('KYTE-DOOLITTLE SCORES', 'RATE SCORES'), prop = fontP)
hB.set_visible(False)
hR.set_visible(False)

plt.show()
########################################################################
						STATISTICS
########################################################################	
							Wiscos test
########################################################################	
import csv
import os
import scipy.stats as stats
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pylab import *
from matplotlib.font_manager import FontProperties
import math

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_aspect("equal")

os.chdir('/home/jacob/Desktop/RESULT/')

with open('DOC sheet.csv','rb') as csvfile:
	reader = csv.reader(csvfile)
	#print reader
	family = []
	conserved_plr = []
	conserved_nonplr = []
	for row in reader:
		if 'TRANSPORT FAMILY' not in row:
			new = row[2].strip()
			family.append(new)
			family = [x for x in family if x!= '']
			new1= row[6].strip()
			conserved_plr.append(new1)
			conserved_plr = [x for x in conserved_plr if x!= '']
			new2= row[10].strip()
			conserved_nonplr.append(new2)
			conserved_nonplr = [x for x in conserved_nonplr if x!= '']


transport_ids = len(family)
index = np.arange(transport_ids)

width = 0.35

con_plr = sorted([float(x) for x in  conserved_plr])
con_nonplr = sorted([float(x) for x in conserved_nonplr])
print len(con_plr), len(con_nonplr)
mean_con_plr = sum(con_plr)/float(len(con_plr))
mean_con_nonplr = sum(con_nonplr)/float(len(con_nonplr))
std_con_plr =  np.std(con_plr)
std_con_nonplr = np.std(con_nonplr)	
#distance_between_means = 0.5/np.std(con_nonplr)

rects1 = ax.bar(index, con_plr, width, color = 'black' , yerr = std_con_plr, error_kw= dict(elinewidth = 2, ecolor = 'red'))
rects2 = ax.bar(index+width, con_nonplr, width, color= 'red', yerr = std_con_nonplr, error_kw = dict(elinewidth = 2, ecolor = 'black'))
		

ax.set_xlim(-width, len(index)+width)
ax.set_ylim(0,5)

plt.title('HYPOTHESIS TESTING ')
plt.ylabel('MEAN DOC for CONSERVED PLRs/nonPLRs', fontsize=14)

xTickmarks = [x for x in  family]
ax.set_xticks(index+width)
xtickNames = ax.set_xticklabels(xTickmarks)
plt.setp(xtickNames, rotation = 90, fontsize = 10)

ax.legend((rects1[0],rects2[0]), ('CONSERVED PLRs','CONSERVED NON-PLRs'))

#plt.show()

ax.set_ylim(0,50)
plt.title('HYPOTHESIS TESTING HISTOGRAM')
plt.ylabel('MEAN DOC', fontsize=14)
plt.xlabel('Bins', fontsize = 14)
hB, = plot([1, 1],'g-')
hR, = plot([1, 1],'b-')
fontP = FontProperties()
fontP.set_size('small')
legend((hB, hR),('NONPLRs', 'PLRs'), prop = fontP)
hB.set_visible(False)
hR.set_visible(False)
#con_fit = stats.norm.pdf(con_plr, np.mean(con_plr), np.std(con_plr))
#noncon_fit = stats.norm.pdf(con_nonplr, np.mean(con_nonplr), np.std(con_nonplr))

con_fit = stats.wilcoxon(con_plr,con_nonplr,zero_method='wilcox', correction=False)[1]
#noncon_fit = stats.wilcoxon(con_nonplr, y=None,zero_method='wilcox', correction=False)
print con_fit
#plt.plot(con_plr,con_fit,'-o')
#plt.plot(con_nonplr,noncon_fit,'-o')
#plt.hist(con_plr, bins=20,normed = True)
#plt.hist(con_nonplr, bins= 20, normed= True)
#plt.scipy.stats.wilcoxon(con_plr, bins=20,normed = True)
#plt.scipy.stats.wilcoxon(con_nonplr, bins= 20, normed= True)

#plt.show()

#print mean_con_plr,mean_con_nonplr
#print std_con_plr, std_con_nonplr


#con_plr_overlap = 0.5 *(1 + math.erf((min(con_plr) - mean_con_nonplr)/(std_con_nonplr * (2**0.5))))
#con_nonplr_overlap = 1- (0.5 *(1 + math.erf((max(con_nonplr) - mean_con_plr)/(std_con_plr * (2**0.5)))))
 
#print con_plr_overlap, con_nonplr_overlap


########################################################################
						HYpergeometric test
########################################################################	
import csv
import os
from scipy.stats.stats import pearsonr 
from scipy.stats import hypergeom
import scipy.stats as stats
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pylab import *
from matplotlib.font_manager import FontProperties
import math
import numpy 
import seaborn as sns
import pandas as pd

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_aspect("equal")

folder_path = '/home/jacob/Desktop/RESULT/'

groups = []

for root, subfolders, files in os.walk(folder_path):
	
	for folders in subfolders:
		groups.append(folders)


group1_larger = folder_path + groups[2]
group2_smaller = folder_path + groups[5]

grp_smaller = []

for root, subfolders, files in os.walk(group2_smaller):
	
	for fils in files:
		
		if 'plr_pos' in fils and fils != 'plr_pos.png':
			
			map_files_smaller = root + '/' + fils
			grp_smaller.append(map_files_smaller)
			


grp_larger = []
for root, subfolders , files in os.walk(group1_larger):
	
	for fils in files:
		
		if 'plr_pos' in fils and fils != 'plr_pos.png' :
			
			map_files_larger = root + '/' + fils
			grp_larger.append(map_files_larger)

whole_set = grp_smaller +grp_larger
print len(grp_smaller),len(grp_larger),len(whole_set)


pvalue = []
hypergom = []
all_names = []

def comp_hyper(pop, plrs, cons, plr_con):
	result = 0
	for i in range(plr_con,pop+1):
		rv = hypergeom(pop, plrs, cons)
		result += rv.pmf(i)
	return result 

for plr in whole_set:
	
	plrfiles = plr
	plr_pos = open(plrfiles).readlines()

	xy = plr_pos[-3:-1]
	pdb  =plr_pos[-1][:-1]

	row_c =[]
	row_p =[]
	len_r = []
	no_c = []
	no_p = []
	no_cp = []
	
	names = plrfiles.split('/')[8]
	all_names.append(names)
	for i,item in enumerate(pdb):
		
		len_r.append(i)
		c,p = xy[0][i], xy[1][i]
		#print c,p,item,plrfiles
		if c == 'c':
			row_c.append(int(1))
			no_c.append(int(1))
	
		elif c == '*':
			row_c.append(int(0))
									
		if p == 'p':
			row_p.append(int(1))
			no_p.append(int(1))
				
		elif p == '*':
			row_p.append(int(0))
		
		if c == 'c' and p == 'p':
			no_cp.append(int(1))
		
		else:
			
			"nothing"
	
	#print len(len_r), len(no_p), len(no_c),len(no_cp),plrfiles
	hypergom.append(comp_hyper(len(len_r), len(no_p), len(no_c), len(no_cp)))
	##print row_c,row_p
	#pearson_corr = pearsonr(row_c,row_p)
	#pearson.append(pearson_corr[0])
	#pvalue.append(pearson_corr[1])
	#print pearson_corr, plrfiles

print(hypergom)
#print(all_names)
#plt.title('Hypergeometric test')
#plt.xlabel('CORRELATION VALUES', fontsize=10)
#plt.ylabel('P-VALUES', fontsize=10)
#ax.axhline(0, color='red', lw=0.5)
sns.set_style("whitegrid")

data=pd.DataFrame({'Protein':all_names, 'P-value':hypergom})
g=sns.FacetGrid(data,aspect=2,size=6)
g.map(sns.barplot, 'Protein','P-value')

plt.show()
	

########################################################################	
		

	
	
	
	

	


	
		
		




	
	
				
				

	
	
	


