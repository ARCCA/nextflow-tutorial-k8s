import argparse
import os, sys, getopt
import pandas as pd
import itertools
import re
import fnmatch
import subprocess

DIR = sys.argv[1]

# DIR = '/scratch/c.mcbsd1/data.nextflow/ds0134'
# RESOURCES = '/scratch/c.mcbsd1/data.nextflow/ds0134/targets.csv'
# sampleName = 'R154-H-001'

RESOURCES = sys.argv[2]

sampleName = sys.argv[3]

a1=[]
fileExtN=[]

f1 = os.listdir(DIR)

for filename in f1:
 if filename.endswith(('.fastq.gz', '.fq.gz', '.filt.fastq.gz', '.filt.fq.gz')):
  a1.append(filename)

#for root, dirs, files in os.walk(DIR):
#    for filename in files:
#        if filename.endswith(('.fastq.gz', '.fq.gz')):
#                a1.append(filename)

if a1[0].endswith(('.fastq.gz')):
        fileExtN.append('.fastq.gz')
elif a1[0].endswith(('.fq.gz')):
        fileExtN.append('.fq.gz')

a1=[x for x in a1 if not x.startswith('Undetermined') ]
targetsfile = str(RESOURCES)
targets = pd.read_csv(targetsfile)

a21 = []

for s in targets['suppliedID']:
	t1 = [i for i in targets['suppliedID'] if i in s]
	a21.append(t1[0])

a21=[]
for a in a1:
	for b in targets['suppliedID']:
		if (a.find(b) == -1):
			continue
		else:
			#print(a)
			a21.append(a)

#print(a21)

r1 = sampleName

#r1='R154-H-001'

y1 = targets[targets['analysisID'] == r1][['suppliedID']]

y2 = list(y1['suppliedID'])

#print(y2)

read1ExtType1 = '_1'+str(fileExtN[0])
read2ExtType1 = '_2'+str(fileExtN[0])
read1ExtType2 = '_R1'+str(fileExtN[0])
read2ExtType2 = '_R2'+str(fileExtN[0])
read1ExtType3 = '_R1_001'+str(fileExtN[0])
read2ExtType3 = '_R2_001'+str(fileExtN[0])
read1ExtType4 = '_1_001'+str(fileExtN[0])
read2ExtType4 = '_2_001'+str(fileExtN[0])
read1ExtType5 = '.R1'+str(fileExtN[0])
read2ExtType5 = '.R2'+str(fileExtN[0])
read1ExtType6 = 'R1_001_combined'+str(fileExtN[0])
read2ExtType6 = 'R2_001_combined'+str(fileExtN[0])

#a31=[]
#for a in a21:
#	if (a.find(y2[0]) == -1):
#		continue
#	else:
#		#print(a)
#		a31.append(a)

a31=[]
for a in a21:
	if (a.find(y2[0]) == -1):
		continue
	else:
		#print(a)
		a31.append(a)

r1filesN = []
r2filesN = []
#print(a31)
mode = 0o666
os.mkdir(r1,mode)
subprocess.call(['chmod', '755', r1])
for j in range(0,len(a31)):
	if read1ExtType1 in a31[j]:
		r1filesN.append(DIR+'/'+a31[j])
		r1filesJ1 = ' '.join(r1filesN)
		fname1 = r1+'_1.fastq.gz'
		fname2 = r1+'_2.fastq.gz'
		r1filesJ2 = 'cat ' + r1filesJ1 + ' > ' +  r1 + '/' + fname1
		print(r1filesJ2)
		os.system(r1filesJ2)
		r2filesJN = []
		for string in r1filesN:
			r2s = string.replace("_1", "_2")
			r2filesJN.append(r2s)
		r2filesJ1 = ' '.join(r2filesJN)
		r2filesJ2 = 'cat ' + r2filesJ1 + ' > ' +  r1 + '/' + fname2
		print(r2filesJ2)
		os.system(r2filesJ2)
	elif read1ExtType2 in a31[j]:
		r1filesN.append(DIR+'/'+a31[j])
		r1filesJ1 = ' '.join(r1filesN)
		fname1 = r1+'_1.fastq.gz'
		fname2 = r1+'_2.fastq.gz'
		r1filesJ2 = 'cat ' + r1filesJ1 + ' > ' +  r1 + '/' + fname1
		print(r1filesJ2)
		os.system(r1filesJ2)
		r2filesJN = []
		for string in r1filesN:
			r2s = string.replace("_R1_", "_R2_")
			r2filesJN.append(r2s)
		r2filesJ1 = ' '.join(r2filesJN)
		r2filesJ2 = 'cat ' + r2filesJ1 + ' > ' +  r1 + '/' + fname2
		print(r2filesJ2)
		os.system(r2filesJ2)
	elif read1ExtType3 in a31[j]:
		r1filesN.append(DIR+'/'+a31[j])
		r1filesJ1 = ' '.join(r1filesN)
		fname1 = r1+'_1.fastq.gz'
		fname2 = r1+'_2.fastq.gz'
		r1filesJ2 = 'cat ' + r1filesJ1 + ' > ' +  r1 + '/' + fname1
		print(r1filesJ2)
		os.system(r1filesJ2)
		r2filesJN = []
		for string in r1filesN:
			r2s = string.replace("_R1_001", "_R2_001")
			r2filesJN.append(r2s)
		r2filesJ1 = ' '.join(r2filesJN)
		r2filesJ2 = 'cat ' + r2filesJ1 + ' > ' +  r1 + '/' + fname2
		print(r2filesJ2)
		os.system(r2filesJ2)
	elif read1ExtType4 in a31[j]:
		r1filesN.append(DIR+'/'+a31[j])
		r1filesJ1 = ' '.join(r1filesN)
		fname1 = r1+'_1.fastq.gz'
		fname2 = r1+'_2.fastq.gz'
		r1filesJ2 = 'cat ' + r1filesJ1 + ' > ' +  r1 + '/' + fname1
		print(r1filesJ2)
		os.system(r1filesJ2)
		r2filesJN = []
		for string in r1filesN:
			r2s = string.replace("_1_001", "_2_001")
			r2filesJN.append(r2s)
		r2filesJ1 = ' '.join(r2filesJN)
		r2filesJ2 = 'cat ' + r2filesJ1 + ' > ' +  r1 + '/' + fname2
		print(r2filesJ2)
		os.system(r2filesJ2)
	elif read1ExtType5 in a31[j]:
		r1filesN.append(DIR+'/'+a31[j])
		r1filesJ1 = ' '.join(r1filesN)
		fname1 = r1+'_1.fastq.gz'
		fname2 = r1+'_2.fastq.gz'
		r1filesJ2 = 'cat ' + r1filesJ1 + ' > ' +  r1 + '/' + fname1
		print(r1filesJ2)
		os.system(r1filesJ2)
		r2filesJN = []
		for string in r1filesN:
			r2s = string.replace(".R1", ".R2")
			r2filesJN.append(r2s)
		r2filesJ1 = ' '.join(r2filesJN)
		r2filesJ2 = 'cat ' + r2filesJ1 + ' > ' +  r1 + '/' + fname2
		print(r2filesJ2)
		os.system(r2filesJ2)
	elif read1ExtType6 in a31[j]:
		r1filesN.append(DIR+'/'+a31[j])
		r1filesJ1 = ' '.join(r1filesN)
		fname1 = r1+'_1.fastq.gz'
		fname2 = r1+'_2.fastq.gz'
		r1filesJ2 = 'cat ' + r1filesJ1 + ' > ' +  r1 + '/' + fname1
		print(r1filesJ2)
		os.system(r1filesJ2)
		r2filesJN = []
		for string in r1filesN:
			r2s = string.replace("R1_001", "R2_001")
			r2filesJN.append(r2s)
		r2filesJ1 = ' '.join(r2filesJN)
		r2filesJ2 = 'cat ' + r2filesJ1 + ' > ' +  r1 + '/' + fname2
		print(r2filesJ2)
		os.system(r2filesJ2)
