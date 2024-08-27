#! /usr/bin/env python2
#encoding=utf-8
import os
import sys
import operator
#import multiprocessing
#from Bio import SeqIO
import math
import random
num_str=''

def desStartCode(codes):
	if(codes in ("ATG","atg")):
		return True
	return False

def desEndCode(codes):
	if(codes in ("tga", "TGA")):
		return True
	return False

def readByThree(string, offset):
	flag = True
	length = len(string)
	start = end = -1
	i = 0
	result = set()
	while i < length-2:
		codes = string[i:i+3]
		if(desStartCode(codes) and flag):
			start = i
			flag = False
		if(desEndCode(codes) and not flag):
			end = i + 2
			flag = True
		if( (end > start) and (start != -1) ):
			result.add((start+offset, end+offset))
		i = i + 3
	return result

def getGC(string):
	GC = string.count("G") + string.count("C") + string.count("g") + string.count("c")
	GC2 = GC/(len(string)+0.0)
	return GC2

def getInfo(string, pos):
	length = pos[1] - pos[0] +1
	gc = getGC(string[pos[0]:pos[1]+1])
	return str(pos[0]), str(pos[1]), str(length), str(gc)

def writeInfo(resultInfo,name,num_str):
	nf = open('ORF'+num_str+'.temp', "a")
	for r in resultInfo:
		nf.write("\t".join(r) + "\n")
	a = name.strip().split()
	nf.write(a[0] + "\n")
	nf.close()

def length_GC(string):
	G = string.count("G") + string.count("g")
	C = string.count("C") + string.count("c")
	L = len(string)
	return str(L) + '\t' + str(G) + '\t' + str(C) + '\n'

def ORF(f1,num_str):
#	f = open(f1,'rU')

#	records = SeqIO.parse(f1, 'fasta')
#	SeqDict = SeqIO.to_dict(SeqIO.parse(f, "fasta"))

#	for record in records:
#		resultInfo = []
#		line = SeqDict[record.id].seq
#		strings = [line, line[1:], line[2:]]
#		for index, string in enumerate(strings):
#			positions = readByThree(string, index)
#			positions = sorted(positions, key = operator.itemgetter(0))
#			for pos in positions:
#				resultInfo.append(getInfo(line, pos))
#		writeInfo(resultInfo,record.id,num_str)
#	f.close()


	f = open(f1)
	name = f.readline().split('>')
	name1 = name[1].split(' ')
	lines = ''
	for line in f:
		if line.startswith(">"):						
			resultInfo = []
			strings = [lines, lines[1:], lines[2:]]
			for index, string in enumerate(strings):
				positions = readByThree(string, index)
				positions = sorted(positions, key = operator.itemgetter(0))
				for pos in positions:
					resultInfo.append(getInfo(lines, pos))
			writeInfo(resultInfo,name1[0],num_str)
			name = line.split('>')
			name1 = name[1].split(' ')
			lines = ''
		else:
			lines = lines + line.strip()
	resultInfo = []
	strings = [lines, lines[1:], lines[2:]]
	for index, string in enumerate(strings):
		positions = readByThree(string, index)
		positions = sorted(positions, key = operator.itemgetter(0))
		for pos in positions:
			resultInfo.append(getInfo(lines, pos))
	writeInfo(resultInfo,name1[0],num_str)
	f.close()


def LenTop3(path,opath,top):
	f = open(path,'r')
	nf = open(opath,'a')
	length = []
	for line in f:
		a = line.strip()
		a = a.split('\t')
		if len(a) == 4:
 			a[2] = int(a[2])
			length.append(a)
			continue
		else:
 			if len(length) != 0:
				length.sort(key=operator.itemgetter(2), reverse = True)
 				for i in range(0,min(top,len(length))):
					length[i][2] = str(length[i][2])
					if int(length[i][2])>100:
						nf.write('\t'.join(length[i]) + '\n')
			nf.write(line.strip() + '>\n')
			length = []
	f.close()
	nf.close()

def main(f2, f3):
	#print "Running LGC with %d process(es)" % processes
	print "Input: " + f2
	#print "param:" + f4
	print "Output: " + f3

	#f=open(f4,'r')
	#s = f.read().split('\n')
	s = ['-0.05215	0.10796	-0.071697	0.01759','0.1014	-0.16872	0.0858349	-0.0083821']
	#f.close()
	c = s[0].split('\t')
	n = s[1].split('\t')
	#print c
	#print n

	num_list = ['0','1','2','3','4','5','6','7','8','9']
	random.shuffle(num_list)
	num_str = ''.join(num_list)
	print 'Scan ORF ...'
	ORF(f2,num_str)
	print 'Done'
	oopath = 'ORF'+num_str+'.temp'
	outpath1 = 'top'+num_str+'.temp'
	outpath2 = f3


	LenTop3(oopath,outpath1,3)
	oo = open(outpath2,'w')
	a1 = eval(c[0])
	a2 = eval(c[1])
	a3 = eval(c[2])
	a4 = eval(c[3])
	b1 = eval(n[0])
	b2 = eval(n[1])
	b3 = eval(n[2])
	b4 = eval(n[3])
	g = open('top'+num_str+'.temp','r')
	para = g.read().split('>\n')
	m = len(para)
	Fc = [1]*m
	Fn = [1]*m
	line = ['',-1,-1,-1,-1]*m
	Jcomp  = [-1]*m
	string = [0]*m
	J      = [0]*m
	Cflag  = [0]*m
	oo.write("# Sequence Name: name of transcript\n")
	oo.write("# ORF Length: length of the longest ORF\n")
	oo.write("# GC Content: GC content of the longest ORF\n")
	oo.write("# Coding Potential Score: coding potential score for a transcript, which is protein-coding RNA if greater than 0 or ncRNA if smaller than 0. '0' indicates that mRNA probability equals lncRNA probability. Also, if the ORF length is shorter than 100nt, '0' is output.\n")
	oo.write("# pc: ORF probability for coding sequence\n")
	oo.write("# pnc: ORF probability for non-coding sequence\n")
	oo.write("# fc: Stop-codon probability for coding sequence\n")
	oo.write("# fnc: Stop-codon probability for non-coding sequence\n")
	oo.write("# Coding Label: 'Coding' represents mRNA and 'Non-coding' represents ncRNA\n")
	oo.write("# \n")
	oo.write("# Sequence Name\t")
	oo.write("ORF Length\t")
	oo.write("GC Content\t")
	oo.write("Coding Potential Score\t")
	oo.write("Coding Label\t")
	oo.write("pc\t")
	oo.write("pnc\t")
	oo.write("fc\t")
	oo.write("fnc\t")
	oo.write("\n")
	for i in range(0,m):

		segment = para[i].split('\n')

		if len(segment) == 1:
			J[i] = 1


		number  = segment.pop()
		string[i] = len(segment)
		fc  = [0]*len(segment)
		fn  = [0]*len(segment)
		C   = [0]*len(segment)
		LL_f = 0
		GC_f = 0

		ffc = 1
		ffn = 1


		for j in range(0,len(segment)):
			seg = segment[j].split('\t')
			C[j]= eval(seg[3])

			GC = eval(seg[3])
			GC2 = GC**2
			GC3 = GC**3

			L  = eval(seg[2])
			N = L/3-1

			pc = 3*(a1*GC3+a2*GC2+a3*GC+a4)
			pn = 3*(b1*GC3+b2*GC2+b3*GC+b4)
			fc[j] = (1-pc)**min(N,5000)*pc
			fn[j] = (1-pn)**min(N,5000)*pn
			if j == 0:
				ffc = pc
				ffn = pn
			if LL_f <= L:
				LL_f = L
				GC_f = GC

		for j in range(0,len(segment)):
			Fc[i] = Fc[i]*fc[j]
			Fn[i] = Fn[i]*fn[j]

		S = Fc[i]/Fn[i]

		compare = [x<y for x,y in zip(C,[0.345]*len(segment))]
		if False not in compare:
			Cflag[i] = True

		if S >= 1:
			Jcomp[i] = 'Coding'
		else:
			Jcomp[i] = 'Non-coding'

		if Cflag[i] is True:
			Jcomp[i] = 'Non-coding'

		if J[i] == 1:
			Jcomp[i] = 'Non-coding'

		if S <= 0:
			SS = -10000
		else:
			SS = math.log(S)

		line[i] = [number,LL_f,GC_f,SS,Jcomp[i]]
		if len(number) != 0:
		#version 2.0 April 2018 zz
			#sequence name
			oo.write(line[i][0]+'\t')
			#orf length
			oo.write(str(line[i][1])+'\t')
			
			#GC content
			a_1 = round(line[i][2],3)
			if a_1==0:
				a_1 = format(line[i][2],'.3e')
			oo.write(str(a_1)+'\t')
			
			#coding potential score
			a_2 = round(line[i][3],3)
			if a_2==0:
				a_2 = format(line[i][3],'.3e')
			oo.write(str(a_2)+'\t')
			
			#coding label
			oo.write(str(line[i][4])+'\t')
			
			#fc
			oo.write(str(format(Fc[i],'.3e'))+'\t')
			
			#fnc
			oo.write(str(format(Fn[i],'.3e'))+'\t')
			
			#pc
			a_3 = round(ffc,3)
			if a_3==0:
				a_3 = format(ffc,'.3e')
			oo.write(str(a_3)+'\t')
			
			#pnc
			a_4 = round(ffn,3)
			if a_4==0:
				a_4 = format(ffn,'.3e')
			oo.write(str(a_4)+'\t')
			
			oo.write('\n')
	Chow = Cflag.count(True)
	oo.close()
	os.remove('ORF'+num_str+'.temp')
	os.remove('top'+num_str+'.temp')


if __name__ == '__main__':
	from optparse import OptionParser
	usage = "LGC version 1.0\nusage: lgc.py [options] INPUT OUTPUT\nPlease cite: \"LGC: Characterization and Identification of Long Non-coding RNAs Based on Feature Relationship\""
	parser = OptionParser(usage=usage)

	parser.add_option("-p", "--processes", dest="processes", type="int",
	                  help="Processes to run LGC", metavar="PROCESSES", default=1)

	(options, args) = parser.parse_args()
	#print len(args)
	if len(args) != 2:
		print usage
		exit()
	input = args[0]
	#print input
	#path = args[0]
	#print path
	output = args[1]
	#print output
	import time
	start_time = time.time()
	#print start_time
	main(input, output)
	print "Computation time %s senconds" % (time.time() - start_time)
