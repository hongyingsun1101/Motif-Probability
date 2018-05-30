#!/usr/bin/python

'''
this code runs rnamotif and gets the motif information
then run ProbStemloop in RNA structure to get the probability for the exterior loops.
the codes takes the genome as an input file, first split the whole genome into each gene file
then process each gene file
if the gene file contains desired motif, convert into ct file
run efn2
convert free energy change to equilibrium constant
multiply equillibrium constant by the output of ProbStemloop in RNAstructure.
sort the lists
done!!!
'''
import os
import sys
import subprocess
from subprocess import Popen, PIPE
import string
import getopt
import glob
import itertools
import re
import ntpath
import math
from itertools import chain
#from glob import glob

new='/home/hsun/rnamotif-3.1.1/hs/new'
cvsReport='/home/hsun/rnamotif-3.1.1/hs/newsummary'
summaryfiles='/home/hsun/rnamotif-3.1.1/hs/newsummary/summaryfiles'
#set up working direcoty
os.getcwd()
#print ('printing the current directory which contains the running python program now:')
#print (os.getcwd()); # Prints the working directory
os.chdir(new) # set up the working directory.
#print ('printing working directory now:')
#print (os.getcwd())



def pretitle(bedFile):
	genomeName=ntpath.basename(bedFile).split('.')
	print('now printing the genomeName')
	print(genomeName)
	pretitle=genomeName[0]+genomeName[2]
	return pretitle

def splitBed(bedFile):
	genomeName=ntpath.basename(bedFile).split('.')
	print('now printing the genomeName')
	print(genomeName)
	pretitle=genomeName[0]+genomeName[2]
	with open(bedFile, 'r') as bed:
		lines=bed.readlines()
		for line in lines:
			if line.startswith('>'):
				header=pretitle+'_'+line.strip()[1:]+'_fasta'
				splitFasta=open(str(header),'w')
				splitFasta.writelines(line.strip()+"\n")
			else:
				splitFasta.writelines(line.strip()+"\n")


def GetFileList(directory_path): #return a list of the files in the given directoty
	files = [f for f in os.listdir(directory_path) if os.path.isfile(os.path.join(directory_path, f))]
	return files

def countLength(f):
	with open(f,'r') as fasta:
		lines=fasta.readlines()
		for line in lines:
			#print('now printing the line of fasta file')
			#print(line+'\n')
			if line.startswith( '>' ):
				pass
			else:
				fastaLen=len(line.strip())
	return fastaLen
def nonEmpymotif(direc):
	for f in glob.glob(os.path.join(direc, '*_fasta')):
		print('now printing f')
		print (f)
		title=ntpath.basename(f).rstrip('_fasta')
		outf=open(str(title)+'_motif','w')
		process=subprocess.Popen(["rnamotif", "-descr", "/home/hsun/rnamotif-3.1.1/hs/hairpin.descr", f], stdout=outf) #takes a descriptor and a sequence file, returns the positions of the motif
		process.communicate()    # execute it, the output goes to the stdout
		exit_code = process.wait()    # when finished, get the exit code
	for motifFile in glob.glob(os.path.join(direc, '*_motif')):
		if (os.stat(motifFile).st_size == 0):
			#nopncontainMotifFasta=motifFile.rstrip('_motif')+'_fasta'
			#process=subprocess.Popen("%s %s %s" % ("rm", motifFile,noncontainMotifFasta), shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE, close_fds=True)
			process=subprocess.Popen(["rm", motifFile])
			process.communicate()    # execute it, the output goes to the stdout
			exit_code = process.wait()    # when finished, get the exit code			

def splitMotif(new):
	for filename in glob.glob(os.path.join(new, '*_motif')):
		title=ntpath.basename(filename).rstrip('_motif').split('_')[0]+'_'
		#print ('now printing tile')
		#print (title)	
		#outf=open(str(title)+'_msplit', 'w')
		fi=open(filename,'r')
		lines=fi.readlines()  #read all lines into list
		for line in lines:
			if line.startswith('N'):
				#headerold=line.strip(" ").replace(" ", "-").rstrip("\n")+ ".new"
				headerold=title+line.strip(" ").rstrip("\n")+ "_msplit"
				#header=' '.join(headerold.split())
				header=re.sub( '\s+', ' ', headerold ).strip().replace(" ", "-")
				#print('now priting header for msplit files')
				#print (header)

				fo=open(str(header), 'w')
				fo.write(lines[0])
				fo.write(lines[1])
				fo.write(lines[2])
				fo.write('>'+str(header).split()[0]+'\n')
				#print count,num_line+count,num_line
				fo.writelines(line.strip()+'\n')
				fo.close()				
def removeFasta(new):
	for fastaFile in glob.glob(os.path.join(new, '*_fasta')):
		containMotif=fastaFile.rstrip('_fasta')+'_motif'
		if not os.path.exists(containMotif):
			process=subprocess.Popen(["rm", fastaFile])
			process.communicate()    # execute it, the output goes to the stdout
			exit_code = process.wait()    # when finished, get the exit code	
			#process=subprocess.Popen("%s %s %s" % ("rm", motifFile,noncontainMotifFasta), shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE, close_fds=True)

def trancateFasta(new):
	for msplitFile in glob.glob(os.path.join(new, '*_msplit')):
		title=ntpath.basename(msplitFile).rstrip('_msplit').split('-')
		print("Now printing title list")
		print (title)
		geneName=title[0]
		freeEnergy = title[1]
		strandInformation=title[2]
		startPostion=title[3]
		motifLength=title[4]
		descr_h5=title[5]
		descr_ss=title[6]
		descr_h3=title[7]
		for fastaFile in glob.glob(os.path.join(new, '*_fasta')):

			fasta=ntpath.basename(fastaFile)
			geneFasta=fasta.rstrip('_fasta')
			seqLen=countLength(fastaFile)
			print('now priting geneName: '+ geneName)
			print('now priting geneFasta: '+ geneFasta)
			print('now printing whether the pretitlename are the same for fasta and msplit files')
			print(geneName==geneFasta)
			trueCount=0
			if geneFasta==geneName:
				trueCount = trueCount+1
				print('now printing trueCount: '+ str(trueCount))
				
				trancatedTitle=ntpath.basename(msplitFile).rstrip('_msplit')+'_notDone'
				
				with open(fastaFile,'r') as oriFasta:
					lines=oriFasta.readlines()
					fo=open(str(trancatedTitle),'w')
					for line in lines:
						#print('now printing lines in the fasta file')
						#print(line)
						if line.startswith('>'):
							fo.writelines(line.strip()+'\n')
						else:
							if int(startPostion) <=400:
								newStartPostion=startPostion
								if seqLen-(int(startPostion)+int(motifLength)) <=400:
									newline1=line[0:seqLen]
									fo.writelines(newline1.strip()+'\n')
								else:
									newline1=line[0:(int(startPostion)+int(motifLength)+400)]
									fo.writelines(newline1.strip()+'\n')		
							else:
								newStartPostion=400
								if seqLen-(int(startPostion)+int(motifLength)) <=400:
									newline1=line[int(startPostion)-400:seqLen]
									fo.writelines(newline1.strip()+'\n')
								else:
									newline1=line[int(startPostion)-400:(int(startPostion)+int(motifLength)+400)]
									fo.writelines(newline1.strip()+'\n')

					print('now printing trueCount: '+ str(trueCount))
					newTitle=geneName+'-'+freeEnergy+'-'+str(strandInformation)+'-'+str(startPostion)+'-'+motifLength+'-'+descr_h5+'-'+descr_ss+'-'+descr_h3+'-'+str(newStartPostion)+'_newsplit'
					#newTitle=geneName+'-'+freeEnergy+'-'+str(strandInformation)+'-'+str(newStartPostion)+'-'+motifLength+'-'+descr_h5+'-'+descr_ss+'-'+descr_h3+'_msplit'				
					print('now printing new msplit titles')
					print(newTitle)
				process=subprocess.Popen("%s %s %s" % ("cp", msplitFile, newTitle), shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE, close_fds=True)
				#process=subprocess.Popen(["rm", motifFile])
				process.communicate()    # execute it, the output goes to the stdout
				exit_code = process.wait()    # when finished, get the e																						



			
			



def addOldPosition():
	for msplitFile in glob.glob(os.path.join(new, '*_newsplit')):
		title=ntpath.basename(msplitFile).rstrip('_newsplit').split('-')
		#print("Now printing newsplit title list")
		#print (title)
		geneName=title[0]
		freeEnergy = title[1]
		strandInformation=title[2]
		startPostion=title[3]
		motifLength=title[4]
		descr_h5=title[5]
		descr_ss=title[6]
		descr_h3=title[7]
		genePosition=title[8]
		com1=geneName+strandInformation+startPostion+motifLength
		for fastaFile in glob.glob(os.path.join(new, '*_notDone')):
			fasta=ntpath.basename(fastaFile).rstrip('_notDone').split('-')
			#print('now printing not done fasta title list\n')
			#print(fasta)
			notDone_geneName=title[0]
			notDone_freeEnergy = title[1]
			notDone_strandInformation=title[2]
			notDone_startPostion=title[3]
			notDone_motifLength=title[4]
			notDone_descr_h5=title[5]
			notDone_descr_ss=title[6]
			notDone_descr_h3=title[7]
			com2=notDone_geneName+notDone_strandInformation+notDone_startPostion+notDone_motifLength
			print('printing com1: '+com1)
			print('printting com2: ' + com2)
			print('print whether com1 and com2 are the same gene')
			print(com1==com2)
			if com1==com2:
				upload=ntpath.basename(msplitFile).rstrip('_newsplit')+'.upload'
		
				process=subprocess.Popen("%s %s %s" % ("cp",fastaFile, upload), shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE, close_fds=True)
				#process=subprocess.Popen(["rm", motifFile])
				process.communicate()    # execute it, the output goes to the stdout
				exit_code = process.wait()    # when finished, get the e	


def toUpper(fileDir):
	for x in glob.glob(os.path.join(fileDir, '*upload')):
		inputFile = open(x, "r")
		content=inputFile.read()
		
		with open(x+'.upper',"w") as outputFile:
			outputFile.write(content.upper())																					

def getMotifCt():
	for filename in glob.glob(os.path.join(new, '*_newsplit')):
		title=ntpath.basename(filename).rstrip('_newsplit')
		with open(filename, 'r') as myfile:
			filelines = myfile.readlines()
			outf=open(str(title+'_ct'),'w')
			process=subprocess.Popen(["rm2ct", filename], stdout=outf)
			process.communicate()    # execute it, the output goes to the stdout
			exit_code = process.wait()    # when finished, get the exit code
	

def efn2onOut():
	#fastaPath = '/home/hsun/rnamotif-3.1.1/hs/human19utr5'
	#print (fastaPath)
	for filename in glob.glob(os.path.join(new, '*_ct')):
		header=ntpath.basename(filename).rstrip('_ct')
		#print ("Now print filename")
		#print (filename)

		fi=open(filename,'r')
		lines=fi.readlines()  #read all lines into list
		#print ("Now printing lines")
		#print (lines)
		for line in lines:
			if (re.search("utr", line) ):
				#header=line.strip(" ").replace(" ", "-").rstrip("\n")+ ".dg" 
				#print ("Now printing header to store dg" + "\n"+ header)
			
				#os.path.join(dir_name, base_filename + "." + filename_suffix)
				#fo=open(str(header),'w')  # string formatting header?
				#print ("Now printing fo" + "\t" + str(header))
			
				process=subprocess.Popen(["efn2", filename, str(header)+'_dg'])
				process.communicate()    # execute it, the output goes to the stdout
				exit_code = process.wait()    # when finished, get the exit code
#modify the code to include the motif sequence information in the output csv file. Instead of splitting ct, split motif files in the previous functions.
def getPfs(uploadDir):
	for filename in glob.glob(os.path.join(uploadDir, '*upper')):
		process=subprocess.Popen(["partition-smp",filename, filename+'.pfs'])
		process.communicate()    # execute it, the output goes to the stdout
		exit_code = process.wait()    # when finished, get the exit code	



def outProb(pfsDir):
	for filename in glob.glob(os.path.join(pfsDir, '*pfs')):
		title=ntpath.basename(filename).rstrip('.upload.upper.pfs')
		titleList = title.split("-") #split title into a list of strings and get information from the list
		#print("Now printing title list")
		#print (titleList)
		geneName = titleList[0]
		strandInformation=titleList[2]
		startPostion=titleList[3]
		motifLength=titleList[4]
		descr_h5=titleList[5]
		descr_ss=titleList[6]
		descr_h3=titleList[7]
		trancatedStartPostion=titleList[8]
		print('now printing file name')
		print(filename)
		process=subprocess.Popen(["ProbStemloop", '-pfs',filename, trancatedStartPostion, motifLength])
		process.communicate()    # execute it, the output goes to the stdout
		exit_code = process.wait()    # when finished, get the exit code

def writedG():  
	#fastaPath = '/home/hsun/rnamotif-3.1.1/hs/human19utr5'
	for filename in glob.glob(os.path.join(new, '*_dg')):
		#print ("Now print filename")
		#print (filename)		
		title=ntpath.basename(filename).rstrip('_dg')
		#print("Now printing file title")
		#print (title)
		#print (type(title))
		titleList = title.split("-") #split title into a list of strings and get information from the list
		#print("Now printing title list")
		#print (titleList)
		geneName = titleList[0]
		strandInformation=titleList[2]
		startPostion=titleList[3]
		motifLength=titleList[4]
		descr_h5=titleList[5]
		descr_ss=titleList[6]
		descr_h3=titleList[7]
		trancatedStartPostion=titleList[8]
		#print('now printing trancatedStartPostion: ')
		#print(trancatedStartPostion)

		#print(geneName+"\n" +strandInformation+"\n"+startPostion+"\n"+motifLength)
		#pfsName=os.path.join(fastaPath, geneName)+".pfs"
		pfsName=title+'.upload.upper.pfs'
		#print ("Now printing genePath")
		pfsPath=os.path.join(new, pfsName)
		#print('now printing pfsPath')
		#print (pfsPath)

		with open(filename, 'r') as myfile:
			data=myfile.read().split()
			#print("Now printing data")
			#print (data)
			dGmotif=float(data[-1]) # this works try whether it still works without using float.
			#dGmotif=data[-1]
			#print("Now printing dGMotif" )
			#print(dGmotif)
			ProbInner=math.exp((-dGmotif)/(310.15*0.00198588))
			#print("Now printing ProbInner" + "\n" )
			#print(ProbInner)
		#print (not os.path.exists(pfsPath))
		if not os.path.exists(pfsPath):
			pass
		else:
			#now pass the fastafile, startPostion and motifLength to ProbStemloop in the RNAstructure to get exterior probablity
			with open(str(title+".txt"),'w') as prob:
				process=subprocess.Popen(["ProbStemloop", '-pfs',pfsPath, trancatedStartPostion, motifLength], stdout=prob)
				process=subprocess.Popen(["ProbStemloop", '-pfs',pfsPath, trancatedStartPostion, motifLength])

				#process=subprocess.Popen(["ProbStemloop", "/rnamotif-3.1.1/hs/human19utr5/NM_006015_utr5", int(1270), int(15)])
				process.communicate()    # execute it, the output goes to the stdout
				exit_code = process.wait()    # when finished, get the exit code
			with open(str(title+".txt"), "r") as exteriorFile:
				lines=exteriorFile.readlines()  
				print('now printing exteriorFiletile')
				print(title+'.txt')
				print('now printing exterior file contents')
				print(lines)
				
				probexteriorline=lines[2]
				probe=probexteriorline.split()
				probExterior=float(probe[0])
				print("Now print ProbExterior")
				print (probExterior)
				probMotif=ProbInner*probExterior
				#probMotif=("{0:.3f}".format(ProbInner*probExterior))

				#float("{0:.2f}".format(13.949999999999999))
				print("Now printing probMotif")
				print (probMotif)
			#write results into summary file.
			with open(str(title+".summary"),'w') as report:
				header = "geneName\tstrandInformation\tstartPostion\tmotifLength\tdescr_h5\tdescr_ss\tdescr_h3\ttrancatedStartPosition\n"
				report.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n".format(geneName,strandInformation,startPostion,motifLength,probMotif,descr_h5,descr_ss,descr_h3,trancatedStartPostion))

#define a function to combine all summary files.
#chang the result file name

#new='/home/hsun/rnamotif-3.1.1/hs/new'
#cvsReport='/home/hsun/rnamotif-3.1.1/hs/newsummary'
#summaryfiles='/home/hsun/rnamotif-3.1.1/hs/newsummary/summaryfiles'
def combineAllReport(genome):
	csvFile=str(cvsReport)+'/'+genome+'.csv'
	header = "geneName\tstrandInformation\tstartPostion\tmotifLength\tmotifProb\tdescr_h5\tdescr_ss\tdescr_h3\ttrancatedStartPosition\n"
	#print('now printing header: ')
	outfile=open(str(csvFile), "w")
	outfile.write(str(header))
	read_files=glob.glob(os.path.join(new, '*.summary'))
	#print('now printing the number of files in the read_files')
	#print(len(read_files))
	#print('now printing read_files: ')
	#print(read_files)
	for file in read_files:
		filetilte=ntpath.basename(file)
		#print('now printing filetile: ')
		#print(filetilte)
		#print ('now printing whether this is the desired genome:')
		#print(filetilte.startswith(genome))
		infile0=open(file, "r")
		infile=infile0.readlines()

		#print('now please printing the input file:')
		#print(infile)
		if (filetilte.startswith(genome)):
			for line in infile:
				#print('now printing lines')
				#print(line)
				outfile.writelines(line)
		infile0.close()
	outfile.close()




def main():
	pretitle(sys.argv[1])
	splitBed(sys.argv[1])
	nonEmpymotif(new)
	splitMotif(new)
	removeFasta(new)
	trancateFasta(new)
	addOldPosition()
	toUpper(new)
	getPfs(new) #only run this on local machine for testing, run partition function on bluehive.
	getMotifCt()
	efn2onOut()
	outProb(new)
	writedG()
	combineAllReport('hg19utr3')



if __name__ == '__main__':
	main()

