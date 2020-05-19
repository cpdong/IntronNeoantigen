# ask input a file 
# ask input a gtf file
# ask input a thread numbers
#
#dir="/N/u/cpdong/Carbonate/Desktop/test"
#gtfFile="gencode.v32.annotation.gtf"
#
import argparse;
import os,shutil,collections,copy,statistics;
from multiprocessing import Pool;
import pandas as pd;
import numpy as np;
#os.chdir(dir)

parser = argparse.ArgumentParser()
parser.add_argument('-gtf', metavar = 'input', dest='gtfFile', help='Give the full path of gtf file');
parser.add_argument('-out', metavar = 'output', dest='output',help='Give output dir for gtf parsing file');
parser.add_argument('-t', dest='thread', type=int, help='number of multiple thread number,>0');

args = parser.parse_args();
gtfFile = args.gtfFile;
dir = args.output;
thread = args.thread;

path, filename = os.path.split(inputFile)
basename, ext = os.path.splitext(filename)

if thread:
    thread=min(os.cpu_count(), thread);
else:
    thread=8; # defualt thread number setting as 8

if dir[-1]=='/':
    dir =dir;
else:
    dir =dir + '/';
if not os.path.exists(dir + '/temp'):
    os.makedirs(dir + '/temp')# create a folder for speeding calculation
#
#https://stackoverflow.com/questions/11968998/remove-lines-that-contain-certain-string
#Filter only protein coding gene for the next step
chr=[];feature=[];start=[];end=[];strand=[];gID=[];txID=[];exonNum=[];
with open(gtfFile, 'r') as gfile:
	for line in gfile:
		if not 'protein_coding' in line:
			continue;
		desc = line.strip().split('\t')[8].split(';');
		#
		gID_temp='NA'
		txID_temp="NA";
		exonNum_temp="NA"
		for attribute in desc: ## for each element of description
			if len(attribute.strip())<2 or len(attribute.strip().split(' '))<2:
				continue; ## probably the last description
			attrName = attribute.strip().split(' ')[0];
			attrVal = attribute.strip().split(' ')[1];
			if attrName.upper() == 'GENE_ID': ## it is a description for gene_id
				gID_temp = attrVal.replace('"', '');
			elif attrName.upper() == 'TRANSCRIPT_ID': ## it is a description for transcript_id
				txID_temp = attrVal.replace('"', '');
			elif attrName.upper() == 'TRANSCRIPT_TYPE': ## it is a description for transcript_id
				txType_temp = attrVal.replace('"', '');
			elif attrName.upper() == 'EXON_NUMBER': ## it is a description for transcript_id
				exonNum_temp = attrVal.replace('"', '');
		if txID_temp =="NA":
			continue; # remove only fully gene row
		if exonNum_temp =="NA":
			continue; # remove only fully transcripts row
		if txType_temp !="protein_coding":
			continue; # remove tx that lncRNA or intron retain or others		
		chr.append(line.strip().split('\t')[0]);
		feature.append (line.strip().split('\t')[2].lower()); ## exon, intron, CDS, start_codon, stop_codon..
		start.append(line.strip().split('\t')[3]); ## start coord, 1-base
		end.append(line.strip().split('\t')[4]); ## end coord, 1-base
		strand.append(line.strip().split('\t')[6]);
		gID.append(gID_temp)
		txID.append(txID_temp)
		exonNum.append(exonNum_temp)
		#
#
df=pd.DataFrame([chr,feature,start,end,strand,gID,txID,exonNum]).T
df.columns =['chr', 'feature', 'start', 'end', 'strand', 'geneid', 'txid', 'exonNum']
#
#remove gFile and other variable to release memory
del gfile,chr,feature,start,end,strand,gID,txID,exonNum;
#uniqGene=df['geneid'].unique()
#
# =============================================================================
# following single gene run function
# =============================================================================
#
def single_gene_run(gene):
	temp= df.loc[(df['geneid'] == gene) & (df['feature'].str.upper() =='CDS')];
	if temp.shape[0]>1: # skip gene with one or less cds
		# stp1. remove duplicated exons
		temp1 = temp.drop_duplicates(subset=['start', 'end'], keep='first');
		# stp2. remove subsetting exon & stp3. get union exon sets
		temp2 = temp1.iloc[0:0]
		a = temp1[['start', 'end']].values.tolist()
		a = [[int(x) for x in list] for list in a]
		b = []
		for begin,end in sorted(a):
			if b and b[-1][1] >= begin - 1:
				b[-1] = [b[-1][0], end]
			else:
				b.append([begin, end])
		for i in range(0, len(b)):# add the new union set
			addRow=temp1.iloc[0];
			addRow['start']=b[i][0]; # update new start
			addRow['end']=b[i][1]; # update new end
			addRow['txid']='cds_new';
			addRow['exonNum']='NA'
			temp2 = temp2.append(addRow, ignore_index = True);
		temp2['exon_endsets']= 'NA'
		temp2['exon_txIds']= 'NA'
		temp2['exonEndLoc']= 'NA'
		temp2['exonStart']= 'NA'
		#
		#temp2['start'] = temp2['start'].astype('int64');
		#temp2=temp2.sort_values('start');
		#exon_new = temp2[['start','end']].values.tolist();
		#check with each transcripts exon end whether match new join exon
		uniqTx=temp['txid'].unique()
		tx_dict = {}
		for tx in uniqTx:
			txTemp=temp[temp['txid'] ==tx];
			txTemp = txTemp.drop_duplicates(subset=['start', 'end'], keep='first');
			#identify the first exon based on cds start points
			if txTemp.shape[0]<=0:
				continue;# for the only one exon transcripts, skip, no intron retain petential
			txTemp['cdsLen']=txTemp['end'].astype(int) - txTemp['start'].astype(int) +1;# calculate the cds-exon length
			#txTemp['exonEndLoc']='NA';#end pharse of 0,1,2 of exon
			tx_dict.update( {tx : txTemp['cdsLen'].sum()} )
			for j in range(0, len(txTemp.index)): #have to run for loop to get reading frame
				#we don't consider the intron before the start_codon and that after the stop codon
				txTemp_1=txTemp[txTemp['exonNum'].astype(int) <= int(txTemp.iloc[j][7])];
				trunc_len=txTemp_1['cdsLen'].sum();
				#txTemp.loc[index,'exonEndLoc']= trunc_len%3
				a_range = range(int(txTemp.iloc[j]['start']), int(txTemp.iloc[j]['end']) +1);
				for i in range(0, len(temp2.index)):
					b_range = range(int(temp2.iloc[i]['start']), int(temp2.iloc[i]['end']) +1);
					overlap_list = [x for x in a_range if x in b_range];#where exon overlap new exon
					if len(overlap_list)>0:#make note of every exon end.
						if temp2.iloc[0]['strand'] =='+':
							exon_extension= int(temp2.iloc[i]['end']) - int(txTemp.iloc[j]['end']);
							temp2.at[i,'exon_endsets']= temp2.at[i,'exon_endsets'] +','+ str((trunc_len%3 + exon_extension)%3);
						else:
							exon_extension= int(txTemp.iloc[j]['start']) - int(temp2.iloc[i]['start']);
							temp2.at[i,'exon_endsets']= temp2.at[i,'exon_endsets'] +','+ str((trunc_len%3 + exon_extension)%3);
						temp2.at[i,'exon_txIds'] = temp2.at[i,'exon_txIds'] + ',' + tx;
		temp2['exon_endsets']= [v.replace('NA,', '') for v in temp2['exon_endsets']]
		temp2['exon_txIds']= [w.replace('NA,', '') for w in temp2['exon_txIds']] 
		# Figure the dominate/primary exon end pharse;
		for k in range(0, len(temp2.index)):# update reasonable reading frame of each exon
			#k=11;
			x1 = temp2.loc[k,'exon_endsets'];# exon end
			x1 = [int(x) for x in x1.split(',')];
			x2 = temp2.loc[k,'exon_txIds'];x2 = x2.split(',');# tx use this exon
			counter=collections.Counter(x1);
			if len(counter) <=1:# if only one tx use this exon
				temp2.loc[k,'exonEndLoc'] = x1[0];
			else:
				cc= counter.most_common(2); # count top 2 most common reading frame
				#common2freq = [x[0] for x in cc]
				if cc[0][1] != cc[1][1]: # use the most common used reading frame
					temp2.loc[k,'exonEndLoc'] = counter.most_common(1)[0][0];
				else:# a little complex, two comon used frame, we choose whose average tx length longer
					top_1 = int(cc[0][0]);top_2 = int(cc[1][0]);#cc[i][0] the endExonLoc
					x1 = np.array(x1);
					x2 = np.array(x2);
					txId_1 = np.unique(x2[x1 == top_1])
					txId_2 = np.unique(x2[x1 == top_2])
					txId_1_len = [tx_dict[x] for x in txId_1]
					txId_2_len = [tx_dict[x] for x in txId_2]
					if statistics.mean(txId_1_len) > statistics.mean(txId_2_len):
						temp2.loc[k,'exonEndLoc'] = cc[0][0];
					else:
						temp2.loc[k,'exonEndLoc'] = cc[1][0];
			if temp2.iloc[0]['strand'] =='+':#figure the first complete codon of exon
				exon_psedu_end = int(temp2.iloc[k]['end']) - int(temp2.iloc[k]['exonEndLoc']);
				exon_codon_num = int((exon_psedu_end - int(temp2.iloc[k]['start']) +1)/3);
				temp2.loc[k,'exonStart'] = exon_psedu_end - exon_codon_num*3 + 1;
			else:
				exon_psedu_end = int(temp2.iloc[k]['start']) + int(temp2.iloc[k]['exonEndLoc']);
				exon_codon_num = int(int(temp2.iloc[k]['end'] - exon_psedu_end +1)/3);
				temp2.loc[k,'exonStart'] = exon_psedu_end + exon_codon_num*3 - 1;
		temp2['start'] = temp2['start'].astype(int);
		temp2['exonEndLoc'].replace(0, 3,inplace=True); # replace 0to3 for easy calculate
		temp2['region'] = 'NA';
		#
		temp2['pep_8mer_start'] = 'NA';
		temp2['pep_8mer_end'] = 'NA';
		temp2['pep_9mer_start'] = 'NA';
		temp2['pep_9mer_end'] = 'NA';
		temp2['pep_10mer_start'] = 'NA';
		temp2['pep_10mer_end'] = 'NA';
		temp2['pep_11mer_start'] = 'NA';
		temp2['pep_11mer_end'] = 'NA';
		temp2['pep_12mer_start'] = 'NA';
		temp2['pep_12mer_end'] = 'NA';
		temp2['pep_13mer_start'] = 'NA';
		temp2['pep_13mer_end'] = 'NA';
		temp2['pep_14mer_start'] = 'NA';
		temp2['pep_14mer_end'] = 'NA';
		# update all information
		if temp2.shape[0]>1:
			#update intron information
			# update peptides information,pre-enlong 7 aa in each end
			temp3 = copy.deepcopy(temp2)
			temp3['start'] = temp3['start'].astype(int)
			temp3=temp3.sort_values('start')
			temp4=copy.deepcopy(temp2)
			temp4=temp4[:-1]
			temp4[['start','end','exonNum','exon_endsets','exon_txIds','exonEndLoc','exonStart']]='NA'
			temp4[['feature']]='intron';
			temp4[['txid']]='intron_region';
			for m in range(temp4.shape[0]):#update the information of each intron region
				temp4.iloc[m]['start'] = int(temp3.iloc[m]['end']) + 1;
				temp4.iloc[m]['end'] = int(temp3.iloc[m+1]['start']) - 1;
				intron_length = int(temp3.iloc[m+1]['start']) - int(temp3.iloc[m]['end']) -1;
				if temp3.iloc[0]['strand']=='+':
					pseudo_end = int(temp3.iloc[m]['end']) - int(temp3.iloc[m]['exonEndLoc']);# here we consider retract 0-based gtf data with pybedtools already!
					temp4.iloc[m]['pep_8mer_start'] = max(pseudo_end -6*3, int(temp3.iloc[m]['exonStart'])-1);
					temp4.iloc[m]['pep_9mer_start'] = max(pseudo_end -7*3, int(temp3.iloc[m]['exonStart'])-1);
					temp4.iloc[m]['pep_10mer_start'] = max(pseudo_end -8*3, int(temp3.iloc[m]['exonStart'])-1);
					temp4.iloc[m]['pep_11mer_start'] = max(pseudo_end -9*3, int(temp3.iloc[m]['exonStart'])-1);
					temp4.iloc[m]['pep_12mer_start'] = max(pseudo_end -10*3, int(temp3.iloc[m]['exonStart'])-1);
					temp4.iloc[m]['pep_13mer_start'] = max(pseudo_end -11*3, int(temp3.iloc[m]['exonStart'])-1);
					temp4.iloc[m]['pep_14mer_start'] = max(pseudo_end -12*3, int(temp3.iloc[m]['exonStart'])-1);
					#
					intron_EndLoc = (int(temp3.iloc[m]['exonEndLoc']) + intron_length)%3;
					temp4.iloc[m]['pep_8mer_end'] = int(temp3.iloc[m+1]['start'])- 1 - intron_EndLoc + 7*3;
					temp4.iloc[m]['pep_9mer_end'] = int(temp3.iloc[m+1]['start'])- 1 - intron_EndLoc + 8*3;
					temp4.iloc[m]['pep_10mer_end'] = int(temp3.iloc[m+1]['start'])- 1 - intron_EndLoc + 9*3;
					temp4.iloc[m]['pep_11mer_end'] = int(temp3.iloc[m+1]['start'])- 1 - intron_EndLoc + 10*3;
					temp4.iloc[m]['pep_12mer_end'] = int(temp3.iloc[m+1]['start'])- 1 - intron_EndLoc + 11*3;
					temp4.iloc[m]['pep_13mer_end'] = int(temp3.iloc[m+1]['start'])- 1 - intron_EndLoc + 12*3;
					temp4.iloc[m]['pep_14mer_end'] = int(temp3.iloc[m+1]['start'])- 1 - intron_EndLoc + 13*3;
				else: # here we consider retract 0-based gtf data with pybedtools already!
					pseudo_start= int(temp3.iloc[m+1]['start']) + int(temp3.iloc[m+1]['exonEndLoc']) -1;
					temp4.iloc[m]['pep_8mer_end'] = min(pseudo_start + 6*3, int(temp3.iloc[m+1]['exonStart']));
					temp4.iloc[m]['pep_9mer_end'] = min(pseudo_start + 7*3, int(temp3.iloc[m+1]['exonStart']));
					temp4.iloc[m]['pep_10mer_end'] = min(pseudo_start + 8*3, int(temp3.iloc[m+1]['exonStart']));
					temp4.iloc[m]['pep_11mer_end'] = min(pseudo_start + 9*3, int(temp3.iloc[m+1]['exonStart']));
					temp4.iloc[m]['pep_12mer_end'] = min(pseudo_start + 10*3, int(temp3.iloc[m+1]['exonStart']));
					temp4.iloc[m]['pep_13mer_end'] = min(pseudo_start + 11*3, int(temp3.iloc[m+1]['exonStart']));
					temp4.iloc[m]['pep_14mer_end'] = min(pseudo_start + 12*3, int(temp3.iloc[m+1]['exonStart']));
					#
					intron_EndLoc = (int(temp3.iloc[m+1]['exonEndLoc']) + intron_length)%3;
					temp4.iloc[m]['pep_8mer_start'] = int(temp3.iloc[m]['end']) + intron_EndLoc - 7*3;
					temp4.iloc[m]['pep_9mer_start'] = int(temp3.iloc[m]['end']) + intron_EndLoc - 8*3;
					temp4.iloc[m]['pep_10mer_start'] = int(temp3.iloc[m]['end']) + intron_EndLoc - 9*3;
					temp4.iloc[m]['pep_11mer_start'] = int(temp3.iloc[m]['end']) + intron_EndLoc - 10*3;
					temp4.iloc[m]['pep_12mer_start'] = int(temp3.iloc[m]['end']) + intron_EndLoc - 11*3;
					temp4.iloc[m]['pep_13mer_start'] = int(temp3.iloc[m]['end']) + intron_EndLoc - 12*3;
					temp4.iloc[m]['pep_14mer_start'] = int(temp3.iloc[m]['end']) + intron_EndLoc - 13*3;
			temp3=temp3.append(temp4);
			temp3=temp3.sort_values('start');
			if temp3.iloc[0]['strand']=='+':
				temp3['region']=range(1,temp3.shape[0]+1);
			else:
				temp3['region']=range(temp3.shape[0], 0, -1);
			temp3['exonEndLoc'].replace(3, 0,inplace=True); # replace 3 to 0 back
			temp3.to_csv(dir + '/temp/' + gene + '_exon_intron_coordinate_final.txt',  index=False, sep='\t')
	#remove gFile and other variable to release memory
	#del a,b,temp,temp2,temp3,temp5 geneTemp,txTemp,txTemp_new; !can not delete inside pool operation
	#print the process percentage of total genes
	if (uniqGene.index(gene) + 1)%1000 ==0:
		print(uniqGene.index(gene)+1) #print(gene) order
	#return geneTemp.values.tolist()
	#return geneTemp.values.tolist() # return a list from processing
#
#
uniqGene=df['geneid'].unique().tolist()
#
if __name__=='__main__':
	pool=Pool(processes=thread); # defualt thread number setting as 8
	pool.map(single_gene_run, uniqGene)
	print('Waiting for all subprocesses done...')
	pool.close()
	pool.join()
	print('All subprocesses done.')
#
#
#combine all the split file with cat function
filenames=os.listdir(dir + '/temp');
with open(dir + '/' + basename + '.re-annotation.gtf', 'w') as outfile:
	for fname in filenames:
		with open(dir + '/temp/'+ fname) as infile:
			next(infile) #skip the first header row
			for line in infile:
				outfile.write(line)
#
shutil.rmtree(dir + '/temp') # remnove the temp folder
#
#
