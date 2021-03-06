import argparse,os,re,json;
import gzip,sys,shutil;
import pandas as pd;
import numpy as np;
import time,subprocess;
import pybedtools;
from Bio import SeqIO;
from multiprocessing import Pool;
from os.path import dirname, abspath
dir = dirname(dirname(abspath(__file__)));# get parent dir of the scripts
#
#print(time.strftime("%Y-%m-%d %H:%M:%S"))
# =============================================================================
# initial setting for future update
# =============================================================================
#
parser = argparse.ArgumentParser()
parser.add_argument('-genotype', metavar = 'inputfile', dest='genotype', help='Give arcasHLA calling result fullname');
parser.add_argument('-hla', metavar = 'inputstring', dest='inputHLA', help='Give arcasHLA calling result fullname');
parser.add_argument('-file', metavar = 'input', dest='inputFile', help='Give the intron calling result files');
parser.add_argument('-fasta', metavar = 'fasta', dest='fastaFile', help='Give the fasta file');
parser.add_argument('-gtfindex', metavar = 'index', dest='gtfIndex', help='Give the index file');
parser.add_argument('-len', metavar = 'pep length', type=str, dest='pepLen', help='Give the peptide files');
parser.add_argument('-thread', dest='thread', type=int, help='number of multiple thread number,>0');
parser.add_argument('-refprotein', dest='refprotein', help='ref protein setting');
parser.add_argument('-outdir', metavar = 'output', dest='outputFile',help='Give output file fullname');
#parser.add_argument('-allelenames', metavar = 'allelenamesCheck', dest='allelenames', help='checking whether the input alleles belong to IEDB');
#parser.add_argument('-hla_pseudo', metavar = 'input_hla_pseudo', dest='hla_pseudo', help='Give MHC pseudo data');

args = parser.parse_args();
#allelenames  = args.allelenames;
#hla_pseudo =  args.hla_pseudo;
genotype  = args.genotype;
inputHLA  = args.inputHLA;
inputFile = args.inputFile; # all the intron derived 8-14 mer strings,total 24216054!
gtfIndex= args.gtfIndex;
fasta= args.fastaFile;
pepLen= args.pepLen;
outdir = args.outputFile;
thread = args.thread;
ref_proteins_file= args.refprotein

allelenames = dir + '/data/allelenames'
hla_pseudo = dir + '/data/MHC_pseudo.dat'

if thread:
    thread=min(os.cpu_count(), thread);
else:
    thread=4; # defualt thread number setting as 4

if ref_proteins_file:
    ref_proteins_file = ref_proteins_file;
else:
    ref_proteins_db_name= 'gencode.v32.pc_translations.fa.gz'
    ref_proteins_file = dir + '/data/' + ref_proteins_db_name; # defualt ref protein from GencodeV32

if outdir[-1]=='/':
    outdir =outdir;
else:
    outdir =outdir + '/';

pepLen_list= pepLen.split(',');
# python version check
pyversion= str(sys.version)[0]
#================================================================
#####
hla_ref = [l.rstrip() for l in open(allelenames)]
if genotype: # get hla results from arcasHLA json file
    hla_json = open(genotype);
    data= json.load(hla_json);
    # genes= list(data.keys());
    # for gene in genes:
        # hla_result= data[gene];
    hla_result= data['A'] + data['B'] + data['C']
    hla_list= ['HLA-' + x.replace('*','') for x in hla_result]
    hla_list= [x.split(':')[0] + ':' + x.split(':')[1] for x in hla_list];
    hla_list= [ x for x in hla_list if x in hla_ref ];
    hla_list = list(set(hla_list)); # combine the identical HLA
    hla_list.sort();
    hla_input = ','.join(hla_list)
    hla_String= [x[:5] + '*' + x[5:] for x in hla_list];
elif inputHLA: #get HLA from user input
    hla_list= inputHLA.split(',')
    hla_list= [ x for x in hla_list if x in hla_ref ];
    hla_list = list(set(hla_list)); # combine the identical HLA
    hla_list.sort();
    hla_input = ','.join(hla_list)
    hla_String= [x[:5] + '*' + x[5:] for x in hla_list]
else:
    hla_input=None;
    hla_String=[];
#================================================================
path, filename = os.path.split(inputFile)
basename, ext = os.path.splitext(filename)

# def split_file(file, lineCut=5000): #split file every 5000 rows
# #http://code.activestate.com/recipes/578045-split-up-text-file-by-line-count/
    # """Split a file based on a number of lines."""
    # path, filename = os.path.split(file)
    # basename, ext = os.path.splitext(filename)
    # # open input file
    # with open(file, 'r') as f_in:
        # try:
            # # open the first output file
            # f_out = open(os.path.join(path, '{}_{}{}'.format(basename, 0, ext)), 'w')
            # # loop over all lines in the input file, and number them
            # for i, line in enumerate(f_in):
                # if i % lineCut == 0:
                    # f_out.close()
                    # f_out = open(os.path.join(path, '{}_{}{}'.format(basename, i//lineCut, ext)), 'w')
                # f_out.write(line)
            # fileNum= i//lineCut + 1
        # finally:
            # f_out.close()
    # return(path, basename,ext,fileNum)

#
def seq2peptide(sequence):
    cds = sequence
    codontable = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M','ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K','AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L','CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q','CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V','GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E','GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S','TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TGC':'C', 'TGG':'W', 'TGT':'C', 'TAA':'_', 'TAG':'_','TGA':'_',
    }
    peptide_seq = '';
    for n in range(0, len(cds), 3):
        current_triple = cds[n:n+3]
        if current_triple.upper() in codontable:
            aa= codontable[current_triple.upper()];
            if aa=='_':
                break
            peptide_seq += aa;
        else:
            peptide_seq += '*'
            print ('unknown codon', current_triple.upper())
            break;
    return peptide_seq

def count_kmers(peptide, k):
    counts = {}
    num_kmers = len(peptide) - k + 1;
    # Loop over the kmer start positions
    for i in range(num_kmers):
        kmer = peptide[i:i+k];
        if kmer not in counts:
            counts[kmer] = 0;
        counts[kmer] += 1;
    return counts

def single_txSeq_kmer(proSeq):
    proSeq= proSeq
    tx_dict={}
    for z in pepLen_list:
        tx_dict.update(count_kmers(proSeq,int(z)));
    return tx_dict;

def single_cmd_run(file):
    # pass hla from main app
    p = subprocess.Popen(['netMHCpan','-hlapseudo', hla_pseudo,'-a', hla_input, '-p ', file], shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()[0];

    if pyversion =='3':
        netMHCpan_out= p.decode('utf-8');
    elif pyversion =='2':
        netMHCpan_out = p;

    netMHCpan_data= netMHCpan_out.splitlines();
    result = [];
    for line in netMHCpan_data:
        if any(hla for hla in hla_String if hla in line) and not any(x for x in ['Protein', 'allele','netMHCpan', 'neighbor'] if x in line):
            line_data= [x.strip() for x in line.split(' ') if x != ''];
            if float(line_data[12])!=0: # some very strong signal with 0.000
                result.append([line_data[1],line_data[2],line_data[12]]);
            else:
                result.append([line_data[1],line_data[2], str(round(1- float(line_data[11]),8))]);
    
    peplist = [line.strip() for line in open(file, 'r')]
    peplist.sort();
    resultS= [peplist]
    for h in hla_String:
        hlaList= [[x[1],x[2]] for x in result if x[0] == h];
        hlaList.sort();
        hla_value= [x[1] for x in hlaList];
        resultS.append(hla_value);
    newResult=np.transpose(resultS).tolist()
    #newResult= [list(x) for x in zip(*resultS)];
    file_suffix=file[-18:]
    file_num=re.search('peplist_(.*).txt', file_suffix).group(1);
    print(time.strftime("%Y-%m-%d %H:%M:%S"), 'Finish with ' + str(int(file_num)+1) +'/' + str(temp_file_num));
    return newResult;

def intron_summary(i):
    global summaryList;
    global data_list;
    global weakpeptide;# count num of weak pep for this IR event
    global strongpeptide;# count num of strong pep for this IR event

    intron= data_list[i]; # processed intron event
    geneid= intron[0].split('@')[0];
    geneid= geneid.split('.')[0];
    position= intron[0].split('@')[1]
    newid= geneid + '@' + position;
    
    intron_pep_sum= [newid];
    subpep_kmer= intron[1]; # all pep in this intron
    ############################################################
    weak_pep_set = [x for x in weakpeptide if x in subpep_kmer]
    weakNum= len(weak_pep_set);
    weak_peps = ','.join(weak_pep_set);
    strong_pep_set = [x for x in strongpeptide if x in subpep_kmer]
    strongNum= len(strong_pep_set);
    if strongNum>0:
        strong_peps = ','.join(strong_pep_set);
    else:
        strong_peps = '';
    ############################################################
    #print(subpep_kmer);print(len(subpep_kmer))
    subpeplist=[sublist for sublist in summaryList if sublist[0] in subpep_kmer];
    if len(subpeplist)>0:
        #print(subpeplist)
        bestScore={}
        for e in range(len(hla_String)):
            valueList= [float(item[e + 1]) for item in subpeplist];
            # print(valueList)
            bestPeptide= [item[0] for item in subpeplist if float(item[e + 1]) == min(valueList)][0];
            #print(bestPeptide)
            bestRank= min(valueList);
            #print(bestRank)
            bestScore[hla_String[e]]= bestRank;
            # stat= [round(np.quantile(valueList, x),2) for x in [0, 0.25, 0.5, 0.75,1]];#mimic describe function
            # meanRankscore = round((sum(valueList)/len(valueList)),2);
            # peptideNum = len(valueList);
            # strongbind_Num = sum(item1 <= 0.5 for item1 in valueList);# count strong binds<=0.5
            # weakbind_Num =  sum((item2>0.5) and (item2<=2) for item2 in valueList);# count weak binds 0.5~2
            # additonInfo = str(peptideNum) + '_' + str(meanRankscore) + '_' + str(strongbind_Num) + '_' + str(weakbind_Num);

            # intron_pep_sum.append(bestPeptide + "_" + "_".join(str(x) for x in stat) + "_" + additonInfo);
            intron_pep_sum.append(bestPeptide + "_" + str(bestRank));
        #print(bestScore)
        minAlleleScore= min([value for key, value in bestScore.items()]);
        intron_pep_sum.append(minAlleleScore);

        # here we assume we genotype A/B/C genes
        # Patient Harmonic-mean Best Rank
        alleleNum=0;
        recip_list=[]
        for g in ['A','B','C']: # 
            gene = [x for x in hla_String if 'HLA-'+ g in x];
            if len(gene)==1:
                geneScore = [bestScore[gene[0]]]*2;
                alleleNum +=2;
            elif len(gene)==2:
                geneScore = [bestScore[ge] for ge in gene];
                alleleNum +=2;
            elif len(gene)==0:
                geneScore = [];
                alleleNum +=0;
            recip_list= recip_list + geneScore;
        recip_list= [1/x for x in recip_list];
        intron_phbr= round(len(recip_list)/sum(recip_list),4)
        intron_pep_sum.append(intron_phbr);
        intron_pep_sum.append(weakNum);
        intron_pep_sum.append(strongNum);
        intron_pep_sum.append(weak_peps);
        intron_pep_sum.append(strong_peps);
        #print(intron_pep_sum)
    else:
        additional= [100]*(len(hla_String) + 2) + [0,0,0,0]
        intron_pep_sum= intron_pep_sum + additional;
    return intron_pep_sum;

if __name__ == '__main__':
    inputFile= pd.read_csv(inputFile, header=0, sep="\t");

    colnames= ['chr', 'feature', 'start', 'end', 'strand', 'geneid', 'txid', 'exonNum', 'exon_endsets','exon_txIds','exonEndLoc','exonStart','region', 'pep_8mer_start','pep_8mer_end', 'pep_9mer_start','pep_9mer_end', 'pep_10mer_start','pep_10mer_end',
    'pep_11mer_start','pep_11mer_end', 'pep_12mer_start','pep_12mer_end', 'pep_13mer_start','pep_13mer_end',
    'pep_14mer_start','pep_14mer_end']
    file=pd.read_csv(gtfIndex, header=None, names=colnames,  sep="\t");
    # select intron from input list
    file= file.loc[(file['geneid'].str[:15] + '-' + file['region'].astype(str) + '@' + file['chr'] + ':' +  file['start'].astype(str) + '-' + file['end'].astype(str)).isin(inputFile['id'])]
    file=file[(file['end'] - file['start'])<=10000];# second filter-skip very large intron

    # translating the intron sequence into amino acids
    print(time.strftime("%Y-%m-%d %H:%M:%S"), 'Translating the intron sequence into amino acids')
    df= pd.DataFrame();
    df['name']= file['geneid']+ '@' + file['chr'].astype(str) + ':' + file['start'].astype(str) + '-' + file['end'].astype(str);
    #pepLen_list= pepLen.split(',');
    for i in pepLen_list:
        #print(i);
        x= file[['chr']]
        x['start']= file['pep_' + str(i) + 'mer_start'].astype(int);
        x['end']= file['pep_' + str(i) + 'mer_end'].astype(int);
        x['name']= file['geneid']+ '@' + file['chr'].astype(str) + ':' + file['start'].astype(str) + '-' + file['end'].astype(str);
        x['score']= 1
        x['strand'] = file['strand']
        #
        x_input = pybedtools.BedTool.from_dataframe(x)
        fasta = pybedtools.example_filename(fasta)
        x_cDNA = x_input.sequence(fi=fasta,s=True); p=[];
        with open(x_cDNA.seqfn) as f:
            for count, line in enumerate(f, start=1):
                if count % 2 == 0:
                    cdna_seq = line.replace("\n", "");
                    aa_line= [seq2peptide(cdna_seq)][0]
                    if len(aa_line)>=int(i):# make sure at least one letter from intron!
                        aa_line = aa_line;
                    else:
                        aa_line = '';
                    p.append(aa_line);
        df['pep_' + str(i)]= p

    maxL= max([int(x) for x in pepLen_list]);
    df = df[df['pep_' + str(maxL)]!=''];
    data_dict={};
    global data_list;
    data_list=[];# for future used with estimate for each events; this may cause additional mem cost
    for j in range(df.shape[0]):
        dict_intron={}
        for z in pepLen_list:
            dict_intron.update(count_kmers(df.iloc[j]['pep_' + str(z)],int(z)));
        #update kmer for each intron    
        data_list.append([df.iloc[j]['name'],list(dict_intron.keys())]);# use for stat every intron
        #update for all kmer for all introns
        data_dict.update(dict_intron);

    # processing the background host peptide sequence hold immune tolerance
    txSeq=[]
    if ref_proteins_file.lower().endswith('.gz'):
        print("You are using ref protein file:", ref_proteins_db_name)
        with gzip.open(ref_proteins_file, "rt") as ref_proteins:
            for record in SeqIO.parse(ref_proteins, "fasta"):
                txSeq.append(str(record.seq))
        #print(len(txSeq), txSeq[:5])
    elif ref_proteins_file.lower().endswith(('.fasta', '.fa')):
        print("You are using ref pretein file:", ref_proteins_file)
        with open(ref_proteins_file, "rt") as ref_proteins:
            for record in SeqIO.parse(ref_proteins, "fasta"):
                txSeq.append(str(record.seq));
    else:
        print("check you reference proteins input!!!")

    host_pep_dict={}
    pool=Pool(processes=thread); # defualt thread number setting as 4
    process = pool.map(single_txSeq_kmer, txSeq); #multiple process
    for tx_dict in process:
        host_pep_dict.update(tx_dict); 
    print(time.strftime("%Y-%m-%d %H:%M:%S"), 'load ref protein database information!')
    print(time.strftime("%Y-%m-%d %H:%M:%S"), len(host_pep_dict))

    intron_novel_dict= {k:v for k,v in data_dict.items() if k not in host_pep_dict}
    print(time.strftime("%Y-%m-%d %H:%M:%S"), 'Filter the intron peptide with reference proteins!')

    list=[];
    for key, value in intron_novel_dict.items():
        list.append(key);

    # split list into every 5000 list;
    print(time.strftime("%Y-%m-%d %H:%M:%S"), 'Split list into small temp peptide files')
    split_list= [list[i:i + 3000] for i in range(0, len(list), 3000)];
    temp_file_num= len(split_list);
    for m in range(temp_file_num):
        with open(outdir + basename + '_peplist_' + str(m) + '.txt', 'w') as f:
            f.write("\n".join(split_list[m]));
    tmpfiles= [outdir + basename + '_peplist_' + str(m) + '.txt' for m in range(len(split_list))];# for multiprocess
    ################################################################################

    print(time.strftime("%Y-%m-%d %H:%M:%S"), 'Doing with Multiprcessing jobs')
    pool=Pool(processes=thread); # defualt thread number setting as 4

    process = pool.map(single_cmd_run, tmpfiles); #multiple process

    global summaryList;
    summaryList = [item for result_list in process for item in result_list];
    print(time.strftime("%Y-%m-%d %H:%M:%S"), 'Summary all the piece file results!')
    #print(len(summaryList))
    pool.close();
    pool.join();
    ################################################################################
    header= ['peptide'] + hla_String;
    # print out the weak and strong binding peptides
    weakList= [sublist for sublist in summaryList if any(x for x in sublist[1:] if float(x) < 2)];
    global weakpeptide;
    weakpeptide=[x[0] for x in weakList]; # for count num of weak pep for each intron
    weakList= [header] + weakList;
    df = pd.DataFrame(weakList[1:],columns=weakList[0]);
    df.to_csv(outdir + basename + '_weak_bind_peptide.txt', index=False, sep='\t');
    
    ################################################################################
    strongList= [sublist for sublist in summaryList if any(x for x in sublist[1:] if float(x) <= 0.5)];
    global strongpeptide;
    strongpeptide=[x[0] for x in strongList]; # for count num of strong pep for each intron
    strongList= [header] + strongList;
    df = pd.DataFrame(strongList[1:],columns=strongList[0]);
    df.to_csv(outdir + basename + '_strong_bind_peptide.txt', index=False, sep='\t');

    # summary intron events presentation with multiprocess
    ################################################################################
    pool=Pool(processes=thread); # defualt thread number setting as 4
    process2 = pool.map(intron_summary, range(len(data_list))); #multiple process
    final_result= [item for item in process2];
    print(time.strftime("%Y-%m-%d %H:%M:%S"), 'Summary all the results for each intron!')

    final_result= [x for x in final_result if x[-6] < 2]; #only get the singnificant introns
    #print(final_result)
    final_result= [['peptide'] + hla_String + ['minRank', 'PHBR','weakBind_Num','strongBind_Num','weakBind_peps','strongBind_peps']] + final_result;
    df = pd.DataFrame(final_result[1:],columns=final_result[0]);
    
    df.to_csv(outdir + basename + '_Final_output.txt', index=False, sep='\t');
    
    for f1 in tmpfiles:
        os.remove(f1);#delete the temp processed file!
    print(time.strftime("%Y-%m-%d %H:%M:%S"), 'All jobs done!')
#
#
