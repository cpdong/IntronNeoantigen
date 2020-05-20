import argparse;
import os,shutil;
import subprocess;
import HTSeq,re,time;
import collections;
import pandas as pd;
from multiprocessing import Pool;
#
parser = argparse.ArgumentParser()
parser.add_argument('-b', metavar = 'input', dest='bamFile', help='Give input file fullname');
parser.add_argument('-g', metavar = 'gtffile', dest='gtfFile', help='Give parsed gtf file fullname');
parser.add_argument('-t', dest='thread', type=int, help='number of multiple thread number,>0');
parser.add_argument('-o', metavar = 'output', dest='outputFile',help='Give output file fullname');
parser.add_argument('-j', metavar = 'anchor length', type=int, dest='junction',help='Give a specific anchor length');

parser.add_argument('-c', metavar = 'intronRead filter', dest='intron_read', type=int,help='Set a intron read filter');
parser.add_argument('-p', metavar = 'intronPSI filter', dest='intron_psi', type=float,help='threshhold of psi value filter');
parser.add_argument('-n', metavar = 'novel intron read read filter', dest='novel_read', type=int,help='Set a novel intron read filter threshhold');

args = parser.parse_args();
bamFile = args.bamFile;
gtf = args.gtfFile;
thread = args.thread;
output = args.outputFile;
junction = args.junction;
intron_read_filter = args.intron_read;
intron_psi_filter = args.intron_psi;
novel_read_filter = args.novel_read;

cpus= os.cpu_count();# get all available cpu counts
if thread:
    thread=thread;
    thread= min(thread, cpus); # make sure the maximum thread not over available cpus!
else:
    thread=1;

if junction: # the anchor length setting
    junction= junction;
else:
    junction= 10;

if intron_read_filter: # the anchor length setting
    intron_read_filter = intron_read_filter;
else:
    intron_read_filter = 50;

if intron_psi_filter: # the anchor length setting
    intron_psi_filter = intron_psi_filter;
else:
    intron_psi_filter = 0.05;

if novel_read_filter: # the anchor length setting
    novel_read_filter = novel_read_filter;
else:
    novel_read_filter= 15;
#
#
def no_parallel_run():
    global features;
    global bamfile;

    counts = collections.Counter()
    bam = HTSeq.BAM_Reader(bamfile)
    for read in bam: 
        #filter1, uniq mapp filter:the option as NM/NH/CC/CP/HI
        NH_list= [x for x in read.optional_fields if x[0]=='NH'];
        if len(NH_list)==0:
            NH_value= 10; # if not NH tag find assign as 10;
        elif len(NH_list)==1:
            NH_value = NH_list[0][1];
        if NH_value == 1 or read.aQual>=30: #filter1, uniq mapp filter
            cigar_list= [cstring.type for cstring in read.cigar]
            cigar_check= [e for e in cigar_list if e in ['D','I','S','H','P','X','=' ]];
            if len(cigar_check)==0:
                counts['uniq_mapped_read'] += 1;#count every unique mapped read for calculate rpkm value
                if cigar_list.count('N') ==0:# call directly
                    gene_ids = []; iv_lens=[]
                    for iv, val in features[ read.iv ].steps():
                        if len(val) ==1:
                            gene_ids.append(list(val)[0]);
                            ivList= re.split('[:[,)/]',str(iv));
                            iv_lens.append(str(int(ivList[3]) - int(ivList[2])))
                    if len(gene_ids)>0 and len(set([x[:15] for x in gene_ids]))==1: # remove ambigous mapped read-map to two diff genes
                        if len(gene_ids)==1:
                            gene_id = gene_ids[0];
                            counts[gene_id] += 1;
                        elif len(gene_ids)>1: #strict the intron at least have length large than anchor
                            for g in range(len(gene_ids)):
                                gene_id = gene_ids[g];
                                iv_len = int(iv_lens[g]);
                                geneInfo = re.split('-|@', gene_id)[1]
                                if(int(geneInfo)%2==1):
                                    counts[gene_id] += 1;
                                elif (int(geneInfo)%2==0) and (iv_len > junction):
                                    counts[gene_id] += 1;
                elif cigar_list.count('N') ==1:
                    gene_ids = []; iv_lens=[]
                    for cigop in read.cigar:
                        if cigop.type == "M":
                            for iv, val in features[cigop.ref_iv].steps():
                                if len(val) ==1:
                                    gene_ids.append(list(val)[0]);
                                    ivList= re.split('[:[,)/]',str(iv));
                                    iv_lens.append(str(int(ivList[3]) - int(ivList[2])))
                    if len(set([x[:15] for x in gene_ids]))==1: # remove ambigous mapped read-map to two diff genes
                        # get two terminal end of skip region
                        if len(gene_ids)>1:
                            skip_iv= re.split('[:[,)/]',str(read.cigar[1].ref_iv));
                            new_query= HTSeq.GenomicInterval(skip_iv[0], int(skip_iv[2])-2,int(skip_iv[3])+3, skip_iv[5]); # make sure the enlarged map window to clarify the mapping region/s!
                            query_ids=set()
                            for iv, val in features[new_query].steps():
                                query_ids |= val;
                            if len(query_ids)>0:
                                query_list= [int(re.split('-|@', x)[1]) for x in list(query_ids)];
                                query_list.sort();# the mapped region list check
                                if query_list[0]%2!=0 and query_list[-1]%2!=0:# if both end at exons
                                    for g in range(len(gene_ids)):
                                        gene_id = gene_ids[g];
                                        iv_len = int(iv_lens[g]);
                                        geneInfo = re.split('-|@', gene_id)[1]
                                        if(int(geneInfo)%2==1):
                                            counts[gene_id] += 1;
                                        elif (int(geneInfo)%2==0) and (iv_len > junction):
                                            counts[gene_id] += 1;
                                else: # all should be  add as novel intron/partial retained events
                                    novel_id= (list(query_ids)[0])[0:15] + '@' + skip_iv[0] + ':u' + skip_iv[2] + '-d' + skip_iv[3];
                                    counts[novel_id] += 1;
                        elif len(gene_ids)==1: # if disruption occurs inner one feature
                            skip_iv= re.split('[:[,)/]',str(read.cigar[1].ref_iv));
                            gene_id = list(gene_ids)[0];
                            geneInfo = re.split('-|@', gene_id)[1]
                            if int(geneInfo)%2==0: # if read mapped within a intron, count it as novel events
                                novel_id= gene_id[0:15] + '@' + skip_iv[0] + ':u' + skip_iv[2] + '-d' + skip_iv[3];
                                counts[novel_id] += 1;
    print(time.strftime("%Y-%m-%d %H:%M:%S") + ": HTSeq count finished!")
    return counts

def do_parallel_run(part_N):

    global thread;
    global features;
    global tmp_prefix;

    bamfile_temp= tmp_prefix + "{0:0=2d}".format(part_N) + '.bam'
    counts = collections.Counter()
    bam = HTSeq.BAM_Reader(bamfile_temp)
    for read in bam: 
        #filter1, uniq mapp filter:the option as NM/NH/CC/CP/HI
        NH_list= [x for x in read.optional_fields if x[0]=='NH'];
        if len(NH_list)==0:
            NH_value= 10; # if not NH tag find assign as 10;
        elif len(NH_list)==1:
            NH_value = NH_list[0][1];
        if NH_value == 1 or read.aQual>=30: #filter1, uniq mapp filter
            cigar_list= [cstring.type for cstring in read.cigar]
            cigar_check= [e for e in cigar_list if e in ['D','I','S','H','P','X','=' ]];
            if len(cigar_check)==0:
                counts['uniq_mapped_read'] += 1;#count every unique mapped read for calculate rpkm value
                if cigar_list.count('N') ==0:# call directly
                    gene_ids = []; iv_lens=[]
                    for iv, val in features[ read.iv ].steps():
                        if len(val) ==1:
                            gene_ids.append(list(val)[0]);
                            ivList= re.split('[:[,)/]',str(iv));
                            iv_lens.append(str(int(ivList[3]) - int(ivList[2])))
                    if len(gene_ids)>0 and len(set([x[:15] for x in gene_ids]))==1: # remove ambigous mapped read-map to two diff genes
                        if len(gene_ids)==1:
                            gene_id = gene_ids[0];
                            counts[gene_id] += 1;
                        elif len(gene_ids)>1: #strict the intron at least have length large than anchor
                            for g in range(len(gene_ids)):
                                gene_id = gene_ids[g];
                                iv_len = int(iv_lens[g]);
                                geneInfo = re.split('-|@', gene_id)[1]
                                if(int(geneInfo)%2==1):
                                    counts[gene_id] += 1;
                                elif (int(geneInfo)%2==0) and (iv_len > junction):
                                    counts[gene_id] += 1;
                elif cigar_list.count('N') ==1:
                    gene_ids = []; iv_lens=[]
                    for cigop in read.cigar:
                        if cigop.type == "M":
                            for iv, val in features[cigop.ref_iv].steps():
                                if len(val) ==1:
                                    gene_ids.append(list(val)[0]);
                                    ivList= re.split('[:[,)/]',str(iv));
                                    iv_lens.append(str(int(ivList[3]) - int(ivList[2])))
                    if len(set([x[:15] for x in gene_ids]))==1: # remove ambigous mapped read-map to two diff genes
                        # get two terminal end of skip region
                        if len(gene_ids)>1:
                            skip_iv= re.split('[:[,)/]',str(read.cigar[1].ref_iv));
                            new_query= HTSeq.GenomicInterval(skip_iv[0], int(skip_iv[2])-2,int(skip_iv[3])+3, skip_iv[5]); # make sure the enlarged map window to clarify the mapping region/s!
                            query_ids=set()
                            for iv, val in features[new_query].steps():
                                query_ids |= val;
                            if len(query_ids)>0:
                                query_list= [int(re.split('-|@', x)[1]) for x in list(query_ids)];
                                query_list.sort();# the mapped region list check
                                if query_list[0]%2!=0 and query_list[-1]%2!=0:# if both end at exons
                                    for g in range(len(gene_ids)):
                                        gene_id = gene_ids[g];
                                        iv_len = int(iv_lens[g]);
                                        geneInfo = re.split('-|@', gene_id)[1]
                                        if(int(geneInfo)%2==1):
                                            counts[gene_id] += 1;
                                        elif (int(geneInfo)%2==0) and (iv_len > junction):
                                            counts[gene_id] += 1;
                                else: # all should be  add as novel intron/partial retained events
                                    novel_id= (list(query_ids)[0])[0:15] + '@' + skip_iv[0] + ':u' + skip_iv[2] + '-d' + skip_iv[3];
                                    counts[novel_id] += 1;
                        elif len(gene_ids)==1: # if disruption occurs inner one feature
                            skip_iv= re.split('[:[,)/]',str(read.cigar[1].ref_iv));
                            gene_id = list(gene_ids)[0];
                            geneInfo = re.split('-|@', gene_id)[1]
                            if int(geneInfo)%2==0: # if read mapped within a intron, count it as novel events
                                novel_id= gene_id[0:15] + '@' + skip_iv[0] + ':u' + skip_iv[2] + '-d' + skip_iv[3];
                                counts[novel_id] += 1;
    print(time.strftime("%Y-%m-%d %H:%M:%S") + ": Process " + str(part_N) + " finished!")
    return counts

def single_thread_RPKM(counter):
    geneSum = collections.Counter()
    countSets = collections.Counter()
    for k, v in counter.items():
        countSets[k] +=v;
        geneSum[k[:15]] +=v;

    del countSets["uniq_mapped_read"];
    keys= [key for key, value in countSets.items()];
    rpkm_list = [];
    for key in keys:
        value = countSets[key];
        list=[key, value] + [None]*10;
        extron_list= re.split('-|@|:|-', key);
        list[3] = geneSum[extron_list[0]];
        if len(extron_list)==5:
            extron_len = int(extron_list[4]) - int(extron_list[3]) + 1;
            list[9] = value*1000/extron_len;
            list[4] = round(value*(10**9)/(geneSum['uniq_mapped_rea'] * extron_len), 4);
            if (int(extron_list[1])%2==0):
                list[2] = 'intron';
                uExon = extron_list[0] + '-' + str(int(extron_list[1]) -1); 
                dExon = extron_list[0] + '-' + str(int(extron_list[1]) + 1); 
                upperExon = [[key, value] for key, value in countSets.items() if uExon in key];
                downExon = [[key, value] for key, value in countSets.items() if dExon in key]
                if len(upperExon) ==1:
                    upperList = re.split('-|@|:|-', upperExon[0][0]);
                    list[6] = upperExon[0][1]; list[10] = int(upperList[4])-int(upperList[3])+1; 
                if len( downExon) == 1:
                    downList = re.split('-|@|:|-', downExon[0][0]);
                    list[7] = downExon[0][1]; list[11] = int(downList[4])-int(downList[3])+1;
            elif (int(extron_list[1])%2==1):
                list[2] = 'exon';
        elif len(extron_list)==4:
            list[2] = 'novel-intron';
        
        rpkm_list.append(list)
    return rpkm_list

def do_RPKM_run(key):
    global geneSum;
    global countSets;

    value = countSets[key];
    list=[key, value] + [None]*10;
    extron_list= re.split('-|@|:|-', key);
    list[3] = geneSum[extron_list[0]];
    if len(extron_list)==5:
        extron_len = int(extron_list[4]) - int(extron_list[3]) + 1;
        list[9] = value*1000/extron_len;
        list[4] = round(value*(10**9)/(geneSum['uniq_mapped_rea'] * extron_len), 4);
        if (int(extron_list[1])%2==0):
            list[2] = 'intron';
            uExon = extron_list[0] + '-' + str(int(extron_list[1]) -1); 
            dExon = extron_list[0] + '-' + str(int(extron_list[1]) + 1); 
            upperExon = [[key, value] for key, value in countSets.items() if uExon in key];
            downExon = [[key, value] for key, value in countSets.items() if dExon in key]
            if len(upperExon) ==1:
                upperList = re.split('-|@|:|-', upperExon[0][0]);
                list[6] = upperExon[0][1]; list[10] = int(upperList[4])-int(upperList[3])+1; 
            if len( downExon) == 1:
                downList = re.split('-|@|:|-', downExon[0][0]);
                list[7] = downExon[0][1]; list[11] = int(downList[4])-int(downList[3])+1;
        elif (int(extron_list[1])%2==1):
            list[2] = 'exon';
    elif len(extron_list)==4:
        list[2] = 'novel-intron';
    return list

def single_thread_TPM(rpkm_list):
    listNum = len(rpkm_list);
    tpmList = [x[9] for x in rpkm_list];
    tpmSum = sum([x for x in tpmList if x != None]);
    tpm_list = [];
    for i in range(listNum):
        list_tpm = rpkm_list[i];
        if list_tpm[2] == 'intron':
            extron_list= re.split('-|@|:|-', list_tpm[0]);
            extron_len = int(extron_list[4]) - int(extron_list[3]) + 1;
            tpm_val = list_tpm[1]*(10**9)/(extron_len * tpmSum);
            list_tpm[5]= round(tpm_val, 4);
            if list_tpm[6] !=None:
                tpm_upper =  list_tpm[6]*(10**9)/(list_tpm[10]* tpmSum);
            else:
                tpm_upper =  0;
            if list_tpm[7] !=None:
                tpm_down =  list_tpm[7]*(10**9)/(list_tpm[11]* tpmSum);
            else:
                tpm_down =  0;
            tpm_flake=[tpm_upper, tpm_down]
            
            if sum(tpm_flake)!=0:
                list_tpm[8] = round(tpm_val*2/sum(tpm_flake), 3);
        else:
            list_tpm = list_tpm;
        tpm_list.append(list_tpm);
    return tpm_list

def do_TPM_run(index):
    list_tpm = rpkm_list[index];
    if list_tpm[2] == 'intron':
        extron_list= re.split('-|@|:|-', list_tpm[0]);
        extron_len = int(extron_list[4]) - int(extron_list[3]) + 1;
        tpm_val = list_tpm[1]*(10**9)/(extron_len * tpmSum);

        list_tpm[5]= round(tpm_val, 4);

        if list_tpm[6] !=None:
            tpm_upper =  list_tpm[6]*(10**9)/(list_tpm[10]* tpmSum);
        else:
            tpm_upper =  0;
        if list_tpm[7] !=None:
            tpm_down =  list_tpm[7]*(10**9)/(list_tpm[11]* tpmSum);
        else:
            tpm_down =  0;
        tpm_flake=[tpm_upper, tpm_down]
        
        if sum(tpm_flake)!=0:
            list_tpm[8] = round(tpm_val*2/sum(tpm_flake), 3);
    else:
        list_tpm = list_tpm
    return list_tpm


if __name__ == '__main__':
    global features;
    features =  HTSeq.GenomicArrayOfSets( "auto", stranded=True )
    for line in open(gtf):
        fields = line.split( "\t" );
        chrom= fields[0]; start=int(fields[2]); end= int(fields[3])+1; strand=fields[4];
        name= fields[5].split('.')[0] + '-' + str(fields[12]) + '@' + fields[0] + ':' + str(fields[2]) + '-' + str(fields[3]);
        iv = HTSeq.GenomicInterval(chrom, start, end, strand)
        features[ iv ] += name;

    global bamfile;
    bamfile = bamFile;

    dir= '/'.join(bamfile.split('/')[0:-1]);
    s1= bamfile.split('/')[-1];
    sample= s1.replace('.bam', '');

    if output:
        if output[-1] =='/':
            output = output;
        else:
            output = output + '/'
    else:
        output = dir + '/'; # defualt located at your input bam file 

    if thread==1:
        print("No parallel mode, thread =",thread)
        results= no_parallel_run();

        with open(output + sample + '_rawCount.txt', 'w') as f: 
            for tag, count in results.items():  
                f.write('{}\t{}\n'.format(tag, count))
        rpkm_list= single_thread_RPKM(results);
        print(time.strftime("%Y-%m-%d %H:%M:%S") + ": RPKM expression calculated!")
        tpm_list = single_thread_TPM(rpkm_list);
        print(time.strftime("%Y-%m-%d %H:%M:%S") + ": TPM and PSI calculated!")

        tpm_list = [x[:9] for x in tpm_list if x[2] !='exon']
        colname = [['id','count','feature','fullGeneCount','rpkm','tpm','uppExonCount','downExonCount','psi']];
        data= colname + tpm_list;
        data = pd.DataFrame(data[1:],columns=data[0]);
        data.to_csv(output + sample + '_intron_calling_Rawresult.txt',  index=False, sep='\t');

        intron_list = [x for x in tpm_list if x[2] =='intron' and x[8]!=None and float(x[1]) >= intron_read_filter and float(x[8]) >= intron_psi_filter]
        intron_list = colname + intron_list;
        intron_data = pd.DataFrame(intron_list[1:], columns=intron_list[0]);
        intron_data.to_csv(output + sample + '_intron_candidates.txt',  index=False, sep='\t');

        novel_list = [x for x in tpm_list if x[2] =='novel-intron' and x[1] >= novel_read_filter and x[1]/x[3] > 0.001]
        novel_list = colname + novel_list;
        novel_data = pd.DataFrame(novel_list[1:], columns=novel_list[0]);
        novel_data.to_csv(output + sample + '_novel-intron_candidates.txt',  index=False, sep='\t');

        print(time.strftime("%Y-%m-%d %H:%M:%S") + ": All intron-calling job done!")
#
    else:
        print("thread =",thread)

        if not os.path.exists(dir + '/tmp_' + sample):
            os.makedirs(dir + '/tmp_' + sample)# create a folder for speeding calculation

        global tmp_prefix;
        tmp_prefix = dir + '/tmp_' + sample + '/tmp_'

        p0=subprocess.Popen("samtools view -@ " + str(thread) + " -c " + bamfile, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()[0];
        readNum= int(p0.decode('utf-8').replace('\n', ''));# get bam total read number
        split_num= int(readNum/thread + 50);# set a split thread to split big bam to samll ones
        print(time.strftime("%Y-%m-%d %H:%M:%S") + ": There are " + str(readNum) + " reads in this BAM file!" )

        p1=subprocess.Popen("samtools view -H " + bamfile +" > " + tmp_prefix + "header; samtools view -@ " + str(thread) +" " + bamfile + " | split - "+ tmp_prefix + " --numeric-suffixes -l " + str(split_num) + " --filter='cat " + tmp_prefix + "header - | samtools view -@ " + str(thread) + " -b - > $FILE.bam'", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate();
        print(time.strftime("%Y-%m-%d %H:%M:%S") + ": Split BAM file into " +  str(thread) + " small BAM files!")
        
        pool=Pool(processes=thread);
        process = pool.map(do_parallel_run, range(thread)); #multiple process

        result_list = [result for result in process];

        results=collections.Counter();
        for result in result_list:
            results += result;
        print(time.strftime("%Y-%m-%d %H:%M:%S") + ": It may take 2~10 minutes to get the final results")

        with open(output + sample + '_rawCount.txt', 'w') as f: 
            for tag, count in results.items():  
                f.write('{}\t{}\n'.format(tag, count))
        shutil.rmtree(dir + '/tmp_' + sample) # remnove the temp folder

        global geneSum, countSets;
        geneSum = collections.Counter()
        countSets = collections.Counter()
        for k, v in results.items():
            countSets[k] +=v;
            geneSum[k[:15]] +=v;

        del countSets["uniq_mapped_read"];
        keys= [key for key, value in countSets.items()]

        pool=Pool(processes=thread);
        rpkm_process = pool.map(do_RPKM_run, keys); #multiple process
        rpkm_list = [list for list in rpkm_process];
        print(time.strftime("%Y-%m-%d %H:%M:%S") + ": RPKM expression calculated!")

        # calculate tpm sum and psi value from flaking exons everate tpm
        listNum = len(rpkm_list);
        tpmList = [x[9] for x in rpkm_list];

        global rpkmList, tpmSum;
        tpmSum = sum([x for x in tpmList if x != None]);
        rpkmList = rpkm_list;
        pool=Pool(processes=thread);
        tpm_process = pool.map(do_TPM_run, range(listNum)); #multiple process
        tpm_list = [list_tpm for list_tpm in tpm_process];
        pool.close()
        pool.join()
        print(time.strftime("%Y-%m-%d %H:%M:%S") + ": TPM and PSI calculated!")

        tpm_list = [x[:9] for x in tpm_list if x[2] !='exon']
        colname = [['id','count','feature','fullGeneCount','rpkm','tpm','uppExonCount','downExonCount','psi']];
        data= colname + tpm_list;
        data = pd.DataFrame(data[1:],columns=data[0]);
        data.to_csv(output + sample + '_intron_calling_Rawresult.txt',  index=False, sep='\t');

        intron_list = [x for x in tpm_list if x[2] =='intron' and x[8]!=None and float(x[1]) >= intron_read_filter and float(x[8]) >= intron_psi_filter and x[6]!=None and x[7]!=None and x[6]>= 10 and x[7]>=10];
        intron_list = colname + intron_list;
        intron_data = pd.DataFrame(intron_list[1:], columns=intron_list[0]);
        intron_data.to_csv(output + sample + '_intron_candidates.txt',  index=False, sep='\t');

        novel_list = [x for x in tpm_list if x[2] =='novel-intron' and x[1] >= novel_read_filter]
        novel_list = colname + novel_list;
        novel_data = pd.DataFrame(novel_list[1:], columns=novel_list[0]);
        novel_data.to_csv(output + sample + '_novel-intron_candidates.txt',  index=False, sep='\t');

        print(time.strftime("%Y-%m-%d %H:%M:%S") + ": All intron-calling job done!")
#
