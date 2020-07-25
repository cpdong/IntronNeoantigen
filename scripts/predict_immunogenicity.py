import argparse;
import subprocess,os,shutil;
import pickle,time;
import pandas as pd;
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from Bio import SeqIO;
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
from multiprocessing import Pool;
from math import log, exp
from os.path import dirname, abspath
dir = dirname(dirname(abspath(__file__)));# get parent dir of the scripts
#
#
parser = argparse.ArgumentParser()
parser.add_argument('-file', metavar = 'input', dest='inputFile', help='Give the intron calling result files');
parser.add_argument('-t', dest='thread', type=int, help='number of multiple thread number,>0');
parser.add_argument('-outdir', metavar = 'output', dest='outputFile',help='Give output file fullname');
#
args = parser.parse_args();
inputFile = args.inputFile; # all the intron derived 8-14 mer strings,total 24216054!
outdir = args.outputFile;
thread = args.thread;
#
#
input_dir= '/'.join(inputFile.split('/')[0:-1]);
tmp_dir= input_dir + '/tmp'
s1= inputFile.split('/')[-1];
sample= s1.replace('.txt', '');
#
#outdir='/N/dc2/scratch/cpdong/randomForest/GSE78220/STARmapping/SRR3184279'
if outdir: # if output dir is give
    if outdir[-1]=='/':
        outdir =outdir;
    else:
        outdir =outdir + '/';
else:
    outdir = input_dir + '/';
#
if thread:
    thread=min(os.cpu_count(), thread);
else:
    thread=4; # defualt thread number setting as 4
#
ref_IEDB= dir + '/data/IEDB_infectious.fasta'
iedb_seq = []
for record in SeqIO.parse(ref_IEDB, "fasta"):
    seqString= str(record.seq).upper()
    AA_list= [char for char in seqString]
    if not any(x for x in AA_list if x not in "ACDEFGHIKLMNPQRSTVWY"):
        iedb_seq.append(seqString)
#
#
def getWTpep(pep_input):
    rid = pep_input.split('@')[0]
    mut_pep = pep_input.split('@')[1]
    pepmatch = '/N/slate/cpdong/software/IntronNeoantigen/pepmatch_db_x86_64'
    ref_seq = '/N/slate/cpdong/software/IntronNeoantigen/data/reference_peptide_{}.txt'.format(len(mut_pep))
    
    pepFile= tmp_dir + '/' + rid + '_' + mut_pep + '.fa'
    with open(pepFile,'w') as fout:
        fout.write('{}\n'.format(mut_pep))
    p= subprocess.Popen([pepmatch,'-thr', '6',pepFile, ref_seq], shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()[0];
    out= p.decode('utf-8').splitlines();
    os.remove(pepFile);
    line=[x for x in out if '#' not in x]
    line_data= line[0].split();
    print(rid, line_data)
    if line_data[5] =='unique':
        match_pep= '---'
        mismatch= '---'
    else:
        match_pep= line_data[3]
        mismatch= line_data[5]
    return [rid,match_pep, mismatch]

def aligner(seq1,seq2):
    matrix = matlist.blosum62
    gap_open = -11
    gap_extend = -1
    almnt = pairwise2.align.localds(seq1.upper(), seq2.upper(), matrix, gap_open, gap_extend)
    return almnt

def logSum(v):
    max_v = max(v)
    return log(sum(map(lambda x: exp(x-max_v),v))) + max_v

def getDAI(neo_aff,wt_aff): # Differential agretopicity index
    if neo_aff and wt_aff:
        DAI=log(neo_aff/wt_aff)
    else:
        print("check you input affinity IC50 value")
    return round(DAI,4)

def getR(neo_seq,iedb_seq):
    align_score = []
    a = 26
    k = 4.86936
    for seq in iedb_seq:
        aln_score = aligner(neo_seq,seq)
        if aln_score:
            localds_core = max([line[2] for line in aln_score])
            align_score.append(localds_core)

    bindingEnergies = list(map(lambda x: -k * (a - x), align_score))
    lZk = logSum(bindingEnergies + [0])
    lGb = logSum(bindingEnergies)
    R=exp(lGb-lZk)
    return round(R,4)

def getNRP(DAI,R):
    NRP_score= DAI*R
    return round(NRP_score,4)

def hydro_score(pep):
    hydro_dict=dict(R=-4.5,K=-3.9,N=-3.5,D=-3.5,Q=-3.5,E=-3.5,H=-3.2,P=-1.6,Y=-1.3,W=-0.9,S=-.8,T=-0.7,G=-0.4,A=1.8,M=1.9,C=2.5,F=2.8,L=3.8,V=4.2,I=4.5)
    hydro_vector = []
    for i in pep:
        hydro_vector.append(hydro_dict[i.upper()])
    hydro_sum= round(sum(hydro_vector),2)
    return hydro_sum

def polarity_score(pep):
    polarity_dict = dict(A=0.00,R=52.00,N=3.38,D=49.70,C=1.48,Q=3.53,E=49.90,G=0.00,H=51.60,I=0.13,L=0.13,K=49.50,M=1.43,F=0.35,P=1.58,S=1.67,T=1.66,W=2.10,Y=1.61,V=0.13)
    polar_vector = []
    for i in pep:
        polar_vector.append(polarity_dict[i.upper()])
    polar_vector_new= polar_vector[1:3] + polar_vector[4:6]; # choose 2,3,5,6 position
    polar_sum= round(sum(polar_vector_new),2)
    return polar_sum

def charge_score(pep):
    netcharge_dict=dict(A=0,R=1,N=0,D=-1,C=0,Q=0,E=-1,G=0,H=0,I=0,L=0,K=1,M=0,F=0,P=0,S=0,T=0,W=0,Y=0,V=0)
    charge_vector = []
    for i in pep:
        charge_vector.append(netcharge_dict[i.upper()])
    charge_vector_new= charge_vector[1:3] + charge_vector[4:6]; # choose 2,3,5,6 position
    charge_sum= round(sum(charge_vector_new),2)
    return charge_sum

def entropy_score(pep):
    entropy_dict=dict(A=154.33,R=341.01,N=207.90,D=194.91,C=219.79,Q=235.51,E=223.16,G=127.90,H=242.54,I=233.21,L=232.30,K=300.46,M=202.65,F=204.74,P=179.93,S=174.06,T=205.80,W=237.01,Y=229.15,V=207.60)
    entropy_vector = []
    for i in pep:
        entropy_vector.append(entropy_dict[i.upper()])
    entropy_sum= round(sum(entropy_vector),2)
    return entropy_sum

def molecule_size(pep):
    molsize_dict=dict(A=89,R=174,N=132,D=133,C=121,E=147,Q=146,G=75,H=155,I=131,L=131,K=146,M=149,F=165,P=115,S=105,T=119,W=204,Y=181,V=117)
    molsize_vector = []
    for i in pep:
        molsize_vector.append(molsize_dict[i.upper()])
    molsize_sum= sum(molsize_vector) - len(pep)*18 + 18; # minus H2O molecule size
    return molsize_sum

def IEDB_immunogenecity(pep): #http://tools.iedb.org/immunogenicity/
    peptide = pep.upper()
    peplen = len(peptide)
    
    #result_list = []
    invalidAA=[]; # report if npeptide exist invalid amino acid out of 20 AA!
    for amino_acid in peptide.strip():
        if not amino_acid.upper() in "ACDEFGHIKLMNPQRSTVWY":
            print("Sequence: '%s' contains an invalid character: '%c' at position %d." %(peptide, amino_acid, peptide.find(amino_acid)))
            invalidAA.append(amino_acid.upper())
    if len(invalidAA)==0:
        # pre-defined scale and weight from IEDB
        immunoscale = {"A":0.127, "C":-0.175, "D":0.072, "E":0.325, "F":0.380, "G":0.110, "H":0.105, "I":0.432, "K":-0.700, "L":-0.036, "M":-0.570, "N":-0.021, "P":-0.036, "Q":-0.376, "R":0.168, "S":-0.537, "T":0.126, "V":0.134, "W":0.719, "Y":-0.012}
        immunoweight = [0.00, 0.00, 0.10, 0.31, 0.30, 0.29, 0.26, 0.18, 0.00]
        #
        if peplen > 9:
            pepweight = immunoweight[:5] + ((peplen - 9) * [0.30]) + immunoweight[5:]
        else:
            pepweight = immunoweight
        #
        cterm = peplen - 1;
        mask_num  = [0, 1, cterm]; # mask at 1,2,9 if peptide len 9
        score = 0;
        count = 0;
        for pos in peptide:
            if count not in mask_num:
                score += pepweight[count] * immunoscale[pos]
                count += 1
            else:
                count += 1
        #result_list.append([peptide, peplen, round(score, 5)])
        score= round(score, 4)
    else:
        score=None
    #
    return score

def netMHCpan_run(id_hla_peps_string): #exp nid-2@HLA-A02:01@AAFDRKSDAK@AAFDYKSDAK
    id = id_hla_peps_string.split('@')[0]
    hla = id_hla_peps_string.split('@')[1]
    mut_pep = id_hla_peps_string.split('@')[2]
    wt_pep = id_hla_peps_string.split('@')[3]
    #
    pepfile= tmp_dir + '/' + id + '_' + mut_pep + '.fa'
    with open(pepfile,'w') as fout:
        fout.write('{}\n{}\n'.format(mut_pep,wt_pep))
    #
    hla_string= hla[:5] + '*' + hla[5:]
    p = subprocess.Popen(['netMHCpan','-BA','-a', hla, '-p ', pepfile], shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()[0];
    #
    netMHCpan_out= p.decode('utf-8');
    netMHCpan_data= netMHCpan_out.splitlines();
    result = [];
    for line in netMHCpan_data:
        if hla_string in line and not 'Protein' in line:
            line_data= [x.strip() for x in line.split(' ') if x != ''];
            if float(line_data[12])!=0: # some very strong signal with 0.000
                result.append([line_data[1],line_data[2],line_data[12],line_data[15]]);
            else:
                result.append([line_data[1],line_data[2],str(round(1- float(line_data[11]),8)),line_data[15]]);
    os.remove(pepfile);
    return result;

def netCTLpan_run(id_hla_pep_string): #exp nid-2@HLA-A02:01@AAFDRKSDAK
    id = id_hla_pep_string.split('@')[0]
    hla = id_hla_pep_string.split('@')[1]
    mut_pep = id_hla_pep_string.split('@')[2]
    #
    pepfile= tmp_dir + '/' + id + '_' + mut_pep + '.fa'
    with open(pepfile,'w') as fout:
        fout.write('{}\n{}\n'.format('>seq1',mut_pep))
    #
    pepLen = len(mut_pep)
    #
    hla_string= hla[:5] + '*' + hla[5:]
    p = subprocess.Popen(['netCTLpan','-a', hla, '-l',str(pepLen),'-f', pepfile], shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()[0];
    #
    netCTLpan_out= p.decode('utf-8');
    netCTLpan_data= netCTLpan_out.splitlines();
    result = [];
    for line in netCTLpan_data:
        if hla_string in line and not 'Protein' in line and not 'NetCTLpan' in line:
            line_data= [x.strip() for x in line.split(' ') if x != ''];
            result.append([line_data[2],line_data[3],line_data[5],line_data[6],line_data[7]]);
    os.remove(pepfile);
    return result;

def getFeatures(id_hla_pep_string): # nid-2@HLA-A02:01@AAFDRKSDAK@AAFDYKSDAK
    string_list= id_hla_pep_string.split('@')
    id = string_list[0]
    hla = string_list[1]
    mut_pep = string_list[2]
    if len(string_list) ==3: # no wt peptide provide, using pepmatch 
        wt_pep = getWTpep(id + '@' + mut_pep)[1]
    elif len(string_list) ==4: # if already have WT peptide information
        wt_pep = string_list[3]
    #
    pep_hydro= hydro_score(mut_pep)
    pep_polar= polarity_score(mut_pep)
    pep_charge= charge_score(mut_pep)
    pep_entropy= entropy_score(mut_pep)
    pep_molsize=  molecule_size(mut_pep)
    #
    netCTL_out= netCTLpan_run(id + '@' + hla + '@' + mut_pep)
    pep_TAP= netCTL_out[0][2];
    pep_Cle= netCTL_out[0][3];
    pep_Comb= netCTL_out[0][4];
    #
    netMHCpan_out= netMHCpan_run(id + '@' + hla + '@' + mut_pep + '@' + wt_pep)
    mut_aff = netMHCpan_out[0][3]
    mut_rank = netMHCpan_out[0][2]
    wt_aff = netMHCpan_out[1][3]
    #
    iedb_immunescore =  IEDB_immunogenecity(mut_pep)
    pep_DAI = getDAI(float(mut_aff),float(wt_aff))
    #
    pep_Rscore = getR(mut_pep,iedb_seq)
    pep_NRP = getNRP(pep_DAI,pep_Rscore)
    
    new_list= [id,hla,mut_pep,wt_pep,pep_hydro,pep_polar,pep_charge,pep_entropy,pep_molsize,pep_TAP,pep_Cle,pep_Comb,mut_aff,mut_rank,iedb_immunescore,pep_DAI,pep_Rscore,pep_NRP]
    #
    current_num=id.split('-')[1];
    if(int(current_num)%100 ==0):
        print(time.strftime("%Y-%m-%d %H:%M:%S"), 'Finish with ' + current_num +' peptides!');
    return new_list

if __name__ == '__main__':
    data = pd.read_csv(inputFile, header=0, sep='\t');
    ################################################the input data format
    #peptide	HLA-A*01:01	HLA-B*44:02	HLA-C*05:01
    #AAMVTGIL	43.667	31	8.349
    #AAVEIRGSV	26.756	12.076	5.51
    #ADSAHPCF	6.891	3.916	1.62
    #AELPPGRPT	17.96	0.829	34.833
    ################################################the input data format
    #
    hla_list= data.columns[1:]
    hla_list= [x.replace('*','') for x in hla_list]
    data_list= data.values.tolist()

    attr_list=[]
    for idx, val in enumerate(data_list):
        id = 'id-' + str(idx + 1)
        pep= val[0]
        pred_val= val[1:]
        hla= hla_list[pred_val.index(min(pred_val))]
        attr_list.append(id + '@' + hla + '@' + pep)
    #
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir);# create a tmp folder
    #
    print(time.strftime("%Y-%m-%d %H:%M:%S") + ": START run!")
    pool=Pool(processes=thread); # defualt thread number setting as 4
    process = pool.map(getFeatures, attr_list); #multiple process
    summaryList = [item for item in process]
    #
    header= ["id","hla","mut_pep","wt_pep","pep_hydro","pep_polar","pep_charge","pep_entropy","pep_molsize","pep_TAP","pep_Cle","pep_Comb","neo_aff","neo_rank","iedb_immunescore","pep_DAI","pep_Rscore","pep_NRP"]
    
    newdata= [header] + summaryList;
    df = pd.DataFrame(newdata[1:],columns=newdata[0]);
    # do prediction step
    #
    ###################################################################
    features=["pep_hydro","pep_polar","pep_charge","pep_entropy","pep_molsize","pep_TAP","pep_Cle","pep_Comb","neo_aff","neo_rank","iedb_immunescore","pep_DAI","pep_Rscore","pep_NRP"];
    with open('/N/slate/cpdong/software/IntronNeoantigen/data/RFmodel.dat', 'rb') as f:
        rfmodel = pickle.load(f)
    ###################################################################
    pred= rfmodel.predict(df[features])
    df["rf_predict"]=pred.tolist()
    df['rf_predict'] = df['rf_predict'].replace({True: 'Positive', False: 'Negative'})
    #
    #
    df.to_csv(outdir + sample + '_predict_result.txt', index=False, sep='\t');
    #
    shutil.rmtree(tmp_dir) # remnove the temp folder
    print(time.strftime("%Y-%m-%d %H:%M:%S") + ": Job complete!")

