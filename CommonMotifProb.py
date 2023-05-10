#!/usr/bin/env python
# coding: utf-8

# In[1]:


from identify_motifs import *
from Bio import Entrez
from Bio import SeqIO
from Bio import AlignIO
import pickle as pickle
from Bio.Align.Applications import ClustalwCommandline
from Bio.Seq import Seq
import pandas as pd
from dna_features_viewer import (GraphicFeature, GraphicRecord,
                                 CircularGraphicRecord)
import matplotlib.pyplot as plt
import numpy as np
import copy


# In[2]:


import subst_matrix as sm


# In[11]:


transmat = sm.GTRSubstMatrix(a=10.628284,
                      b=1.436060,
                      c=1.000000,
                      d=2.351774,
                      e=0.364864,
                      f=9.600893,
                      freqs=[0.214756, 0.219200, 0.274362, 0.291682])
P = transmat.getP(0.07)

#Gives resulting PDICT


# In[4]:


PolyproteinDICT={'C':(1,122),'PrM':(123,290),'E':(291,790), 'NS1':(791,1142),'NS2A':(1143,1368),'NS2B':(1369,1498),'NS3':(1499,2115),'NS4A':(2116,2265),'NS4B':(2266,2516),'NS5':(2517,3419)}


# In[17]:


def importELMs2():
    with open ('elm_classes.tsv', 'r',encoding="utf8", errors='ignore') as doc:
        for i in range(6):#just getting rid of some header lines
            next(doc)

        ELMlist = csv.reader(doc,delimiter='\t')
        return [i for i in ELMlist]


# In[5]:


ELMs=importELMs()
ELMsDICT={}
for i in ELMs:
    if '{' not in i[4] and '$' not in i[4] and '(' not in i[4] and i[4][0] != '^':
        ELMsDICT[i[1]]=i[4]
ELMsDICT['LIG_14-3-3_CanoR_1_C']='R[^DEPG][ST][^PRIKGN]....[VILMFWYP]'
ELMsDICT['DOC_PIKK_1_B']='[DEN][DEN]..[ILMVA][DEN][DEN]L'
ELMsDICT['DEG_SCF_TRCP1_1_B']='DSG..[ST]'


# In[ ]:


PDICT={}
for i in ['A','G','U','C']:
    for j in ['A','G','U','C']:
        PDICT[i+j]=P[i,j]


# In[2]:


Codons = ['AUU', 'AUC', 'AUA','CUU', 'CUC', 'CUA', 'CUG', 'UUA', 'UUG','GUU', 'GUC', 'GUA', 'GUG','AUG','GCU', 'GCC', 'GCA', 'GCG','GGU', 'GGC', 'GGA', 'GGG','CCU', 'CCC', 'CCA', 'CCG','UUU', 'UUC','UGG','UAU', 'UAC','AAU', 'AAC','CAA', 'CAG','UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC','UGU', 'UGC','ACU', 'ACC', 'ACA', 'ACG','AAA', 'AAG','CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG','GAA', 'GAG','GAU', 'GAC','CAU', 'CAC','UAA','UGA','UAG']
CodonsDF = pd.DataFrame(index=Codons, columns=Codons)


# In[4]:


fullPDICT={}
for i in Codons:
    for j in Codons:
        fullPDICT[i+j]=PDICT[i[0]+j[0]]*PDICT[i[1]+j[1]]*PDICT[i[2]+j[2]]


# In[3]:


CodontoAADICT={'AUU': 'ILE','AUC': 'ILE','AUA': 'ILE','CUU': 'LEU','CUC': 'LEU','CUA': 'LEU','CUG': 'LEU','UUA': 'LEU','UUG': 'LEU','GUU': 'VAL','GUC': 'VAL','GUA': 'VAL','GUG': 'VAL','AUG': 'MET','GCU': 'ALA','GCC': 'ALA','GCA': 'ALA','GCG': 'ALA','GGU': 'GLY','GGC': 'GLY','GGA': 'GLY','GGG': 'GLY','CCU': 'PRO','CCC': 'PRO','CCA': 'PRO','CCG': 'PRO','UUU': 'PHE','UUC': 'PHE','UGG': 'TRP','UAU': 'TYR','UAC': 'TYR','AAU': 'ASN','AAC': 'ASN','CAA': 'GLN','CAG': 'GLN','UCU': 'SER','UCC': 'SER','UCA': 'SER','UCG': 'SER','AGU': 'SER','AGC': 'SER','UGU': 'CYS','UGC': 'CYS','ACU': 'THR','ACC': 'THR','ACA': 'THR','ACG': 'THR','AAA': 'LYS','AAG': 'LYS','CGU': 'ARG','CGC': 'ARG','CGA': 'ARG','CGG': 'ARG','AGA': 'ARG','AGG': 'ARG','GAA': 'GLU','GAG': 'GLU','GAU': 'ASP','GAC': 'ASP','CAU': 'HIS','CAC': 'HIS','UAG':'*','UAA':'*','UGA':'*'
}


# In[4]:


DNACodontoAADICT={'ATT': 'ILE','ATC': 'ILE','ATA': 'ILE','CTT': 'LEU','CTC': 'LEU','CTA': 'LEU','CTG': 'LEU','TTA': 'LEU','TTG': 'LEU','GTT': 'VAL','GTC': 'VAL','GTA': 'VAL','GTG': 'VAL','ATG': 'MET','GCT': 'ALA','GCC': 'ALA','GCA': 'ALA','GCG': 'ALA','GGT': 'GLY','GGC': 'GLY','GGA': 'GLY','GGG': 'GLY','CCT': 'PRO','CCC': 'PRO','CCA': 'PRO','CCG': 'PRO','TTT': 'PHE','TTC': 'PHE','TGG': 'TRP','TAT': 'TYR','TAC': 'TYR','AAT': 'ASN','AAC': 'ASN','CAA': 'GLN','CAG': 'GLN','TCT': 'SER','TCC': 'SER','TCA': 'SER','TCG': 'SER','AGT': 'SER','AGC': 'SER','TGT': 'CYS','TGC': 'CYS','ACT': 'THR','ACC': 'THR','ACA': 'THR','ACG': 'THR','AAA': 'LYS','AAG': 'LYS','CGT': 'ARG','CGC': 'ARG','CGA': 'ARG','CGG': 'ARG','AGA': 'ARG','AGG': 'ARG','GAA': 'GLU','GAG': 'GLU','GAT': 'ASP','GAC': 'ASP','CAT': 'HIS','CAC': 'HIS','TAG': '*','TAA': '*','TGA': '*','---': '-'
}


# In[ ]:





# In[5]:


AAtoCodonDICT={
'ILE': ['AUU','AUC', 'AUA'],
'LEU': ['CUU', 'CUC', 'CUA', 'CUG', 'UUA', 'UUG'],
'VAL': ['GUU', 'GUC', 'GUA', 'GUG'],
'MET': ['AUG'],
'ALA': ['GCU', 'GCC', 'GCA', 'GCG'],
'GLY': ['GGU', 'GGC', 'GGA', 'GGG'],
'PRO': ['CCU', 'CCC', 'CCA', 'CCG'],
'PHE': ['UUU', 'UUC'],
'TRP': ['UGG'],
'TYR': ['UAU', 'UAC'],
'ASN': ['AAU', 'AAC'],
'GLN': ['CAA', 'CAG'],
'SER': ['UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC'],
'CYS': ['UGU', 'UGC'],
'THR': ['ACU', 'ACC', 'ACA', 'ACG'],
'LYS': ['AAA', 'AAG'],
'ARG': ['CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
'GLU': ['GAA', 'GAG'],
'ASP': ['GAU', 'GAC'],
'HIS': ['CAU', 'CAC'],
'*': ['UAA','UAG','UGA']}


# In[6]:


SingleLetterDICT={
'I': 'ILE',\
'L': 'LEU',\
'V': 'VAL',\
'M': 'MET',\
'A': 'ALA',\
'G': 'GLY',\
'P': 'PRO',\
'F': 'PHE',\
'W': 'TRP',\
'Y': 'TYR',\
'N': 'ASN',\
'Q': 'GLN',\
'S': 'SER',\
'C': 'CYS',\
'T': 'THR',\
'K': 'LYS',\
'R': 'ARG',\
'E': 'GLU',\
'D': 'ASP',\
'H': 'HIS',\
'*': '*'\
}


# In[7]:


SingleLetterREVDICT={
'ILE': 'I',\
'LEU': 'L',\
'VAL': 'V',\
'MET': 'M',\
'ALA': 'A',\
'GLY': 'G',\
'PRO': 'P',\
'PHE': 'F',\
'TRP': 'W',\
'TYR': 'Y',\
'ASN': 'N',\
'GLN': 'Q',\
'SER': 'S',\
'CYS': 'C',\
'THR': 'T',\
'LYS': 'K',\
'ARG': 'R',\
'GLU': 'E',\
'ASP': 'D',\
'HIS': 'H',\
'*': '*',\
'-': '-'
}


# In[8]:


AminoAcid1LTR={'F','L','I','M','V','S','P','T','A','Y','H','Q','N','K','D','E','C','W','R','S','G','*'}


# In[12]:


stopCodons=['UAG','UGA','UAA']


# In[11]:


def translateRNA(seq):
    prot=''.join([SingleLetterREVDICT[DNACodontoAADICT[seq[i:i+3]]] if '-' not in seq[i:i+3] and 'N' not in seq[i:i+3] else '-' for i in range(0,len(seq),3) if len(seq)-i>2])
    return prot


# In[13]:


def compareCodon(Codon1, Codon2):
    distance = 0
    if Codon1[0]!=Codon2[0]:
        distance+=1
    if Codon1[1]!=Codon2[1]:
        distance+=1
    if Codon1[2]!=Codon2[2]:
        distance+=1
    return distance


# In[ ]:


for i in Codons:
    for j in Codons:
        

        CodonsDF.loc[i,j]=compareCodon(i,j)


# In[14]:


def MotifRegex(motif,Loss=False):#Creates a list of all amino acids for each position from regex eg [[I,J,K],[R,L,I],[F,C,Y,T]]. Limited to fixed length motifs. Can't deal with terminal motifs
    finalMotif=[]
    i=0
    while i<len(motif):
        if motif[i]=='[':
            j=i
            while j!=len(motif):
                if motif[j]==']':
                    currentGroup=motif[i+1:j]
                    if not Loss:
                        if currentGroup[0]=='^':
                            finalMotif.append(list(AminoAcid1LTR-set(currentGroup)))
                        else:
                            finalMotif.append(currentGroup)
                    else:
                        if currentGroup[0]=='^':
                            finalMotif.append(currentGroup)
                        else:
                            finalMotif.append(list(AminoAcid1LTR-set(currentGroup)))
                    i=j
                    break
                j+=1
        elif motif[i]=='.':
            if not Loss:
                finalMotif.append(list(AminoAcid1LTR))
            else:
                finalMotif.append([])
        else:
            if not Loss:
                finalMotif.append(motif[i])
            else:
                finalMotif.append(list(AminoAcid1LTR-set(motif[i])))
        i+=1

    return finalMotif


# In[ ]:


def CountDefinedPos(motif):#Counts number of positions that accept single amino acids
    finalMotif=[]
    i=0
    definedpositioncount=0
    while i<len(motif):
        
        if motif[i]=='[':
            j=i
            while j!=len(motif):
                if motif[j]==']':
                    currentGroup=motif[i+1:j]
                    if currentGroup[0]=='^':
                        finalMotif.append(list(AminoAcid1LTR-set(currentGroup)))
                    elif len(currentGroup)<=4:
                        definedpositioncount+=1
                    i=j
                    break
                j+=1
        elif motif[i]=='.':
            finalMotif.append(list(AminoAcid1LTR))
        else:
            definedpositioncount+=1
        i+=1

    return definedpositioncount


# In[15]:


def MotifProbability(seq,motif,branchDist):#Motifprobability gain
    
    #probs for influenza
    transmat = sm.GTRSubstMatrix(a=10.628284,
                          b=1.436060,
                          c=1.000000,
                          d=2.351774,
                          e=0.364864,
                          f=9.600893,
                          freqs=[0.214756, 0.219200, 0.274362, 0.291682])
    P = transmat.getP(branchDist)

 #probs for zika first instance   
#    transmat = sm.GTRSubstMatrix(a=29.406208,
#                          b=1.543007,
#                          c=1.000000,
#                          d=1.195367,
#                          e=0.433305,
#                          f=10.008729,
#                          freqs=[0.235079, 0.212114, 0.306874, 0.245933])
#    P = transmat.getP(branchDist)
  
    motif = MotifRegex(motif)
    AllCodons=[]
    for i in motif:
        CodonsForAminoacids = [AAtoCodonDICT[SingleLetterDICT[f]] for f in i] #list of all codons for each amino acid at each position
        CodonsForAminoacids = [item for sublist in CodonsForAminoacids for item in sublist]#retarded way of making it one list instead of list of lists
        AllCodons.append(CodonsForAminoacids)
    
    compoundProbability=[0 for i in range(len(motif))]
    probabilityList=[]

    i=0
    motifsize=len(motif)*3
    while i<len(seq)-motifsize+1:

        currentMotif= seq[i:i+motifsize]
        i+=3
        j=0
        count=0
        while j<len(currentMotif)-3+1:
            seqcodon= currentMotif[j]+currentMotif[j+1]+currentMotif[j+2]
            #currentcodonAA = CodontoAADICT[codon]

            if seqcodon in stopCodons:
                compoundProbability[count]=0
                #print codon, currentcodonAA, motif[count], compoundProbability[count]
                
                
                if count==len(motif)-1:
                    #print '|||||||||||||',compoundProbability
                    probabilityList.append(np.prod(compoundProbability))
                
                count+=1
                j+=3
                
                continue


            #elif SingleLetterREVDICT[CodontoAADICT[seqcodon]] in motif[count]:
            #    compoundProbability[count]=1
            #    
            #   # print codon, currentcodonAA, motif[count], compoundProbability[count]
            #    
            #    if count==len(motif)-1:
            #        #print '|||||||||||||',compoundProbability
            #        probabilityList.append(np.prod(compoundProbability))
            #
            #    count+=1
            #    j+=3
            #    continue

            else:
                relevantcodons=AllCodons[count]
                tempCurrentPosProbList=[]
                for motifCodon in relevantcodons:
                    codonprobability=P[seqcodon[0],motifCodon[0]]*P[seqcodon[1],motifCodon[1]]*P[seqcodon[2],motifCodon[2]]
                    tempCurrentPosProbList.append(codonprobability)
                compoundProbability[count]=sum(tempCurrentPosProbList)
                
                #countones=len([cdon for cdon in relevantcodons if CodonsDF.loc[seqcodon,cdon]==1 ])
                #counttwos=len([cdon for cdon in relevantcodons if CodonsDF.loc[seqcodon,cdon]==2 ])
                #countthrees=len([cdon for cdon in relevantcodons if CodonsDF.loc[seqcodon,cdon]==3 ])
                #compoundProbability[count] = (mutationRate*countones/3.0+mutationRate**2*counttwos/9.0+mutationRate**3*countthrees/27.0) #3 depends on for each instance which number of positions yield a favourable mutations, could be 1-3
                #print codon, SingleLetterREVDICT[currentcodonAA], motif[count], compoundProbability[count], relevantcodons,countones,counttwos

            if count==len(motif)-1:
                #print '|||||||||||||',compoundProbability
                probabilityList.append(np.prod(compoundProbability))
            
            #print codon, currentcodonAA, motif[count], compoundProbability[count]
            count+=1
            j+=3
        
    return probabilityList


# In[15]:


def MotifProbabilityCodonList(seq,motif,P):#Motifprobability gain
    
    #probs for influenza (takes ages running this every function call) use P as argument instead?
    '''
    transmat = sm.GTRSubstMatrix(a=10.628284,
                          b=1.436060,
                          c=1.000000,
                          d=2.351774,
                          e=0.364864,
                          f=9.600893,
                          freqs=[0.214756, 0.219200, 0.274362, 0.291682])
    P = transmat.getP(branchDist)
    '''

 #probs for zika first instance   
#    transmat = sm.GTRSubstMatrix(a=29.406208,
#                          b=1.543007,
#                          c=1.000000,
#                          d=1.195367,
#                          e=0.433305,
#                          f=10.008729,
#                          freqs=[0.235079, 0.212114, 0.306874, 0.245933])
#    P = transmat.getP(branchDist)
  
    motif = MotifRegex(motif)
    AllCodons=[]
    for i in motif:
        CodonsForAminoacids = [AAtoCodonDICT[SingleLetterDICT[f]] for f in i] #list of all codons for each amino acid at each position
        CodonsForAminoacids = [item for sublist in CodonsForAminoacids for item in sublist]#retarded way of making it one list instead of list of lists
        AllCodons.append(CodonsForAminoacids)
    
    compoundProbability=[0 for i in range(len(motif))]
    probabilityList=[]

    i=0
    motifsize=len(motif)*3
    while i<len(seq)-motifsize+1:

        currentMotif= seq[i:i+motifsize]
        i+=3
        j=0
        count=0
        while j<len(currentMotif)-3+1:
            seqcodon= currentMotif[j]+currentMotif[j+1]+currentMotif[j+2]
            #currentcodonAA = CodontoAADICT[codon]

            if seqcodon in stopCodons:
                compoundProbability[count]=0
                #print codon, currentcodonAA, motif[count], compoundProbability[count]
                
                
                if count==len(motif)-1:
                    #print '|||||||||||||',compoundProbability
                    probabilityList.append(np.prod(compoundProbability))
                
                count+=1
                j+=3
                
                continue


            #elif SingleLetterREVDICT[CodontoAADICT[seqcodon]] in motif[count]:
            #    compoundProbability[count]=1
            #    
            #   # print codon, currentcodonAA, motif[count], compoundProbability[count]
            #    
            #    if count==len(motif)-1:
            #        #print '|||||||||||||',compoundProbability
            #        probabilityList.append(np.prod(compoundProbability))
            #
            #    count+=1
            #    j+=3
            #    continue

            else:
                relevantcodons=AllCodons[count]
                tempCurrentPosProbList=[]
                for motifCodon in relevantcodons:
                    codonprobability=P[seqcodon[0],motifCodon[0]]*P[seqcodon[1],motifCodon[1]]*P[seqcodon[2],motifCodon[2]]
                    tempCurrentPosProbList.append(codonprobability)
                    
                compoundProbability[count]=sum(tempCurrentPosProbList)
                
                #countones=len([cdon for cdon in relevantcodons if CodonsDF.loc[seqcodon,cdon]==1 ])
                #counttwos=len([cdon for cdon in relevantcodons if CodonsDF.loc[seqcodon,cdon]==2 ])
                #countthrees=len([cdon for cdon in relevantcodons if CodonsDF.loc[seqcodon,cdon]==3 ])
                #compoundProbability[count] = (mutationRate*countones/3.0+mutationRate**2*counttwos/9.0+mutationRate**3*countthrees/27.0) #3 depends on for each instance which number of positions yield a favourable mutations, could be 1-3
                #print codon, SingleLetterREVDICT[currentcodonAA], motif[count], compoundProbability[count], relevantcodons,countones,counttwos

            if count==len(motif)-1:
                #print '|||||||||||||',compoundProbability
                assert len(compoundProbability)==len(AllCodons)
                probabilityList.append(compoundProbability)
            
            #print codon, currentcodonAA, motif[count], compoundProbability[count]
            count+=1
            j+=3
        
    return probabilityList


# In[15]:


def MotifProbabilityCodonListV2(seq,motif,P):#Motifprobability gain
    
    #probs for influenza (takes ages running this every function call) use P as argument instead?

    currentMotif=seq
    motif = MotifRegex(motif)

    compoundProbability=[0 for i in range(len(motif))]
    probabilityList=[]
    for i in range(len(motif)):
        CodonsForAminoacids = [AAtoCodonDICT[SingleLetterDICT[f]] for f in motif[i]] #list of all codons for each amino acid at each position
        CodonsForAminoacids = [item for sublist in CodonsForAminoacids for item in sublist]#retarded way of making it one list instead of list of lists


        seqcodon= currentMotif[i*3]+currentMotif[(i*3)+1]+currentMotif[(i*3)+2]
        
        tempCurrentPosProb=0
        for motifCodon in CodonsForAminoacids:
            codonprobability=fullPDICT[seqcodon+motifCodon]
            tempCurrentPosProb+=codonprobability

        compoundProbability[i]=tempCurrentPosProb
        
    return compoundProbability


# In[17]:


#currentseq= SeqIO.read('NCBI_Zika/'+'AY632535.2', "genbank")#AY632535.2
#ORF=readingframe(currentseq)
#problem with indexing since I am using PolyproteinDICT for ranges. WIll give me ranges not based on actual refseq index
#but actually uniprot index

#RNASeq=str(currentseq.seq[ORF[0]:ORF[1]].transcribe())

#probbs=MotifProbability(RNASeq,ELMsDICT['LIG_LIR_Nem_3_C'])


# In[ ]:


'''
            elif len(motifregex[count])==0: #if current position can accept any AA, probabiliity of loss is 0
                compoundProbability[count]=0

               # print codon, currentcodonAA, motif[count], compoundProbability[count]

                if count==len(motifregex)-1:
                    #print '|||||||||||||',compoundProbability
                    entry.probabilityLoss=(np.sum(compoundProbability))

                count+=1
                j+=3
                continue

            else:
                relevantcodons=AllCodons[count]
                countones=len([cdon for cdon in relevantcodons if CodonsDF.loc[codon,cdon]==1 ])
                counttwos=len([cdon for cdon in relevantcodons if CodonsDF.loc[codon,cdon]==2 ])
                countthrees=len([cdon for cdon in relevantcodons if CodonsDF.loc[codon,cdon]==3 ])
                compoundProbability[count] = (mutationRate*countones/3.0+mutationRate**2*counttwos/9.0+mutationRate**3*countthrees/27.0) #3 depends on for each instance which number of positions yield a favourable mutations, could be 1-3
                #print codon, SingleLetterREVDICT[currentcodonAA], motif[count], compoundProbability[count], relevantcodons,countones,counttwos
                
'''    

def MotifProbabilityLoss(RefStrain, motiflist, ELMsDICT,branchDist):
    #best way is to just take the list of motifs lost (motif objects) and add attribute loss probability and return a list of 
    #motif objects with this added attribute. Only thing I need to fix is the motif regex by importing ELMDICT as well.
    
    #motiflist must be based on RefStrain sequence amino acid seq
    if isinstance(RefStrain, Seq):
        currentseq=RefStrain
    else:
        currentseq= SeqIO.read('NCBI_Zika/'+RefStrain, "genbank")#AY632535.2
        ORF=readingframe(currentseq)
        currentseq=currentseq.seq[ORF[0]:ORF[1]]
#problem with indexing since I am using PolyproteinDICT for ranges. WIll give me ranges not based on actual refseq index
#but actually uniprot index
    #RNASeq=str(currentseq.seq[ORF[0]+(i.position[0]*3):ORF[0]+(i.position[1]*3)].transcribe())

    transmat = sm.GTRSubstMatrix(a=29.406208,
                          b=1.543007,
                          c=1.000000,
                          d=1.195367,
                          e=0.433305,
                          f=9.600893,
                          freqs=[0.214756, 0.219200, 0.274362, 0.291682])
    P = transmat.getP(branchDist)
    
    
    
    
    for entry in motiflist:
    #motif = MotifRegexLoss(motif)#have to make special motifRegex or add choice to normalregex
        AllCodons=[]
        #for i in motif:
        try:
            motifregex=MotifRegex(ELMsDICT[entry.motifID])
        except KeyError:
            print(entry.motifID, 'Not in ELMsDICT, probability cant be determined')
            continue
        for i in motifregex:
            CodonsForAminoacids = [AAtoCodonDICT[SingleLetterDICT[f]] for f in i]
            CodonsForAminoacids = [item for sublist in CodonsForAminoacids for item in sublist]
            AllCodons.append(CodonsForAminoacids)


        compoundProbability=[0 for i in range(len(entry.motifSeq))]
        probabilityList=[]
        entry.AllCodons=AllCodons

        #i=0
        #motifsize=len(entry.motifSeq)*3
        #while i<len(seq)-motifsize+1:
        currentMotif=str(currentseq[entry.position[0]*3:entry.position[1]*3].transcribe())
        #currentMotif= seq[i:i+motifsize]
        #i+=3
        j=0
        count=0
        assert len(AllCodons)==len(entry.motifSeq)
        while j<len(currentMotif)-3+1:
            seqcodon= currentMotif[j]+currentMotif[j+1]+currentMotif[j+2]
            #currentcodonAA = CodontoAADICT[codon]

            if seqcodon in stopCodons:
                compoundProbability[count]=0
                #print codon, currentcodonAA, motif[count], compoundProbability[count]


                if count==len(motifregex)-1:
                    #print '|||||||||||||',compoundProbability
                    entry.probabilityLoss=(np.sum(compoundProbability))

                count+=1
                j+=3

                continue
          
            relevantcodons=AllCodons[count]
            
            tempCurrentPosProbList=[]
            
            for motifCodon in relevantcodons:
                codonprobability=P[seqcodon[0],motifCodon[0]]*P[seqcodon[1],motifCodon[1]]*P[seqcodon[2],motifCodon[2]]
                tempCurrentPosProbList.append(codonprobability)
            compoundProbability[count]=sum(tempCurrentPosProbList)

                #countones=len([cdon for cdon in relevantcodons if CodonsDF.loc[seqcodon,cdon]==1 ])
                #counttwos=len([cdon for cdon in relevantcodons if CodonsDF.loc[seqcodon,cdon]==2 ])
                #countthrees=len([cdon for cdon in relevantcodons if CodonsDF.loc[seqcodon,cdon]==3 ])
                #compoundProbability[count] = (mutationRate*countones/3.0+mutationRate**2*counttwos/9.0+mutationRate**3*countthrees/27.0) #3 depends on for each instance which number of positions yield a favourable mutations, could be 1-3
                #print codon, SingleLetterREVDICT[currentcodonAA], motif[count], compoundProbability[count], relevantcodons,countones,counttwos

            if count==len(motifregex)-1:
                print('|||||||||||||',compoundProbability)
                entry.probabilityLoss=(1-np.prod(compoundProbability))
            
            #print codon, currentcodonAA, motif[count], compoundProbability[count]
            count+=1
            j+=3
            
            
                
        
    return motiflist


# In[16]:


#adding len feature check makes it work for both zika and flu for now
def readingframe(NCBIrecord):
    for feature in NCBIrecord.features:
        if feature.type=='CDS' and len(feature.location.parts)!=2:
            return (feature.location.start,feature.location.end)


# In[ ]:


def SeqMatch(seq1,seq2):
    distance=0
    try:
        for i in range(len(seq1)):
            if seq1[i]!=seq2[i] and seq1[i] not in ['X','N'] and seq2[i] not in ['X','N']:
                distance+=1
    except IndexError:
        print('Sequences different length, returning match value up until len seq1= %i' % (len(seq1)))
        return distance
    
    return distance


# ## Two functions to identify new motifs that accept different input, either NCBI records or a list of motifs and NCBI strain IDs. They compare motifs by specific sequence, exact position and codon sequence similarity to determine if the motifs are likely to be completely new. Only accurate for relatively similar sequences as is the case in zika. Deals with framshifts easily but would misidentify motifs as new if big deletions or fusions have occurred, shifting the relative motif position.

# In[ ]:


########24/05/2018##########
#wont work with ancestral reconstructed sequences as of now, only NCBI records. But this is no longer in use, see NewMotifFind3
#for current function
def NewMotifFind(ancestralRecord, derivedRecord):
    motifList = []   
    identifier = 'Zika_Polyprotein,H8XX12_ZIKV'
    #define reading frame and sequences
    frameAncest = readingframe(ancestralRecord)
    sequence1=str(ancestralRecord.seq[frameAncest[0]:frameAncest[1]].translate())
    frameDerived = readingframe(derivedRecord)
    sequence2=str(derivedRecord.seq[frameDerived[0]:frameDerived[1]].translate())

    
    #find motifs in the sequences
    ancestralSeq= AllMotifs(identifier=identifier,sequence=sequence1[:-1],textfile='zika2polyprotein',iupredOp=None)


    derivedSeq = AllMotifs(identifier=identifier,sequence=sequence2[:-1],textfile='zika2polyprotein',iupredOp=None)


    for i in derivedSeq.motifs:
        FoundMatch=False

        for j in ancestralSeq.motifs:
            if i.motifID==j.motifID:
                if i.motifSeq==j.motifSeq:
                    if i.position==j.position:
                        #Definitely same motif
                        FoundMatch=True
                        break
                    elif abs(i.position[0]-j.position[0])<10:
                        #Very likely same motif
                        FoundMatch=True
                        break

                    elif abs(i.position[0]-j.position[0])>40:
                        a=1
                        #no match

                    elif len(ancestralRecord.seq[frameAncest[0]+i.position[0]*3:frameAncest[0]+i.position[1]*3])>9 and SeqMatch(str(ancestralRecord.seq[frameAncest[0]+i.position[0]*3:frameAncest[0]+i.position[1]*3]), str(derivedRecord.seq[frameDerived[0]+j.position[0]*3:frameDerived[0]+j.position[1]*3]))<3:

                        #Very likely true
                        FoundMatch=True
                        break

                    #else:
                        #not the same instance
                elif len(i.motifSeq)!=len(j.motifSeq):
                    FoundMatch = False
                elif abs(i.position[0]-j.position[0])<10:
                    #Very likely same motif
                    FoundMatch=True
                    break

                elif abs(i.position[0]-j.position[0])>40:
                    a=1
                    #no match
                elif SeqMatch(str(ancestralRecord.seq[frameAncest[0]+i.position[0]*3:frameAncest[0]+i.position[1]*3]), str(derivedRecord.seq[frameDerived[0]+j.position[0]*3:frameDerived[0]+j.position[1]*3]))<3 and len(ancestralRecord.seq[frameAncest[0]+i.position[0]*3:frameAncest[0]+i.position[1]*3])>9:
                    #Very likely true
                    FoundMatch=True
                    break
                #else:
                    #not the same instance
            #else:
                #not the same for sure

        if FoundMatch==False:
            motifList.append(i)
            
    return motifList


# In[ ]:


########24/05/2018##########
#wont work with ancestral reconstructed sequences as of now, only NCBI records. But this is no longer in use, see NewMotifFind3
#for current function

def NewMotifFind2(motiflist1,motiflist2,NCBIId1,NCBIId2):
    motifList = [] 
    
    NCBISeq1=SeqIO.read('NCBI_Zika/'+NCBIId1, "genbank")
    NCBISeq2=SeqIO.read('NCBI_Zika/'+NCBIId2, "genbank")
    frameId1 = readingframe(NCBISeq1)
    frameId2 = readingframe(NCBISeq2)
    for i in motiflist2.motifs:
        FoundMatch=False

        for j in motiflist1.motifs:
            if i.motifID==j.motifID:
                if i.motifSeq==j.motifSeq:
                    if i.position==j.position:
                        #Definitely same motif
                        FoundMatch=True
                        break
                    elif abs(i.position[0]-j.position[0])<10:
                        #Very likely same motif
                        FoundMatch=True
                        break

                    elif abs(i.position[0]-j.position[0])>40:
                        a=1
                        #no match

                    elif len(NCBISeq1.seq[frameId1[0]+i.position[0]*3:frameId1[0]+i.position[1]*3])>9 and SeqMatch(str(NCBISeq1.seq[frameId1[0]+i.position[0]*3:frameId1[0]+i.position[1]*3]), str(NCBISeq2.seq[frameId2[0]+j.position[0]*3:frameId2[0]+j.position[1]*3]))<3:

                        #Very likely true
                        FoundMatch=True
                        break

                    #else:
                        #not the same instance
                elif len(i.motifSeq)!=len(j.motifSeq):
                    FoundMatch = False
                elif abs(i.position[0]-j.position[0])<10:
                    #Very likely same motif
                    FoundMatch=True
                    break

                elif abs(i.position[0]-j.position[0])>40:
                    a=1
                    #no match
                elif SeqMatch(str(NCBISeq1.seq[frameId1[0]+i.position[0]*3:frameId1[0]+i.position[1]*3]), str(NCBISeq2.seq[frameId2[0]+j.position[0]*3:frameId2[0]+j.position[1]*3]))<3 and len(NCBISeq1.seq[frameId1[0]+i.position[0]*3:frameId1[0]+i.position[1]*3])>9:
                    #Very likely true
                    FoundMatch=True
                    break
                #else:
                    #not the same instance
            #else:
                #not the same for sure

        if FoundMatch==False:
            motifList.append(i)
            
    return motifList


# In[ ]:


#improved function to determine if motif A and motif B from strain X and strain Y are the same motif or if motif B is a new motif
#evolved in strain Y. Incorporated the findAlignPostn function to compare the aligned sequences to say for certain what the
#relevant position in the strain is.
def NewMotifFind3(motiflist1,motiflist2,NCBIId1,NCBIId2,alignment):
    motifList = [] 

    #NCBISeq1=SeqIO.read('NCBI_Zika/'+NCBIId1, "genbank")
    #NCBISeq2=SeqIO.read('NCBI_Zika/'+NCBIId2, "genbank")
    #frameId1 = readingframe(NCBISeq1)
    #frameId2 = readingframe(NCBISeq2)
    for querymotif in motiflist2:
        
        FoundMatch=False
        closeNeig=0
        for refmotif in motiflist1:
            if querymotif.motifID==refmotif.motifID and len(querymotif.motifSeq)==len(refmotif.motifSeq) and SeqMatch(querymotif.motifSeq,refmotif.motifSeq)<=1:
                #find the actual index in the MSA of a given (real) sequence position
                try:
                    AdjustedIndexqueryStart=findAlignPostn(NCBIId2,querymotif.position[0]+1,alignment)
                    AdjustedIndexrefStart=findAlignPostn(NCBIId1,refmotif.position[0]+1,alignment)
                except:
                    print('#######ERROR#####',querymotif.motifID,'Cant be read')
                

                #AdjustedIndex1End=findAlignPostn(NCBIId1,j.position[1],alignment)
                #in the MSA all indexes refer to the same position per definition, so we convert the common index from above
                #and find out in our second sequence which real position that index refers to, and we can then compare
                #with the position of the motif in that sequence.
                Position2alignedwith1=findAlignPostn(NCBIId2,AdjustedIndexrefStart, alignment, refIndex=False)
                #if querymotif.motifID=='CLV_C14_Caspase3-7' :
                    #print 'pos new',querymotif.position,querymotif.motifSeq,'pos anc recon',refmotif.position, refmotif.motifSeq,'adjusted new',Position2alignedwith1
                    #if abs(Position2alignedwith1-querymotif.position[0]+1)<2:
                        #print 'Adjusted: ',Position2alignedwith1, 'LIG_PDZ_Class_2_B position start: ',querymotif.position[0]+1, 'RefMotif:', refmotif.position
                if Position2alignedwith1==querymotif.position[0]+1:#then the motifs are actually aligned and in the same position in the protein.
                    #Same motif, same sequence and same (aligned) position if we get to this point. Means definitely the same motif.
                    #what what margins should I allow to say that somehing is a new motif?
                    FoundMatch=True
                    #if querymotif.motifID=='LIG_PDZ_Class_2_B':
                        #print 'matches', refmotif.motifID,'at',refmotif.position[0],refmotif.position[1], refmotif.motifSeq
                    break
                    #maybe omit motifs if aligned sequence position is different but if sequence position is very near, write
                    #out the motif instance to allow manual asignment.

                #if the instance amino acid sequence is not the same (but it is the same type of motif) it could be a silent mutation
                #(which we want to call the same motif) or it could be another instance. Which we want to call not a match.
                #So first thing to ask might be the aligned position. Second thing might be the similarity of the NT seq.
                #Should not first check the length. A silent mutation could lead to a shortened form of a motif (rare).
                
                #Need to consider two extreme cases. Eg. a short prevalent motif in many different copies. And a new motif evolves
                #overlapping an old motif of the same type. How can I make sure this motif gets classified as a new motif.
                #elif abs(int(Position2alignedwith1)-int(j.position[0]+1))<len(j.motifSeq):
                #    print '\n','\t',i.motifID,'at', i.position[0],i.position[1], i.motifSeq
                #    print '###','Ambiguous case', j.motifSeq, j.position[0], j.position[1]
                #else:
                #    if closeNeig:
                #        if closeNeig>abs(int(Position2alignedwith1)-int(j.position[0]+1)):
                #            closeNeig=abs(int(Position2alignedwith1)-int(j.position[0]+1))
                #            closeN=j
                #    else:
                #        closeNeig=abs(int(Position2alignedwith1)-int(j.position[0]+1))
                #        closeN=j

        if not FoundMatch:
            #print '\n','\t',i.motifID,'at', i.position[0],i.position[1], i.motifSeq
            #print 'found no match'
            #print 'Closest neighbour',closeN.motifID, closeN.position[0],closeN.position[1],closeN.motifSeq
            motifList.append(querymotif)
            #print 'saved', querymotif.motifID, querymotif.position[0],querymotif.motifSeq
            
    return motifList


# In[ ]:


#improved function to determine if motif A and motif B from strain X and strain Y are the same motif or if motif B is a new motif
#evolved in strain Y. Incorporated the findAlignPostn function to compare the aligned sequences to say for certain what the
#relevant position in the strain is.
'''
def NewMotifFind3(motiflist1,motiflist2,NCBIId1,NCBIId2,alignment):
    motifList = [] 
    
    NCBISeq1=SeqIO.read('NCBI_Zika/'+NCBIId1, "genbank")
    NCBISeq2=SeqIO.read('NCBI_Zika/'+NCBIId2, "genbank")
    frameId1 = readingframe(NCBISeq1)
    frameId2 = readingframe(NCBISeq2)
    for i in motiflist2.motifs:
        
        FoundMatch=False
        closeNeig=0
        for j in motiflist1.motifs:
            if i.motifID==j.motifID and i.motifSeq==j.motifSeq:
                #find the actual index in the MSA of a given (real) sequence position
                AdjustedIndex1Start=findAlignPostn(NCBIId1,i.position[0]+1,alignment)
                AdjustedIndex1End=findAlignPostn(NCBIId1,i.position[1],alignment)
                #in the MSA all indexes refer to the same position per definition, so we convert the common index from above
                #and find out in our second sequence which real position that index refers to, and we can then compare
                #with the position of the motif in that sequence.
                Position2alignedwith1=findAlignPostn(NCBIId2,AdjustedIndex1Start, alignment, refIndex=False)
                if Position2alignedwith1==j.position[0]+1:#then the motifs are actually aligned and in the same position in the protein.
                    #Same motif, same sequence and same (aligned) position if we get to this point. Means definitely the same motif.
                    #what what margins should I allow to say that somehing is a new motif?
                    FoundMatch=True
                    #print 'matches', j.motifID,'at',j.position[0],j.position[1], j.motifSeq
                    break
                    #maybe omit motifs if aligned sequence position is different but if sequence position is very near, write
                    #out the motif instance to allow manual asignment.

                #if the instance amino acid sequence is not the same (but it is the same type of motif) it could be a silent mutation
                #(which we want to call the same motif) or it could be another instance. Which we want to call not a match.
                #So first thing to ask might be the aligned position. Second thing might be the similarity of the NT seq.
                #Should not first check the length. A silent mutation could lead to a shortened form of a motif (rare).
                
                #Need to consider two extreme cases. Eg. a short prevalent motif in many different copies. And a new motif evolves
                #overlapping an old motif of the same type. How can I make sure this motif gets classified as a new motif.
                elif abs(int(Position2alignedwith1)-int(j.position[0]+1))<len(j.motifSeq):
                    print '\n','\t',i.motifID,'at', i.position[0],i.position[1], i.motifSeq
                    print '###','Ambiguous case', j.motifSeq, j.position[0], j.position[1]
                else:
                    if closeNeig:
                        if closeNeig>abs(int(Position2alignedwith1)-int(j.position[0]+1)):
                            closeNeig=abs(int(Position2alignedwith1)-int(j.position[0]+1))
                            closeN=j
                    else:
                        closeNeig=abs(int(Position2alignedwith1)-int(j.position[0]+1))
                        closeN=j

        if not FoundMatch:
            print '\n','\t',i.motifID,'at', i.position[0],i.position[1], i.motifSeq
            print 'found no match'
            print 'Closest neighbour',closeN.motifID, closeN.position[0],closeN.position[1],closeN.motifSeq
            motifList.append(i)
            
    return motifList
'''


# In[ ]:


#not generalised, can only be used for zika in current state
def FilterNCBIUnknown(NCBI_IDList,directory):
    filteredStrainList=[]
    for i in NCBI_IDList:
        hasReadingframe=False
        if i!='KR872956.1': #mislaigned readingframe
            currentseq= SeqIO.read(directory+'/'+i, "genbank")
            for feature in currentseq.features:
                if feature.type=='CDS':
                    hasReadingframe=True
            if hasReadingframe:
                currentreadingframe=readingframe(currentseq)
                if abs(len(currentseq.seq)-10779)<500: #compare to len of old strain to filter out non-zika
                    countNs=0
                    try:
                        for j in currentseq.seq[currentreadingframe[0]:currentreadingframe[1]]:
                            if j=='N':
                                countNs+=1
                        if countNs < len(currentseq.seq[currentreadingframe[0]:currentreadingframe[1]])/20: #unknowns can be max 5% of seq len
                            filteredStrainList.append((i,countNs))
                    except:
                        print('Error', i)
    return filteredStrainList


# In[ ]:


def FilterNCBIUnknownFLU(NCBI_IDList,directory):
    filteredStrainList=[]
    for i in NCBI_IDList:

        currentseq= SeqIO.read(directory+'/'+i, "genbank")
        for feature in currentseq.features:
            if feature.type=='CDS' and len(feature.location.parts)!=2:

                currentreadingframe=(feature.location.start,feature.location.end)

                countNs=0
                try:
                    for j in currentseq.seq[currentreadingframe[0]:currentreadingframe[1]]:
                        if j=='N':
                            countNs+=1
                    if countNs < len(currentseq.seq[currentreadingframe[0]:currentreadingframe[1]])/20: #unknowns can be max 5% of seq len
                        filteredStrainList.append((i,countNs))
                except:
                    print('Error', i)
    return filteredStrainList


# In[ ]:


#problematic when filtering with thresholds above identical ie removing all sequences that have 
#fewer than 4 substitutions. This is because I will group all sequences within this range in 
#one tempsamelist and add to checkedlist and then only chose one to return, however the outcome
#depends very much on starting sequence in these cases eg if sequence 1 is within 4 subs with
# seq 2, 3 and 4 only 1 will be kept, but seq 4 might only be within 4 of seq 1 so if 4 starts 2, 3 
#4 is kept instead inreasing the sequence variety in the final sample.
def FilterNCBISeqs(NCBI_IDList, directory='NCBI_Zika',seqMatchThreshold=4):
    if directory=='NCBI_Zika':
        filteredStrainList=FilterNCBIUnknown(NCBI_IDList,directory)
    else:
        filteredStrainList=FilterNCBIUnknownFLU(NCBI_IDList,directory)
    tempfilter=[]
    checkedlist=[]
    for i in filteredStrainList:
        tempsamelist=[]
        tempsamelist.append(i)
        if i not in checkedlist:

            primaryseq= SeqIO.read(directory+'/'+i[0], "genbank")
            primaryreadingframe=readingframe(primaryseq)
            sequence1=str(primaryseq.seq[primaryreadingframe[0]:primaryreadingframe[1]].translate())
            

            checkedlist.append(i)
            for j in filteredStrainList:
                if j not in checkedlist:



                    secondaryseq= SeqIO.read(directory+'/'+j[0], "genbank")
                    secondaryreadingframe=readingframe(secondaryseq)
                    sequence2=str(secondaryseq.seq[secondaryreadingframe[0]:secondaryreadingframe[1]].translate())
                    
                    if len(sequence1)==len(sequence2):
                        if SeqMatch(sequence1,sequence2)<seqMatchThreshold:
                            tempsamelist.append(j)
                            checkedlist.append(j)
            for k in tempsamelist:
                if k[1]==0:
                    minN=k
                    break
                elif k == tempsamelist[0]:
                    minN=k
                elif k[1] < minN[1]:
                    minN=k
            tempfilter.append(minN[0])
    return tempfilter


# In[ ]:


#read in TM predictions from file (ouput from TMHMM) into dictinary with keys being NCBI entries and for each a list of features + positions
def TMHMMRead(filename):
    TMMpredDICT={}
    newEntry=True
    with open (filename, 'r') as doc:
        TMMpred = csv.reader(doc,delimiter='\t')

        for i in TMMpred:

            if '#' not in i[0]:
                loc=i[2]
                start=i[3].split()[0]
                stop=i[3].split()[1]
                if newEntry:
                    TMMpredDICT[i[0]]=[]
                    TMMpredDICT[i[0]].append([loc,start,stop])
                    newEntry=False
                else:
                    TMMpredDICT[i[0]].append([loc,start,stop])

            else:
                newEntry=True
    return TMMpredDICT


# In[ ]:


#Takes an (string) NCBI code entry or a Bio.Seq.Seq sequence
#changed folder for flu purposes temporarily (for zika use NCBI_Zika/)
def MotifFind(NCBIEntry):
    if isinstance(NCBIEntry,Seq):
        sequence1=str(NCBIEntry.translate())
    else:
        NCBISeq=SeqIO.read('NCBI_NS1_H7N7/'+NCBIEntry, "genbank")
        frameQ = readingframe(NCBISeq)
        sequence1=str(NCBISeq.seq[frameQ[0]:frameQ[1]].translate())
        
    identifier = 'Zika_Polyprotein,H8XX12_ZIKV'
    #define reading frame and sequences

    MotifList= AllMotifs(identifier=identifier,sequence=sequence1[:-1],textfile='zika2polyprotein',iupredOp=None)
    return MotifList


# In[ ]:


#Strainlist can be either list of NCBI codes or ancestral reconstructed sequences of Bio.Seqrecord.SeqRecord format or both mixed
def NCBIstrainMotifLoad(StrainList,filename):
    Nofile=False
    try:
        f = open(filename+'.pkl','rb')
        print('Loaded '+filename+'.pkl',' successfully')
    except IOError:
        print('No such file, creating new file')
        Nofile=True
    except:
        print('Unknown error. Scary')
    else:
        MotifDICT=pickle.load(f)
        for i in StrainList:
            if isinstance(i,SeqIO.SeqRecord):
                if i.id not in list(MotifDICT.keys()):
                    MotifDICT[i.id]=MotifFind(i.seq)
                    print('Added %s to DICT in loaded file %s' % (i.id,filename))
            else:
                if i not in list(MotifDICT.keys()):
                    MotifDICT[i]=MotifFind(i)
                    print('Added %s to DICT in loaded file %s' % (i,filename))
                
        f.close()
        return MotifDICT
    
    if Nofile:
        MotifDICT={}
        for i in StrainList:
            if isinstance(i,SeqIO.SeqRecord):
                MotifDICT[i.id]=MotifFind(i.seq)
            else:
                MotifDICT[i]=MotifFind(i)
        with open(filename+'.pkl','wb') as f:
            pickle.dump(MotifDICT, f, protocol=pickle.HIGHEST_PROTOCOL)
        return MotifDICT


# In[ ]:


#remove all motifs from motifDICT that fall within predicted TM domains as predicted for each strain 
###ERROR: CUrrently adding each motif to the filtered DICT for every TM helix it is not in, need to fix so I only add motifs once
##Error should be fixed, but might need to verify quickly that it is running well.
def TMMotifFilter(TMHMMfile,MotifDICT):
    TMMpredDICT=TMHMMRead(TMHMMfile)
    TMfilteredMotifDICT={}
    for i in TMMpredDICT:
        
        try:
            TMfilteredMotifDICT[i]=[]
            for motif in MotifDICT[i].motifs:
                MotifInHelix=False
                motifStart=motif.position[0]
                motifEnd=motif.position[1]

                for TMhelix in TMMpredDICT[i]:
                    if TMhelix[0]=='TMhelix':

                        TMstart=int(TMhelix[1])
                        TMend=int(TMhelix[2])

                        if not (motifStart > TMend or motifEnd < TMstart):
                            #if True we are in TM helix

                            MotifInHelix=True
                            break
                if not MotifInHelix:
                    TMfilteredMotifDICT[i].append(motif)
        except:
            print('%s Not Found' % (i))
    return TMfilteredMotifDICT


# In[ ]:


def TMMotifFilter2(TMHMMfile,Motifs, NCBIID='all'):
    TMMpredDICT=TMHMMRead(TMHMMfile)
    TMfilteredMotifDICT={}
    if NCBIID=='all':
        for i in Motifs:


            #try:
            print('################',i)
            TMfilteredMotifDICT[i]=[]
            for motif in Motifs[i]:
                MotifInHelix=False
                motifStart=motif.position[0]
                motifEnd=motif.position[1]

                for TMhelix in TMMpredDICT[i]:
                    if TMhelix[0]=='TMhelix':

                        TMstart=int(TMhelix[1])
                        TMend=int(TMhelix[2])

                        if not (motifStart > TMend or motifEnd < TMstart):
                            #if True we are in TM helix

                            MotifInHelix=True
                            break
                if not MotifInHelix:
                    TMfilteredMotifDICT[i].append(motif)
            #except:
                #print '%s Not Found' % (i)
        return TMfilteredMotifDICT
    else:#only looking at a filtering a single NCBI entry for a list of motifs, not DICT

        TMfilteredMotifDICT[NCBIID]=[]
        for motif in Motifs:
            #print motif.motifID
            MotifInHelix=False
            motifStart=motif.position[0]
            motifEnd=motif.position[1]

            for TMhelix in TMMpredDICT[NCBIID]:
                if TMhelix[0]=='TMhelix':

                    TMstart=int(TMhelix[1])
                    TMend=int(TMhelix[2])

                    if not (motifStart > TMend or motifEnd < TMstart):
                        #print motifStart,motifEnd,TMstart,TMend
                        #if True we are in TM helix

                        MotifInHelix=True
                        break
            if not MotifInHelix:
                TMfilteredMotifDICT[NCBIID].append(motif)
        return TMfilteredMotifDICT


# In[2]:


#using clustalw in location /home/agunnar/clustalw-2.1-linux-x86_64-libcppstatic/clustalw2
#make an alignment file of a fasta of sequences or if a second fasta is included with an already existing alignment file 
#the second is aligned with the already existing one. Filenames to be called without .fasta or .aln
def StrainAlign(FASTAfile1, FASTAfile2add=None, OW=False):
    if not FASTAfile2add:
        try:
            with open(FASTAfile1+'.aln', 'r') as f:
                if OW:
                    print('%s.aln already exists, overwriting...' % (FASTAfile1))
                else:
                    print('%s.aln already exists, not overwriting (set OW argument to True to change). Returning old alignment' % (FASTAfile1))
                    alignment = AlignIO.read(FASTAfile1+'.aln','clustal')
                    return alignment
        except IOError:
            print('Creating new align file %s.aln' % (FASTAfile1))
            
    clustal_run = r"~/clustalw-2.1-linux-x86_64-libcppstatic/clustalw2"
    in_file = FASTAfile1 #r"ZikaNewMotifStrainsPreAlign.fasta"
    if FASTAfile2add:
        print('Adding %s.fasta to already aligned file %s.aln' % (FASTAfile2add,FASTAfile1))
        in_file = FASTAfile1+'.aln'
        new_seqs=FASTAfile2add+'.fasta'
        clustalw_cline=ClustalwCommandline(clustal_run, quicktree=True, profile1=in_file, profile2=new_seqs, sequences=True)
        stdout, stderr = clustalw_cline()
        alignment = AlignIO.read(FASTAfile2add+'.aln','clustal')
        return alignment
    else:
        in_file = FASTAfile1+'.fasta'
        clustalw_cline = ClustalwCommandline(clustal_run, infile=in_file, quicktree=True)
        stdout, stderr = clustalw_cline()
        alignment = AlignIO.read(FASTAfile1+'.aln','clustal')
        return alignment


# In[ ]:


def true2index(seq, truepos):
    for i in range(len(seq)):
        if seq[i][0]==truepos:
            return i


# In[ ]:


#Use this function to convert positions from eg uniprot for polyprotein features to the relevant index positions ina custom
#alignment (refIndex=True), and then again to convert from alignment index to a true protein position in another strain
#that corresponds to the location in the uniprot reference strain. Use this to then split polyprotein into individual proteins

#This function will also be able to align and adjust positions for the reconstructed sequences without alterations. Alignment 
#will be including reconstructed, and I just send the name of the ancestral seq (N4 etc)
def findAlignPostn(refSeqID,position,alignment, refIndex=True):
    alignmentids=[i.id for i in alignment]
    assert refSeqID in alignmentids
    for record in alignment:
        
        if record.id == refSeqID:
            trueCount=0
            for n in range(len(record.seq)):
                if record.seq[n]=='-':
                    continue
                else:
                    trueCount+=1
                if refIndex:
                    if trueCount==position:
                        return n
                else:
                    if n==position:
                        return trueCount


# In[ ]:


########24/05/2018##########
#Needs to work for reconstructed sequences as well.

def disorderpolyproteinsplit(querySeq, PolyproteinDICT, alignment):
    disordersplit={}

    for i in PolyproteinDICT:

        identifier = 'Zika_Polyprotein,H8XX12_ZIKV'
        ancestralseq='NC_012532.1'#this is the global permanent ancestral seq for things related to the positions of subproteins
        if isinstance(querySeq,SeqIO.SeqRecord):
            queryseq=querySeq.id
        else:
            queryseq=querySeq#'KY785475.1'
        
        AdjustedIndex1Start=findAlignPostn(ancestralseq,PolyproteinDICT[i][0],alignment)
        AdjustedIndex1End=findAlignPostn(ancestralseq,PolyproteinDICT[i][1],alignment)

        queryseqTrueStart=findAlignPostn(queryseq,AdjustedIndex1Start, alignment, refIndex=False)
        queryseqTrueEnd=findAlignPostn(queryseq,AdjustedIndex1End, alignment, refIndex=False)

        if isinstance(querySeq,SeqIO.SeqRecord):
            currentseq=copy.deepcopy(querySeq)
            currentseq.seq=currentseq.seq.translate()[queryseqTrueStart-1:queryseqTrueEnd]
        else:
            currentseq= SeqIO.read('NCBI_Zika/'+queryseq, "genbank")
            ORF=readingframe(currentseq)
            currentseq.seq = currentseq.seq[ORF[0]:ORF[1]].translate()[queryseqTrueStart-1:queryseqTrueEnd]
        
        disordersplit[i]=DisorderConsensus(currentseq,'iupred',queryseqTrueStart-1)
    return disordersplit


# In[ ]:


#seq must be a bioseq record with the sequence of interest
#think about adding the index to adjust amino acid positions to fit with reference polyprotein. Will make it easier to compare 
#motifs
def DisorderConsensus(seq, predictor,index=0):
    if predictor=='iupred':
        iupred_loc_path='/home/agunnar/iupred/'
        with open(iupred_loc_path+seq.id+".fasta", "w") as output_handle:
            SeqIO.write(seq, output_handle, "fasta")
            
        run_iupred_cmd = [iupred_loc_path+'program', iupred_loc_path+seq.id+".fasta",'short']
        iupred_output = subprocess.Popen(run_iupred_cmd, stdout=subprocess.PIPE)
        iupred_data=[]
        for line in iupred_output.stdout:
        
            entry = line.split()
            if entry[0] !='#':
                entry[0]=index+int(entry[0])
                iupred_data.append(entry)
        iupred_output.wait()
        assert iupred_output.returncode==0
        os.remove(iupred_loc_path+seq.id+".fasta")
        return iupred_data


# In[ ]:


#one shortcoming atm is that i plot disorder based on RefStrain but I use disorder information for new motifs from each relevant
#new strain and so should plot the most representative disorder value for the motifs I use. For Zika plotting disorder for
#RefStrain is gonna be completely fine since sequence vary so little.

# relies on #polyporteinDIct which is a global imported DICT defined in CommonMotifProb.py


def plotMotifDis(RefStrain,polyproteinchoice, motiflist, alignment, colourlist, distance,MotifDICT, ELMlist=[]):
    #make disordersplit based on current ancestral sequence
    disordersplit=disorderpolyproteinsplit(RefStrain,PolyproteinDICT, alignment)
    
    lowestSeq=polyproteinchoice[0]
    highestSeq=polyproteinchoice[-1]
    
    Seq=[]
    SeqIndex=[]
    count=0
    for seq in polyproteinchoice:
        if count==0:
            SeqIndex.append((0,len(disordersplit[seq])))
        else:
            SeqIndex.append((SeqIndex[count-1][1],SeqIndex[count-1][1]+len(disordersplit[seq])))
        count+=1
        Seq=Seq+disordersplit[seq]

    rangeStart=Seq[0][0] #Uniprot True position
    rangeEnd=Seq[-1][0] #Uniprot True position
    print('Seq', Seq)
    #if motifGain:
    
    relevantProbs=motifprobs(RefStrain,motiflist,rangeStart,rangeEnd, Seq, distance, ELMlist)
    
    #else:
    #    relevantProbs
        
        
    height_ratios=[2]+[1.5]*len(relevantProbs)+[0.5]
    ax_handles= plt.subplots(2+len(relevantProbs), 1, gridspec_kw = {'height_ratios':height_ratios}, figsize=(80, 4*len(height_ratios)), sharex=True)
    features=[]
    count=0
    for i in SeqIndex:
        
        print(i)
        features.append(GraphicFeature(start=i[0], end=i[1], strand=+1, color = '#ffcccc', label=polyproteinchoice[count]))
        count+=1
        
    for motif in motiflist:
        if motif.position[0]>rangeStart-1 and motif.position[1]<rangeEnd:
            adjustedStart=true2index(Seq, motif.position[0]+1)
            adjustedEnd=true2index(Seq, motif.position[1]+1)
            features.append(GraphicFeature(start=adjustedStart, end=adjustedEnd, strand=+1, color = colourlist[motif.motifID], label=motif.motifID))


    def _create_highlights(motiflist):
        highlights={}
        for motif in motiflist:
            if motif.position[0]>rangeStart-1 and motif.position[1]<rangeEnd:#polyporteinDIct is global imported DICT defined in CommonMotifProb.py
                adjustedStart=true2index(Seq, motif.position[0]+1)
                print(motif.motifID,motif.motifSeq,Seq[adjustedStart:adjustedStart+len(motif.motifSeq)])
                try:
                    highlights[motif.motifID].append(adjustedStart)
                except KeyError:
                    highlights[motif.motifID]=[adjustedStart]
        return highlights
    
    def _motifs_already_present(motifID, motifDICT=MotifDICT):
        highlightPrior=[]
        if isinstance(RefStrain, SeqIO.SeqRecord):
            for motif in motifDICT[RefStrain.id]:#Only currently works with ancestral seqs need to make general case
                if motif.motifID==motifID and motif.position[0]>rangeStart-1 and motif.position[1]<rangeEnd:
                    adjustedMotif=true2index(Seq, motif.position[0]+1)
                    highlightPrior.append(adjustedMotif)
        else:
            for motif in motifDICT[RefStrain]:#Only currently works with ancestral seqs need to make general case
                if motif.motifID==motifID and motif.position[0]>rangeStart-1 and motif.position[1]<rangeEnd:
                    adjustedMotif=true2index(Seq, motif.position[0]+1)
                    highlightPrior.append(adjustedMotif)            
        return highlightPrior
            
    
    highlights=_create_highlights(motiflist)
    
    def plot_iupred(ax):
        yy = [float(i[2]) for i in Seq]
        xx = np.arange(len(Seq))
        ax.fill_between(xx, yy, alpha=0.3)
        ax.set_ylabel("Iupred disorder")
    def plot_motifprob(ax,prob,title):
        yyy=np.asarray(prob)
        xxx = np.arange(len(prob))
        ax.set_ylabel("Probability(Log)")
        ax.set_yscale('linear')
        ax.plot(xxx,yyy,'.-',c=colourlist[title], label=title)
        if title not in ELMlist:
            ax.plot(highlights[title],yyy[highlights[title]], linestyle='none', color='r', marker='o')
        priorMotifs=_motifs_already_present(title)
        print(title)
        print(priorMotifs)
        ax.plot(priorMotifs,yyy[priorMotifs], linestyle='none', color='b', marker='o')
        ax.legend(loc='lower right')

    record = GraphicRecord(sequence_length=len(Seq), features=features)
    #ax, levels = record.plot()
    plot_iupred(ax=ax_handles[1][-1])#ax4)
    axindex=1
    for i in relevantProbs:
        plot_motifprob(ax_handles[1][axindex],relevantProbs[i],i)
        axindex+=1
    record.plot(ax=ax_handles[1][0], with_ruler=True)
    if isinstance(RefStrain, SeqIO.SeqRecord):
        ax_handles[0].savefig(RefStrain.id+lowestSeq+'-'+highestSeq+"DisPlot.png")
    else:
        ax_handles[0].savefig(RefStrain+lowestSeq+'-'+highestSeq+"DisPlot.png")


# In[ ]:


#add choice of probability filter option. If option is true filter each probability as a motif list and change irrelevant probability values to lowest value for now.

##############24/05/2018#################
#added ability to handle Bio.SeqRecord.SeqRecord sequences which will contain ancestral reconstructions, as well as 
#retaining normal handling of NCBI ids. If it bugs it will likely be because of differences in ORF ranges and input ranges.


def motifprobs(RefStrain,motiflist,rangeStart,rangeEnd, disSeq, distance,ELMlist):
    if isinstance(RefStrain,SeqIO.SeqRecord):
        currentseq=RefStrain
    else:
        currentseq= SeqIO.read('NCBI_Zika/'+RefStrain, "genbank")#AY632535.2
        ORF=readingframe(currentseq)
        currentseq.seq=currentseq.seq[ORF[0]:ORF[1]]
    #problem with indexing since I am using PolyproteinDICT for ranges. WIll give me ranges not based on actual refseq index
    #but actually uniprot index
    rangestart=rangeStart-1 #Converted to python Index
    rangeend=rangeEnd #Stays Uniprot True because when slicing this position is excluded meaning final indluded is rangeEnd-1
    RNASeq=str(currentseq.seq[rangestart*3:rangeend*3].transcribe())
    ####make sure it is the same sequence
    for k in disSeq: 
        try:
            disseqpepseq+=k[1]
        except NameError:
            disseqpepseq=k[1]
    #print str(currentseq.seq[ORF[0]+(rangestart*3):ORF[0]+(rangeend*3)].translate())
    #print str(disseqpepseq)
    assert str(currentseq.seq[rangestart*3:rangeend*3].translate())==str(disseqpepseq)
    ####
    relevantProbs={}
    for i in motiflist:
        if i.position[0]>rangestart and i.motifID not in list(relevantProbs.keys()) and i.motifID in list(ELMsDICT.keys()) and i.position[1]<rangeend:
            relevantProbs[i.motifID]=[]
            probs=MotifProbability(RNASeq,ELMsDICT[i.motifID], distance)
            
            ### Filter probs basde on disorder (need to add TM filter)

            minprobs=min(probs)
            for j in range(len(probs)):
                currentposition=(j,j+len(i.motifSeq))
                if probs[j]!=1 and isAvgDisorder(disSeq,currentposition[0], currentposition[1]):
                    relevantProbs[i.motifID].append(probs[j])
                else:
                    relevantProbs[i.motifID].append(minprobs)
    if len(ELMlist)>0:
        for p in ELMlist:
            if p not in list(relevantProbs.keys()) and p in list(ELMsDICT.keys()):
                relevantProbs[p]=[]
                probs=MotifProbability(RNASeq,ELMsDICT[p], distance)

                ### Filter probs basde on disorder (need to add TM filter)

                minprobs=min(probs)
                for j in range(len(probs)):
                    currentposition=(j,j+len(ELMsDICT[p]))
                    if probs[j]!=1 and isAvgDisorder(disSeq,currentposition[0], currentposition[1]):
                        relevantProbs[p].append(probs[j])
                    else:
                        relevantProbs[p].append(minprobs)
            
    return relevantProbs


# In[ ]:


#think about difference in doing it from and to true
def adjust_motif_positions(motifSeq, RefSeq, motifposStart,motifposEnd, alignment):
    #using indexed motif position here directly rather than converting it to protein True position since I correct for that
    #later in the plot function so I want to keep it being motif position index rather than True Uniprot numbering True ie
    #difference between 0,1,2... 1,2,3...
    AdjustedIndex1Start=findAlignPostn(motifSeq,motifposStart,alignment)
    AdjustedIndex1End=findAlignPostn(motifSeq,motifposEnd,alignment)

    motifStartinRef=findAlignPostn(RefSeq,AdjustedIndex1Start, alignment, refIndex=False)
    motifEndinRef=findAlignPostn(RefSeq,AdjustedIndex1End, alignment, refIndex=False)
    
    return (motifStartinRef,motifEndinRef)


# In[ ]:


#filter motifs based on disorder score (iupred). SHould also include others/consensus maybe? iupred short is extremely stringent
#Criterium should be average value over motif sequence and no stretch of values under cutoff longer than 20ish% of motif length
#or maybe over 3 residues.
#takes iupred ouput from disorder consensus function (seq and disorder data) and a list of motifs (currently only works for motifs
#generated specifically for input sequence, otherwise there will be out of index issues) and returns a filtered list of motifs
#only retaining the ones that fall above a cutoff value. The default cutoff is fairly low to remove motifs in regions that are 
#highly unlikely to be disordered but keeps them in ambiguous regions and regions highly likely to be disordered.
def disorderFilter(disSeq,motiflist,cutoff=0.3):
    motifsdisfilter=[]
    seqrange=(disSeq[0][0],disSeq[-1][0])
    print(seqrange)
    for motif in motiflist:
        if motif.motifID=='LIG_PDZ_Class_1':
            print(motif.motifID, motif.position)
        if motif.position[0]+1>=seqrange[0] and motif.position[1]<=seqrange[1]:#make sure i correctly identify terminal motifs
            contiguousLow=0
            motifStart=true2index(disSeq, motif.position[0]+1)
            motifEnd=true2index(disSeq, motif.position[1]+1)
            currentMotif=[float(i[2]) for i in disSeq[motifStart:motifEnd]]#cant be index since i need it based on full seq position
            disorderAvg=sum(currentMotif)/(len(currentMotif))
            print(motif.motifID,motif.motifSeq, end=' ')
            print(disorderAvg, end=' ')
            if disorderAvg>=cutoff:
                maxStretch=0
                for AAdisorder in currentMotif: #if there is a stretch longer than 3 within motif of lower than cutoff disorder,exclude motif even if average is higher
                    if AAdisorder<cutoff:
                        contiguousLow+=1
                        maxStretch=contiguousLow
                    else:
                        contiguousLow=0
                if maxStretch<3:
                    print('saved')
                    motifsdisfilter.append(motif)
                else:
                    print(maxStretch)
            print(' ')
    return motifsdisfilter


# In[ ]:


def isAvgDisorder(disSeq,RangeStart, RangeEnd,cutoff=0.3, isIndex=True):
    if isIndex:
        motifStart=RangeStart
        motifEnd=RangeEnd
    elif not isIndex:
        motifStart=true2index(disSeq, RangeStart+1)
        motifEnd=true2index(disSeq, RangeEnd+1)
    contiguousLow=0    
    currentMotif=[float(i[2]) for i in disSeq[motifStart:motifEnd]]#cant be index since i need it based on full seq position
    disorderAvg=sum(currentMotif)/(len(currentMotif))
    if disorderAvg>=cutoff:
        maxStretch=0
        for AAdisorder in currentMotif: #if there is a stretch longer than 3 within motif of lower than cutoff disorder,exclude motif even if average is higher
            if AAdisorder<cutoff:
                contiguousLow+=1
                maxStretch=contiguousLow
            else:
                contiguousLow=0
        if maxStretch<3:
            return True
        else:
            return False
    else:
        return False

