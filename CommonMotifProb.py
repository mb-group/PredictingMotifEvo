#!/usr/bin/env python
# coding: utf-8

# In[1]:

from Bio import SeqIO
import pandas as pd




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





#adding len feature check makes it work for both zika and flu for now
def readingframe(NCBIrecord):
    for feature in NCBIrecord.features:
        if feature.type=='CDS' and len(feature.location.parts)!=2:
            return (feature.location.start,feature.location.end)


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
