
# coding: utf-8

# In[41]:

import sys, getopt
import numpy as np
import random
import pandas as pd


# In[42]:





# In[46]:




# In[26]:

AminoAcid1LTR={'F','L','I','M','V','S','P','T','A','Y','H','Q','N','K','D','E','C','W','R','S','G','*'}


# In[27]:

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


# In[28]:

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


# In[29]:

def createmotif(motif, prot=True):
    finalmotif=[]
    
    for i in motif:
        if i[0]=='.':
            finalmotif.append(list(AminoAcid1LTR))
        elif i[0]=='^':
            aas=i[1:]
            finalmotif.append(list(AminoAcid1LTR-set(aas)))
        else:
            finalmotif.append(i)
    if prot:
        return finalmotif
    else:
        codonfinalmotif=[]
        for k in finalmotif:
            temp=[]
            for aa in k:
                temp+=[k.replace('U','T') for k in AAtoCodonDICT[SingleLetterDICT[aa]]]
            codonfinalmotif.append(temp)
        return codonfinalmotif
        


# In[30]:

def randomNucleotide(nt,weights={'A':0.25,'T':0.25,'C':0.25,'G':0.25},randnum=0.0):#A,T,C,G

    if nt=='A':
        if randnum<weights['A']+weights['T']: #Aweight+Tweight
            return 'T'
        elif randnum<weights['A']+weights['T']+weights['C']: #Aweight+Tweight+Cweight
            return 'C'
        else:
            return 'G'
    elif nt=='T':
        if randnum<weights['A']+weights['T']: #Aweight+Tweight
            return 'A'
        elif randnum<weights['A']+weights['T']+weights['C']: #Aweight+Tweight+Cweight
            return 'C'
        else:
            return 'G'
    elif nt=='C':
        if randnum<weights['A']+weights['C']: #Aweight+Tweight
            return 'A'
        elif randnum<weights['A']+weights['T']+weights['C']: #Aweight+Tweight+Cweight
            return 'T'
        else:
            return 'G'
    elif nt=='G':
        if randnum<weights['A']+weights['G']: #Aweight+Tweight
            return 'A'
        elif randnum<weights['A']+weights['G']+weights['C']: #Aweight+Tweight+Cweight
            return 'C'
        else:
            return 'T'
    
    
    


# In[31]:

ALTNTevolveprobDICT={'A':[1-(1.3*10**-5+3.4*10**-5+3.0*10**-4), 1.3*10**-5, 3.4*10**-5, 3.0*10**-4],
                  'T':[4.5*10**-6, 1-(4.5*10**-6+3.1*10**-4+3.6*10**-5), 3.1*10**-4, 3.6*10**-5],
                  'C':[1.7*10**-5, 4.6*10**-5, 1-(1.7*10**-5+4.6*10**-5+9.7*10**-6), 9.7*10**-6],
                  'G':[7.2*10**-5, 6.0*10**-5, 2.8*10**-5, 1-(7.2*10**-5+6.0*10**-5+2.8*10**-5)]}


# In[32]:

#define some commonly calculated variables used by all instances to save speed

#motif1=createmotif([['R','K'],['R','K'],['S','T'],['^','P']],prot=False)
motif1=createmotif([['R'],['R'],['T','S','N'],['L'],['R']],prot=False)
motif2=createmotif([['R'],['R'],['M'],['L'],['R']],prot=False)
motif3=createmotif([['R'],['R'],['I'],['L'],['R']],prot=False)
#motifpositions=range(len(motif1))


# In[33]:

def ismotifmatch(seq,motif):
    motifpositions=range(len(motif))
    for i in motifpositions:

        if seq[i*3:(i*3)+3] not in motif[i]:
            return False
    return True
    


# In[34]:


#NTevolveprobDICT={'A':[1-(1.8*10**-5+1.5*10**-5+2.0*10**-4), 1.8*10**-5, 1.5*10**-5, 2.0*10**-4],
#                  'T':[1.4*10**-5, 1-(2.3*10**-4+3.5*10**-5+1.4*10**-5), 2.3*10**-4, 3.5*10**-5],
#                  'C':[7.7*10**-6, 2.7*10**-5, 1-(7.7*10**-6+2.7*10**-5+5.1*10**-6), 5.1*10**-6],
#                  'G':[3.1*10**-5, 3.5*10**-5, 5.4*10**-5, 1-(3.1*10**-5+3.5*10**-5+5.4*10**-5)]}

NTevolveprobDICTDICT={'A':{'A':1-(1.8*10**-5+1.5*10**-5+2.0*10**-4), 'T': 1.8*10**-5, 'C': 1.5*10**-5, 'G': 2.0*10**-4},
                  'T':{'A': 1.4*10**-5, 'T': 1-(2.3*10**-4+3.5*10**-5+1.4*10**-5), 'C': 2.3*10**-4, 'G': 3.5*10**-5},
                  'C':{'A': 7.7*10**-6, 'T': 2.7*10**-5, 'C': 1-(7.7*10**-6+2.7*10**-5+5.1*10**-6), 'G': 5.1*10**-6},
                  'G':{'A': 3.1*10**-5, 'T': 3.5*10**-5, 'C': 5.4*10**-5, 'G': 1-(3.1*10**-5+3.5*10**-5+5.4*10**-5)}}

#NTevolveprobDICT={'A':[1-(1.5*10**-5+3.5*10**-5+2.5*10**-4), 1.5*10**-5, 3.5*10**-5, 2.5*10**-4],
#                  'T':[1.4*10**-5, 1-(2.3*10**-4+3.5*10**-5+1.4*10**-5), 2.3*10**-4, 3.5*10**-5],
#                  'C':[7.7*10**-6, 2.7*10**-5, 1-(7.7*10**-6+2.7*10**-5+1.0*10**-5), 1.0*10**-5],
#                  'G':[3.1*10**-5, 3.5*10**-5, 5.4*10**-5, 1-(3.1*10**-5+3.5*10**-5+5.4*10**-5)]}

class QuasiSeq():
    def __init__(self,parentseq,randnums, seqlen,parentfitness=1.0):#, parentid):
        #self.parent=parentid
        self.seq= self._evolve(parentseq,seqlen,randnums)
        if self.seq==parentseq:
            self.fitness=parentfitness
        else:
            self.fitness=self._fitness(self.seq)
        #else:
        #    self.fitness=1
        
        
    def _evolve(self,seq,seqlen,randnums):
        #newseq=''.join([randomNucleotide(nt,NTevolveprobDICT[nt],randnums[i]) for i,nt in enumerate(seq)])
        
        #for i,nt in enumerate(seq):
            #print randnums[i],NTevolveprobDICTDICT[nt][nt]
            #print randomNucleotide(NTevolveprobDICTDICT[nt],randnums[i])
        newseq=''.join([nt if randnums[i]<NTevolveprobDICTDICT[nt][nt] else randomNucleotide(nt,NTevolveprobDICTDICT[nt],randnums[i]) for i,nt in enumerate(seq)]) 
        
        
        return newseq
    
    def _fitness(self,seq):
        ismotif3=ismotifmatch(self.seq,motif3)
        ismotif1=ismotifmatch(self.seq,motif1)
        ismotif2=ismotifmatch(self.seq,motif2)
        if ismotif1:
            return 1.2
        elif ismotif2:
            return 1.4
        elif ismotif3:
            return 1.0
        else:
            return 0.5
            
        
    #def hasdied(self):
    #    #deathfunction(len(self.children))
    #    if len(self.children)>8:
    #        self.isDead=True
    #    else:
    #        self.isDead=False
    
    #def _generate_ID(self):
    #    seqs=['1','2','3','4','5','6','7','8','9','A','B','C','D','E','F','G','H','I','J','K','L','M','N']
    #    ID=[0]*12
    #    for i in range(12):        
    #        randomIndex=random.randint(0,len(seqs)-1)
    #        ID[i]=seqs[randomIndex]
    #    return ID
    
    
    #def add_children(self,child):
    #    self.children.append(child)
        
        
            
            


# In[35]:

def selection(seqpop,averagefit=1.0,carrymod=1.0):
    averagefit=averagefit/len(seqpop)

    nextfounders=[]
    for particle in seqpop:
        randnum=random.random()
        if randnum<averagefit*particle.fitness*carrymod:
            nextfounders.append(particle)

    return nextfounders


# In[36]:

#add maximum carrying capacity for the population
def EvolveFix(founderseq,seqlen,popcarry):

    if type(founderseq)==str:
        segmentEnd=seqlen
        vRNApop=runSim(founderseq,segmentEnd)

        poptotal=[]
        poptotal+=vRNApop
        lenpoptotal=len(poptotal)
        while lenpoptotal<5:
            if type(poptotal)==list:
                newround=list(poptotal)
            else:
                newround=list([poptotal])
            poptotal=[]

            for founder in newround: 
                
                vRNApop=[founder]
                cRNApop=[]

                while len(cRNApop)<300:
                    cRNApop+=runSim(vRNApop,segmentEnd)
      
                    #print len(cRNApop)
                    vRNApop+=runSim(cRNApop,segmentEnd)

                #print 'End of cRNA production'
                #print 'vRNApop size=',len(vRNApop)
                while len(vRNApop)<11000:
                    vRNApop+=runSim(cRNApop,segmentEnd)


                founders=selection(vRNApop, averagefit=10.0)


                poptotal+=founders
                lenpoptotal=len(poptotal)
    else:
        poptotal=[]
        poptotal+=founderseq
        segmentEnd=seqlen
        lenpoptotal=len(poptotal)
        #if lenpoptotal<popcarry:
        #    modi=1.0+((popcarry-lenpoptotal)/float(popcarry))
#
        #    carrymodifier=modi
#        else:
        carrymodifier=1.0
        if type(poptotal)==list:
            newround=list(poptotal)
        else:
            newround=list([poptotal])
        poptotal=[]
        for founder in newround: 

            vRNApop=[founder]
            cRNApop=[]

            while len(cRNApop)<300:
                cRNApop+=runSim(vRNApop,segmentEnd)

                #print len(cRNApop)
                vRNApop+=runSim(cRNApop,segmentEnd)

            #print 'End of cRNA production'
            #print 'vRNApop size=',len(vRNApop)
            while len(vRNApop)<11000:
                vRNApop+=runSim(cRNApop,segmentEnd)

            #countMotif1=0
            #countMotif2=0
            #countNoMotif=0
            #for k in vRNApop:
            #    if k.fitness==1.4:
            #        countMotif1+=1
            #    elif k.fitness==1.2:
            #        countMotif2+=1
            #    else:
            #        countNoMotif+=1
            #print 'Motif1:',countMotif1,'\n','Motif2:',countMotif2,'\n','NoMotif:',countNoMotif,'\n'
            founders=selection(vRNApop,averagefit=2.0,carrymod=carrymodifier,)


            poptotal+=founders
            
            #if len(poptotal)>popcarry:
            #    print 'Poptotal:',len(poptotal)
            #    survivors=[]
            #    deathprob=(len(poptotal)-float(popcarry))/float(len(poptotal))
            #    for k in poptotal:
            #        if random.random()>(deathprob/k.fitness):
            #            survivors.append(k)
            #    poptotal=list(survivors)
            #    print 'Survivors:',len(poptotal)
    return poptotal


# In[37]:

#if called for a string to found new population use randpool=a list of random numbers of len=seqlen
def runSim(founderseq,seqlen):
    
    if type(founderseq)==str:
        randpool=np.random.random(size=seqlen)
        Founder=QuasiSeq(founderseq,randpool,seqlen)#,None)
        VirusPool=[Founder]
    else:
        VirusPool=founderseq
    randpool=np.random.random(size=len(founderseq)*seqlen)
    newGeneration=[]
    randIndex=0
    for virus in VirusPool:
        randnums=randpool[randIndex:randIndex+seqlen]
        #if not virus.isDead:
        viralChild=QuasiSeq(virus.seq, randnums,seqlen,parentfitness=virus.fitness)#,virus.id)
        newGeneration+=[viralChild]
        randIndex+=seqlen
            #virus.add_children(viralChild.id)
            #if len(VirusPool)>3000000:
            #    randnum=random.random()
            #    if randnum<1.0/len(VirusPool):
            #        virus.isDead=True
    VirusPool=newGeneration
        
    return VirusPool
#for virus in VirusPool:
    #print len(VirusPool), VirusPool[0].__dict__,'\n',VirusPool[1].__dict__,'\n',VirusPool[2].__dict__,'\n',VirusPool[3].__dict__,'\n',VirusPool[4].__dict__


# In[39]:

#AGG AGG GGT ATC CTT GGT GCA (non robust codons apart from S/T position which still needs to evolve)

def main(argv):

    outputfile = ''
    try:
        opts, args = getopt.getopt(argv,"hi:o:",["ofile="])
    except getopt.GetoptError:
        print ('test.py -o <outputfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print( 'test.py -o <outputfile>')
            sys.exit()
        elif opt in ("-o", "--ofile"):
            outputfile = arg
    return outputfile

if __name__ == "__main__":
    outputfile=main(sys.argv[1:])
    #print outputfile
    columns=['Founder','Gen','PopSize','Motif1','Motif2','Selected Offspring']
    alldata=[]
    lenNewGen=0
    l=0
    while lenNewGen<5000:
        if l>=16:
            break
        gen=l
        #print gen
        if l==0:     
            poptotal='CGACGAATCCTCCGA'
            seqlen=len(poptotal)#seq1='CGGCGGATCAAC'
            nextgen=['start']
            newnextgen=[]
        else:
            poptotal=[]
            newnextgen=[]


        fixationdata=[]

        Fixation=False
        #print nextgen
        for seq in range(len(nextgen)):

            #print seq
            if type(nextgen[seq])!=str:     
                poptotal=[nextgen[seq]]
                founder=nextgen[seq].seq
            else:
                founder=poptotal
            lentotal=len(nextgen)
            while lentotal<3000: 
                poptotal=EvolveFix(poptotal,seqlen=seqlen, popcarry=10000)

                lentotal=len(poptotal)
                #print lentotal
                if lentotal==0:
                    alldata.append([founder,gen,lentotal,None,None,None])
                    break
                motifcount1=0

                motifcount2=0


                for m in poptotal:
                    #genotypes.append(m.seq)
                    if m.fitness==1.2:
                        motifcount1+=1
                    elif m.fitness==1.4:
                        motifcount2+=1

            if lentotal>0:
                addnewgen=selection(poptotal,averagefit=2.0)
                newnextgen+=addnewgen
                if len(newnextgen)==0:
                    while len(newnextgen)==0:
                        addnewgen=selection(poptotal,averagefit=2.0)
                        newnextgen+=addnewgen
                currentseqlist=[]
                for t in addnewgen:
                    currentseqlist.append(t.seq)
                alldata.append([founder,gen,lentotal,motifcount1,motifcount2,currentseqlist])
        nextgen=newnextgen
        lenNewGen=len(nextgen)
        l+=1

        #with open(outputfile+'-gen'+str(l)+'.csv', 'wb') as myfile:
        #    wr = csv.writer(myfile)
        #    wr.writerow(alldata)
        FixDataDF=pd.DataFrame(alldata,columns=columns)
        FixDataDF.to_csv(outputfile+str(l)+'.csv')


