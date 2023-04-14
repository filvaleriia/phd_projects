import rdkit
from rdkit import Chem
from rdkit.Chem import Draw,AllChem,DataStructs
import selfies as sf
from group_selfies import fragment_mols, Group, MolecularGraph, GroupGrammar, group_encoder
import deepsmiles
import csv
import numpy as np
import itertools as it
import random
import copy
import time
import pickle
import lzma
from rdkit import RDLogger                                                                                                                                                               
RDLogger.DisableLog('rdApp.*') 

converter = deepsmiles.Converter(rings=True, branches=True)
grammar = GroupGrammar.from_file('data/nfa_test_grammar.txt') #load

#with loaded grammar
def get_gsf(m):
    extracted = grammar.extract_groups(m)
    encoded = grammar.encoder(m, extracted)
    return encoded

def decode_gsf(gsf):
    mol = grammar.decoder(gsf)
    return mol

#print(get_gsf(Chem.MolFromSmiles("CCCc1ccccc1CCC")))

def bredt_violation(mol):
    """filter to check if there's a violation of bredts rule
    (no bridgehead double bonds). in this case also no aromatic bonds are
    allowed so it also filters cyclophane type compounds (intendedly)
    """
    problem=False
    sssr_idx = [set(x) for x in list(Chem.GetSymmSSSR(mol))]
    for i,ring1 in enumerate(sssr_idx):
        for j,ring2 in enumerate(sssr_idx):
            if i>j:
                intersect_idx=ring1.intersection(ring2)
                if len(intersect_idx)>2:
                    for idx in intersect_idx:
                        bh=False
                        neighbors = [a.GetIdx() for a
                                     in mol.GetAtomWithIdx(idx).GetNeighbors()]
                        for nidx in neighbors:
                            if nidx not in intersect_idx:
                                bh=True
                        if bh==True:
                            bondtypes = [b.GetBondType() for b in
                                         mol.GetAtomWithIdx(idx).GetBonds()]
                            if set(bondtypes) != {Chem.rdchem.BondType.SINGLE}:
                                problem=True
    return problem
    
def filter_problematic_frags(mol):
    reason=""
    problem=False
    if mol==None:
        problem=True
        reason = "the input mol being None"
    patts=["[*]1[*]=[*][*]=1","[*]1[*]=[*]1","[*]1~[*]~[*]~[*]~1",
           "[!#6]-[F,Cl,Br,I]","[!a;#6]=[!a;#6]-[!a;#6]=[!a;#6]",
           "[N,O,S]-[CX4]-[N,O,S]","[SX3,SX4,SX5,SX6;h1,h2,h3]",
          "[C;!r]=[N;!r]-[C;!r]=[C;!r]","[C;!r](=[O;!r])-[N;!r]-[C;!r](=[O;!r])-[N;!r]"]
    if problem==False and len([a for a in mol.GetAtoms() if a.GetAtomicNum()!=1])<10:#think about this number
        problem=True
        reason="an atom count that is under 10"
    if problem==False:
        for patt in patts:
            if mol.HasSubstructMatch(Chem.MolFromSmarts(patt))==True:
                problem=True
                reason="a problematic pattern:"+patt
    if problem==False:
        ringsizes = [len(i) for i in list(Chem.GetSymmSSSR(mol))]
        if len(ringsizes)>0:
            if max(ringsizes)>7: #filtered 8membered rings and bigger
                problem=True
                reason="a ringsize over 7"
    if problem==False:
        problem=bredt_violation(mol)
        if problem==True:
            reason="a Bredt violation or cyclophane type structure"
    return problem,reason

class Generator:
    def __init__(self, mode="groupselfies", depth=7,max_n=1000000,verbose=False):
        self.mode=mode
        self.depth=depth
        self.max_n=max_n
        self.verbose=verbose
        self.all_chars_sfs = [] # array of all possible chars or tokens
        self.all_chars_gsfs = []
        self.all_chars_smi = []
        self.all_chars_dsmi = []
        self.all_list_sfs = [] # array of all smiles/selfies in the set, as arrays split per token
        self.all_list_gsfs = []
        self.all_list_smi = []
        self.all_list_dsmi = []
        self.inchikeys=[] # fps formatching later
        self.ecfps=[]
        self.markov_dict={} #contains transition probabilities
        self.starters=[] #starting token sequence
        self.bad_mols=[]
        self.good_mols=[]

        
    def load(self, filepath='data/chembl100K.csv'): 
        """ load a csv of SMILES. if the file ends with xz it will 
        decompress using lzma.
        it also stores ecfps and inchi of the corpus for usage in 
        metrics and duplication detection later.
        selfies,groupselfies,deepsmiles and smiles get calculated
        for the set
        max_n can be used to load only max_n first SMILES from the
        input set
        """
        curr_n=0
        if filepath[-2:]=="xz":
            f = lzma.open(filepath, mode='rt',encoding="utf-8")
        else:
            f = open(filepath, newline='')
        reader = csv.reader(f, delimiter=',', quotechar='"')
        for row in reader:
            if row[0]!="smiles" and curr_n<self.max_n:
                try:
                    mol=Chem.MolFromSmiles(row[0])
                    selfie=sf.encoder(row[0])
                    gselfie=get_gsf(mol)
                    if selfie!=None:
                        self.inchikeys.append(Chem.MolToInchiKey(mol))
                        self.ecfps.append(AllChem.GetMorganFingerprintAsBitVect(mol,2,nBits=1024))
                        sflist = list(sf.split_selfies(selfie))
                        self.all_chars_sfs += sflist
                        self.all_list_sfs.append(sflist)
                        gsflist = list(sf.split_selfies(gselfie))
                        self.all_chars_gsfs += gsflist
                        self.all_list_gsfs.append(gsflist)
                        smilist=[x for x in row[0]]
                        self.all_chars_smi += smilist
                        self.all_list_smi.append(smilist)
                        dsmilist=[x for x in converter.encode(row[0])]
                        self.all_chars_dsmi += dsmilist
                        self.all_list_dsmi.append(dsmilist)
                        curr_n+=1
                except Exception as e:
                    pass
                    print(e)
        self.all_chars_sfs = tuple(list(set(self.all_chars_sfs))+['t'])
        self.all_chars_gsfs = tuple(list(set(self.all_chars_gsfs))+['t'])
        self.all_chars_dsmi = tuple(list(set(self.all_chars_dsmi))+['t'])
        self.all_chars_smi = tuple(list(set(self.all_chars_smi))+['t'])

        self.inchikeys=set(self.inchikeys)
        
    def populate_matrix(self):
        """ read the corpus and collect the transition probabilities
        """
        if self.mode=="SMILES":
            all_chars=self.all_chars_smi
            all_list=self.all_list_smi
        elif self.mode=="DeepSMILES":
            all_chars=self.all_chars_dsmi
            all_list=self.all_list_dsmi
        elif self.mode=="selfies":
            all_chars=self.all_chars_sfs
            all_list=self.all_list_sfs
        elif self.mode=="groupselfies":
            all_chars=self.all_chars_gsfs
            all_list=self.all_list_gsfs
        else:
            print("invalid mode selected")
            
        weights_dict={}
        self.markov_dict={}
        self.starters = [] #starting options
        for key in all_chars:
            weights_dict[key]=0
        for csf in all_list:
            for i in range(len(csf)-self.depth):
                curr_key="".join(csf[i:i+self.depth])
                if curr_key in self.markov_dict:
                    self.markov_dict[curr_key][csf[i+self.depth]]+=1
                else:
                    self.markov_dict[curr_key]=copy.deepcopy(weights_dict)
                    self.markov_dict[curr_key][csf[i+self.depth]]+=1
            curr_key="".join(csf[-1*self.depth:])
            if curr_key in self.markov_dict:
                self.markov_dict[curr_key]['t']+=1 #termination token
            else:
                self.markov_dict[curr_key]=copy.deepcopy(weights_dict)
                self.markov_dict[curr_key]['t']+=1
            self.starters.append(csf[:self.depth])
            
    def gen_mol(self):
        """ generate one molecule using the markov chain
        """
        newsf=random.choice(self.starters)
        prob_dict = self.markov_dict["".join(newsf[-1*self.depth:])]
        nextchar = random.choices(list(prob_dict.keys()),
                                  weights=prob_dict.values())
        while nextchar != ['t']:
            newsf+=nextchar
            prob_dict = self.markov_dict["".join(newsf[-1*self.depth:])]
            nextchar = random.choices(list(prob_dict.keys()),
                                      weights=prob_dict.values())
        if self.mode=="SMILES":
            mol=Chem.MolFromSmiles("".join(newsf))
        elif self.mode=="DeepSMILES":
            try:
                decoded = converter.decode("".join(newsf))
                mol=Chem.MolFromSmiles(decoded)
            except:
                mol = None
        elif self.mode=="selfies":
            mol = Chem.MolFromSmiles(sf.decoder("".join(newsf)))
        elif self.mode=="groupselfies":
            mol = decode_gsf("".join(newsf))
        return mol,newsf
        
    def generate_and_filter_mols(self,nmols=1000):
        t1=time.time()
        mols=[]
        sf_efficiency=[]
        #filter them
        self.good_mols=[]
        self.bad_mols=[]
        filt_inchis=[]
        duplicount=0
        trys=0
        while len(self.good_mols)<nmols:
            mol,newsf=self.gen_mol()
            if self.mode == "selfies":
                currsf="".join(newsf)
                roundtripsf = list(sf.split_selfies(sf.encoder(sf.decoder(currsf))))
                sf_efficiency.append(len(roundtripsf)/len(newsf))
            elif self.mode == "groupselfies":
                smi="".join(newsf)
                roundtripsf = list(sf.split_selfies(get_gsf(decode_gsf(smi))))
                sf_efficiency.append(len(roundtripsf)/len(newsf))
            else:
                sf_efficiency.append(1)
            try:
                Chem.SanitizeMol(mol) #especially to clean chiral flags
                mol=Chem.MolFromSmiles(Chem.MolToSmiles(mol).replace("[C]","C").replace("[CH]","C").replace("[CH2]","C")) #remove radical artefacts
                if filter_problematic_frags(mol)[0]==False:# and sf_efficiency[i]>0.85:#85:
                    self.good_mols.append(mol)
                    filt_inchis.append(Chem.MolToInchiKey(mol))
                    if filt_inchis[-1] in inchikeys:
                        duplicount += 1
                else:
                    self.bad_mols.append(mol)
            except Exception as e:
                pass#print("sanitization failure:")
                #print(e)
            trys+=1
        if self.verbose==True:
            print("amount of strings parsed is:",len(self.good_mols)+len(self.bad_mols))
            print("efficiency is",len(self.good_mols)/(len(self.good_mols)+len(self.bad_mols)))
            if len(self.good_mols)>0:
                print("duplicate rate of filt compounds is",duplicount/len(self.good_mols))
                print("rate of uniques is",len(set(filt_inchis))/len(filt_inchis))
                dur=time.time()-t1
                print("generation took {} s, thats {} sec per unique molecule after filtering".format(round(dur,2),dur/len(set(filt_inchis))))
                print("total amount of trys was {}. this means a validity of {} %".format(trys,100*(len(self.good_mols)+len(self.bad_mols))/trys))
            else:
                print("0 mols remain!!!!")
                
    def save_state(self,path="data/weights.pkl.gz"):
        with lzma.open(path, 'wb') as f:
            pickle.dump(self.__dict__,f, pickle.HIGHEST_PROTOCOL)
            
    def load_state(self,path="data/weights.pkl.gz"):
        with lzma.open(path, 'rb') as f:
            self.__dict__ = pickle.load(f)

        
    
