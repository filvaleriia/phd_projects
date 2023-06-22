import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import AllChem
from rdkit.Chem.Scaffolds.MurckoScaffold import GetScaffoldForMol
from rdkit.Chem.Scaffolds.MurckoScaffold import MurckoScaffoldSmiles
from rdkit.Chem.Scaffolds.MurckoScaffold import MakeScaffoldGeneric
import seaborn as sns
import matplotlib.pyplot as plt
import sys
from rdkit.Chem import Lipinski
from rdkit.Chem.Lipinski import NumAromaticHeterocycles
from rdkit.Chem.Lipinski import NumAliphaticRings
from rdkit.Chem.Lipinski import NumAromaticHeterocycles
from rdkit.Chem.Lipinski import NumAromaticRings
from rdkit.Chem.SaltRemover import SaltRemover
from rdkit.Chem import Draw
import lzma
import csv




def convert_to_scaffold_(df):
    a = []
    print("Convert_")
    for x in range(len(df)):
        try:
            a.append(MurckoScaffoldSmiles(\
                                     Chem.MolToSmiles(MakeScaffoldGeneric\
                                   (Chem.MolFromSmiles(df.loc[x][0])))))
        except:
            print("Faild to create scaffold_csk")
            print("Index",x)
            print(df.loc[x][0])

    dff = pd.DataFrame(data = a)
    return dff

def add_columns_same_like_input_function(df_generated, inputt):
    print("add_columns_same_like_input_function")
    df_generated = pd.DataFrame(df_generated)   
    df = pd.DataFrame()
    df[0] = inputt
    df[1] = int(0)
    for x in range(len(inputt)):
        y = df[0][x]
        df[1][x] = [df_generated[0].value_counts()[y] if y in df_generated[0].unique() else 0][0]       
    
    df[2] = df[1].apply(lambda x : 1 if x > 0 else 0)
    #print(df)
    return df

class Metrics:
    def __init__(self,var):
        self.generated_compounds = []
        self.test_set = []
        self.var = var
        self.unique_compounds = []
        self.unique_input = []
    
    def load(self,filepath_generated_com, filepath_test_set):
        #nacitani vstupu
        generated_compounds = []
        test_set = []

        if filepath_generated_com[-2:]=="xz":
            f = lzma.open(filepath_generated_com, mode='rt',encoding="utf-8")
            generated_compounds = csv.reader(f, delimiter=',', quotechar='"')
        else:
            with open(filepath_generated_com, 'r') as f:
                for line in f.readlines():
                    generated_compounds.append(line)
        self.generated_compounds = pd.DataFrame(generated_compounds)

        if filepath_test_set[-2:] == "xz":
            f = lzma.open(filepath_test_set, mode = 'rt', encoding = "utf-8")
            test_set = csv.reader(f,delimiter=',', quotechar = '"')
        else:
            with open(filepath_test_set, 'r') as f:
                for line in f.readlines():
                    test_set.append(line)
        #self.test_set = pd.DataFrame(test_set)
        print(test_set)
        #print(len(self.generated_compounds))
        #print(len(self.test_set))

    def scaffolds_unique(self,df_generated,df_input):
        
        return df_generated, df_input
    
    

    def main_function_return(self,df_generated, df_input):

        df_input = convert_to_scaffold_(df_input)
        df_generated = convert_to_scaffold_(df_generated)

        self.unique_compounds, self.unique_input = self.scaffolds_unique(df_generated, df_input)



        unique_scaffolds_in_whole_generated_set = df_generated[0].unique()
        unique_input = df_input[0].unique()

        #df obsahuje tri sloupce a unique_input pocet radku. kde sloupec 1 obsahuje count tRS a 2 obsahuje uniqe RS
        df = add_columns_same_like_input_function(df_generated, unique_input)

        #df1 obsahuje tri sloupce a df_input/test_set pocet radku. kde sloupec 1 obsahuje count tRS a 2 obsahuje uniqe RS
        df1 = add_columns_same_like_input_function(df_generated, df_input)

        ns = len(unique_scaffolds_in_whole_generated_set)
        ss = len(df_generated)
        sescy = ns/ss
        sescr = df[2].sum()/ns
        sescry = sescy*sescr
        asescr = df[1].sum()/ss
        
        try:
            tpra = df[2].value_counts()[1]/len(df)
            tpra_text = f"{df[2].value_counts()[1]}/{len(df)}"
            #print(f"TPRA:  {df[2].value_counts()[1]}/{len(df)} = {tpra}", )
        except:
            #print(f"TPRA:  {0}/{len(df)} = {0/len(df)}", )
            tpra = 0
            tpra_text = f"{0}/{len(df)}"

        try:
            tprar = df1[2].value_counts()[1]/len(df1)
            tprar_text = f"{df1[2].value_counts()[1]}/{len(df1)}"
            #print(f"TPRAR:  {df1[2].value_counts()[1]}/{len(df1)} = {tprar}", )
        except:
            #print(f"TPRAR:  {0}/{len(df1)} = {0/len(df1)}", )
            tprar = 0
            tprar_text = f"{0}/{len(df1)}"

    
        return self.var , ns,ss,tpra_text,tpra, tprar_text,tprar,sescy,sescr, sescry,asescr


    def calculate_metrics_with_return(self):
        res = self.main_function_return(self.generated_compounds, self.test_set)

        return res
    

    def unique(self):
        return self.unique_compounds, self.unique_input

    
       



















