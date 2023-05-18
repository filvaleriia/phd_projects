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

#Funkce pro opravu problemu z SF5- prevede na CF3
#funkce od Wima
rxn = AllChem.ReactionFromSmarts\
('[*:0][S:1]([F:2])([F:3])([F:4])([F:5])[F:6]>>[*:0]-[C](-[F])(-[F])-[F].[*:1]([*:2])([*:3])([*:4])([*:5])[*:6]')
def MakeScaffoldGeneric_fixed(mol):
    for i in range(len(mol.GetSubstructMatches(Chem.MolFromSmiles("S(F)(F)(F)(F)F")))):
        products = rxn.RunReactants((mol,)) # tuple
        if len(products)>0:
            mol = products[0][0]
    return MakeScaffoldGeneric(mol)


def convert_to_scaffold(df):
    print("Convert_")
    a = []
    delete_element = 0

    for x in range(len(df)):

        try:
            Chem.MolFromSmiles(df.loc[x][0])
            remover = SaltRemover()
            res = remover.StripMol(Chem.MolFromSmiles(df.loc[x][0]))
            df.loc[x][0] = Chem.MolToSmiles(res)
            if NumAromaticRings(Chem.MolFromSmiles(df.loc[x][0])) == 0 and\
                NumAliphaticRings(Chem.MolFromSmiles(df.loc[x][0])) == 0:
                dff = df.drop([x])
                delete_element = 1

            if delete_element != 1:

                try:
                    a.append(MurckoScaffoldSmiles(\
                                             Chem.MolToSmiles(MakeScaffoldGeneric\
                                           (Chem.MolFromSmiles(df.loc[x][0])))))
                except:
                    print("Faild to create scaffold_csk")
                    print("Index",x)
                    print(df.loc[x][0])                
                    try:
                        mol = MakeScaffoldGeneric_fixed(Chem.MolFromSmiles(df.loc[x][0]))
                        a.append(MurckoScaffoldSmiles(Chem.MolToSmiles(mol)))
                    except:
                        df = df.drop([x])
            delete_element = 0
        except:
            print("Nepovedlo se")

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
    
    def load(self,filepath_generated_com, filepath_test_set):
        #nacitani vstupu
        generated_compounds = []
        test_set = []
        with open(filepath_generated_com, 'r') as f:
            for line in f.readlines():
                generated_compounds.append(line)
        self.generated_compounds = pd.DataFrame(generated_compounds)

        with open(filepath_test_set, 'r') as f:
            for line in f.readlines():
                test_set.append(line)
        self.test_set = pd.DataFrame(test_set)

        #print(len(self.generated_compounds))
        #print(len(self.test_set))

    def main_function(self,df_generated, df_input):

        df_input = convert_to_scaffold(df_input)
        df_generated = convert_to_scaffold(df_generated)
        unique_compounds_in_whole_generated_set = df_generated[0].unique()
        unique_input = df_input[0].unique()

        df = add_columns_same_like_input_function(df_generated, unique_input)

        df1 = add_columns_same_like_input_function(df_generated, df_input)

        sescy = len(unique_compounds_in_whole_generated_set)/len(df_generated)
        sescr = df[2].sum()/len(unique_compounds_in_whole_generated_set)
        sescry = sescy*sescr
        asescr = df[1].sum()/len(df_generated)
        print("Pocet_unikatnich_scaffoldu", len(unique_compounds_in_whole_generated_set))
        print("Size of dataset", len(df_generated))
        print("SeScY",sescy)
        try:
            tpra = df[2].value_counts()[1]/len(df)
            print(f"TPRA:  {df[2].value_counts()[1]}/{len(df)} = {tpra}", )
        except:
            print(f"TPRA:  {0}/{len(df)} = {0/len(df)}", )

        try:
            tprar = df1[2].value_counts()[1]/len(df)
            print(f"TPRAR:  {df1[2].value_counts()[1]}/{len(df1)} = {tprar}", )
        except:
            print(f"TPRAR:  {0}/{len(df1)} = {0/len(df1)}", )
        
        print("SeScR",sescr)
        print("aSeScR",asescr)
        print("SeScY * SeScR: ", sescry)

    def main_function_return(self,df_generated, df_input):

        df_input = convert_to_scaffold(df_input)
        df_generated = convert_to_scaffold(df_generated)

        unique_compounds_in_whole_generated_set = df_generated[0].unique()
        unique_input = df_input[0].unique()

        df = add_columns_same_like_input_function(df_generated, unique_input)

        df1 = add_columns_same_like_input_function(df_generated, df_input)

        uniq_compounds = len(unique_compounds_in_whole_generated_set)
        set_size = len(df_generated)
        sescy = len(unique_compounds_in_whole_generated_set)/len(df_generated)
        sescr = df[2].sum()/len(unique_compounds_in_whole_generated_set)
        sescry = sescy*sescr
        asescr = df[1].sum()/len(df_generated)
        
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

    
        return self.var , uniq_compounds,set_size,tpra_text,tpra, tprar_text,tprar,sescy,sescr, sescry,asescr


    

    def calculate_metrics(self):
        #print(len(self.generated_compounds))
        #print(len(self.test_set))
        self.main_function(self.generated_compounds, self.test_set)

    def calculate_metrics_with_return(self):
        res = self.main_function_return(self.generated_compounds, self.test_set)

        return res

    
       



















