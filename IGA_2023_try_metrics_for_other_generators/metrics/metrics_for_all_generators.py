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
        #print(x)
        #print(df.loc[x][0])
        remover = SaltRemover()
        res = remover.StripMol(Chem.MolFromSmiles(df.loc[x][0]))
        df.loc[x][0] = Chem.MolToSmiles(res)
        #print(df.loc[x][0])
        if NumAromaticRings(Chem.MolFromSmiles(df.loc[x][0])) == 0 and\
            NumAliphaticRings(Chem.MolFromSmiles(df.loc[x][0])) == 0:
            print("Not contain ring")
            print(df.loc[x][0])
            
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
    dff = pd.DataFrame(data = a)
    #print(dff)
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
    print(df)
    return df


def main_function(df_generated, df_input):
    df_input = convert_to_scaffold(df_input)
    df_generated = convert_to_scaffold(df_generated)
    unique_compounds_in_whole_generated_set = df_generated[0].unique()
    unique_input = df_input[0].unique()
    print(len(unique_compounds_in_whole_generated_set))
    print(len(unique_input))
    df = add_columns_same_like_input_function(df_generated, unique_input)
    
    df1 = add_columns_same_like_input_function(unique_compounds_in_whole_generated_set, unique_input)

    print("Pocet_unikatnich_scaffoldu", len(unique_compounds_in_whole_generated_set))
    print("Size of dataset", len(df_generated))
    print("SeScY",len(unique_compounds_in_whole_generated_set)/len(df_generated))
    try:
        print(f"TPRA:  {df[2].value_counts()[1]}/{len(df)} = {df[2].value_counts()[1]/len(df)}", )
    except:
        print(f"TPRA:  {0}/{len(df)} = {0/len(df)}", )
    print("SeScR",df1[2].sum()/len(unique_compounds_in_whole_generated_set))
    print("aSeScR",df[1].sum()/len(df_generated))
    
    


if __name__ == "__main__":
    generated_compounds = []
    input_coumpounds = []
    
    #nacitani vstupu
    with open(sys.argv[1], 'r') as f:
        for line in f.readlines():
            generated_compounds.append(line)
    df_generated = pd.DataFrame(generated_compounds)

    with open(sys.argv[2], 'r') as f:
        for line in f.readlines():
            input_coumpounds.append(line)
    df_input = pd.DataFrame(input_coumpounds)

    print(type(df_generated))
    print(type(df_input))
    
    main_function(df_generated, df_input)

    
       



















