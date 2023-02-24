"""
Run this script to get all Chembl assays and associated compounds that meet the
criteria, and store them in the chembl_data folder. change DBNAME and USERNAME
global variable to the un and version of chembl you are using.
"""
import csv,os,pickle
import psycopg2
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem


HOSTNAME = "localhost"
DBNAME = "chembl_31"
USERNAME = "filvaler"
PASSWORD = "lera0609"


def get_targets(cursor,organisms=("Homo sapiens","Mus musculus", "Rattus norvegicus")):
    """this function gets chembl assays that meet the critera and exports their
    names, sequence, categorization, chemblid etc"""
    cursor.execute('''
    SELECT 
        DISTINCT ON (target_dictionary.tid) target_dictionary.tid,  pref_name,  
        target_dictionary.organism, 
        protein_family_classification.l1, 
        protein_family_classification.l2 
        
        FROM 
        target_dictionary JOIN target_components ON target_components.tid = target_dictionary.tid 
        JOIN component_sequences ON component_sequences.component_id = target_components.component_id 
        JOIN component_class on component_class.component_id = component_sequences.component_id 
        JOIN protein_family_classification on protein_family_classification.protein_class_id = component_class.protein_class_id 
        
        WHERE target_dictionary.organism IN ('{}');'''.format("\',\'".join(organisms)))
    targets = cursor.fetchall()
    return targets

def get_target(cursor, organisms = ("Homo sapiens")):
    print(organisms)
    target_type = "SINGLE PROTEIN"
    l1 = "Other nuclear protein"
    l2 = "Nuclear receptor"
    #l2 = "Protease"
    cursor.execute ('''
    SELECT 
    DISTINCT ON (target_dictionary.tid) target_dictionary.tid, target_dictionary.chembl_id, 
    target_dictionary.pref_name, 
    protein_family_classification.l1, 
    protein_family_classification.l2 
    
    FROM 
    target_dictionary JOIN target_components ON target_components.tid = target_dictionary.tid 
    JOIN component_class on component_class.component_id = target_components.component_id 
    JOIN protein_family_classification on protein_family_classification.protein_class_id = component_class.protein_class_id
    
    WHERE target_dictionary.organism = \'%s\' AND target_dictionary.target_type = \'%s\' AND (protein_family_classification.l1 = \'%s\' OR protein_family_classification.l2 = \'%s\')''' % (organisms, target_type, l1, l2 ))
    target = cursor.fetchall()
    #print(target)
    return target

def get_target_only_with_l2(cursor, organisms = ("Homo sapiens")):
    print(organisms)
    target_type = "SINGLE PROTEIN"
    #l1 = "Other nuclear protein"
    #l2 = "Nuclear receptor"
    l2 = "Protease"
    cursor.execute ('''
    SELECT 
    DISTINCT ON (target_dictionary.tid) target_dictionary.tid, target_dictionary.chembl_id, 
    target_dictionary.pref_name, 
    protein_family_classification.l1, 
    protein_family_classification.l2 
    
    FROM 
    target_dictionary JOIN target_components ON target_components.tid = target_dictionary.tid 
    JOIN component_class on component_class.component_id = target_components.component_id 
    JOIN protein_family_classification on protein_family_classification.protein_class_id = component_class.protein_class_id
    
    WHERE target_dictionary.organism = \'%s\' AND target_dictionary.target_type = \'%s\' AND protein_family_classification.l2 = \'%s\'''' % (organisms, target_type, l2 ))
    target = cursor.fetchall()
    #print(target)
    return target
    
def extract_assay(cursor,tid):
    """ this function extracts the compounds and their activities (expressed as
    pchembl values), it also grabs their SMILES, it also generates their
    morgan fingerprints."""
    EC50, IC50,Ki,Kd = None,None,None,None
    cursor.execute('''
    SELECT 
    activities.molregno,
    activities.standard_type,
    activities.pchembl_value,
    activities.standard_value,
    compound_structures.canonical_smiles,
    compound_structures.standard_inchi_key,
    molecule_dictionary.chembl_id AS cmpd_chembl_id,
    assays.tid,
    target_dictionary.pref_name 
    
    FROM 
    activities,
    assays,
    target_dictionary,
    target_components,
    molecule_dictionary,
    compound_properties,
    compound_structures 
    
    WHERE 
    activities.potential_duplicate = '0' AND activities.data_validity_comment is null AND activities.assay_id = assays.assay_id AND assays.tid = target_dictionary.tid AND target_dictionary.tid = target_components.tid AND molecule_dictionary.molregno = activities.molregno AND molecule_dictionary.molregno = compound_properties.molregno AND molecule_dictionary.molregno = compound_structures.molregno AND activities.standard_relation = '=' AND activities.standard_units = 'nM' AND activities.standard_value <= 1000 AND assays.confidence_score IN ('7','8', '9') AND activities.standard_type IN ('EC50', 'IC50','Ki','Kd') AND compound_properties.full_mwt>100 AND compound_properties.full_mwt < 1000 AND activities.pchembl_value > 0 AND activities.pchembl_value < 100 AND assays.tid = {}'''.format(tid))
    bioactivities = cursor.fetchall()
    if len(bioactivities)>50:
    	EC50, IC50,Ki,Kd = split_and_clean(bioactivities)
        #print(bioactivities[0])
    return EC50, IC50,Ki,Kd

def extract_assay_(cursor,tid):
    """ this function extracts the compounds and their activities (expressed as
    pchembl values), it also grabs their SMILES, it also generates their
    morgan fingerprints."""
    EC50, IC50,Ki,Kd = None,None,None,None
    cursor.execute('''
    SELECT 
    activities.molregno,
    activities.standard_type,
    activities.pchembl_value,
    activities.standard_value,
    compound_structures.canonical_smiles,
    compound_structures.standard_inchi_key,
    molecule_dictionary.chembl_id AS cmpd_chembl_id,
    assays.tid,
    target_dictionary.pref_name 
    
    FROM 
    activities,
    assays,
    target_dictionary,
    target_components,
    molecule_dictionary,
    compound_properties,
    compound_structures 
    
    WHERE 
    activities.potential_duplicate = '0' AND activities.data_validity_comment is null AND activities.assay_id = assays.assay_id AND assays.tid = target_dictionary.tid AND target_dictionary.tid = target_components.tid AND molecule_dictionary.molregno = activities.molregno AND molecule_dictionary.molregno = compound_properties.molregno AND molecule_dictionary.molregno = compound_structures.molregno AND activities.standard_relation = '=' AND activities.standard_units = 'nM' AND activities.standard_value <= 1000 AND assays.confidence_score IN ('7','8', '9') AND activities.standard_type IN ('EC50', 'IC50','Ki','Kd') AND compound_properties.full_mwt>100 AND compound_properties.full_mwt < 1000 AND assays.tid = {}'''.format(tid))
    bioactivities = cursor.fetchall()
    
    #print(bioactivities)
    EC50, IC50,Ki,Kd = split_and_clean(bioactivities)
    
    return EC50, IC50,Ki,Kd
    
def split_and_clean(bioactivities,mincompounds=50):
    """ this function takes the bioactivities output from extract_assay func and
    splits it into IC50,EC50,Ki and Kd sets. It also removes duplicates. If the
    sets contain less than mincompounds, it outputs None becayse the set is too
    small for subsequent modelbuilding """
    activities=[None,None,None,None]
    for i,assaytype in enumerate(['EC50', 'IC50','Ki','Kd']):
        compounds = [x for x in bioactivities if x[1]==assaytype]
        inchis = set([x[4] for x in compounds])
        clean_compounds = []
        for inchi in inchis:
            compound = [x for x in compounds if x[4]==inchi]
            if len(compound)>1:
                average_score = sum([x[2] for x in compound])/len(compound)
                compound = list(compound[0])
                compound[2] = average_score
                compound = tuple(compound) #convert back to tuple for consistency
            else:
                compound = compound[0]
            clean_compounds.append(compound) #convert back to tuple for consistency
        if len(clean_compounds) >= mincompounds:
            activities[i]=clean_compounds
    return tuple(activities)
    
def save_sets(compounds, tid,assaytype,name):
    """
    save sets for later usage using pickle
    """
    
    #with open("chembl_data/sets/{}/{}_{}.csv".format(name,tid,assaytype), 'w') as fp:
    #    pickle.dump(compounds, fp)
    
    out = open("chembl_data/sets/{}/{}_{}.csv".format(name,tid,assaytype),"a")
    for x in compounds:
        out.write("{0},{1},{2},{3},{4},{5},{6},{7},{8} \n".format(x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8]))
    out.close()

    #fpdict = make_fp_dict(compounds)
    #with open("chembl_data/fingerprints/{}_{}.pkl".format(tid,assaytype), 'wb') as fp:
    #    pickle.dump(fpdict, fp)
    return
    
def make_fp_dict(compounds,radius=2,fp_length=1024):
    """ make a dict (key = molregno) that returns the morgan fp with the given
    specifications"""
    fpdict = {}
    
    for compound in compounds:
        mol = Chem.MolFromSmiles(compound[4])
        if mol!=None:
            fpdict[compound[0]]=AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=fp_length)
        else:
            print("Mol was None. This should never happen. Check what's wrong with this smiles: {}".format(compound[3]))
    return fpdict
    
    
def targets2fasta(targets):
    """
    function to export a fasta file which has all the targets that occur at
    least once in the assay sets.
    """
    fp = open("chembl_data/all_targets.fasta", "w")
    for target in targets:
        fp.write(">{}|{}\n{}\n".format(target[0],target[1],target[5]))
    fp.close()
    return
    
def get_target_class(target):
    """
    input target info, output the class (gpcr, ion channel etc)
    """
    l1 = target[3]
    l2 = target[4]
    if l1 != "Membrane receptor" and l1 != "Ion channel" and l1 != "Transcription factor" and l1 != "Enzyme" and l1 != "Transporter":
        l1 = "Other non-enzyme target"
    if l1 == "Enzyme":
        l1 = l2
        if l1 != "Kinase" and l1 != "Protease":
            l1 = "Other enzyme"
    if l1 == "Membrane receptor":
        l1 = l2
        if l1 != None:
            if l1.find("G protein") != -1:
                l1 = "GPCR"
            else:
                l1 = "Other non-enzyme target"
        else:
            if target[1].find("G-protein") != -1:
                l1 = "GPCR"
            else: 
                l1 = "Other non-enzyme target"
    if l1 == "Transcription factor":
        l1 = "Nuclear receptor" #check if other nuclear protein should be included or not
    return l1



if __name__ == "__main__":
    conn = psycopg2.connect("host={} dbname={} user={} password={}".format(HOSTNAME,DBNAME,USERNAME,PASSWORD))
    cur = conn.cursor()
    #targets = get_target(cur) #with L1,L2
    targets = get_target_only_with_l2(cur)
    sets = []
    fastatargets = []
    #name = "nuclear_threshold_100nM"  
    name = "protease_threshold_1000nM"
    for p,target in enumerate(targets):
       
        EC50,IC50,Ki,Kd = extract_assay_(cur,target[0])
        #IC50,Ki,Kd = extract_assay(cur,target[0])
        
        if EC50 != None:
            sets.append(["_".join([str(target[0]),"pec50"]),target[1], target[2]])
            save_sets(EC50,target[0],"pec50", name)
        if IC50 != None:
            sets.append(["_".join([str(target[0]),"pic50"]),target[1], target[2]])
            save_sets(IC50,target[0],"pic50", name)
        if Ki != None:
            sets.append(["_".join([str(target[0]),"pki"]),target[1],target[2]])
            save_sets(Ki,target[0],"pki", name)
        if Kd != None:
            sets.append(["_".join([str(target[0]),"pkd"]),target[1],target[2]])
            save_sets(Kd,target[0],"pkd", name)
        

    #targets2fasta(fastatargets)   
    
    with open(f'chembl_data/{name}_all_sets.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=',',quotechar='"', quoting=csv.QUOTE_MINIMAL)
        #writer.writerow(["tid+assay","name","class","len","sequence"])
        writer.writerow(["tid+assay","name","class"])
        for row in sets:
            writer.writerow(row)      	
   
         	
