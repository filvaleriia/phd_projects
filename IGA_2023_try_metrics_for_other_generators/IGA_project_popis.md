31.1.2023-24.2.2023

## Important informations:
* ChEMBL 31: download 18.1.2023
* conda environment is rdkit-env(conda activate rdkit-env)
* data are download in 22.02.2023


# Nazev: Vývoj metodiky pro vyhodnocování kvality molekulových generátorů

## Vyber vhodneho biologickeho cile:
Je potreba vybrat dva biologicke cile pro tento projet, bylo zvoleno ze budou se vybirat z jadernych receptoru a proteaz. 

Pomoci scriptu Wima (*get_chembl_data_own_scripts.py*) byli vybrane nejprve targety pro jaderny receptor a pak nasledne byli stazene ligandy, ktere byli anotovane EC50,IC50,Ki,Kd a pro ketre byla jejich standrardni typ '=' a standardni hodnota
pro jaderny receptor  <0-100> nM.
pro proteazy <0-1000> nM 
**podle clanku https://druggablegenome.net/ProteinFam byli zvolene thresholdy aktivit** 

### statistika jednotlivych zastoupeni EC50,IC50,Ki,Kd
Znela otazka jestli ma cenu michat dohromady ruzne hodnoty namerenych assay (EC50,IC50,Ki,Kd)
Proto znel navrh pripravit statistiku a se podivat jakych hodnot bylo vypocitano nejvic pro jendotlive receptory
![alt text](img/statistika_zastoupeni_EC50,IC50,Ki,Kd_pro_jaderne_receptory.png)

![alt text](img/statistika_zastoupeni_EC50,IC50,Ki,Kd_pro_proteazy.png)



Pro prevod na fingerprinty byli odstranene slouceniny ktere neobsahuji ring system, tak slouceniny ktere obsahuji SF6 byli prevedene na CF3, plus slouceniny ktere obsahuji sul tak sul byla odstranena
Pak byli byli jednotlive smiles prevedene na morganuv fingerprint r=3,  nBit = 2048, pro ktere byla spocitana tanimotova distance(distsimilarity) takze cimn vetsi hodnota tim je to lepsi. Nejprve byla spocitana tanimotova prumerna hodnota a zastoupeni jednotlivych trid a taky byla provedene shulkova anlyza pro vizualizaci jednotlivych ligandu.
![alt text](img/sort_tanimotova_na_jaderny_receptory.png)



Podle me je nejlepsi volba LXR-alpha ma nejvic diverznejsi scaffodly a plus ma pomerny pocet ligandu.
![alt text](img/t-sne_nuclear_receptor.png)
l1 = "Other nuclear protein"
l2 = "Nuclear receptor"

Dalsim krokem bude prozkoumani dalsiho rodiny biologickych cilu a to jsou napriklad proteozomy.

- 1.2.2023
podle klasifikace je to 
l1 = 'enzym'
l2 = 'protease'

takze zahrneme jenom l2, protoze chceme jenom protease.
Podle vypoceitane prumerne hodnotu Tanimotove vdalenosti a poctu jednotlivych zastupcu, byl vybran Cathepsin K
![alt text](img/sort_tanimotova_na_proteaze.png)


![alt text](img/t-sne_proteaze.png)

EC50:

IC50:



## Rozdeleni na testovaci a trenovaci sady


## Beh jednotlivych generatoru


## Vypocet jednotlivych metrik

## Statisticka analyza dat a vizualizace