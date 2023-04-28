31.1.2023-24.2.2023

## Important informations:
* ChEMBL 31: download 18.1.2023
* conda environment is rdkit-env(conda activate rdkit-env)
* data are download in 22.02.2023


# Nazev: Vývoj metodiky pro vyhodnocování kvality molekulových generátorů

## Vyber vhodneho biologickeho cile:
Je potreba vybrat dva biologicke cile pro tento projet, bylo zvoleno ze budou se vybirat z jadernych receptoru a proteaz. 

Pomoci scriptu Wima (*get_chembl_data_own_scripts.py*) byli vybrane nejprve targety pro jaderny receptor a pak nasledne byli stazene ligandy, ktere byli anotovane EC50,IC50,Ki,Kd a pro ketre byla jejich standrardni typ '=' a standardni hodnota
* pro jaderny receptor  <0-100> nM.
* pro proteazy <0-1000> nM 
* **podle clanku https://druggablegenome.net/ProteinFam byli zvolene thresholdy aktivit** 

* pro jaderny receptor: l1 = "Other nuclear protein",l2 = "Nuclear receptor"
* pro protease: l1 = 'enzym',l2 = 'protease'

* obecne parametry: 
    * activities.potential_duplicate = '0' 
    * activities.data_validity_comment is null  
    * activities.standard_relation = '=' 
    * activities.standard_units = 'nM' 
    * assays.confidence_score IN ('7','8', '9') 
    * activities.standard_type IN ('EC50', 'IC50','Ki','Kd') 
    * compound_properties.full_mwt>100 AND compound_properties.full_mwt < 1000 
    * activities.pchembl_value > 0 AND activities.pchembl_value < 100 

### statistika jednotlivych zastoupeni EC50,IC50,Ki,Kd
Znela otazka jestli ma cenu michat dohromady ruzne hodnoty namerenych assay (EC50,IC50,Ki,Kd)
Proto znel navrh pripravit statistiku a se podivat jakych hodnot bylo vypocitano nejvic pro jendotlive receptory
![alt text](img/statistika_zastoupeni_EC50,IC50,Ki,Kd_pro_jaderne_receptory.png)

![alt text](img/statistika_zastoupeni_EC50,IC50,Ki,Kd_pro_proteazy.png)


Na obrazku jsou zobrazene sumarni hodnoty pro jaderny receptory a proteazy.
Spocitana data pro jendotlive targety jsou ulozene v slozce */results/results_data_for_chosen_perfect_target/*
* statistics_for_nuclear_receptors.csv
* statistics_for_proteases.csv

Podle statistiky muzeme videt ze vetsinou jsou spocitane hodnoty IC50 v obou pripadech a proto jsme rozhodli misto toho abychom michali hodnoty budou vybrane jenom slouceniny s definovanou IC50 a na z nich budou vytvorene datove mnoziny.

Z ChEMBL byli stazene slouceniny s definovoanou IC50 a pak na nich byla provedene scaffoldova analyza, kterou potrebujeme pro vyber vhodneho biologickeho targetu. potrebujeme aby scaffoldova diverzita byla co nejvic. Takze aby ty scaffoldy co nejvic byli nepodobne jeden druhemu.

Byli vzate slouceniny SMILES a byli provedene upravy datasetu/sloucenin:
* odstranene slouceniny ktere neobsahuji ring system, 
* slouceniny ktere obsahuji SF6 byli prevedene na CF3, 
* slouceniny ktere obsahuji sul tak sul byla odstranena


Pro vyber typu scaffoldu ktery bude pouzit, jsme zkusili tri ruzne scaffoldy a podivali jsme se na Tanimotovou similarity a Pocet unikatnich scaffoldu. Scaffoldy ktere byli pouzite:
* Murcko scaffold - podle definice RDKit
* Generic scaffold - podle definice RDKit
* CSK (cyclic skeleton)(https://pubs.acs.org/doi/pdf/10.1021/ci200179y) - pro vytvoreni tohoto scaffldu byl pouzit nastroj RDKit - nejprve sloucenina byla prevedena na Generic_scaffold a pak nasledne na Murcko scaffold podle definice RDKit.

Nasledne jednotlive scaffoldy byli prevedene na Morganuv fingerprint r=3,  nBit = 2048, a byla spocitana Tanimotova similarita, pak vysledky byli serazene podle Tanimitova similarity od nejnizsi k nejvyssi:


Pak byli byli jednotlive smiles prevedene na morganuv fingerprint r=3,  nBit = 2048, pro ktere byla spocitana tanimotova distance(distsimilarity) takze cimn vetsi hodnota tim je to lepsi. Nejprve byla spocitana tanimotova prumerna hodnota a zastoupeni jednotlivych trid a taky byla provedene shulkova anlyza pro vizualizaci jednotlivych ligandu.
![alt text](img/nuclear_Tanim_sim_and_unic_scaf_for_murck,gener,csk.png)

![alt text](img/proteaz_Tanim_sim_and_unic_scaf_for_murck,gener,csk.png)

### Zobrazeni jednotlivych scaffoldu murcko a csk, abychom se podivali jak vypadaji
![alt text](img/nuclear_draw_scaffolds_murco_Glucocorticoid_r.png)

![alt text](img/nuclear_draw_scaffolds_csk_Glucocorticoid_r.png)

Podle nas nejlepsi vyber scaffoldu budou CSK, nejsou moc specificke jako murco a generic.
Podle vypocitaneho Tanimotova distance a poctu sloucenin a poctu uniatnich caffoldu csk zvolili jsme:
* pro jaderny receptor: **Glucocorticoid receptor** 25_pic50.csv
* pro proteazu: **Leukocyte elastase**  235_pic50.csv

Data sety z vypocitanou hodnotou Tanimota similarity a unikatnich scaffoldu jsou v slozce *results/results_data_for_chosen_perfect_target/*
* nuclear_IC50_value_for_Tanimot_similarity.csv
* protease_IC50_value_for_Tanimot_similarity.csv

Data sety ktere byli pouzite a jsou ulozene v slozce *results/results_data_for_chosen_perfect_target/*:
* nuclear_all_assay_together_EC50,IC50,Ki,Kd.csv
* protease_all_assay_together_EC50,IC50,Ki,Kd.csv
* nuclear_IC50_all_targets.csv
* protease_IC50_all_targets.csv
tyhle ctyri obsahuji vypocitane jednotlive scaffoldy (murco, generic, csk)


IC50:


# Uprava 14.4.23

## Rozdeleni na testovaci a trenovaci sady
Rozdeleni dat je napsano v jupyter *priprava_dat.ipynb*
Pro rozdeleni na jednotlive shluky byla pouzita metoda KMedoids, je to podobna metoda jako KMeans ale da se v KMedoids zmenit distancni funkci a misto mean pouziva median
Protoze chceme rozdelit dataset na jednotlive shliky s pouzitiv Jaccardove distance (1-Tanimotova similarita). A protoze to neslo zmenit v KMeans, kde je pouzita Euclidovska vzdalenost a ji nejde zmenit.
Taky jsem zkousela pouzit jinou klastrovaci metodu napriklad Rdkit metodu klastrovani, ktera bohuzel rozdelila shluky na 5, ale ten prvni shluk obsahoval 95% vsech dat a my tedy potrebujeme podobne shluky
Dalsi bode byla zkousena balicek *pyclustering*, protoze tam jde napsat vlastni distancni funkci, ale bohuzel jsem objevila ze to nefunguje protoze u prepoctu venter ty bitove vektory zacinaji byt odlisne delky, coz znamna ze mozna maji problem ve zdrojovem kode
Takze bylo rozhodnuto pouziti KMedoids pro rozdeleni na 5 shluku

Protoze potrebujeme pomerne podobnou velikost shluku, udelali jsme trochu cheating. Takze tim ze velikost shluku zavisi na random_state takze zkusila jsem 1500 ruznych random state a nasla jsem ten random state kde velikost shluku je co nejvic podobna
Je to pro gluccocorticoidni receptor:
* Random_state = 1160
* a Velikost jednotlivych shluku je 42,36,42,40,46
Rozdeleni na jednotlive shluky je ulozeno v *data/input_data/glucocor_recp_split_to_clusters_using_KMedoids.csv*

Takze rozdeleni na trenovaci a testovaci shluky:
1. Dissimilarity (exploracni vlastnost):
Takze data byli rozdeleny na 5 shluku a pokazde byl jeden shluk dan bokem jako testovaci data
jednotlive data jsou ulozene v slozce *data/input_data/*
Vzdycky u trenovaci a testovaci mnozine je napsano dis a cislo 0-4(shluk ktery byl dan bokem)
* Pro data pro Molpher *data/input_data* trenovaci
    - input_Molpher_glucocor_dis_0.csv
    - input_Molpher_glucocor_dis_1.csv
    - input_Molpher_glucocor_dis_2.csv
    - input_Molpher_glucocor_dis_3.csv
    - input_Molpher_glucocor_dis_4.csv
* Pro generatory bez ID *data/input_data* trenovaci
    - input_gener_glucocor_dis_0.csv
    - input_gener_glucocor_dis_1.csv
    - input_gener_glucocor_dis_2.csv
    - input_gener_glucocor_dis_3.csv
    - input_gener_glucocor_dis_4.csv
* Pro generatory s ID trenovaci
    - input_gener_glucocor_with_ID_dis_0.csv
    - input_gener_glucocor_with_ID_dis_1.csv
    - input_gener_glucocor_with_ID_dis_2.csv
    - input_gener_glucocor_with_ID_dis_3.csv
    - input_gener_glucocor_with_ID_dis_4.csv
* Pro generatory bez ID testovaci:
    - input_test_glucocor_dis_0.csv
    - input_test_glucocor_dis_1.csv
    - input_test_glucocor_dis_2.csv
    - input_test_glucocor_dis_3.csv
    - input_test_glucocor_dis_4.csv
* Pro generatory s ID testovaci:
    - input_test_glucocor_with_ID_dis_0.csv
    - input_test_glucocor_with_ID_dis_1.csv
    - input_test_glucocor_with_ID_dis_2.csv
    - input_test_glucocor_with_ID_dis_3.csv
    - input_test_glucocor_with_ID_dis_4.csv

2. Simmilarity (exploatacni vlastnosti):
Abychom udelali tenhle test takze byli vzate s kazdeho shluku 20% do testovaci mnoziny
* Pro data pro Molpher *data/input_data* trenovaci
    - input_Molpher_glucocor_sim_0.csv
    - input_Molpher_glucocor_sim_1.csv
    - input_Molpher_glucocor_sim_2.csv
    - input_Molpher_glucocor_sim_3.csv
    - input_Molpher_glucocor_sim_4.csv
* Pro generatory bez ID *data/input_data* trenovaci
    - input_gener_glucocor_sim_0.csv
    - input_gener_glucocor_sim_1.csv
    - input_gener_glucocor_sim_2.csv
    - input_gener_glucocor_sim_3.csv
    - input_gener_glucocor_sim_4.csv
* Pro generatory s ID trenovaci
    - input_gener_glucocor_with_ID_sim_0.csv
    - input_gener_glucocor_with_ID_sim_1.csv
    - input_gener_glucocor_with_ID_sim_2.csv
    - input_gener_glucocor_with_ID_sim_3.csv
    - input_gener_glucocor_with_ID_sim_4.csv
* Pro generatory bez ID testovaci:
    - input_test_glucocor_sim_0.csv
    - input_test_glucocor_sim_1.csv
    - input_test_glucocor_sim_2.csv
    - input_test_glucocor_sim_3.csv
    - input_test_glucocor_sim_4.csv
* Pro generatory s ID testovaci:
    - input_test_glucocor_with_ID_sim_0.csv
    - input_test_glucocor_with_ID_sim_1.csv
    - input_test_glucocor_with_ID_sim_2.csv
    - input_test_glucocor_with_ID_sim_3.csv
    - input_test_glucocor_with_ID_sim_4.csv

## Beh jednotlivych generatoru
- Tak byl pousten vypocet Molpher na lich-compute a na metacentru
s CPU=8, MEM = 14GB, wall_time = 48
- Byl zkousen Wim generators ale funguje hodne spatne

## Vypocet jednotlivych metrik
Pro vypocet jednotlivych metric byl napsan module je to v slozce *IGA_2023_try_metrics_for_other_generators/metrics_.py*
A pak byl importovan do jupyter notebook *calculate_the_metrics.ipynb*

V DrugEx je pouzivan model QSAR pro transfer learnign, takze je tam pouzita hodnota p_chebl_value pro korekci navrzenych sloucenin, takze Martin pouziva Median hodnotu kdyz jsou duplicitni hodnoty v datech a to me zarazilo takze jsem zacala zkoumat co je vlastne lepsi pouzit, takze jsem spocitala hodnoty mean a median pro svoje data a objevila jsem ze Max se lisi na 0.3 a Min -0.66 takze ta zmena neni ani na 1 lisi od sebe. Takze puvodne ve skriptu Wima byl pouzity mean, takze ho necham protoze v danem pripade ta zmena je zanedbatelna
Takze u stahovani dat u funkce split_and_clean byl pouzity mean u duplicitnich hodnotach.
Jeste nasledovne vzniklo pri zpracovani dat k duplicitnim hodnotam protoze treba data byli stazene se soli a po odstraneni ty soli ty data byli duplicitni podle SMILES ale se lisili podle INCHI, takze u tech dat byl vzat taky mean a jedna hodnota byla odstranena.
## Statisticka analyza dat a vizualizace