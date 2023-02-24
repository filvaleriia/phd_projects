18.1.2023-24.2.2023-
# Postgres databaze

## Version of database
* ChEMBL 31: download 18.1.2023

# Content

## Folders
1. chembl_data
2. chembl_db_postgres

## Scripts
1. get_chembl_data_Wim.py - original script by Wim, didn't rewrite
2. get_chembl_data_own_scripts.py - rewrite Wim scripts
3. all_set_connect.sh - connect different sets to one


### chembl_data/
in this folder download all sets after run scripts *get_chembl_data_own_scripts.py*
contains folders: sets,fingerprints


### chembl_db_postgres

chembl_db_postgres/chembl_31/chembl_31_postgresql

'''
for run DB:
pg_ctl start -D chembl_31


for end DB:
pg_ctl stop -D chembl_31

'''
