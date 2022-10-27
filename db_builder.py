#!/usr/bin/env python3
import os
import sqlite3
import tarfile
import subprocess


def get_files():
    # Populate DB
    with os.scandir('proteomes') as it:
        for entry in it:
            if entry.name.endswith('.tar'):
                taxonomy_id = entry.name.replace('proteome-tax_id-', '').replace('-0_v3.tar', '')
                print(f'Working on... {taxonomy_id}')
                with tarfile.TarFile(entry.path) as tf:
                    for file in tf:
                        if file.name.endswith('-F1-model_v3.cif.gz'):
                            yield taxonomy_id, file.name.replace('-F1-model_v3.cif.gz', '')


def get_pdb_mappings():
    if not os.path.exists('idmapping.pdb'):
        url = 'https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.gz'
        cmd = f'curl {url} | gzip -dc | grep PDB > idmapping.pdb'
        subprocess.run(cmd, shell=True, check=True)
    with open('idmapping.pdb', 'r') as pdb_records:
        for line in pdb_records:
            datum = line.split()
            try:
                yield datum[2], datum[0]
            except IndexError:
                break


with sqlite3.connect('test.sqlite') as sqlite_conn:
    cursor = sqlite_conn.cursor()

    # Set up taxonomy<->uniprot DB
    cursor.execute('DROP TABLE IF EXISTS taxonomy_tmp;')
    cursor.execute('CREATE TABLE taxonomy_tmp (taxonomy_id text, uniprot_id text);')
    cursor.executemany("INSERT INTO taxonomy_tmp(taxonomy_id, uniprot_id) values (?,?)", get_files())
    cursor.execute('DROP INDEX IF EXISTS uni_index;')
    cursor.execute('CREATE INDEX uni_index ON taxonomy_tmp(uniprot_id, taxonomy_id);')
    cursor.execute('DROP TABLE IF EXISTS taxonomy;')
    cursor.execute('ALTER TABLE taxonomy_tmp RENAME TO taxonomy;')
    sqlite_conn.commit()

    # Set up the PDB<->uniprot DB
    cursor.execute('DROP TABLE IF EXISTS pdb_tmp;')
    cursor.execute('CREATE TABLE pdb_tmp (pdb_id text, uniprot_id text);')
    cursor.executemany("INSERT INTO pdb_tmp(pdb_id, uniprot_id) values (?,?)", get_pdb_mappings())
    cursor.execute('DROP INDEX IF EXISTS pdb_index;')
    cursor.execute('CREATE INDEX pdb_index ON pdb_tmp (pdb_id, uniprot_id);')
    cursor.execute('DROP TABLE IF EXISTS pdb;')
    cursor.execute('ALTER TABLE pdb_tmp RENAME TO pdb;')
    sqlite_conn.commit()
