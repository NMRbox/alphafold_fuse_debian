#!/usr/bin/env python3
import argparse
import os
import sqlite3
import subprocess
import multiprocessing


def get_files(root_dir: str):
    # Populate DB
    with os.scandir(root_dir) as it:
        for entry in it:
            if entry.name.endswith('.tar'):
                taxonomy_id = entry.name.split('-')[2]
                print(f'Working on... {taxonomy_id}')

                data = subprocess.check_output(['/usr/bin/tar', '--list', '-f', entry.path])
                for file in data.decode().split('\n'):
                    if file.endswith('-F1-model_v3.cif.gz'):
                        yield taxonomy_id, file.replace('-F1-model_v3.cif.gz', '')

                # with tarfile.open(entry.path) as tf:
                #     for file in tf:
                #         if file.name.endswith('-F1-model_v3.cif.gz'):
                #             yield taxonomy_id, file.name.replace('-F1-model_v3.cif.gz', '')


def get_pdb_mappings(download=False):
    if not os.path.exists('idmapping.pdb') or download:
        print("Downloading Uniprot<->PDB id mapping file...")
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


def create_or_update_sqlite(args: argparse.Namespace):
    with sqlite3.connect('alphafold.sqlite') as sqlite_conn:
        cursor = sqlite_conn.cursor()

        # Set up taxonomy<->uniprot DB
        cursor.execute('DROP TABLE IF EXISTS taxonomy_tmp;')
        cursor.execute('CREATE TABLE taxonomy_tmp (taxonomy_id text, uniprot_id text);')
        cursor.executemany("INSERT INTO taxonomy_tmp(taxonomy_id, uniprot_id) values (?,?)",
                           get_files(args.alphafold_path))
        cursor.execute('DROP INDEX IF EXISTS uni_index;')
        cursor.execute('CREATE INDEX uni_index ON taxonomy_tmp(uniprot_id, taxonomy_id);')
        cursor.execute('DROP TABLE IF EXISTS taxonomy;')
        cursor.execute('ALTER TABLE taxonomy_tmp RENAME TO taxonomy;')
        sqlite_conn.commit()

        # Set up the PDB<->uniprot DB
        cursor.execute('DROP TABLE IF EXISTS pdb_tmp;')
        cursor.execute('CREATE TABLE pdb_tmp (pdb_id text, uniprot_id text);')
        cursor.executemany("INSERT INTO pdb_tmp(pdb_id, uniprot_id) values (?,?)", get_pdb_mappings(args.download_pdb))
        cursor.execute('DROP INDEX IF EXISTS pdb_index;')
        cursor.execute('CREATE INDEX pdb_index ON pdb_tmp (pdb_id, uniprot_id);')
        cursor.execute('DROP TABLE IF EXISTS pdb;')
        cursor.execute('ALTER TABLE pdb_tmp RENAME TO pdb;')
        sqlite_conn.commit()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('alphafold_path',
                        action='store',
                        default='/extra/alphafoldorig/proteomes/',
                        help='Where the source AlphaFold proteomes folder is.')
    parser.add_argument('-d', '--download',
                        action='store_true',
                        dest='download_pdb',
                        help='Force re-download the PDB index before processing')
    args = parser.parse_args()

    create_or_update_sqlite(args)
