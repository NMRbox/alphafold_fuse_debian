#!/usr/bin/env python3
import argparse
import gzip
import multiprocessing
import os
import sqlite3
import subprocess
import sys
import tarfile
from typing import Union, Literal


def round_to_512(number):
    if number == 0:
        return number

    remainder = number % 512

    if remainder == 0:
        return number

    return number + 512 - remainder


def get_files_from_tar(argument):
    name, path = argument
    offset = 0
    split = name.split('-')
    taxonomy_id = split[2]
    chunk = split[3].split("_")[0]
    print(f'Working on... {taxonomy_id}-{chunk}')

    files = []

    data = subprocess.check_output(['/usr/bin/tar', '--list', '--verbose', '-f', path])
    for line in data.decode().split('\n'):
        parts = line.split()
        if parts:
            size = int(parts[2])
            if parts[5].endswith('-F1-model_v3.cif.gz'):
                files.append([taxonomy_id, chunk, parts[5].replace('-F1-model_v3.cif.gz', ''), offset, size])
            offset += size + 512
            offset = round_to_512(offset)

    return files

    # with tarfile.open(entry.path) as tf:
    #     for file in tf.getmembers():
    #         if file.name.endswith('-F1-model_v3.cif.gz'):
    #             if file.name.replace('-F1-model_v3.cif.gz', '') == 'AF-A0A1Q1MKJ4':
    #                 yield taxonomy_id, chunk, file.name.replace('-F1-model_v3.cif.gz', ''), file.offset, file.size


def index_files(root_dir: str):
    # Populate DB
    def get_files_as_iterator():
        with os.scandir(root_dir) as it:
            for entry in it:
                if entry.name.endswith('.tar'):
                    yield entry.name, entry.path

    with multiprocessing.Pool(processes=250) as p:
        map = p.imap_unordered(get_files_from_tar, get_files_as_iterator(), 500)
        for result in map:
            for row in result:
                yield row


def get_id_mappings(download=False, action: Union[Literal['pdb'], Literal['uniprot']] = 'uniprot'):
    if not os.path.exists('idmapping_selected.tab.gz') or download:
        print("Downloading Uniprot<->PDB id mapping file...")
        url = 'https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz'
        cmd = f'curl -z idmapping_selected.tab.gz -o idmapping_selected.tab.gz {url}'
        subprocess.run(cmd, shell=True, check=True)
    with gzip.open('idmapping_selected.tab.gz', 'r') as id_mapping:
        for line in id_mapping:
            datum = line.decode().split('\t')
            try:
                if action == 'uniprot':
                    yield datum[0], datum[12]
                elif action == 'pdb':
                    for pdb in set([_.split(":")[0] for _ in datum[5].split('; ')]):
                        yield pdb, datum[0]
            except IndexError:
                break


def create_or_update_sqlite(args: argparse.Namespace):
    with sqlite3.connect(args.sqlite_location) as sqlite_conn:
        cursor = sqlite_conn.cursor()

        # Set up taxonomy<->uniprot DB
        print("Doing Uniprot<->Taxonomy ID")
        cursor.execute('DROP TABLE IF EXISTS taxonomy_tmp;')
        cursor.execute('CREATE TABLE taxonomy_tmp (taxonomy_id text, chunk number, uniprot_id text, offset numeric, '
                       'size numeric);')
        cursor.executemany("INSERT INTO taxonomy_tmp(taxonomy_id, chunk, uniprot_id, offset, size) values (?,?,?,?,?)",
                           index_files(args.alphafold_path))
        sqlite_conn.commit()
        print('Building index...')
        cursor.execute('DROP INDEX IF EXISTS uni_index;')
        cursor.execute('CREATE UNIQUE INDEX uni_index ON taxonomy_tmp(uniprot_id, taxonomy_id, chunk);')
        cursor.execute('DROP INDEX IF EXISTS taxon_index;')
        cursor.execute('CREATE UNIQUE INDEX taxon_index ON taxonomy_tmp(taxonomy_id, uniprot_id);')
        cursor.execute('DROP TABLE IF EXISTS taxonomy;')
        cursor.execute('ALTER TABLE taxonomy_tmp RENAME TO taxonomy;')
        sqlite_conn.commit()

        # Set up the PDB<->uniprot DB
        print("Doing PDB<->UniProt")
        cursor.execute('DROP TABLE IF EXISTS pdb_tmp;')
        cursor.execute('CREATE TABLE pdb_tmp (pdb_id text, uniprot_id text);')
        cursor.executemany("INSERT INTO pdb_tmp(pdb_id, uniprot_id) values (?,?)",
                           get_id_mappings(args.download_pdb, 'pdb'))
        cursor.execute('DROP INDEX IF EXISTS pdb_index;')
        cursor.execute('CREATE UNIQUE INDEX pdb_index ON pdb_tmp (pdb_id, uniprot_id);')
        cursor.execute('DROP TABLE IF EXISTS pdb;')
        cursor.execute('ALTER TABLE pdb_tmp RENAME TO pdb;')
        sqlite_conn.commit()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--alphafold-path',
                        dest='alphafold_path',
                        action='store',
                        default='/extra/alphafoldorig/proteomes/',
                        help='Where the source AlphaFold proteomes folder is.')
    parser.add_argument('-s', '--sql-file',
                        action='store',
                        dest='sqlite_location',
                        default='alphafold.sqlite',
                        help='Where to store the sqlite file.')
    parser.add_argument('-d', '--download',
                        action='store_true',
                        dest='download_pdb',
                        help='Force re-download the PDB index before processing')
    args = parser.parse_args()

    create_or_update_sqlite(args)
