#!/usr/bin/env python3
import argparse
import gzip
import multiprocessing
import os
import sqlite3
import struct
import subprocess
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
    split = name.split('-')
    taxonomy_id = split[2]
    chunk = split[3].split("_")[0]
    files = []

    print(f"Processing {taxonomy_id}-{chunk}...")

    # # This is faster when only getting the file names and offsets, but it can't get the uncompressed sizes
    # offset = 0
    # data = subprocess.check_output(['/usr/bin/tar', '--list', '--verbose', '-f', path])
    # for line in data.decode().split('\n'):
    #     parts = line.split()
    #     if parts:
    #         size = int(parts[2])
    #         if parts[5].endswith('-F1-model_v3.cif.gz'):
    #             files.append([taxonomy_id, chunk, parts[5].split('-')[1], offset, size])
    #         offset += size + 512
    #         offset = round_to_512(offset)

    with tarfile.open(path) as tf, open(path, 'rb') as raw:
        for file in tf:
            if file.name.endswith('-F1-model_v3.cif.gz'):
                # Note - this only works as long as the biggest extracted file is <4gb. If the compressed data is >
                #  (1/1024)*gzip_size, we assume it may expand to be too big and use the thorough size calculation,
                #   but otherwise use the lazy uncompressed file size check.
                #  When written (10/31/22) the largest uncompressed file was only 2.6MB so this logic shouldn't trigger.

                if file.size > 4194304:
                    with gzip.open(tf.extractfile(file)) as gzip_file:
                        gzip_file.seek(0)
                        files.append([taxonomy_id, chunk, file.name.split('-')[1], file.offset,
                                      file.size, len(gzip_file.read())])
                else:
                    raw.seek((file.offset + 512) + (file.size - 4))
                    files.append([taxonomy_id, chunk, file.name.split('-')[1], file.offset,
                                  file.size, struct.unpack("<I", raw.read(4))[0]])
    return files


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


def get_id_mappings(download=False):
    if not os.path.exists('idmapping_selected.tab.gz') or download:
        print("Downloading Uniprot<->PDB id mapping file...")
        url = 'https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz'
        cmd = f'curl -z idmapping_selected.tab.gz -o idmapping_selected.tab.gz {url}'
        subprocess.run(cmd, shell=True, check=True)
    with gzip.open('idmapping_selected.tab.gz', 'r') as id_mapping:
        for line in id_mapping:
            datum = line.decode().split('\t')
            try:
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
                       'size numeric, expanded_size numeric);')
        cursor.executemany("INSERT INTO taxonomy_tmp(taxonomy_id, chunk, uniprot_id, offset, size, expanded_size) "
                           "VALUES (?,?,?,?,?,?)",
                           index_files(args.alphafold_path))
        sqlite_conn.commit()
        print('Building UniProt location index...')
        cursor.execute('DROP INDEX IF EXISTS uni_index;')
        cursor.execute('CREATE UNIQUE INDEX uni_index ON taxonomy_tmp(uniprot_id);')
        print('Building taxonomy ID index...')
        cursor.execute('DROP INDEX IF EXISTS taxon_index;')
        cursor.execute('CREATE UNIQUE INDEX taxon_index ON taxonomy_tmp(taxonomy_id, uniprot_id);')
        print('Building substring index on taxonomy...')
        cursor.execute('DROP INDEX IF EXISTS taxon_substr;')
        cursor.execute('CREATE INDEX taxon_substr ON taxonomy_tmp(substr(taxonomy_id, 1, 2));')
        print('Building substring index on UniProt...')
        cursor.execute('DROP INDEX IF EXISTS uniprot_substr;')
        cursor.execute('CREATE INDEX uniprot_substr ON taxonomy_tmp(substr(uniprot_id, 1, 2));')
        cursor.execute('DROP TABLE IF EXISTS taxonomy;')
        cursor.execute('ALTER TABLE taxonomy_tmp RENAME TO taxonomy;')
        sqlite_conn.commit()

        # Set up the PDB<->uniprot DB
        if args.rebuild_pdb:
            print("Doing PDB<->UniProt")
            cursor.execute('DROP TABLE IF EXISTS pdb_tmp;')
            cursor.execute('CREATE TABLE pdb_tmp (pdb_id text, uniprot_id text);')
            cursor.executemany("INSERT INTO pdb_tmp(pdb_id, uniprot_id) values (?,?)",
                               get_id_mappings(args.download_pdb))
            print('Building index on PDB IDs (pdb_id -> uniprot_id)...')
            cursor.execute('DROP INDEX IF EXISTS pdb_uniprot_index;')
            cursor.execute('CREATE UNIQUE INDEX pdb_uniprot_index ON pdb_tmp (pdb_id, uniprot_id);')
            print('Building index on PDB IDs (uniprot_id -> pdb_id)...')
            cursor.execute('DROP INDEX IF EXISTS uniprot_pdb_index;')
            cursor.execute('CREATE UNIQUE INDEX uniprot_pdb_index ON pdb_tmp (uniprot_id);')
            print('Building index on PDB ID substrings...')
            cursor.execute('DROP INDEX IF EXISTS pdb_substr;')
            cursor.execute('CREATE INDEX pdb_substr ON pdb_tmp (substr(pdb_id, 1, 2));')
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
    parser.add_argument('-n', '--no-pdb',
                        action='store_false',
                        default=True,
                        dest='rebuild_pdb',
                        help='Reload the PDB data. Default turned off.')
    args = parser.parse_args()

    create_or_update_sqlite(args)
