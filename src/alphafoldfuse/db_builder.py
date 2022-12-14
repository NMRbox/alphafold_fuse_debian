#!/usr/bin/env python3
import argparse
import gzip
import multiprocessing
import os
import sqlite3
import struct
import subprocess
import time
from typing import Tuple, Generator, List, Union, Literal

import ciso8601


def round_to_512(number):
    if number == 0:
        return number

    remainder = number % 512

    if remainder == 0:
        return number

    return number + 512 - remainder


def get_files_from_tar(args: Tuple[str, str]) -> List[Tuple[str, str, str, int, int, int, float]]:
    """ Returns a list of lists (rows) of records from one single tar file of data. Called by the multiprocessing
    code."""
    relpath, full_path = args
    files = []

    print(f"Processing {relpath}...")

    # We manually read bytes from the tar file to determine the size of the gzipped files inside
    with open(full_path, 'rb') as raw:
        offset = 0
        # Use the tar binary to get the offsets of files inside the tar file
        data = subprocess.check_output(['/usr/bin/tar', '--list', '--verbose', '-f', full_path])
        for line in data.decode().split('\n'):
            parts = line.split()
            if parts:
                size = int(parts[2])
                if parts[5].endswith('.cif.gz') and "F1-model" in parts[5]:
                    ns = parts[5].split('-')
                    version = ns[3].split('_')[1].split('.')[0]
                    ctime = time.mktime(ciso8601.parse_datetime(f'{parts[3]} {parts[4]}').timetuple())

                    # Note - this only works as long as the biggest extracted file is <4gb. If the compressed data is >
                    #  (1/1024)*gzip_size, we assume it may expand to be too big and use the thorough size calculation,
                    #   but otherwise use the lazy uncompressed file size check.
                    #  When written (10/31/22) the largest uncompressed file was only 2.6MB so this logic shouldn't
                    #   trigger, though it has been tested.

                    if size > 4194304:
                        raw.seek(offset + 512)
                        uncompressed_size = len(gzip.decompress(raw.read(size)))
                        files.append((relpath, version.replace('v', ''), ns[1], offset, size, uncompressed_size, ctime))
                    else:
                        raw.seek((offset + 512) + (size - 4))
                        files.append((relpath, version.replace('v', ''), ns[1], offset, size,
                                      struct.unpack("<I", raw.read(4))[0], ctime))
                offset += size + 512
                offset = round_to_512(offset)

    return files


def index_files(args: argparse.Namespace) -> Generator[Tuple[str, str, str, int, int, int], None, None]:
    """ Returns a generator which spits out lists (rows) with information about structure records in the AlphaFold tar
    files which need to be inserted into the database. Uses multiprocessing to ensure that IO is the bottleneck rather
    than CPU.

    Current record format:

    [taxonomy_id, chunk_id, UniProt_id, offset (into tar file, prior to tar metadata), size (of structure inside of tar
    file), uncompressed_file_size]
    """

    # Populate DB
    def get_files_as_iterator(source_path: str):
        """Iterates through the AlphaFold directory using scandir, which is a generator which spits out file
        names as needed. os.listdir() and similar get the full list of file names before any work starts being done.

        For now assumes that the top level folder contains one or more subfolders containing versioned files."""

        if not source_path.endswith('/'):
            source_path = source_path + '/'

        with os.scandir(source_path) as root_scanner:
            for version_folder in root_scanner:
                if version_folder.is_dir():
                    with os.scandir(version_folder.path) as inner_scanner:
                        for tar in inner_scanner:
                            if tar.name.endswith('.tar'):
                                yield tar.path.replace(source_path, ''), tar.path

    with multiprocessing.Pool(processes=250) as p:
        iterator = get_files_as_iterator(args.alphafold_path)
        map = p.imap_unordered(get_files_from_tar, iterator, 500)
        for result in map:
            for row in result:
                yield row


def get_id_mappings(download=False, mode: Union[Literal['pdb'], Literal['taxonomy']] = 'pdb') -> \
        Generator[Tuple[str, str], None, None]:
    """ Returns a generator which spits out PDB_id,UniProt_id tuples. """

    if not os.path.exists('idmapping_selected.tab.gz') or download:
        print("Downloading Uniprot<->PDB id mapping file...")

        # This will redownload if the file on the server is newer
        url = 'https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz'
        cmd = f'curl -z idmapping_selected.tab.gz -o idmapping_selected.tab.gz {url}'
        subprocess.run(cmd, shell=True, check=True)

    with gzip.open('idmapping_selected.tab.gz', 'r') as id_mapping:
        for line in id_mapping:
            datum = line.decode().split('\t')
            if mode == 'taxonomy':
                yield datum[0], datum[12]
            else:
                try:
                    pdb_ids = set([_.split(":")[0] for _ in datum[5].split('; ')])
                    if pdb_ids == {''}:
                        continue
                    for pdb in pdb_ids:
                        yield datum[0], pdb
                except IndexError:
                    break


def create_or_update_sqlite(args: argparse.Namespace) -> None:
    """ Creates (if necessary) and updates the SQLite data file. """

    with sqlite3.connect(args.sqlite_location) as sqlite_conn:
        cursor = sqlite_conn.cursor()

        if args.rebuild_entries:
            # Set up taxonomy<->uniprot DB
            print("Doing Uniprot<->Taxonomy ID")
            cursor.execute('DROP TABLE IF EXISTS files_tmp;')
            cursor.execute('CREATE TABLE files_tmp (relpath text, version int, uniprot_id text,'
                           'offset numeric, size numeric, expanded_size numeric, modification_time numeric,'
                           'PRIMARY KEY(uniprot_id, version)) WITHOUT ROWID;')
            cursor.executemany("INSERT INTO files_tmp(relpath, version, uniprot_id, offset, size, "
                               "expanded_size, modification_time) VALUES (?,?,?,?,?,?,?)",
                               index_files(args))
            sqlite_conn.commit()
            print('Building substring index on UniProt...')
            cursor.execute('DROP INDEX IF EXISTS uniprot_substr;')
            cursor.execute('CREATE INDEX uniprot_substr ON files_tmp(substr(uniprot_id, -3, 2));')
            cursor.execute('DROP TABLE IF EXISTS files;')
            cursor.execute('ALTER TABLE files_tmp RENAME TO files;')
            # Now prepare the versions table
            print('Preparing versions table...')
            cursor.execute('CREATE TABLE IF NOT EXISTS versions (version int);')
            cursor.execute('DELETE FROM versions;')
            cursor.execute('INSERT INTO versions (version) SELECT DISTINCT(version) FROM files;')
            sqlite_conn.commit()

        # Set up the PDB<->uniprot DB
        if args.rebuild_pdb:
            print("Doing PDB ID mapping table.")
            cursor.execute('DROP TABLE IF EXISTS pdb_tmp;')
            cursor.execute('CREATE TABLE pdb_tmp (uniprot_id text, pdb_id text,'
                           'PRIMARY KEY (uniprot_id, pdb_id)) WITHOUT ROWID;')
            cursor.executemany("INSERT INTO pdb_tmp(uniprot_id, pdb_id) values (?,?)",
                               get_id_mappings(args.download_pdb, mode='pdb'))
            sqlite_conn.commit()

            print("Doing taxonomy ID mapping table.")
            cursor.execute('DROP TABLE IF EXISTS taxonomy_tmp;')
            cursor.execute('CREATE TABLE taxonomy_tmp (uniprot_id text, taxonomy_id text,'
                           'PRIMARY KEY (uniprot_id, taxonomy_id)) WITHOUT ROWID;')
            cursor.executemany("INSERT INTO taxonomy_tmp(uniprot_id, taxonomy_id) values (?,?)",
                               get_id_mappings(args.download_pdb, mode='taxonomy'))
            sqlite_conn.commit()

            # PDB table indexes first
            print('CREATE UNIQUE INDEX pdb_index ON pdb_tmp(pdb_id);')
            cursor.execute('DROP INDEX IF EXISTS pdb_index;')
            cursor.execute('CREATE INDEX pdb_index ON pdb_tmp(pdb_id);')
            print('CREATE INDEX pdb_substr ON pdb_tmp(substr(pdb_id, -3, 2));')
            cursor.execute('DROP INDEX IF EXISTS pdb_substr;')
            cursor.execute('CREATE INDEX pdb_substr ON pdb_tmp(substr(pdb_id, -3, 2));')
            # This is used for the second directory folder name lookup
            print('CREATE INDEX pdb_2level ON pdb(substr(pdb_id, -3, 1), pdb_id);')
            cursor.execute('DROP INDEX IF EXISTS pdb_2level;')
            cursor.execute('CREATE INDEX pdb_2level ON pdb(substr(pdb_id, -3, 1));')

            print("Creating taxonomy ID lookup table...")
            cursor.execute('DROP TABLE IF EXISTS taxonomy_unique_tmp;')
            cursor.execute('CREATE TABLE taxonomy_unique_tmp(taxonomy_id text PRIMARY KEY) WITHOUT ROWID;')
            cursor.execute('INSERT INTO taxonomy_unique_tmp(taxonomy_id) SELECT DISTINCT(taxonomy_id) FROM taxonomy;')
            print('CREATE INDEX taxon_substr ON taxonomy_unique_tmp(substr(taxonomy_id, -3, 2));')
            cursor.execute('DROP INDEX IF EXISTS taxon_substr;')
            cursor.execute('CREATE INDEX taxon_substr ON taxonomy_unique_tmp(substr(taxonomy_id, -3, 2));')

            # Taxon table indexes
            print('CREATE INDEX taxon_index ON taxonomy_tmp(taxonomy_id);')
            cursor.execute('DROP INDEX IF EXISTS taxon_index;')
            cursor.execute('CREATE INDEX taxon_index ON taxonomy_tmp(taxonomy_id);')

            print('Moving tables into position...')
            cursor.execute('DROP TABLE IF EXISTS pdb;')
            cursor.execute('ALTER TABLE pdb_tmp RENAME TO pdb;')
            cursor.execute('DROP TABLE IF EXISTS taxonomy;')
            cursor.execute('ALTER TABLE taxonomy_tmp RENAME TO taxonomy;')
            cursor.execute('DROP TABLE IF EXISTS taxonomy_unique;')
            cursor.execute('ALTER TABLE taxonomy_unique_tmp RENAME TO taxonomy_unique;')
            sqlite_conn.commit()
    print("Done!")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--alphafold-path',
                        dest='alphafold_path',
                        action='store',
                        default='/extra/alphafold/',
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
    parser.add_argument('--no-pdb',
                        action='store_false',
                        default=True,
                        dest='rebuild_pdb',
                        help='Don\'t reload the PDB ID mapping data.')
    parser.add_argument('--no-entry',
                        action='store_false',
                        default=True,
                        dest='rebuild_entries',
                        help='Don\'t reload the entry location data.')

    args = parser.parse_args()

    if not args.rebuild_entries and not args.rebuild_pdb:
        raise ValueError('You have asked to do nothing. Please specify one or neither of --no-entry and --no-pdb.')

    create_or_update_sqlite(args)
