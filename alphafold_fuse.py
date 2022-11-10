#!/usr/bin/env python
import dataclasses
import errno
import functools
import gzip
import logging
import os
import sqlite3
import stat
import tarfile
import typing
from string import ascii_uppercase, digits
from typing import Union, Literal, Optional, Generator, List

import fuse

# No clean way to access "open" from within the class otherwise
# __builtins__.__getattribute__('open')() works but is worse than this
fs_open = open

# Check and set up fuse
if not hasattr(fuse, '__version__'):
    raise RuntimeError("You don't appear to be using the version of fuse-python specified in the requirements.")
fuse.fuse_python_api = (0, 2)


def get_alphanum() -> Generator[str, None, None]:
    for letter in ascii_uppercase + digits:
        yield fuse.Direntry(letter)


def get_numeric() -> Generator[str, None, None]:
    for number in "123456789":
        yield fuse.Direntry(number)


path_config = {
    'pdb': get_alphanum,
    'uniprot': get_alphanum,
    'taxonomy': get_numeric
}


@dataclasses.dataclass
class LocationAwareStat(fuse.Stat):
    """ Implements a Stat object, but adds some additional functionality.
    Specifically, if it references a structure in AlphaFold, it keeps track
    of where to actually get the file so that SQLite doesn't need to be
    queried again. """
    tar_size: Optional[int] = None
    relpath: Optional[str] = None
    offset: Optional[int] = None
    uniprot_id: Optional[str] = None
    version: Optional[str] = None

    def __init__(self, **kw):
        super().__init__(**kw)
        self.st_mode = 0
        self.st_ino = 0
        self.st_dev = 0
        self.st_atime = 0
        self.st_mtime = 0
        self.st_ctime = 0
        self.st_nlink = 1

        if 'modification_time' in kw:
            self.st_mtime = kw['modification_time']
            self.st_ctime = kw['modification_time']

        if 'st_size' in kw:
            self.st_size = kw['st_size']
        else:
            self.st_size = 0

        if 'st_mode' in kw:
            self.st_mode = kw['st_mode']
            if kw['st_mode'] == (stat.S_IFDIR | 0o555):
                self.st_nlink = 2
        else:
            self.st_mode = stat.S_IFREG | 0o444

        self.st_gid = os.getgid()
        self.st_uid = os.getuid()

    def __hash__(self):
        """ This and __eq__ are implemented so that this class can be used as an input
        to a function which has functools.lru_cache applied on it. Basically, if the
        AlphaFold data is the same, then two LocationAwareStat are the same. """
        return hash(self.uniprot_id)

    def __eq__(self, other):
        return self.uniprot_id == other.uniprot_id


def dirent_gen_from_list(original_list: List[str]) -> Generator[fuse.Direntry, None, None]:
    """ Takes a list of strings and returns a generator which yields Direntries. """
    for record in original_list:
        yield fuse.Direntry(record)


def dirent_gen_from_result(result: List[sqlite3.Row]) -> Generator[fuse.Direntry, None, None]:
    """ Presents the result of the query like UniProt mmCIF files, including version.
    The input should be a cursor which when iterated through provides ['uniprot_id']
    and ['version'] on each item returned."""

    for record in result:
        yield fuse.Direntry(f'{record["uniprot_id"]}_v{record["version"]}.cif')


class SQLReader:
    """ A class to help get data out of the SQLite database. """

    def __init__(self, sql_file_path):
        self.sql_file_path = sql_file_path

    def __enter__(self):
        # Open the file read-only (uri=True) and use the row factory
        self.sql_connection = sqlite3.connect(f'file:{self.sql_file_path}', uri=True)
        self.sql_connection.row_factory = sqlite3.Row
        self.cursor = self.sql_connection.cursor()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.cursor.close()
        self.sql_connection.close()

    def get_versions(self):
        self.cursor.execute('SELECT DISTINCT(version) as version FROM versions;')
        return [_['version'] for _ in self.cursor.fetchall()]

    def get_valid_pdb_dirnames_l2(self, level_1: str, version: str):
        self.cursor.execute('''
SELECT DISTINCT(substr(pdb_id,-2,1)) AS pdb_id FROM pdb
                                               LEFT JOIN files f ON pdb.uniprot_id = f.uniprot_id
                                               WHERE substr(pdb_id , -3, 1) = ?
                                               AND version <= ?;''',
                            [level_1.upper(), version])
        return dirent_gen_from_list([_['pdb_id'] for _ in self.cursor.fetchall()])

    def get_pdb_from_pdb_substring(self, pdb_substring: str, version: str):
        self.cursor.execute('''
SELECT DISTINCT (pdb_id) AS pdb_id
FROM pdb
         INNER JOIN files f ON pdb.uniprot_id = f.uniprot_id
WHERE substr(pdb.pdb_id, -3, 2) = ?
  AND f.version <= ?;''',
                            [pdb_substring.upper(), version])
        return dirent_gen_from_list([_['pdb_id'] for _ in self.cursor.fetchall()])

    def get_uniprot_from_uniprot_substring(self, uniprot_substring: str, version: str):
        self.cursor.execute('''
SELECT uniprot_id, MAX(version) AS version
FROM files
WHERE substr(uniprot_id, -3, 2) = ?
  AND version <= ?
GROUP BY uniprot_id''', [uniprot_substring, version])
        return dirent_gen_from_result(self.cursor.fetchall())

    def get_taxonomy_from_taxonomy_substring(self, taxonomy_substring: str, version: str):
        self.cursor.execute('''
SELECT DISTINCT(taxonomy_id) AS taxonomy_id FROM taxonomy
                                    LEFT JOIN files f on taxonomy.uniprot_id = f.uniprot_id
                                    WHERE substr(taxonomy_id, -3, 2) = ?
AND f.version <= ?;''',
                            [taxonomy_substring, version])
        return dirent_gen_from_list([_['taxonomy_id'] for _ in self.cursor.fetchall()])

    def get_uniprot_from_taxonomy(self, taxonomy: str, version: str):
        self.cursor.execute('''
SELECT taxonomy.uniprot_id, MAX(files.version) AS version
FROM taxonomy
         LEFT JOIN files ON taxonomy.uniprot_id = files.uniprot_id
WHERE taxonomy_id = ?
  AND files.version <= ?
GROUP BY taxonomy.uniprot_id;''', [taxonomy, version])
        return dirent_gen_from_result(self.cursor.fetchall())

    def get_uniprot_from_pdb(self, pdb: str, version: str):
        self.cursor.execute('''
SELECT pdb.uniprot_id, MAX(files.version) AS version
FROM pdb
         LEFT JOIN files ON pdb.uniprot_id = files.uniprot_id
WHERE pdb.pdb_id = ?
  AND files.version <= ?
GROUP BY pdb.uniprot_id;''', [pdb.upper(), version])

        return dirent_gen_from_result(self.cursor.fetchall())

    @functools.lru_cache(10000)
    def get_uniprot_info(self, uniprot_id, max_version: Optional[int] = None) -> Union[LocationAwareStat, Literal[-2]]:
        """ Load info for one particular UniProt ID from SQLite. Cache recent results. """

        uniprot_id = uniprot_id.replace('.cif', '')
        if '_' in uniprot_id:
            version = uniprot_id.split('_')[1].replace('v', '')
            # They can't request a newer version than the FS allows
            if version > max_version:
                return -2
        else:
            version = None
        uniprot_id = uniprot_id.split('_')[0]

        sql = 'SELECT relpath, offset, size, expanded_size,modification_time, max(version) AS version FROM files ' \
              'WHERE uniprot_id = ?'
        args = [uniprot_id]
        if version:
            sql += ' AND version = ?'
            args.append(version)
        if max_version:
            sql += ' AND version <= ?'
            args.append(max_version)

        self.cursor.execute(sql, args)
        data = self.cursor.fetchone()
        if not data:
            return -2

        return LocationAwareStat(st_size=data['expanded_size'],
                                 tar_size=data['size'],
                                 relpath=data['relpath'],
                                 offset=data['offset'],
                                 uniprot_id=uniprot_id,
                                 modification_time=data['modification_time'],
                                 version=data['version'])


def _send_from_buffer(buffer: bytes, size: int, offset: int) -> bytes:
    slen = len(buffer)
    if offset < slen:
        if offset + size > slen:
            size = slen - offset
        buf = buffer[offset:offset + size]
    else:
        buf = b''
    return buf


class AlphaFoldFS(fuse.Fuse):
    versions: List[str]

    def __init__(self, *args, **kw):
        fuse.Fuse.__init__(self, *args, **kw)
        # This gets overwritten (inside the Fuse code) if specified in the arguments
        self.versions = []
        self.sqlpath = '/extras/alphafold/'
        self.sqlite = None

    def prepare_sqlite(self):
        self.sqlite = SQLReader(self.sqlpath)
        with self.sqlite as sql:
            self.versions = sql.get_versions()

    @functools.lru_cache(50)
    def _read_uniprot_contents(self, stat_info: LocationAwareStat) -> bytes:
        logging.debug(f'Getting data for {stat_info.uniprot_id} from {stat_info.relpath} '
                      f'{stat_info.offset} offset {stat_info.tar_size} bytes')

        with fs_open(os.path.join(self.alphafold_dir, stat_info.relpath), 'rb') as tf:
            tf.seek(stat_info.offset+512)
            compressed_bytes = tf.read(stat_info.tar_size)
            return gzip.decompress(compressed_bytes)

    def _fake_filesystem_logging(self, path: str,
                                 action: Union[Literal['readdir'], Literal['getattr'], Literal['read']],
                                 size: Optional[int] = None, offset: Optional[int] = None) -> \
            Union[fuse.Stat, Literal[-2], Generator[fuse.Direntry, None, None], bytes]:
        result = self._fake_filesystem(path, action, size, offset)
        if isinstance(result, typing.Generator):
            logging.debug(f'{action}: {path, size, offset, list(result)}')
        else:
            logging.debug(f'{action}: {path, size, offset, result}')
        return result

    def _fake_filesystem(self, path: str,
                         action: Union[Literal['readdir'], Literal['getattr'], Literal['read']],
                         size: Optional[int] = None, offset: Optional[int] = None) -> \
            Union[fuse.Stat, Literal[-2], Generator[fuse.Direntry, None, None], bytes]:

        # Check for bugs in implementation
        if action == 'read':
            assert size is not None and offset is not None

        # First determine what they want
        pc = path.split("/")[1:]

        # Handle the root directory getattr
        if pc[0] == '':
            if action == 'getattr':
                return LocationAwareStat(st_mode=stat.S_IFDIR | 0o555)
            elif action == 'readdir':
                return dirent_gen_from_list(self.versions + ['README.md'])
        if pc[0] not in self.versions + ['README.md']:
            return -2

        # They want the readme
        if pc[0] == 'README.md':
            if action == 'getattr':
                with open('README.md', 'r') as readme:
                    return LocationAwareStat(st_size=len(readme.read()))
            if action == 'read':
                with open('README.md', 'rb') as readme:
                    return _send_from_buffer(readme.read(), size, offset)

        # Handle the version part of the FS path
        if len(pc) == 1:
            if action == 'getattr':
                return LocationAwareStat(st_mode=stat.S_IFDIR | 0o555)
            elif action == 'readdir':
                return dirent_gen_from_list(['pdb', 'uniprot', 'taxonomy'])

        # Figure out the version, and shift the path
        version = pc[0].replace('v', '')
        pc = pc[1:]

        # First level ('/uniprot')
        if len(pc) == 1:
            if action == 'getattr':
                return LocationAwareStat(st_mode=stat.S_IFDIR | 0o555)
            elif action == 'readdir':
                try:
                    return path_config[pc[0]]()
                except KeyError:
                    return -2
        # Of the form /pdb/2 or /uniprot/UNIPROT_ID or /taxonomy/taxonomy_ID
        elif len(pc) == 2:
            # First, do the directories since those are easy
            if len(pc[1]) == 1:
                if action == 'getattr':
                    return LocationAwareStat(st_mode=stat.S_IFDIR | 0o555)
                elif action == 'readdir':
                    if pc[0] == 'pdb':
                        with self.sqlite as sql:
                            return sql.get_valid_pdb_dirnames_l2(pc[1], version)
                    else:
                        return path_config[pc[0]]()
            # Now handle actual data requests of one sort or another
            #  These are the direct reference short-cuts, bypassing the directory slices by character
            else:
                if pc[0] == 'uniprot':
                    with self.sqlite as sql:
                        stat_info = sql.get_uniprot_info(uniprot_id=pc[1], max_version=version)
                        if action == 'getattr':
                            return stat_info
                        elif action == 'read':
                            return _send_from_buffer(self._read_uniprot_contents(stat_info), size, offset)
                elif pc[0] == 'taxonomy':
                    if action == 'readdir':
                        with self.sqlite as sql:
                            return sql.get_uniprot_from_taxonomy(pc[1], version=version)
                    elif action == 'getattr':
                        return LocationAwareStat(st_mode=stat.S_IFDIR | 0o555)
                elif pc[0] == 'pdb':
                    if action == 'readdir':
                        with self.sqlite as sql:
                            return sql.get_uniprot_from_pdb(pc[1], version=version)
                    elif action == 'getattr':
                        return LocationAwareStat(st_mode=stat.S_IFDIR | 0o555)
        # Of the form /taxonomy/taxonomy_id/uniprot or /pdb/2/D
        elif len(pc) == 3:
            # First, do the directories since those are easy
            if len(pc[1]) == 1:
                if action == 'getattr':
                    return LocationAwareStat(st_mode=stat.S_IFDIR | 0o555)
                elif action == 'readdir':
                    if pc[0] == 'taxonomy':
                        with self.sqlite as sql:
                            return sql.get_taxonomy_from_taxonomy_substring(f'{pc[1]}{pc[2]}', version=version)
                    elif pc[0] == 'uniprot':
                        with self.sqlite as sql:
                            return sql.get_uniprot_from_uniprot_substring(f'{pc[1]}{pc[2]}', version=version)
                    elif pc[0] == 'pdb':
                        with self.sqlite as sql:
                            return sql.get_pdb_from_pdb_substring(f'{pc[1]}{pc[2]}', version=version)
            if pc[0] == 'taxonomy':
                with self.sqlite as sql:
                    stat_info = sql.get_uniprot_info(uniprot_id=pc[2], max_version=version)
                if action == 'read':
                    return _send_from_buffer(self._read_uniprot_contents(stat_info), size, offset)
                elif action == 'getattr':
                    return stat_info
        # At this level, they are looking inside a PDB ID folder, Taxonomy ID folder, or directly at a uniprot
        # Of the form /pdb/2/D/2DOG or /uniprot/C/4/C4K3Z3
        elif len(pc) == 4:
            if pc[0] == 'uniprot':
                # For uniprot, this is the file level
                with self.sqlite as sql:
                    stat_info = sql.get_uniprot_info(uniprot_id=pc[3], max_version=version)
                if action == 'getattr':
                    return stat_info
                elif action == 'read':
                    return _send_from_buffer(self._read_uniprot_contents(stat_info), size, offset)

            # For taxonomy and PDB, it's a directory
            if action == 'getattr':
                return LocationAwareStat(st_mode=stat.S_IFDIR | 0o555)

            if action == 'readdir':
                if pc[0] == 'taxonomy':
                    with self.sqlite as sql:
                        return sql.get_uniprot_from_taxonomy(pc[3], version=version)
                elif pc[0] == 'pdb':
                    with self.sqlite as sql:
                        return sql.get_uniprot_from_pdb(pc[3], version=version)
        # Of the form /pdb/2/D/2DOG/C4K3Z3
        elif len(pc) == 5:
            # At this level, it's always a uniprot id
            with self.sqlite as sql:
                stat_info = sql.get_uniprot_info(uniprot_id=pc[4], max_version=version)
            if action == 'getattr':
                return stat_info
            elif action == 'read':
                return _send_from_buffer(self._read_uniprot_contents(stat_info), size, offset)
        else:
            return -2

    def getattr(self, path):
        return self._fake_filesystem_logging(path, 'getattr')

    def readdir(self, path, offset):
        return self._fake_filesystem_logging(path, 'readdir')

    def open(self, path, flags):
        logging.debug(f'open {path} {flags}')
        # Only allow reading
        accmode = os.O_RDONLY | os.O_WRONLY | os.O_RDWR
        if (flags & accmode) != os.O_RDONLY:
            return -errno.EACCES

    def read(self, path, size, offset):
        return self._fake_filesystem_logging(path, 'read', size, offset)


def main():
    usage = """
    Userspace AlphaFold archive decompressing file system.
    """ + fuse.Fuse.fusage

    server = AlphaFoldFS(version="%prog " + fuse.__version__,
                         usage=usage,
                         dash_s_do='setsingle')

    server.parser.add_option(mountopt="alphafold_dir",
                             metavar="ALPHAFOLD_PATH",
                             default='/extra/alphafold/',
                             help="Source of AlphaFold tar files [default: %default]")
    server.parser.add_option(mountopt="sqlpath", metavar="SQL_PATH",
                             default='/extra/alphafold/alphafold.sqlite',
                             help="Where to load metadata from [default: %default]")
    server.parse(values=server, errex=1)
    server.prepare_sqlite()
    server.main()


if __name__ == '__main__':
    logging.getLogger().setLevel(logging.DEBUG)
    main()
