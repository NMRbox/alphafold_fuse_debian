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
from string import ascii_uppercase, digits
from typing import Union, Literal, Optional, Generator, List

import fuse
from fuse import Fuse

if not hasattr(fuse, '__version__'):
    raise RuntimeError("your fuse-py doesn't know of fuse.__version__, probably it's too old.")

fuse.fuse_python_api = (0, 2)

fs_open = open


def get_alphanum():
    for letter in ascii_uppercase + digits:
        yield fuse.Direntry(letter)


def get_numeric():
    for number in digits:
        yield fuse.Direntry(number)


path_config = {
    'pdb': get_alphanum,
    'uniprot': get_alphanum,
    'taxonomy': get_numeric
}


@functools.lru_cache(1024)
def get_uniprot(alphafold_path: str, uniprot_id: str, taxonomy_id: str):
    logging.debug(f'Getting data for {alphafold_path} {uniprot_id} {taxonomy_id}')
    chunk = 0
    tar_path = os.path.join(alphafold_path, f'proteome-tax_id-{taxonomy_id}-{chunk}_v3.tar')
    while os.path.isfile(tar_path):
        print("Reading ", tar_path)
        with tarfile.open(tar_path) as tf:
            for member in tf:
                if member.name == f'AF-{uniprot_id}-F1-model_v3.cif.gz':
                    data = tf.extractfile(member)
                    with gzip.open(data) as uncompressed:
                        return member, uncompressed.read()
        chunk += 1
        tar_path = os.path.join(alphafold_path, f'proteome-tax_id-{taxonomy_id}-{chunk}_v3.tar')

    # Didn't find the file - had to seek the entire archive
    return None, None


@dataclasses.dataclass
class LocationAwareStat(fuse.Stat):
    tar_size: Optional[int] = None
    taxonomy_id: Optional[str] = None
    chunk: Optional[str] = None
    offset: Optional[int] = None
    uniprot_id: Optional[str] = None

    def __init__(self, **kw):
        super().__init__(**kw)
        self.st_mode = 0
        self.st_ino = 0
        self.st_dev = 0
        self.st_atime = 0
        self.st_mtime = 0
        self.st_ctime = 0
        self.st_nlink = 1

        if 'st_size' in kw:
            self.st_size = kw['st_size']
        else:
            self.st_size = 0

        if 'st_mode' in kw:
            self.st_mode = kw['st_mode']
            if kw['st_mode'] == (stat.S_IFDIR | 0o555):
                self.st_nlink = 2

        self.st_gid = os.getgid()
        self.st_uid = os.getuid()

    def __hash__(self):
        return hash(self.uniprot_id)

    def __eq__(self, other):
        return self.uniprot_id == other.uniprot_id


def dirent_gen_from_list(original_list: List) -> Generator[fuse.Direntry, None, None]:
    for record in original_list:
        yield fuse.Direntry(record)


class SQLReader:
    def __init__(self, sql_file_path):
        self.sql_file_path = sql_file_path

    def __enter__(self):
        self.sql_connection = sqlite3.connect(f'file:{self.sql_file_path}', uri=True)
        self.sql_connection.row_factory = sqlite3.Row
        self.cursor = self.sql_connection.cursor()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.cursor.close()
        self.sql_connection.close()

    def get_uniprot_from_taxonomy_substring(self, taxonomy_substring):
        self.cursor.execute('SELECT uniprot_id FROM taxonomy WHERE substr(taxonomy_id, 1, 2) = ?', [taxonomy_substring])
        return dirent_gen_from_list([_['uniprot_id'] for _ in self.cursor.fetchall()])

    @functools.lru_cache(10)
    def get_uniprot_from_taxonomy(self, taxonomy):
        self.cursor.execute('SELECT uniprot_id FROM taxonomy WHERE taxonomy_id = ?', [taxonomy])
        return dirent_gen_from_list([_['uniprot_id'] for _ in self.cursor.fetchall()])

    @functools.lru_cache(10000)
    def get_uniprot_info(self, uniprot_id) -> Union[LocationAwareStat, Literal[-2]]:
        self.cursor.execute('SELECT taxonomy_id, chunk, offset, size, expanded_size FROM taxonomy WHERE uniprot_id = ?',
                            [uniprot_id])
        data = self.cursor.fetchone()
        if not data:
            return -2
        else:
            return LocationAwareStat(st_size=data['expanded_size'],
                                     st_mode=stat.S_IFREG | 0o444,
                                     tar_size=data['size'],
                                     taxonomy_id=data['taxonomy_id'],
                                     chunk=data['chunk'],
                                     offset=data['offset'],
                                     uniprot_id=uniprot_id)


def _send_from_buffer(buffer: bytes, size: int, offset: int) -> bytes:
    slen = len(buffer)
    if offset < slen:
        if offset + size > slen:
            size = slen - offset
        buf = buffer[offset:offset + size]
    else:
        buf = b''
    return buf


class AlphaFoldFS(Fuse):

    def __init__(self, *args, **kw):
        Fuse.__init__(self, *args, **kw)
        self.sqlpath = '/extras/alphafoldorig/proteomes/'
        self.sqlite = None

    def prepare_sqlite(self):
        self.sqlite = SQLReader(self.sqlpath)

    @functools.lru_cache(50)
    def _read_uniprot_contents(self, stat_info: LocationAwareStat) -> bytes:
        logging.debug(f'Getting data for {stat_info.uniprot_id} from {stat_info.taxonomy_id} {stat_info.chunk} '
                      f'{stat_info.offset} offset {stat_info.tar_size} bytes')

        tar_path = os.path.join(self.alphafold_dir, f'proteome-tax_id-{stat_info.taxonomy_id}-{stat_info.chunk}_v3.tar')
        with fs_open(tar_path, 'rb') as tf:
            tf.seek(stat_info.offset+512)
            compressed_bytes = tf.read(stat_info.tar_size)
            return gzip.decompress(compressed_bytes)

    def _fake_filesystem_logging(self, path: str,
                                 action: Union[Literal['readdir'], Literal['getattr'], Literal['read']],
                                 size: Optional[int] = None, offset: Optional[int] = None) -> \
            Union[fuse.Stat, Literal[-2], Generator[fuse.Direntry, None, None], bytes]:
        result = self._fake_filesystem(path, action, size, offset)
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
                return dirent_gen_from_list(['.', '..', 'uniprot', 'pdb', 'taxonomy'])

        if pc[0] not in ['uniprot', 'pdb', 'taxonomy', '..', '.']:
            return -2

        # First level ('/uniprot')
        if len(pc) == 1:
            if action == 'getattr':
                return LocationAwareStat(st_mode=stat.S_IFDIR | 0o555)
            elif action == 'readdir':
                try:
                    return path_config[pc[0]]()
                except KeyError:
                    return -2
        elif len(pc) == 2:
            # First, do the directories since those are easy
            if len(pc[1]) == 1:
                if action == 'getattr':
                    return LocationAwareStat(st_mode=stat.S_IFDIR | 0o555)
                elif action == 'readdir':
                    return path_config[pc[0]]()
            # Now handle actual data requests of one sort or another
            else:
                if pc[0] == 'uniprot':
                    with self.sqlite as sql:
                        stat_info = sql.get_uniprot_info(uniprot_id=pc[1])
                        if action == 'getattr':
                            return stat_info
                        elif action == 'read':
                            return _send_from_buffer(self._read_uniprot_contents(stat_info), size, offset)
                if pc[0] == 'taxonomy':
                    if action == 'readdir':
                        with self.sqlite as sql:
                            return sql.get_uniprot_from_taxonomy(pc[1])
                    elif action == 'getattr':
                        return LocationAwareStat(st_mode=stat.S_IFDIR | 0o555)

        # Of the form /taxonomy/taxonomy_id/uniprot
        elif len(pc) == 3:
            # First, do the directories since those are easy
            if len(pc[1]) == 1:
                if action == 'getattr':
                    return LocationAwareStat(st_mode=stat.S_IFDIR | 0o555)
                elif action == 'readdir':
                    if pc[0] == 'taxonomy':
                        with self.sqlite as sql:
                            return sql.get_uniprot_from_taxonomy_substring(f'{pc[1]}{pc[2]}')
                    elif pc[0] == 'uniprot':
                        return dirent_gen_from_list(['not', 'yet', 'implemented'])

            if pc[0] == 'taxonomy':
                with self.sqlite as sql:
                    stat_info = sql.get_uniprot_info(uniprot_id=pc[2])
                if action == 'read':
                    return _send_from_buffer(self._read_uniprot_contents(stat_info), size, offset)
                elif action == 'getattr':
                    return stat_info
        elif len(pc) == 4:
            # At this level, it's always a uniprot id
            with self.sqlite as sql:
                stat_info = sql.get_uniprot_info(uniprot_id=pc[3])
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

        # TODO: Implement check that path isn't one we don't know
        # just in case something buggy calls open before calling getent. Return:
        #     return -errno.ENOENT

        # pc = path.split('/')[1:]
        # if pc[0] == 'uniprot' and len(pc) == 2:
        #     with self.sqlite as sql:
        #         taxonomy = sql.get_taxonomy_from_uniprot(pc[1])
        #         if taxonomy:
        #             metadata, data = get_uniprot(self.alphafold_dir, pc[1], taxonomy)
        #             return _send_from_buffer(data, size, offset)
        #         else:
        #             return -errno.ENOENT
        # if pc[0] == 'taxonomy':
        #     if len(pc) == 3:
        #         metadata, data = get_uniprot(self.alphafold_dir, pc[2], pc[1])
        #         return _send_from_buffer(data, size, offset)
        #     else:
        #         return -errno.ENOENT
        #
        # return -errno.ENOENT


def main():
    usage = """
    Userspace AlphaFold archive decompressing file system.
    """ + Fuse.fusage

    server = AlphaFoldFS(version="%prog " + fuse.__version__,
                         usage=usage,
                         dash_s_do='setsingle')

    server.parser.add_option(mountopt="alphafold_dir", metavar="ALPHAFOLD_PATH",
                             default='/extra/alphafoldorig/proteomes/',
                             help="Source of AlphaFold tar files [default: %default]")
    server.parser.add_option(mountopt="sqlpath", metavar="SQL_PATH",
                             default='/extra/alphafoldorig/proteomes/alphafold.sqlite',
                             help="Where to load metadata from [default: %default]")
    server.parse(values=server, errex=1)
    server.prepare_sqlite()
    server.main()


if __name__ == '__main__':
    logging.getLogger().setLevel(logging.DEBUG)
    main()
