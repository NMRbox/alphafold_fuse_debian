#!/usr/bin/env python

#    Copyright (C) 2006  Andrew Straw  <strawman@astraw.com>
#
#    This program can be distributed under the terms of the GNU LGPL.
#    See the file COPYING.
#

import errno
import functools
import gzip
import os
import sqlite3
import stat
import sys
import tarfile


# pull in some spaghetti to make this stuff work without fuse-py being installed
try:
    import _find_fuse_parts
except ImportError:
    pass
import fuse
from fuse import Fuse

if not hasattr(fuse, '__version__'):
    raise RuntimeError("your fuse-py doesn't know of fuse.__version__, probably it's too old.")

fuse.fuse_python_api = (0, 2)

hello_path = '/hello'
hello_str = b'Hello World!\n'


@functools.lru_cache(1024)
def get_uniprot(alphafold_path: str, uniprot_id: str, taxonomy_id: str):
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
    return None


class MyStat(fuse.Stat):
    def __init__(self, **kw):
        super().__init__(**kw)
        self.st_mode = 0
        self.st_ino = 0
        self.st_dev = 0
        self.st_nlink = 0
        self.st_uid = 0
        self.st_gid = 0
        self.st_size = 0
        self.st_atime = 0
        self.st_mtime = 0
        self.st_ctime = 0


class SQLReader:
    def __init__(self, sql_file_path):
        self.sql_file_path = sql_file_path

    def __enter__(self):
        self.sql_connection = sqlite3.connect(self.sql_file_path)
        self.cursor = self.sql_connection.cursor()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.cursor.close()
        self.sql_connection.close()

    @functools.lru_cache(1024)
    def get_uniprot_from_taxonomy(self, taxonomy):
        self.cursor.execute('SELECT uniprot_id FROM taxonomy WHERE taxonomy_id = ?', [taxonomy])
        return [_[0] for _ in self.cursor.fetchall()]
    @functools.lru_cache(1024)
    def get_taxonomy_from_uniprot(self, uniprot):
        self.cursor.execute('SELECT taxonomy_id FROM taxonomy WHERE uniprot_id = ?', [uniprot])
        result = self.cursor.fetchone()
        if result:
            return result[0]
        else:
            return None


class AlphaFoldFS(Fuse):

    def __init__(self, *args, **kw):
        Fuse.__init__(self, *args, **kw)

    def prepare_sqlite(self):
        self.sqlite = SQLReader(self.sqlpath)

    # def _get_expected_from_path(self, path: str, met):


    def getattr(self, path):


        print('getattr', path)
        if path == hello_path:
            st = MyStat()
            st.st_mode = stat.S_IFREG | 0o444
            st.st_nlink = 1
            st.st_size = len(hello_str)
            return st

        if path in ['/uniprot', '/pdb', '/taxonomy']:
            st = MyStat()
            st.st_mode = stat.S_IFDIR | 0o555
            st.st_nlink = 2
            st.st_gid = os.getgid()
            st.st_uid = os.getuid()
            return st
        pc = path.split('/')[1:]
        print(pc)
        if pc[0] == 'uniprot' and len(pc) == 2:
            with self.sqlite as sql:
                taxonomy = sql.get_taxonomy_from_uniprot(pc[1])
                if taxonomy:
                    st = MyStat()
                    st.st_mode = stat.S_IFREG | 0o444
                    st.st_nlink = 1
                    st.st_size = 1  # TODO: Look this up?
                    return st
                else:
                    return -errno.ENOENT
        return os.lstat("." + path)

    def readdir(self, path, offset):
        if path == "/":
            for r in '.', '..', 'uniprot', 'pdb', 'taxonomy':
                yield fuse.Direntry(r)
            return
        # for r in '.', '..':
        #     yield fuse.Direntry(r)
        #return

    def open(self, path, flags):
        print('open', path, flags)
        # Only allow reading
        accmode = os.O_RDONLY | os.O_WRONLY | os.O_RDWR
        if (flags & accmode) != os.O_RDONLY:
            return -errno.EACCES

    def read(self, path, size, offset):
        print('read', path, size, offset)
        # TODO: Implement check that path isn't one we don't know
        # just in case something buggy calls open before calling getent. Return:
        #     return -errno.ENOENT

        pc = path.split('/')[1:]
        if pc[0] == 'uniprot' and len(pc) == 2:
            print('READING uniprot!')
            with self.sqlite as sql:
                taxonomy = sql.get_taxonomy_from_uniprot(pc[1])
                if taxonomy:
                    metadata, data = get_uniprot(self.alphafold_dir, pc[1], taxonomy)

                    slen = len(data)
                    if offset < slen:
                        if offset + size > slen:
                            size = slen - offset
                        print(f'Getting data from offset {offset} of size {size}. Data length: {slen}')
                        buf = data[offset:offset + size]
                        print('Sending...', buf)
                    else:
                        buf = b''
                    return buf
                else:
                    return -errno.ENOENT

        return -errno.ENOENT


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
    main()
