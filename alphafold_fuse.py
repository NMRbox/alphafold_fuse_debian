#!/usr/bin/env python

#    Copyright (C) 2006  Andrew Straw  <strawman@astraw.com>
#
#    This program can be distributed under the terms of the GNU LGPL.
#    See the file COPYING.
#

import errno
import os
import sqlite3
import stat
import sys

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

    def get_uniprot_from_taxonomy(self, taxonomy):
        self.cursor.execute('SELECT uniprot_id FROM taxonomy WHERE taxonomy_id = %s', [taxonomy])
        return [_[0] for _ in self.cursor.fetchall()]

    def get_taxonomy_from_uniprot(self, uniprot):
        self.cursor.execute('SELECT taxonomy_id FROM taxonomy WHERE uniprot_id = %s', [uniprot])
        return self.cursor.fetchone()[0]


class AlphaFoldFS(Fuse):

    def __init__(self, *args, **kw):
        Fuse.__init__(self, *args, **kw)

    # def getattr(self, path):
    #     st = MyStat()
    #     if path == '/':
    #         st.st_mode = stat.S_IFDIR | 0o755
    #         st.st_nlink = 2
    #     elif path == hello_path:
    #         st.st_mode = stat.S_IFREG | 0o444
    #         st.st_nlink = 1
    #         st.st_size = len(hello_str)
    #     else:
    #         return -errno.ENOENT
    #     return st

    def getattr(self, path):
        print('getattr', path)
        if path in ['/uniprot', '/pdb', '/taxonomy']:
            st = MyStat()
            st.st_mode = stat.S_IFDIR | 0o755
            st.st_nlink = 2
            return st
        return os.lstat("." + path)

    def readdir(self, path, offset):
        if path == "/":
            for r in '.', '..', 'uniprot', 'pdb', 'taxonomy':
                yield fuse.Direntry(r)
            return
        for r in '.', '..':
            yield fuse.Direntry(r)

    def open(self, path, flags):
        if path != hello_path:
            return -errno.ENOENT
        accmode = os.O_RDONLY | os.O_WRONLY | os.O_RDWR
        if (flags & accmode) != os.O_RDONLY:
            return -errno.EACCES

    def read(self, path, size, offset):
        print('read', path, size, offset)
        if path != hello_path:
            return -errno.ENOENT
        slen = len(hello_str)
        if offset < slen:
            if offset + size > slen:
                size = slen - offset
            buf = hello_str[offset:offset + size]
        else:
            buf = b''
        return buf


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
    server.main()


if __name__ == '__main__':
    main()
