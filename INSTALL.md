# AlphaFold FUSE

## Reading AlphaFold on-the-fly

This library implements a FUSE filesystem and associated SQLite database
generation code to allow pulling individual UniProt IDs out of the full
AlphaFold data (22TB+ at this writing) without needing to untar the million
plus tar files first.

It creates an index of which tar files contain which data files (and at what
offsets in the file) and then presents a filesystem allowing for either pulling
out a file directly by UniProt ID, or to view a list of UniProt IDs associated
with a given taxonomy ID or PDB ID.

Before running, you must use `db_builder.py` to generate a SQLite database. This will take
some time. See the help for information on options.

To run it, you must specify the `-o alphafold_dir=PATH` and `-o sqlpath=PATH` options at mount time.
For debugging, you can also specify `-f` to keep the fuse filesystem in the foreground.