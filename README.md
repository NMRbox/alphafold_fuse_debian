# AlphaFold FUSE filesystem

You are viewing a FUSE (Filesystem in userspace: 
https://en.wikipedia.org/wiki/Filesystem_in_Userspace) interface for accessing 
the AlphaFold database. This filesystem is designed to ease access to AlphaFold
generated structures, as the 280+ million files would be difficult to organize 
an operate through on a traditional
filesystem. This file describes the basics of using this filesystem, organized 
by the relative paths inside this file system.

Ultimately, each AlphaFold structure corresponds to a single UniProt record. 
This tool provides different ways to look up and access those records. Many 
users will find it quickest to look through the filesystem themselves, but
for those who would like to read the full documentation, it is available below.
Note that for large folders, it will be much quicker to run 
`ls -f --color=never` which will not sort the files or get their attributes, 
but just print their names.

A few quick notes on data access:

1. For those who intend to read every file in the archive, iterating through 
the taxonomy ID folder (and subfolders) will be the quickest due to how the 
AlphaFold source data files are structured. Alternatively, use the Python
helper module to access the structure files in the fastest possible manner.
2. For those who want a single structure (or some number of structures) the 
files can be accessed using a shortcut path, bypassing the sudirectories: 
`./uniprot/$uniprot_ID_$version.cif` where `$uniprot_ID` is replaced with a 
UniProt ID, and `$version` is replaced with a version string. ('v3' or 'v4' 
at the time this documentation was written.)

### List of files and directories in this filesystem

* ./README.md - This file
* .v3/ OR .v4/ - The top level directories will indicate which versions of the
  AlphaFold files are available. If a given file is not available in the
  requested version, the next highest version of the file will be returned
  instead. For example, version 4 consists of just 0.5% of the files in
  version 3. When requesting version 4, if a version 4 file exists it will be
  returned, and if not version 3 will be returned. On the other hand, `v3` will
  not return `v4` files. (If `v2` files existed where no `v3` file existed it
  would be returned, but version 3 is the first version available.) In future
  paths, this part of the path will be referred to as `$version`.
#### By PDB
* .`$version`/pdb/ - a directory which can be traversed to look up UniProt IDs
  associated with the PDB ID. You must first go through two subdirectories
  which correspond to the third from last and second from last character in
  the PDB ID.
    * .`$version`/pdb/2/A/ - a directory which lists all PDB IDs (which have
      associated UniProt IDs) that have `2A` as the second and third from last
      characters in the PDB ID. For example, this includes the ID `12AS`.
        * .`$version`/pdb/2/A/12AS - a directory which contains all UniProt
          IDs associated with the specified PDB ID.
            * .`$version`/pdb/2/A/12AS/P00963_v3.cif - an example of a
              structure file found in the directory above
#### By taxonomy ID
* .`$version`/taxonomy/ - a directory which can be traversed to look up
  UniProt IDs associated with the given taxonomy ID. You must first go through
  two subdirectories which
  correspond to the third from last and second from last character in the
  taxonomy ID.
    * .`$version`/taxonomy/2/3 - a directory which lists all taxonomy IDs
      (which have associated UniProt IDs) that have `23` as the second and third
      from last characters in the taxonomy ID. For example, this includes the
      ID `2206232`.
        * .`$version`/taxonomy/2/3/2022238  - a directory which contains all
          UniProt IDs associated with the specified taxonomy ID.
            *  .`$version`/taxonomy/2/3/2022238/A0A221S5L8_v3.cif - a structure
               for the taxonomy ID above
#### By UniProt ID
* .`$version`/uniprot/ - a directory which can be traversed to look up UniProt
  IDs associated with the given UniProt ID substrings. You must first go through
  two subdirectories which correspond to the third from last and second from
  last character in the UniProt IDs.
    * .`$version`/uniprot/P/F - a directory which lists all UniProt IDs that
      have `PF` as the second and third from last characters in the UniProt ID.
      For example, this includes the ID `A0A6I3SPF8`.
        * .`$version`/uniprot/P/F/A0A6I3SPF8_v3.cif  - an example structure
          file found in the above directory

#### Shortcut path
* .`$version`/uniprot/A0A496APF7 or .`$version`/uniprot/A0A496APF7.cif - a
  direct reference to the structure for a given UniProt ID - no need to go
  through the subdirectories as shown above.