# AlphaFold FUSE filesystem

You are viewing a FUSE (Filesystem in userspace: https://en.wikipedia.org/wiki/Filesystem_in_Userspace)
interface for accessing the AlphaFold database. This filesystem is designed to ease access to AlphaFold
generated structures, as the 280+ million files would be difficult to organize an operate through on a traditional
filesystem. This file describes the basics of using this filesystem, organized by the relative paths inside this file system.

Ultimately, each AlphaFold structure corresponds to a single UniProt record. This tool provides different ways
to look up and access those records. Many users will find it quickest to look through the filesystem themselves, but
for those who would like to read the full documentation, it is available below. A few quick notes on data access:

1. For those who intend to read every file in the archive, iterating through the taxonomy ID folder (and subfolders) will
be the quickest due to how the AlphaFold source data files are structured. Alternatively, use the Python helper module
to access the structure files in the fastest possible manner.
2. For those who want a single structure (or some number of structures) the files can be accessed using a shortcut path, bypassing
the sudirectories: `./uniprot/$uniprot_ID_$version.cif` where `$uniprot_ID` is replaced with a UniProt ID, and `$version` is replaced with
a version string. ('v3' or 'v4' at the time this documentation was written.)

### List of files and directories in this filesystem

* ./README.md - This file
* .v3/ OR .v4/ - The top level directories will indicate which versions of the AlphaFold files are available. If a given
file is not available in the requested version, the next highest version of the file will be returned instead. For example,
version 4 consists of just 0.5% of the files in version 3. When requesting version 4, if a version 4 file exists it will be returned,
and if not version 3 will be returned. On the other hand, `v3` will not return `v4` files. (If `v2` files existed where no `v3` file 
existed it would be returned, but version 3 is the first version available.) In future paths, this part of the path will be referred to as `$version`.
* .`$version`/pdb/ - a directory which can be traversed to look up UniProt IDs associated with the PDB ID. You must first go through two subdirectories which
correspond to the third from last and second from last character in the PDB ID.
  * .`$version`/pdb/2/D/ - a directory which lists all PDB IDs (which have associated UniProt IDs) that have '2D' as the second and third from last
characters in the PDB ID. For example, this includes the ID 