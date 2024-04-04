
## The FASTA2.x package - protein and DNA sequence similarity searching and alignment programs

4-April-2024

This repository contains a slightly updated version of the fasta2.0
programs, which were last updated in May of 2006.

This code should never be used for similarity searching research -- it
has been entirely superceded by **FASTA3**, which is available from
[github](https://github.com/wrpearson/fasta36).

The code has been made available for historical reasons, and to
provide some very legacy non-FASTA programs:

1. grease/psgrease (implementations of the Kyte-Doolittle algorithm)
2. chofas (an implementation of the Chou-Fasman secondary structure prediction algorithm)
3. garnier (an implementation of the Garnier secondary structure prediction algorithm)

The almost untouched fasta21u1d1 code from 2006 is in the initial version of
this repository.  This code runs, but does not
produce the correct E()-values with gcc8.5, probably because of some
compiler changes in the past 20 years.

The current version of the code does compile and run, with many
warnings.  I have tested `fasta`, `fastx`,`lalign`,`grease`,`chofas`,
and `garnier` and they seem to be working properly.

**FASTA2** is a direct descendent of the first `fasta` described in
1988, when most personal computing was done with floppy disks.  As a
result, it is very interactive, asking whether results should be saved
to a file, and how many results (and alignments) should be displayed.
To get a more *modern* experience, run it with the `-q` (quiet)
option.  It will then run from command line arguments like the current
version of fasta.

Except for additing function declarations to reduce the number of
warnings, this code has not been changed since 2006.

Bill Pearson
wrp@virginia.edu
