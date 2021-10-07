# Introduction

MolEnc: a molecular encoder using rdkit and OCaml.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3546675.svg)](https://doi.org/10.5281/zenodo.3546675)

The implemented fingerprint is J-L Faulon's "Signature Molecular Descriptor"
(SMD [1]).
This is an unfolded-counted chemical fingerprint.
Such fingerprints are less lossy than famous chemical fingerprints like ECFP4.
SMD encoding doesn't introduce feature collisions upon encoding.
Also, a feature dictionary is created at encoding time.
This dictionary can be used later on to map a given feature index to an
atom environment.
Molenc also implements unfolded-counted atom pairs [2].

For SMD, we recommend using a radius of zero to one (molenc.sh -r 0:1 ...) or
zero to two.

Currently, the atom typing scheme being used is:
(#pi-electrons, element symbol, #HA neighbors, formal charge).

In the future, we might add pharmacophore feature points[3]
(Donor, Acceptor, PosIonizable, NegIonizable, Aromatic, Hydrophobe),
to allow a fuzzier description of molecules.

# How to install the software

For beginners/non opam users:
download and execute the latest self-installer
shell script from (https://github.com/UnixJunkie/molenc/releases).

Then execute:
```
./molenc-5.0.1.sh ~/usr/molenc-5.0.1
```

This will create ~/usr/molenc-5.0.1/bin/molenc.sh, among other things
inside the same directory.

For opam users:
```
opam install molenc
```

Do not hesitate to contact the author in case you have problems installing
or using the software or if you have any question.

# Usage

```
molenc.sh -i input.smi -o output.txt
         [-d encoding.dix]: reuse existing feature dictionary
         [-r i:j]: fingerprint radius (default=0:1)
         [--pairs]: use atom pairs instead of Faulon's FP
         [-m <int>]: maximum allowed atom-pair distance
                     (default: no limit)
         [--seq]: sequential mode (disable parallelization)
         [-v]: debug mode; keep temp files
         [-n <int>]: max jobs in parallel
         [-c <int>]: chunk size
         [--no-std]: don't standardize input file molecules
                     ONLY USE IF THEY HAVE ALREADY BEEN STANDARDIZED
```

How to encode a database of molecules:

```
molenc.sh -i molecules.smi -o molecules.txt

```

How to encode another database of molecules, but reusing the feature
dictionary from another database:

```
molenc.sh -i other_molecules.smi -o other_molecules.txt -d molecules.txt.dix
```

# Bibliography

[1] Faulon, J. L., Visco, D. P., & Pophale, R. S. (2003). The signature molecular descriptor. 1. Using extended valence sequences in QSAR and QSPR studies. Journal of chemical information and computer sciences, 43(3), 707-720.

[2] Carhart, R. E., Smith, D. H., & Venkataraghavan, R. (1985). Atom pairs as molecular features in structure-activity studies: definition and applications. Journal of Chemical Information and Computer Sciences, 25(2), 64-73.

[3] Kearsley, S. K., Sallamack, S., Fluder, E. M., Andose, J. D., Mosley, R. T., & Sheridan, R. P. (1996). Chemical similarity using physiochemical property descriptors. Journal of Chemical Information and Computer Sciences, 36(1), 118-127.
