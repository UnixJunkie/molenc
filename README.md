# molenc

Molecular encoder using rdkit and OCaml.

The implemented fingerprint is J-L Faulon's "Signature Molecular Descriptor"
(SFP[Faulon2003][1]).
This is an unfolded-counted chemical fingerprint.
Such fingerprints are less lossy than famous chemical fingerprints like ECFP4.
SFP encoding doesn't introduce feature collisions upon encoding.
Also, a feature dictionary is created at encoding time.
This disctionary can be used later on to map a given feature index to an
atom environment.

We recommend using a radius of zero to one (molenc_d -r 0:1 ...) or
zero to two.

Currently, the fingerprint can be run using atom types
(#pi-electrons, element symbol, #HA neighbors, formal charge).

In the future, we might add pharmacophore feature points[Kearsley1996][3]
(Donor, Acceptor, PosIonizable, NegIonizable, Aromatic, Hydrophobe),
to allow a fuzzier description of molecules.
It is also planned to support atom pairs[Carhart1985][2] in addition
to or in combination with SFP.

# Bibliography

[1] Faulon, J. L., Visco, D. P., & Pophale, R. S. (2003). The signature molecular descriptor. 1. Using extended valence sequences in QSAR and QSPR studies. Journal of chemical information and computer sciences, 43(3), 707-720.

[2] Carhart, R. E., Smith, D. H., & Venkataraghavan, R. (1985). Atom pairs as molecular features in structure-activity studies: definition and applications. Journal of Chemical Information and Computer Sciences, 25(2), 64-73.

[3] Kearsley, S. K., Sallamack, S., Fluder, E. M., Andose, J. D., Mosley, R. T., & Sheridan, R. P. (1996). Chemical similarity using physiochemical property descriptors. Journal of Chemical Information and Computer Sciences, 36(1), 118-127.
