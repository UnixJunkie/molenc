# molenc

Molecular encoder using rdkit and OCaml.

The implemented fingerprint is J-L Faulon's "Signature Molecular Descriptor".
This is an unfolded-counted chemical fingerprint.
Such fingerprints are less lossy than most well-known chemical fingerprints
like ECFP4.

Some advantages: such fingerprints don't create feature collisions upon encoding.
Also, upon encoding, a feature dictionary is created. It can be used
later on to map a given feature index to an atom environment.

We recommend using a radius of zero to one (molenc_d -r 0:1 ...).

The fingerprint can be run using atom types
(#pi-electrons, element symbol, #HA neighbors, formal charge) or rdkit
pharmacophore features (TODO)
(Donor, Acceptor, PosIonizable, NegIonizable, Aromatic, Hydrophobe),
if you want a fuzzier description of your molecules.

# Bibliography

Carhart, R. E., Smith, D. H., & Venkataraghavan, R. (1985). Atom pairs as molecular features in structure-activity studies: definition and applications. Journal of Chemical Information and Computer Sciences, 25(2), 64-73.

Kearsley, S. K., Sallamack, S., Fluder, E. M., Andose, J. D., Mosley, R. T., & Sheridan, R. P. (1996). Chemical similarity using physiochemical property descriptors. Journal of Chemical Information and Computer Sciences, 36(1), 118-127.

Faulon, J. L., Visco, D. P., & Pophale, R. S. (2003). The signature molecular descriptor. 1. Using extended valence sequences in QSAR and QSPR studies. Journal of chemical information and computer sciences, 43(3), 707-720.
