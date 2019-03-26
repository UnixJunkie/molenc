# molenc

Molecular encoder using rdkit and OCaml.

OUTDATED DESCRIPTION
The implemented fingerprint is J-L Faulon's "Signature Molecular Descriptor".
This is a counted, unfolded fingerprint of molecules.

The fingerprint can be run using atom types
(#pi-electrons, element symbol, #HA neighbors, formal charge) or rdkit
pharmacophore features (TODO)
(Donor, Acceptor, PosIonizable, NegIonizable, Aromatic, Hydrophobe),
if you want a fuzzier description of your molecules.

# Bibliography

Carhart, R. E., Smith, D. H., & Venkataraghavan, R. (1985). Atom pairs as molecular features in structure-activity studies: definition and applications. Journal of Chemical Information and Computer Sciences, 25(2), 64-73.

Kearsley, S. K., Sallamack, S., Fluder, E. M., Andose, J. D., Mosley, R. T., & Sheridan, R. P. (1996). Chemical similarity using physiochemical property descriptors. Journal of Chemical Information and Computer Sciences, 36(1), 118-127.

Faulon, J. L., Visco, D. P., & Pophale, R. S. (2003). The signature molecular descriptor. 1. Using extended valence sequences in QSAR and QSPR studies. Journal of chemical information and computer sciences, 43(3), 707-720.
