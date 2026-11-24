# Fine-structure line data

Read by `jaco.models.starforge.fine_structure`.

`glover_jappsen_2007_table5.txt`: Table 5 of Glover & Jappsen 2007, ApJ 666, 1 (2007ApJ...666....1G), "Atomic data
for the fine structure transitions included in our thermal model", p. 50 of the arXiv preprint (0705.0182v2): level
statistical weights, wavelengths, transition energies and Einstein A coefficients of C, O, Si, C+ and Si+. Transcribed
as printed.

`glover_jappsen_2007_table6.txt`: Table 6 of the same paper, "Collisional de-excitation rates for atomic fine-structure
coolants", pp. 51-53 of the preprint, with its reference numbers (note at the end of the table). Transcribed as printed,
each fit written as a named functional form (see the file header), except for three rows marked `corrected`:

- C + e-, 1->0, T > 1000 K: the constant in the exponent is printed -4.44600e2. With +4.44600e2 the fit joins the
  T <= 1000 K fit at 1000 K (as the 2->0 and 2->1 fits do); as printed the rate is ~e^-890 smaller.
- C+ + ortho- and para-H2, T > 250 K: printed 5.85e-10 T^0.07 and 4.85e-10 T^0.07. Both prefactors are the T <= 250 K
  fits at 250 K, and the table's note 9 says the rate keeps the H rate's scaling (T^0.07) above 250 K, so T is taken in
  units of 250 K; as printed the rates jump by 250^0.07 = 1.47 at 250 K.

`hollenbach_mckee_1989_fe_ii_atomic.txt`, `hollenbach_mckee_1989_fe_ii_rates.txt`: Fe+ 26 um (a6D 7/2 -> 9/2), which
Glover & Jappsen do not cover, from Hollenbach & McKee 1989, ApJ 342, 306 (1989ApJ...342..306H) as tabulated by Maio et
al. 2007, MNRAS 379, 963 (Appendix B; arXiv:0704.2182) and Grassi et al. 2014, MNRAS 439, 2386 (Table B4;
arXiv:1311.1070), which agree. Same formats as the Glover & Jappsen files.

Colliders with no rate in these sources, and so left out: H2 for Si+ and Fe+, H+ for C+, Si+ and Fe+ (for C+ Glover &
Jappsen call it negligible: Coulomb repulsion), He for all.
