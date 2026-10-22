# Dust opacity tables

`kext_albedo_WD_MW_3.1_60_D03.all`: B. T. Draine's table for the Milky Way R_V = 3.1 carbonaceous-silicate dust
(Weingartner & Draine 2001 case A, renormalized per Draine 2003 ARA&A 41, 241; PAH properties of Li & Draine 2001;
graphite and silicate properties of Draine 2003 ApJ 598, 1017), computed 2009 Oct 4, from
https://www.astro.princeton.edu/~draine/dust/extcurvs/kext_albedo_WD_MW_3.1_60_D03.all . Stored unmodified; read by
`jaco.bands.draine_mw31()`.

`hensley_draine_2022_astrodust_opacity.dat`: the astrodust + PAH model of Hensley & Draine (2023, ApJ 948, 55):
per-H absorption, scattering and extinction cross sections of the PAH and astrodust populations and the albedo from
0.1 micron to 3 cm, with the header's dust masses per H (M_Ad/M_H = 0.00642, M_PAH/M_H = 0.000659). The copy
distributed with meshoid (radiation/opacity_tables/), stored unmodified; read by `jaco.bands.hensley_draine_astrodust()`.
It does not extend above 12.4 eV.
