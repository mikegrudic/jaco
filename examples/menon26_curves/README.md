# Menon et al. (2026), Fig. 2: T_gas(n_H) of the one-zone free-fall collapse

Menon, Sharda, Burkhart & Grassi 2026, arXiv:2610.17577 (v1), Fig. 2: the solid curves (no FUV field, no cosmic rays)
at its seven metallicities, which `../one_zone_collapse.ipynb` compares with.

- `Tgas_nH_<Z>_chiFUV0.csv`, `<Z>` = `Z0`, `Z1e-5`, `Z1e-4`, `Z1e-3`, `Z1e-2`, `Z1e-1`, `Zsun` (Z in solar units):
  columns `log10_nH_cm-3, log10_Tgas_K`, n_H from 0.1 to 1e18 cm^-3.
- The figure is vector graphics in the paper's PDF, so these are the plotted vertices themselves (read with
  `pdftocairo -svg` and calibrated on the gridlines), not a pixel tracing: calibration residuals 3e-5 dex; matplotlib
  drops vertices within ~0.11 pt of a straight line, about 0.003 dex in T and 0.01 dex in n here.
