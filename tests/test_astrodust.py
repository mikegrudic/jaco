"""The Hensley & Draine astrodust table reads into a DustModel consistent with itself and with Draine's MW model."""
import numpy as np

from jaco.bands import GIZMO_STARFORGE, band_means, draine_mw31, hensley_draine_astrodust


def test_astrodust_reads_and_is_consistent():
    d = hensley_draine_astrodust()
    assert len(d.E) == 1000 and np.all(np.diff(d.E) > 0) and np.all(d.kappa_abs > 0) and np.all(d.kappa_sca >= 0)
    assert abs(d.E.max() - 12.398) < 0.01 and d.E.min() < 1e-4
    assert abs(d.dust_to_gas - 0.00708 / 1.4) < 1e-6
    # the extinction per H rebuilt from the model against the table's own C_ext column, at every tabulated wavelength
    table = np.loadtxt(d.source)
    cext_per_H = (d.kappa_abs + d.kappa_sca) * 0.00708 * 1.6735575e-24
    assert np.allclose(cext_per_H, table[::-1, 7], rtol=2e-4)


def test_astrodust_band_means_against_draine():
    """In the non-ionizing bands the two Milky Way models agree to a factor of 2; the ionizing band is out of the
    astrodust table's range and is not compared."""
    a = band_means(hensley_draine_astrodust(), GIZMO_STARFORGE)
    dr = band_means(draine_mw31(), GIZMO_STARFORGE)
    for b in GIZMO_STARFORGE.names:
        ratio = a[b][0] / dr[b][0] if np.isscalar(a[b][0]) else None
        print(f"{b:>8s}: astrodust/Draine chi_E = {ratio}")
        if ratio is not None and b != "EUV":
            assert 0.5 < ratio < 2.0, (b, ratio)
