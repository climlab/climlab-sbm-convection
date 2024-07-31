import numpy as np
import pytest
from climlab_sbm_convection import betts_miller, escomp, capecalc


def test_saturation():
    temperature = 280.  # in Kelvin
    es = escomp(temperature)/100.  # convert to hPa
    #  relative tolerance for these tests ...
    tol = 1E-5
    assert es == pytest.approx(9.91189, rel=tol)
