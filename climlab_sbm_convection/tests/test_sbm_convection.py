import numpy as np
import pytest
from climlab_sbm_convection import betts_miller, escomp, capecalc


### Define a sounding for testing purposes
#  Physical constants
Cp_air = 1004.0
rdgas = 287.0
rvgas = 461.5
HLv = 2500000.0
kappa = 0.2858565737051793
Grav = 9.8
#  Grid size
num_lev = 75
#  Pressure at mid-levels (Pa)
pfull = np.array([ 7000.,  7070.,  7400.,  7780.,  7860.,  8070.,  8170.,  8290.,
        8590.,  8640.,  9030.,  9540., 10000., 10300., 10800., 11030.,
       11100., 12170., 12790., 13800., 14000., 14110., 14800., 15000.,
       15500., 16380., 16800., 17210., 19600., 19960., 20000., 21970.,
       25000., 25350., 26700., 27780., 30000., 31750., 34650., 35000.,
       39200., 39350., 40000., 41000., 41020., 44530., 46300., 48290.,
       50000., 50250., 52270., 52300., 54340., 55400., 55900., 56500.,
       58730., 58900., 60980., 65730., 70000., 73460., 74400., 76160.,
       76300., 78920., 79200., 80700., 81790., 82300., 84400., 85000.,
       87830., 90300., 92300.])
#  Pressure at interfaces (Pa)
phalf = np.array([     0.,   7035.,   7235.,   7590.,   7820.,   7965.,   8120.,
         8230.,   8440.,   8615.,   8835.,   9285.,   9770.,  10150.,
        10550.,  10915.,  11065.,  11635.,  12480.,  13295.,  13900.,
        14055.,  14455.,  14900.,  15250.,  15940.,  16590.,  17005.,
        18405.,  19780.,  19980.,  20985.,  23485.,  25175.,  26025.,
        27240.,  28890.,  30875.,  33200.,  34825.,  37100.,  39275.,
        39675.,  40500.,  41010.,  42775.,  45415.,  47295.,  49145.,
        50125.,  51260.,  52285.,  53320.,  54870.,  55650.,  56200.,
        57615.,  58815.,  59940.,  63355.,  67865.,  71730.,  73930.,
        75280.,  76230.,  77610.,  79060.,  79950.,  81245.,  82045.,
        83350.,  84700.,  86415.,  89065.,  91300., 100000.])
#  Air temperature (K)
tin = np.array([208.25, 208.05, 209.35, 210.75, 211.05, 210.65, 209.05, 207.25,
       206.25, 206.05, 207.35, 208.85, 208.65, 208.85, 210.65, 210.15,
       210.05, 211.15, 211.75, 212.65, 212.05, 211.75, 209.85, 209.85,
       210.25, 208.45, 207.65, 208.35, 212.25, 212.95, 213.05, 217.35,
       223.25, 223.95, 226.45, 228.85, 233.45, 236.95, 242.25, 242.85,
       249.05, 249.25, 249.85, 250.25, 250.25, 255.15, 257.45, 260.55,
       263.05, 263.35, 266.05, 266.05, 267.15, 267.65, 267.05, 267.45,
       270.25, 270.45, 273.05, 278.65, 283.35, 287.15, 288.15, 289.05,
       289.15, 291.35, 291.55, 290.15, 290.45, 290.55, 289.75, 290.35,
       292.85, 294.95, 297.55])
#  Mixing ratio
rin = np.array([2.51036657e-06, 2.39861183e-06, 2.46047788e-06, 2.55646085e-06,
       2.57537379e-06, 2.78643746e-06, 2.43442613e-06, 2.11972918e-06,
       2.31539671e-06, 2.34287389e-06, 2.40466005e-06, 2.44071303e-06,
       2.24865309e-06, 2.26062093e-06, 2.93907018e-06, 2.64231448e-06,
       2.58103689e-06, 3.77341738e-06, 4.59225623e-06, 6.15337056e-06,
       6.35942270e-06, 6.40993055e-06, 6.92627123e-06, 6.83391993e-06,
       7.03759057e-06, 5.60908726e-06, 5.05454545e-06, 5.50875002e-06,
       8.78820526e-06, 9.57819504e-06, 9.70182772e-06, 1.21851555e-05,
       1.63939345e-05, 1.73324413e-05, 2.07866655e-05, 1.99785187e-05,
       1.85000643e-05, 1.97461366e-05, 2.15300279e-05, 2.18877891e-05,
       3.81802066e-05, 4.36143675e-05, 7.56535070e-05, 1.90317030e-04,
       1.90224209e-04, 2.73274292e-04, 3.22707912e-04, 3.18968156e-04,
       3.17549663e-04, 2.97328696e-04, 1.82927447e-04, 1.80893822e-04,
       3.13515979e-04, 4.09766859e-04, 1.16770878e-03, 1.42359202e-03,
       2.01905507e-03, 2.08084097e-03, 2.25401288e-03, 2.66127096e-03,
       3.04095053e-03, 3.87636203e-03, 4.12588208e-03, 5.29536808e-03,
       5.40208449e-03, 4.67974268e-03, 4.62895382e-03, 1.02771488e-02,
       1.03451661e-02, 1.03494641e-02, 1.13785111e-02, 1.14478408e-02,
       1.16744206e-02, 1.18065917e-02, 1.36767697e-02])
#  Specific humidity
qin = rin / (1+rin)
#  A different profile with moister surface (unstable)
qin_moist = np.array([5.93771129e-06, 5.68577791e-06, 5.80706944e-06, 6.00110638e-06,
       6.03901392e-06, 6.49224559e-06, 5.71463886e-06, 5.01374616e-06,
       5.43522603e-06, 5.49382334e-06, 5.61464845e-06, 5.67469660e-06,
       5.23923820e-06, 5.25597979e-06, 6.70673839e-06, 6.06059720e-06,
       5.92616209e-06, 8.42302362e-06, 1.01030739e-05, 1.32494765e-05,
       1.36556725e-05, 1.37516052e-05, 1.47523832e-05, 1.45556830e-05,
       1.49360133e-05, 1.20225823e-05, 1.08834829e-05, 1.17861379e-05,
       1.81724573e-05, 1.96899848e-05, 1.99275199e-05, 2.45823173e-05,
       3.23065317e-05, 3.40266435e-05, 4.02943267e-05, 3.87277717e-05,
       3.58618542e-05, 3.80270204e-05, 4.10783720e-05, 4.17019517e-05,
       7.02271908e-05, 7.96567986e-05, 1.34236641e-04, 3.22264744e-04,
       3.22107587e-04, 4.52989922e-04, 5.29771424e-04, 5.22887110e-04,
       5.19823623e-04, 4.88115344e-04, 3.06741575e-04, 3.03485816e-04,
       5.11536466e-04, 6.59719040e-04, 1.79220809e-03, 2.16484229e-03,
       3.01856662e-03, 3.10646106e-03, 3.34835210e-03, 3.91258577e-03,
       4.43348176e-03, 5.58148861e-03, 5.92172553e-03, 7.51144610e-03,
       7.65561374e-03, 6.66435178e-03, 6.59419730e-03, 1.41233421e-02,
       1.42051974e-02, 1.42074340e-02, 1.55371996e-02, 1.56232526e-02,
       1.58981688e-02, 1.60528379e-02, 1.84526661e-02])

def test_saturation():
    temperature = 280.  # in Kelvin
    es = escomp(temperature)/100.  # convert to hPa
    #  relative tolerance for these tests ...
    tol = 1E-5
    assert es == pytest.approx(9.91189, rel=tol)

def test_capecalc():
    es0 = 1.0
    avgbl = False  # If true, the parcel is averaged in theta and r up to its LCL -- not actually implemented
    cape, cin, tp,rp,klzb = capecalc(num_lev,pfull,phalf,
                        Cp_air,rdgas,rvgas,HLv,kappa,es0,tin,rin,avgbl)
    #  relative tolerance for these tests ...
    tol = 1E-6
    assert cape == pytest.approx(2605.31640625, rel=tol)
    assert cin == pytest.approx(122.73005676269531, rel=tol)
    assert klzb == 28

def test_betts_miller():
    seconds_per_hour = 3600
    dt = seconds_per_hour * 3
    tau_bm=7200.
    rhbm=0.7
    do_simp=False
    do_shallower=True
    do_changeqref=True
    do_envsat=True
    do_taucape=False
    capetaubm=900.  # only used if do_taucape == True
    tau_min=2400.   # only used if do_taucape == True
    ix = 1; jx = 1; kx = num_lev
    es0 = 1.0

    #  Convert input arrays to the expected dimensionality [ix,jx,kx]
    tin_grid = tin[np.newaxis, np.newaxis, :]
    qin_grid = qin[np.newaxis, np.newaxis, :]
    qin_grid_moist = qin_moist[np.newaxis, np.newaxis, :]
    pfull_grid = pfull[np.newaxis, np.newaxis, :]
    phalf_grid = phalf[np.newaxis, np.newaxis, :]
    rhbm_grid = rhbm * np.ones_like(tin_grid)

    rain, tdel, qdel, q_ref, bmflag, klzbs, cape, cin, t_ref, invtau_bm_t, invtau_bm_q, capeflag = \
        betts_miller(dt, tin_grid, qin_grid, pfull_grid, phalf_grid, 
                                          HLv, Cp_air, Grav,
                                          rdgas,rvgas,kappa, es0, tau_bm, rhbm_grid, 
                                          do_simp, do_shallower, do_changeqref, 
                                          do_envsat, do_taucape, capetaubm, tau_min, 
                                          ix, jx, kx,)
    tol = 1E-6
    assert cape == pytest.approx(2605.31640625, rel=tol)
    assert cin == pytest.approx(122.73005676269531, rel=tol)
    assert rain[0,0] == 0.
    assert bmflag[0,0] == 1.
    assert capeflag[0,0] == 0.
    assert np.all(tdel) == 0.
    assert np.all(qdel) == 0.

    # Now destabilize the surface with added moisture
    rain, tdel, qdel, q_ref, bmflag, klzbs, cape, cin, t_ref, invtau_bm_t, invtau_bm_q, capeflag = \
    betts_miller(dt, tin_grid, qin_grid_moist, pfull_grid, phalf_grid, 
                                          HLv, Cp_air, Grav,
                                          rdgas,rvgas,kappa, es0, tau_bm, rhbm_grid, 
                                          do_simp, do_shallower, do_changeqref, 
                                          do_envsat, do_taucape, capetaubm, tau_min, 
                                          ix, jx, kx,)
    assert cape == pytest.approx(73.06134, rel=tol)
    assert cin == pytest.approx(-5.2824445, rel=tol)
    assert rain[0,0] == pytest.approx(1.4249854, rel=tol)
    assert bmflag[0,0] == 2

    # Can we call the Betts-Miller routine with a non-unform profile of relative humidity?
    rhbm_grid[:] = np.linspace(0.6, 0.9, num_lev)
    rain, tdel, qdel, q_ref, bmflag, klzbs, cape, cin, t_ref, invtau_bm_t, invtau_bm_q, capeflag = \
    betts_miller(dt, tin_grid, qin_grid_moist, pfull_grid, phalf_grid, 
                                          HLv, Cp_air, Grav,
                                          rdgas,rvgas,kappa, es0, tau_bm, rhbm_grid, 
                                          do_simp, do_shallower, do_changeqref, 
                                          do_envsat, do_taucape, capetaubm, tau_min, 
                                          ix, jx, kx,)

    # Make sure we can call the Betts-Miller routine with different option flags
    do_simp = True  # alternative way to conserve energy by adjusting the relaxation timescale on temperature rather than calculating an adjusted reference temperature. 
    rain, tdel, qdel, q_ref, bmflag, klzbs, cape, cin, t_ref, invtau_bm_t, invtau_bm_q, capeflag = \
    betts_miller(dt, tin_grid, qin_grid_moist, pfull_grid, phalf_grid, 
                                          HLv, Cp_air, Grav,
                                          rdgas,rvgas,kappa, es0, tau_bm, rhbm_grid, 
                                          do_simp, do_shallower, do_changeqref, 
                                          do_envsat, do_taucape, capetaubm, tau_min, 
                                          ix, jx, kx,)
    do_shallower = False  # Don't use the scheme that lowers the depth of shallow convection
    #  In this case, since we set do_changeqref = True, we use the scheme that involves
    #  changing the reference profiles for both temperature and humidity so the precipitation is zero
    rain, tdel, qdel, q_ref, bmflag, klzbs, cape, cin, t_ref, invtau_bm_t, invtau_bm_q, capeflag = \
    betts_miller(dt, tin_grid, qin_grid_moist, pfull_grid, phalf_grid, 
                                          HLv, Cp_air, Grav,
                                          rdgas,rvgas,kappa, es0, tau_bm, rhbm_grid, 
                                          do_simp, do_shallower, do_changeqref, 
                                          do_envsat, do_taucape, capetaubm, tau_min, 
                                          ix, jx, kx,)
    do_changeqref = False
    rain, tdel, qdel, q_ref, bmflag, klzbs, cape, cin, t_ref, invtau_bm_t, invtau_bm_q, capeflag = \
    betts_miller(dt, tin_grid, qin_grid_moist, pfull_grid, phalf_grid, 
                                          HLv, Cp_air, Grav,
                                          rdgas,rvgas,kappa, es0, tau_bm, rhbm_grid, 
                                          do_simp, do_shallower, do_changeqref, 
                                          do_envsat, do_taucape, capetaubm, tau_min, 
                                          ix, jx, kx,)
    do_envsat = False  
    # The reference profile of temperature is always the pseudoadiabat computed based on 
    # temperature and humidity at the surface. 
    # The parameter do_envsat determines how the reference profile of humidity is computed.
    # Setting do_envsat=False gives the situation described in Section 2d of the Frierson paper, 
    # namely that the same pseudoadiabat is used to define the target humidity profile, 
    # just multiplying by a specified relative humidity. 
    # In that sense, the target of the moisture adjustment due to convection is dictated by surface conditions.
    rain, tdel, qdel, q_ref, bmflag, klzbs, cape, cin, t_ref, invtau_bm_t, invtau_bm_q, capeflag = \
    betts_miller(dt, tin_grid, qin_grid_moist, pfull_grid, phalf_grid, 
                                          HLv, Cp_air, Grav,
                                          rdgas,rvgas,kappa, es0, tau_bm, rhbm_grid, 
                                          do_simp, do_shallower, do_changeqref, 
                                          do_envsat, do_taucape, capetaubm, tau_min, 
                                          ix, jx, kx,)

    do_taucape = True
    #  this makes the relaxation time proportional to 1/sqrt(CAPE), and so makes the scheme more aggressive when / where the CAPE is larger.
    rain, tdel, qdel, q_ref, bmflag, klzbs, cape, cin, t_ref, invtau_bm_t, invtau_bm_q, capeflag = \
    betts_miller(dt, tin_grid, qin_grid_moist, pfull_grid, phalf_grid, 
                                          HLv, Cp_air, Grav,
                                          rdgas,rvgas,kappa, es0, tau_bm, rhbm_grid, 
                                          do_simp, do_shallower, do_changeqref, 
                                          do_envsat, do_taucape, capetaubm, tau_min, 
                                          ix, jx, kx,)
