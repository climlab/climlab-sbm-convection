'''Translate Dargan Frierson's code to Python... work in progress'''

from climlab.utils.thermo import clausius_clapeyron as escomp
from climlab.utils.thermo import potential_temperature, temperature_from_potential
import numpy as np

def capecalc(p,phalf,cp_air,rdgas,rvgas,
             hlv,kappa,tin,rin,avgbl,cape,cin,tp,rp,klzb):
'''
!
!    Input:
!
!    kx          number of levels
!    p           pressure (index 1 refers to TOA, index kx refers to surface)
!    phalf       pressure at half levels
!    cp_air      specific heat of dry air
!    rdgas       gas constant for dry air
!    rvgas       gas constant for water vapor (used in Clausius-Clapeyron,
!                not for virtual temperature effects, which are not considered)
!    hlv         latent heat of vaporization
!    kappa       the constant kappa
!    tin         temperature of the environment
!    rin         specific humidity of the environment
!    avgbl       if true, the parcel is averaged in theta and r up to its LCL
!
!    Output:
!    cape        Convective available potential energy
!    cin         Convective inhibition (if there's no LFC, then this is set
!                to zero)
!    tp          Parcel temperature (set to the environmental temperature
!                where no adjustment)
!    rp          Parcel specific humidity (set to the environmental humidity
!                where no adjustment, and set to the saturation humidity at
!                the parcel temperature below the LCL)
!    klzb        Level of zero buoyancy
!
!    Algorithm:
!    Start with surface parcel.
!    Calculate the lifting condensation level (uses an analytic formula and a
!       lookup table).
!    Average under the LCL if desired, if this is done, then a new LCL must
!       be calculated.
!    Calculate parcel ascent up to LZB.
!    Calculate CAPE and CIN.
'''

    # integer, intent(in)                    :: kx
    # logical, intent(in)                    :: avgbl
    # real, intent(in), dimension(:)         :: p, phalf, tin, rin
    # real, intent(in)                       :: rdgas, rvgas, hlv, kappa, cp_air
    # integer, intent(out)                   :: klzb
    # real, intent(out), dimension(:)        :: tp, rp
    # real, intent(out)                      :: cape, cin
    #
    # integer            :: k, klcl, klfc, klcl2
    # logical            :: nocape
    # real, dimension(kx)   :: theta
    # real                  :: t0, r0, es, rs, theta0, pstar, value, tlcl, &
    #                          a, b, dtdlnp, d2tdlnp2, thetam, rm, tlcl2, &
    #                          plcl2, plcl, plzb, small

    kx = p.size
    pstar = 1.e5
    # ! so we can run dry limit (one expression involves 1/hlv)
    small = 1.e-10

    nocape = True
    cape = 0.
    cin = 0.
    plcl = 0.
    plzb = 0.
    klfc = 0
    klcl = 0
    klzb = 0
    tp = np.zeros_like(tin); tp[:] = tin #tp(1:kx) = tin(1:kx)
    rp = np.zeros_like(rin); rp[:] = rin #rp(1:kx) = rin(1:kx)

    # ! start with surface parcel
    t0 = tin[-1] #t0 = tin(kx)
    r0 = rin[-1] #r0 = rin(kx)
    # ! calculate the lifting condensation level by the following:
    # ! are you saturated to begin with?
    es = escomp(t0) #call escomp(t0,es)
    rs = rdgas/rvgas*es/p[-1] #rs = rdgas/rvgas*es/p(kx)
    if (r0 >= rs): #if (r0.ge.rs) then
    # ! if you�re already saturated, set lcl to be the surface value.
        plcl = p[-1] #plcl = p(kx)
    # ! the first level where you�re completely saturated.
        klcl = kx
    # ! saturate out to get the parcel temp and humidity at this level
    # ! first order (in delta T) accurate expression for change in temp
       #tp(kx) = t0 + (r0 - rs)/(cp_air/(hlv+small) + hlv*rs/rvgas/t0**2.)
       #call escomp(tp(kx),es)
       #rp(kx) = rdgas/rvgas*es/p(kx)
        tp[-1] = t0 + (r0 - rs)/(cp_air/(hlv+small) + hlv*rs/rvgas/t0**2.)
        rp[-1] = rdgas/rvgas*escomp(tp[-1])/p[-1]
    else:
   #  # ! if not saturated to begin with, use the analytic expression to calculate the
   #  # ! exact pressure and temperature where you�re saturated.
   #     theta0 = tin(kx)*(pstar/p(kx))**kappa
        theta0 = potential_temperature(tin[-1],p[-1])
   #  # ! the expression that we utilize is
   #  # ! log(r/theta**(1/kappa)*pstar*rvgas/rdgas/es00) = log(es/T**(1/kappa))
   #  # ! (the division by es00 is necessary because the RHS values are tabulated
   #  # ! for control moisture content)
   #  # ! The right hand side of this is only a function of temperature, therefore
   #  # ! this is put into a lookup table to solve for temperature.
   #     if (r0.gt.0.) then
        if r0>0:
            value = np.log(theta0**(-1/kappa)*r0*pstar*rvgas/rdgas/es0)
            tlcl = lcltabl(value)
            plcl = pstar*(tlcl/theta0)**(1/kappa)
   #        value = log(theta0**(-1/kappa)*r0*pstar*rvgas/rdgas/es0)
   #        call lcltabl(value,tlcl)
   #        plcl = pstar*(tlcl/theta0)**(1/kappa)
   #  # ! just in case plcl is very high up
            if plcl < p[0]:
                plcl = p[0]
                tlcl = temperature_from_potential(theta0, plcl)
            k = kx
   #        if (plcl.lt.p(1)) then
   #           plcl = p(1)
   #           tlcl = theta0*(plcl/pstar)**kappa
   #           write (*,*) 'hi lcl'
   #        end if
   #        k = kx
        else:
   #  # ! if the parcel sp hum is zero or negative, set lcl to top level
            plcl = p[0]
            tlcl = temperature_from_potential(theta0,plcl)
            tp[:] = temperature_from_potential(theta0,p)
            rp[:] = 0.
            cin += np.sum(rdgas*(tin - tp)*np.log(phalf[1:]/phalf[:-1]))
   #        plcl = p(1)
   #        tlcl = theta0*(plcl/pstar)**kappa
   #  # !            write (*,*) 'zero r0', r0
   #        do k=1,kx
   #           tp(k) = theta0*(p(k)/pstar)**kappa
   #           rp(k) = 0.
   #  # ! this definition of CIN contains everything below the LCL
   #           cin = cin + rdgas*(tin(k) - tp(k))*log(phalf(k+1)/phalf(k))
   #        end do
   #        go to 11
   #     end if
   #  # ! calculate the parcel temperature (adiabatic ascent) below the LCL.
   #  # ! the mixing ratio stays the same
   #     do while (p(k).gt.plcl)
   #        tp(k) = theta0*(p(k)/pstar)**kappa
   #        call escomp(tp(k),es)
   #        rp(k) = rdgas/rvgas*es/p(k)
   #  # ! this definition of CIN contains everything below the LCL
   #        cin = cin + rdgas*(tin(k) - tp(k))*log(phalf(k+1)/phalf(k))
   #        k = k-1
   #     end do
   #  # ! first level where you're saturated at the level
   #     klcl = k
   #  if (klcl.eq.1) klcl = 2
   #  # ! do a saturated ascent to get the parcel temp at the LCL.
   #  # ! use your 2nd order equation up to the pressure above.
   #  # ! moist adaibat derivatives: (use the lcl values for temp, humid, and
   #  # ! pressure)
   #     a = kappa*tlcl + hlv/cp_air*r0
   #     b = hlv**2.*r0/cp_air/rvgas/tlcl**2.
   #     dtdlnp = a/(1. + b)
   #  # ! first order in p
   #  # !         tp(klcl) = tlcl + dtdlnp*log(p(klcl)/plcl)
   #  # ! second order in p (RK2)
   #  # ! first get temp halfway up
   #     tp(klcl) = tlcl + dtdlnp*log(p(klcl)/plcl)/2.
   #     if ((tp(klcl).lt.173.16).and.nocape) go to 11
   #     call escomp(tp(klcl),es)
   #     rp(klcl) = rdgas/rvgas*es/(p(klcl) + plcl)*2.
   #     a = kappa*tp(klcl) + hlv/cp_air*rp(klcl)
   #     b = hlv**2./cp_air/rvgas*rp(klcl)/tp(klcl)**2.
   #     dtdlnp = a/(1. + b)
   #  # ! second half of RK2
   #     tp(klcl) = tlcl + dtdlnp*log(p(klcl)/plcl)
   #  # !         d2tdlnp2 = (kappa + b - 1. - b/tlcl*(hlv/rvgas/tlcl - &
   #  # !                   2.)*dtdlnp)/ (1. + b)*dtdlnp - hlv*r0/cp_air/ &
   #  # !                   (1. + b)
   #  # ! second order in p
   #  # !         tp(klcl) = tlcl + dtdlnp*log(p(klcl)/plcl) + .5*d2tdlnp2*(log(&
   #  # !             p(klcl)/plcl))**2.
   #     if ((tp(klcl).lt.173.16).and.nocape) go to 11
   #     call escomp(tp(klcl),es)
   #     rp(klcl) = rdgas/rvgas*es/p(klcl)
   #  # !         write (*,*) 'tp, rp klcl:kx, new', tp(klcl:kx), rp(klcl:kx)
   #  # ! CAPE/CIN stuff
   #     if ((tp(klcl).lt.tin(klcl)).and.nocape) then
   #  # ! if you�re not yet buoyant, then add to the CIN and continue
   #        cin = cin + rdgas*(tin(klcl) - &
   #             tp(klcl))*log(phalf(klcl+1)/phalf(klcl))
   #     else
   #  # ! if you�re buoyant, then add to cape
   #        cape = cape + rdgas*(tp(klcl) - &
   #              tin(klcl))*log(phalf(klcl+1)/phalf(klcl))
   #  # ! if it�s the first time buoyant, then set the level of free convection to k
   #        if (nocape) then
   #           nocape = .false.
   #           klfc = klcl
   #        endif
   #     end if
   #  end if
   #  # ! then average the properties over the boundary layer if so desired.  to give
   #  # ! a new "parcel".  this may not be saturated at the LCL, so make sure you get
   #  # ! to a level where it is before moist adiabatic ascent!
   #  # !!!! take out all the below (between the exclamation points) if no avgbl !!!!
   #  if (avgbl) then
   #     theta(klcl:kx) = tin(klcl:kx)*(pstar/p(klcl:kx))**kappa
   #     thetam = 0.
   #     rm = 0.
   #     do k=klcl,kx
   #        thetam = thetam + theta(k)*(phalf(k+1) - phalf(k))
   #        rm = rm + rin(k)*(phalf(k+1) - phalf(k))
   #     end do
   #     thetam = thetam/(phalf(kx+1) - phalf(klcl))
   #     rm = rm/(phalf(kx+1) - phalf(klcl))
   #  # ! check if you're saturated at the top level.  if not, then get a new LCL
   #     tp(klcl) = thetam*(p(klcl)/pstar)**kappa
   #     call escomp(tp(klcl),es)
   #     rs = rdgas/rvgas*es/p(klcl)
   #  # ! if you're not saturated, get a new LCL
   #     if (rm.lt.rs) then
   #  # ! reset CIN to zero.
   #        cin = 0.
   #  # ! again, use the analytic expression to calculate the exact pressure and
   #  # ! temperature where you�re saturated.
   #  # ! the expression that we utilize is
   #  # ! log(r/theta**(1/kappa)*pstar*rvgas/rdgas/es00)= log(es/T**(1/kappa))
   #  # ! (the division by es00 is necessary because the RHS values are tabulated
   #  # ! for control moisture content)
   #  # ! The right hand side of this is only a function of temperature, therefore
   #  # ! this is put into a lookup table to solve for temperature.
   #        value = log(thetam**(-1/kappa)*rm*pstar*rvgas/rdgas/es0)
   #        call lcltabl(value,tlcl2)
   #        plcl2 = pstar*(tlcl2/thetam)**(1/kappa)
   #  # ! just in case plcl is very high up
   #        if (plcl2.lt.p(1)) then
   #           plcl2 = p(1)
   #        end if
   #        k = kx
   #  # ! calculate the parcel temperature (adiabatic ascent) below the LCL.
   #  # ! the mixing ratio stays the same
   #        do while (p(k).gt.plcl2)
   #           tp(k) = thetam*(p(k)/pstar)**kappa
   #           call escomp(tp(k),es)
   #           rp(k) = rdgas/rvgas*es/p(k)
   #  # ! this definition of CIN contains everything below the LCL
   #           cin = cin + rdgas*(tin(k) - tp(k))*log(phalf(k+1)/phalf(k))
   #           k = k-1
   #        end do
   #  # ! first level where you�re saturated at the level
   #        klcl2 = k
   #  if (klcl2.eq.1) klcl2 = 2
   #  # ! do a saturated ascent to get the parcel temp at the LCL.
   #  # ! use your 2nd order equation up to the pressure above.
   #  # ! moist adaibat derivatives: (use the lcl values for temp, humid, and
   #  # ! pressure)
   #        a = kappa*tlcl2 + hlv/cp_air*rm
   #        b = hlv**2.*rm/cp_air/rvgas/tlcl2**2.
   #        dtdlnp = a/(1. + b)
   #  # ! first order in p
   #  # !            tp(klcl2) = tlcl2 + dtdlnp*log(p(klcl2)/plcl2)
   #  # ! second order in p (RK2)
   #  # ! first get temp halfway up
   #     tp(klcl2) = tlcl2 + dtdlnp*log(p(klcl2)/plcl2)/2.
   #     if ((tp(klcl2).lt.173.16).and.nocape) go to 11
   #     call escomp(tp(klcl2),es)
   #     rp(klcl2) = rdgas/rvgas*es/(p(klcl2) + plcl2)*2.
   #     a = kappa*tp(klcl2) + hlv/cp_air*rp(klcl2)
   #     b = hlv**2./cp_air/rvgas*rp(klcl2)/tp(klcl2)**2.
   #     dtdlnp = a/(1. + b)
   #  # ! second half of RK2
   #     tp(klcl2) = tlcl2 + dtdlnp*log(p(klcl2)/plcl2)
   #  # !            d2tdlnp2 = (kappa + b - 1. - b/tlcl2*(hlv/rvgas/tlcl2 - &
   #  # !                          2.)*dtdlnp)/ (1. + b)*dtdlnp - hlv*rm/cp_air/ &
   #  # !                          (1. + b)
   #  # ! second order in p
   #  # !            tp(klcl2) = tlcl2 + dtdlnp*log(p(klcl2)/plcl2) + &
   #  # !               .5*d2tdlnp2*(log(p(klcl2)/plcl2))**2.
   #        call escomp(tp(klcl2),es)
   #        rp(klcl2) = rdgas/rvgas*es/p(klcl2)
   #  # ! CAPE/CIN stuff
   #        if ((tp(klcl2).lt.tin(klcl2)).and.nocape) then
   #  # ! if you�re not yet buoyant, then add to the CIN and continue
   #           cin = cin + rdgas*(tin(klcl2) - &
   #                tp(klcl2))*log(phalf(klcl2+1)/phalf(klcl2))
   #        else
   #  # ! if you�re buoyant, then add to cape
   #           cape = cape + rdgas*(tp(klcl) - &
   #                 tin(klcl))*log(phalf(klcl+1)/phalf(klcl))
   #  # ! if it�s the first time buoyant, then set the level of free convection to k
   #           if (nocape) then
   #              nocape = .false.
   #              klfc = klcl2
   #           endif
   #        end if
   #     end if
   #  end if
   #  # !!!! take out all of the above (within the exclamations) if no avgbl !!!!
   #  # ! then, start at the LCL, and do moist adiabatic ascent by the first order
   #  # ! scheme -- 2nd order as well
   #  do k=klcl-1,1,-1
   #     a = kappa*tp(k+1) + hlv/cp_air*rp(k+1)
   #     b = hlv**2./cp_air/rvgas*rp(k+1)/tp(k+1)**2.
   #     dtdlnp = a/(1. + b)
   #  # ! first order in p
   #  # !         tp(k) = tp(k+1) + dtdlnp*log(p(k)/p(k+1))
   #  # ! second order in p (RK2)
   #  # ! first get temp halfway up
   #     tp(k) = tp(k+1) + dtdlnp*log(p(k)/p(k+1))/2.
   #     if ((tp(k).lt.173.16).and.nocape) go to 11
   #     call escomp(tp(k),es)
   #     rp(k) = rdgas/rvgas*es/(p(k) + p(k+1))*2.
   #     a = kappa*tp(k) + hlv/cp_air*rp(k)
   #     b = hlv**2./cp_air/rvgas*rp(k)/tp(k)**2.
   #     dtdlnp = a/(1. + b)
   #  # ! second half of RK2
   #     tp(k) = tp(k+1) + dtdlnp*log(p(k)/p(k+1))
   #  # !         d2tdlnp2 = (kappa + b - 1. - b/tp(k+1)*(hlv/rvgas/tp(k+1) - &
   #  # !               2.)*dtdlnp)/(1. + b)*dtdlnp - hlv/cp_air*rp(k+1)/(1. + b)
   #  # ! second order in p
   #  # !         tp(k) = tp(k+1) + dtdlnp*log(p(k)/p(k+1)) + .5*d2tdlnp2*(log( &
   #  # !             p(k)/p(k+1)))**2.
   #  # ! if you're below the lookup table value, just presume that there's no way
   #  # ! you could have cape and call it quits
   #     if ((tp(k).lt.173.16).and.nocape) go to 11
   #     call escomp(tp(k),es)
   #     rp(k) = rdgas/rvgas*es/p(k)
   #     if ((tp(k).lt.tin(k)).and.nocape) then
   #  # ! if you�re not yet buoyant, then add to the CIN and continue
   #        cin = cin + rdgas*(tin(k) - tp(k))*log(phalf(k+1)/phalf(k))
   #     elseif((tp(k).lt.tin(k)).and.(.not.nocape)) then
   #  # ! if you have CAPE, and it�s your first time being negatively buoyant,
   #  # ! then set the level of zero buoyancy to k+1, and stop the moist ascent
   #        klzb = k+1
   #        go to 11
   #     else
   #  # ! if you�re buoyant, then add to cape
   #        cape = cape + rdgas*(tp(k) - tin(k))*log(phalf(k+1)/phalf(k))
   #  # ! if it�s the first time buoyant, then set the level of free convection to k
   #        if (nocape) then
   #           nocape = .false.
   #           klfc = k
   #        endif
   #     end if
   #  end do
   #  11   if(nocape) then
   #  # ! this is if you made it through without having a LZB
   #  # ! set LZB to be the top level.
   #     plzb = p(1)
   #     klzb = 0
   #     klfc = 0
   #     cin = 0.
   #     tp(1:kx) = tin(1:kx)
   #     rp(1:kx) = rin(1:kx)
   #  #end if
   #  # !      write (*,*) 'plcl, klcl, tlcl, r0 new', plcl, klcl, tlcl, r0
   #  # !      write (*,*) 'tp, rp new', tp, rp
   #  # !       write (*,*) 'tp, new', tp
   #  # !       write (*,*) 'tin new', tin
   #  # !       write (*,*) 'klcl, klfc, klzb new', klcl, klfc, klzb
    return cape, cin, tp, rp, klzb


def lcltabl(value):
'''
! lookup table for the analytic evaluation of LCL
      subroutine lcltabl(value,tlcl)
!
! Table of values used to compute the temperature of the lifting condensation
! level.
!
! the expression that we utilize is
! log(r/theta**(1/kappa)*pstar*rvgas/rdgas/es00) = log(es/T**(1/kappa))
! the RHS is tabulated for control moisture content, hence the division
! by es00 on the LHS
!
! Gives the values of the temperature for the following range:
!   starts with -23, is uniformly distributed up to -10.4.  There are a
! total of 127 values, and the increment is .1.
!
'''
      # implicit none
      # real, intent(in)     :: value
      # real, intent(out)    :: tlcl
      #
      # integer              :: ival
      # real, dimension(127) :: lcltable
      # real                 :: v1, v2

      lcltable = np.array([  1.7364512e+02,   1.7427449e+02,   1.7490874e+02,
      1.7554791e+02,   1.7619208e+02,   1.7684130e+02,   1.7749563e+02,
      1.7815514e+02,   1.7881989e+02,   1.7948995e+02,   1.8016539e+02,
      1.8084626e+02,   1.8153265e+02,   1.8222461e+02,   1.8292223e+02,
      1.8362557e+02,   1.8433471e+02,   1.8504972e+02,   1.8577068e+02,
      1.8649767e+02,   1.8723077e+02,   1.8797006e+02,   1.8871561e+02,
      1.8946752e+02,   1.9022587e+02,   1.9099074e+02,   1.9176222e+02,
      1.9254042e+02,   1.9332540e+02,   1.9411728e+02,   1.9491614e+02,
      1.9572209e+02,   1.9653521e+02,   1.9735562e+02,   1.9818341e+02,
      1.9901870e+02,   1.9986158e+02,   2.0071216e+02,   2.0157057e+02,
      2.0243690e+02,   2.0331128e+02,   2.0419383e+02,   2.0508466e+02,
      2.0598391e+02,   2.0689168e+02,   2.0780812e+02,   2.0873335e+02,
      2.0966751e+02,   2.1061074e+02,   2.1156316e+02,   2.1252493e+02,
      2.1349619e+02,   2.1447709e+02,   2.1546778e+02,   2.1646842e+02,
      2.1747916e+02,   2.1850016e+02,   2.1953160e+02,   2.2057364e+02,
      2.2162645e+02,   2.2269022e+02,   2.2376511e+02,   2.2485133e+02,
      2.2594905e+02,   2.2705847e+02,   2.2817979e+02,   2.2931322e+02,
      2.3045895e+02,   2.3161721e+02,   2.3278821e+02,   2.3397218e+02,
      2.3516935e+02,   2.3637994e+02,   2.3760420e+02,   2.3884238e+02,
      2.4009473e+02,   2.4136150e+02,   2.4264297e+02,   2.4393941e+02,
      2.4525110e+02,   2.4657831e+02,   2.4792136e+02,   2.4928053e+02,
      2.5065615e+02,   2.5204853e+02,   2.5345799e+02,   2.5488487e+02,
      2.5632953e+02,   2.5779231e+02,   2.5927358e+02,   2.6077372e+02,
      2.6229310e+02,   2.6383214e+02,   2.6539124e+02,   2.6697081e+02,
      2.6857130e+02,   2.7019315e+02,   2.7183682e+02,   2.7350278e+02,
      2.7519152e+02,   2.7690354e+02,   2.7863937e+02,   2.8039954e+02,
      2.8218459e+02,   2.8399511e+02,   2.8583167e+02,   2.8769489e+02,
      2.8958539e+02,   2.9150383e+02,   2.9345086e+02,   2.9542719e+02,
      2.9743353e+02,   2.9947061e+02,   3.0153922e+02,   3.0364014e+02,
      3.0577420e+02,   3.0794224e+02,   3.1014515e+02,   3.1238386e+02,
      3.1465930e+02,   3.1697246e+02,   3.1932437e+02,   3.2171609e+02,
      3.2414873e+02,   3.2662343e+02,   3.2914139e+02,   3.3170385e+02])

      # v1 = value
      # if (value.lt.-23.0) v1 = -23.0
      # if (value.gt.-10.4) v1 = -10.4
      # ival = floor(10.*(v1 + 23.0))
      # v2 = -230. + ival
      # v1 = 10.*v1
      # tlcl = (v2 + 1.0 - v1)*lcltable(ival+1) + (v1 - v2)*lcltable(ival+2)
      v1 = np.maximum(value, -23.0)
      v1 = np.minimum(v1, -10.4)
      ival = int(10.*(v1+23.0))
      v2 = -230. + ival
      v1 = 10.*v1
      tlcl = (v2 + 1.0 - v1)*lcltable[ival] + (v1 - v2)*lcltable[ival+1]
      return tlcl
