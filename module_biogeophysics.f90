module module_biogeophysics

  use shr_const_mod, ONLY: SHR_CONST_TKFRZ
  use shr_const_mod, only : r8 => shr_kind_r8, r4 => SHR_KIND_R4

contains

subroutine Biogeophysics1 (clm,cgrnd,cgrndl,cgrnds,tg,emg,htvp, dlrad,ulrad,tssbef)

!-----------------------------------------------------------------------
!
!  CLMCLMCLMCLMCLMCLMCLMCLMCL  A community developed and sponsored, freely   
!  L                        M  available land surface process model.  
!  M --COMMON LAND MODEL--  C  
!  C                        L  CLM WEB INFO: http://clm.gsfc.nasa.gov
!  LMCLMCLMCLMCLMCLMCLMCLMCLM  CLM ListServ/Mailing List: 
!
!-----------------------------------------------------------------------
! Purpose:
! This is the main subroutine to execute the calculation of leaf temperature
! and surface fluxes. Biogeophysics2.F90 then determines soil/snow and ground
! temperatures and updates the surface fluxes for the new ground
! temperature.
!
! Calling sequence is:
!  Biogeophysics1:                   surface biogeophysics driver
!    -> QSatclm:                        saturated vapor pressure, specific humidity,
!                                     and derivatives at ground surface
!    -> SurfaceRadiation:            surface solar radiation
!    -> BareGroundFluxes:            surface fluxes for bare soil or
!                                     snow-covered vegetation patches
!          -> MoninObukIni:          first-guess Monin-Obukhov length and
!                                     wind speed
!          -> FrictionVelocity:      friction velocity and potential 
!                                     temperature and humidity profiles
!    -> CanopyFluxes:                leaf temperature and surface fluxes
!                                     for vegetated patches 
!          -> QSatclm                   saturated vapor pressure, specific humidity,
!                                     and derivatives at leaf surface
!          -> MoninObukIni           first-guess Monin-Obukhov length and 
!                                     wind speed
!          -> FrictionVelocity       friction velocity and potential
!                                     temperature and humidity profiles
!          -> Stomata                stomatal resistance and photosynthesis
!                                     for sunlit leaves
!          -> Stomata                stomatal resistance and photosynthesis
!                                     for shaded leaves
!          -> SensibleHCond          sensible heat conductance for air, leaf,
!                                     and ground
!          -> LatentHCond            latent heat conductance for ground and
!                                     leaf
!          -> QSatclm                   recalculation of saturated vapor pressure,
!                                     specific humidity, and derivatives at
!                                     leaf surface using updated leaf temperature
!  Leaf temperature
!   Foliage energy conservation is given by the foliage energy budget 
!   equation:
!                  Rnet - Hf - LEf = 0 
!   The equation is solved by Newton-Raphson iteration, in which this 
!   iteration includes the calculation of the photosynthesis and 
!   stomatal resistance, and the integration of turbulent flux profiles. 
!   The sensible and latent heat transfer between foliage and atmosphere 
!   and ground is linked by the equations:  
!                  Ha = Hf + Hg and Ea = Ef + Eg
!
! Author:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision 
!
!-----------------------------------------------------------------------
! $Id: Biogeophysics1.F90,v 1.8 2004/11/24 22:56:17 jim Exp $
!-----------------------------------------------------------------------

  use clmtype
  use clm_varcon, only : denh2o, denice, roverg, hvap, hsub, istice, istwet 
  use clm_varpar, only : nlevsoi 

  implicit none

!----Arguments----------------------------------------------------------

  type (clm1d), intent(inout)  :: clm !CLM 1-D Module

!----Local Variables----------------------------------------------------

  integer i      ! loop indices
  real(r8) qred    ! soil surface relative humidity
  real(r8) avmuir  ! ir inverse optical depth per unit leaf area
  real(r8) eg      ! water vapor pressure at temperature T [pa]
  real(r8) qsatg   ! saturated humidity [kg/kg]
  real(r8) degdT   ! d(eg)/dT
  real(r8) qsatgdT ! d(qsatg)/dT
  real(r8) fac     ! soil wetness of surface layer
  real(r8) psit    ! negative potential of soil
  real(r8) hr      ! relative humidity
  real(r8) wx      ! partial volume of ice and water of surface layer
  real(r8) zlnd
  real(r8) zsno
  real(r8) smpmin
  real(r8) qg
  real(r8) dqgdT
  real(r8) emv
  real(r8) z0mg
  real(r8) z0hg
  real(r8) z0qg
  real(r8) z0mv
  real(r8) z0hv
  real(r8) z0qv
  real(r8) beta
  real(r8) zii
  real(r8) thm
  real(r8) thv
  real(r8) ur

  real(r8) :: cgrnd
  real(r8) :: cgrndl
  real(r8) :: cgrnds
  real(r8) :: tg
  real(r8) :: emg
  real(r8) :: htvp
  real(r8) :: dlrad
  real(r8) :: ulrad
  real(r8) :: tssbef(-5:10)
!----End Variable List--------------------------------------------------

!
! Initial set 
!
  zlnd = 0.01
  zsno = 0.0024
  smpmin = -1.e8
  clm%eflx_sh_tot    = 0.
  clm%qflx_evap_tot  = 0.
  clm%eflx_lh_tot    = 0.
  clm%eflx_sh_veg    = 0.  
  clm%qflx_evap_veg  = 0.  
  clm%qflx_tran_veg  = 0.  
  cgrnd          = 0._r8
  cgrnds         = 0._r8
  cgrndl         = 0._r8
  clm%t_ref2m        = 0.
!
! Ground and soil temperatures from previous time step
!


  tg = clm%t_soisno(clm%snl+1)
  do i = clm%snl+1, nlevsoi
     tssbef(i) = clm%t_soisno(i)
  enddo

!
! Saturated vapor pressure, specific humidity and their derivatives
! at ground surface
!
  qred = 1.
  if (clm%itypwat/=istwet .AND. clm%itypwat/=istice) then
     wx   = (clm%h2osoi_liq(1)/denh2o+clm%h2osoi_ice(1)/denice)/clm%dz(1)
     fac  = min(1._r4, wx/clm%watsat(1))
     fac  = max( fac, 0.01_r4 )
     psit = -clm%sucsat(1) * fac ** (- clm%bsw(1))
     psit = max(smpmin, psit)
     hr   = exp(psit/roverg/tg)
     qred = (1.-clm%frac_sno)*hr + clm%frac_sno
  endif

  call QSatclm(tg, clm%forc_pbot, eg, degdT, qsatg, &
            qsatgdT)
  qg = qred*qsatg  
  dqgdT = qred*qsatgdT

  if (qsatg > clm%forc_q .AND. clm%forc_q > qred*qsatg) then
     qg = clm%forc_q
     dqgdT = 0.
  endif

!
! Emissivity
!
  if (clm%h2osno>0. .OR.clm%itypwat==istice) then
     emg = 0.97
  else
     emg = 0.96
  endif
  avmuir=1.
  emv=1.-exp(-(clm%elai+clm%esai)/avmuir)

!
! Latent heat. We arbitrarily assume that the sublimation occurs 
! only as h2osoi_liq = 0
!

  htvp = hvap
  if (clm%h2osoi_liq(clm%snl+1) <= 0. .AND. clm%h2osoi_ice(clm%snl+1) > 0.) htvp = hsub

!
! Switch between vaporization and sublimation causes rapid solution
! separation in perturbation growth test
!

#if (defined PERGRO)
  htvp = hvap
#endif

!
! Roughness lengths
!

  if (clm%frac_sno > 0.) then
     z0mg = zsno
     z0hg = z0mg            ! initial set only
     z0qg = z0mg            ! initial set only
  else
     z0mg = zlnd
     z0hg = z0mg            ! initial set only
     z0qg = z0mg            ! initial set only
  endif

  clm%z0m = clm%z0mr*clm%htop
  clm%displa = clm%displar*clm%htop
!print*, 'displa, displar, htop:', clm%displa, clm%displar, clm%htop 
  z0mv = clm%z0m
  z0hv = z0mv
  z0qv = z0mv

!=== LDAS modifications****
!=== New code added for LDAS framework
!=== Added specification of forcing heights here instead of in the ldasdrv.f
!=== since displacement heights set here as opposed to being read in from a
!=== text file

!=== This change was necessary because the forcing heights would be less than the
!=== canopy heights and cause model crash 
#if ( defined COUPLED )
!!clm%forc_hgt  =clm%forc_hgt+clm%displa+clm%z0m
  clm%forc_hgt_u=clm%forc_hgt+clm%displa+clm%z0m
  clm%forc_hgt_t=clm%forc_hgt+clm%displa+clm%z0m
  clm%forc_hgt_q=clm%forc_hgt+clm%displa+clm%z0m
!!print*, 'forc_hgt', clm%forc_hgt, clm%forc_hgt_u, clm%forc_hgt_t, clm%forc_hgt_q
#else
  clm%forc_hgt  =10.0+clm%displa+clm%z0m
  clm%forc_hgt_u=10.0+clm%displa+clm%z0m
  clm%forc_hgt_t=2.0+clm%displa+clm%z0m
  clm%forc_hgt_q=2.0+clm%displa+clm%z0m
#endif
!    print*, 'in bio..',clm%t_veg

!
! Potential, virtual potential temperature, and wind speed at the 
! reference height
!

  beta=1.
  zii = 1000.
  thm = clm%forc_t + 0.0098*clm%forc_hgt_t              
  thv = clm%forc_th*(1.+0.61*clm%forc_q)
  ur = max(1.0_r8,sqrt(clm%forc_u*clm%forc_u+clm%forc_v*clm%forc_v))
!  print*, ur

!
! Surface Radiation
!
  call SurfaceRadiation (clm)
!
! Surface Temperature and Fluxes
!

!
! BARE SOIL OR SNOW-COVERED VEGETATION
! Ground fluxes
! NOTE: in the current scheme clm%frac_veg_nosno is EITHER 1 or 0
!

! print*, clm%itypveg,clm%t_veg,clm%t_grnd,clm%displa,clm%z0m
!  write(*,33) clm%itypveg,clm%t_veg,clm%t_grnd,clm%displa,clm%z0m,thm,&
!& clm%forc_th,thm,clm%tg,clm%qg,clm%dqgdT,emv,clm%emg,clm%dlrad,clm%ulrad,&
!& clm%cgrnds,clm%cgrndl,clm%cgrnd
! 33 format(i3,1x,f8.4,1x,f8.4,1x,f8.4,1x,f8.4,1x,3(f8.4,1x),3(f8.4,1x),1x,7(f8.4,1x))
  if (clm%frac_veg_nosno == 0) then

!     call BareGroundFluxes (tg,     thm,   qg,    thv,   z0mg,   &
!                            z0hg,   z0qg,  dqgdT, htvp,  beta,   &
!                            zii,    ur,    dlrad, ulrad, cgrnds, &
!                            cgrndl, cgrnd, clm    )
     clm%psnsun = 0.
     clm%psnsha = 0. !put these lines here to avoid psn = NaN

!         print*, 'in bio..',clm%t_veg

!
! VEGETATION
! Calculate canopy temperature, latent and sensible fluxes from the canopy,
! and leaf water change by evapotranspiration
!

  else

     call CanopyFluxes (z0mv,   z0hv,  z0qv,  thm,   clm%forc_th, &
                        thv,    tg,    qg,    dqgdT, htvp,        &
                        emv,    emg,   dlrad, ulrad, cgrnds,      &
                        cgrndl, cgrnd, clm    )

!         print*, clm%t_veg, clm%t_grnd, clm%displa, clm%z0m, clm%frac_veg_nosno

  endif
  
!  print*, clm%itypveg, clm%t_veg, clm%t_grnd, clm%displa, clm%z0m, clm%frac_veg_nosno
  

end subroutine Biogeophysics1

!.........................................
subroutine QSatclm (T, p, es, esdT, qs, &
                 qsdT)

!-----------------------------------------------------------------------
!
!  CLMCLMCLMCLMCLMCLMCLMCLMCL  A community developed and sponsored, freely   
!  L                        M  available land surface process model.  
!  M --COMMON LAND MODEL--  C  
!  C                        L  CLM WEB INFO: http://clm.gsfc.nasa.gov
!  LMCLMCLMCLMCLMCLMCLMCLMCLM  CLM ListServ/Mailing List: 
!
!-----------------------------------------------------------------------
! Purpose:
! Computes saturation mixing ratio and the change in saturation
! mixing ratio with respect to temperature.

! Method:
! Reference:  Polynomial approximations from:
!             Piotr J. Flatau, et al.,1992:  Polynomial fits to saturation 
!             vapor pressure.  Journal of Applied Meteorology, 31, 1507-1513.
!
! Author:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision 
!
!-----------------------------------------------------------------------
! $Id: QSatclm.F90,v 1.6 2004/11/24 22:56:33 jim Exp $
!-----------------------------------------------------------------------

  implicit none

!----Arguments----------------------------------------------------------

  real(r8), intent(in)  :: T        ! temperature (K)
  real(r8), intent(in)  :: p        ! surface atmospheric pressure (pa)

  real(r8), intent(out) :: es       ! vapor pressure (pa)
  real(r8), intent(out) :: esdT     ! d(es)/d(T)
  real(r8), intent(out) :: qs       ! humidity (kg/kg)
  real(r8), intent(out) :: qsdT     ! d(qs)/d(T)

!----Local Variables----------------------------------------------------

  real(r8) T_limit

  real(r8) td,vp,vp1,vp2
  real(r8) a0,a1,a2,a3,a4,a5,a6,a7,a8
  real(r8) b0,b1,b2,b3,b4,b5,b6,b7,b8

  real(r8) c0,c1,c2,c3,c4,c5,c6,c7,c8
  real(r8) d0,d1,d2,d3,d4,d5,d6,d7,d8

!
! For water vapor (temperature range 0C-100C)
!

  data a0/6.11213476   /,a1/ 0.444007856    /,a2/0.143064234e-01/  &
      ,a3/0.264461437e-03/,a4/ 0.305903558e-05/,a5/0.196237241e-07/  &
      ,a6/0.892344772e-10/,a7/-0.373208410e-12/,a8/0.209339997e-15/

!
! For derivative:water vapor
!

  data b0/0.444017302  /,b1/ 0.286064092e-01/,b2/ 0.794683137e-03/ &
      ,b3/ 0.121211669e-04/,b4/ 0.103354611e-06/,b5/ 0.404125005e-09/ &
      ,b6/-0.788037859e-12/,b7/-0.114596802e-13/,b8/ 0.381294516e-16/

!
! For ice (temperature range -75C-0C) 
!

  data c0/6.11123516     /,c1/0.503109514    /,c2/0.188369801e-01/ &
      ,c3/0.420547422e-03/,c4/0.614396778e-05/,c5/0.602780717e-07/ &
      ,c6/0.387940929e-09/,c7/0.149436277e-11/,c8/0.262655803e-14/

!
! For derivative:ice  
!

  data d0/0.503277922    /,d1/0.377289173e-01/,d2/0.126801703e-02/ &
      ,d3/0.249468427e-04/,d4/0.313703411e-06/,d5/0.257180651e-08/ &
      ,d6/0.133268878e-10/,d7/0.394116744e-13/,d8/0.498070196e-16/

!----End Variable List--------------------------------------------------

  T_limit = T - SHR_CONST_TKFRZ
  if (T_limit > 100.0) T_limit=100.0
  if (T_limit < -75.0) T_limit=-75.0

  td       = T_limit
  if (td >= 0.0) then
     es   = a0 + td*(a1 + td*(a2 + td*(a3 + td*(a4 &
          + td*(a5 + td*(a6 + td*(a7 + td*a8)))))))
     esdT = b0 + td*(b1 + td*(b2 + td*(b3 + td*(b4 &
          + td*(b5 + td*(b6 + td*(b7 + td*b8)))))))
  else
     es   = c0 + td*(c1 + td*(c2 + td*(c3 + td*(c4 &
          + td*(c5 + td*(c6 + td*(c7 + td*c8)))))))
     esdT = d0 + td*(d1 + td*(d2 + td*(d3 + td*(d4 &
          + td*(d5 + td*(d6 + td*(d7 + td*d8)))))))
  endif

  es    = es    * 100.            ! pa
  esdT  = esdT  * 100.            ! pa/K

  vp    = 1.0   / (p - 0.378*es)
  vp1   = 0.622 * vp
  vp2   = vp1   * vp

  qs    = es    * vp1             ! kg/kg
  qsdT  = esdT  * vp2 * p         ! 1 / K
end subroutine QSatclm

!...............................................
subroutine SurfaceRadiation (clm)

!-----------------------------------------------------------------------
!
!  CLMCLMCLMCLMCLMCLMCLMCLMCL  A community developed and sponsored, freely
!  L                        M  available land surface process model.
!  M --COMMON LAND MODEL--  C
!  C                        L  CLM WEB INFO: http://clm.gsfc.nasa.gov
!  LMCLMCLMCLMCLMCLMCLMCLMCLM  CLM ListServ/Mailing List:
!
!-----------------------------------------------------------------------
! Purpose: 
! Solar fluxes absorbed by vegetation and ground surface
! 
! Method: 
! Note possible problem when land is on different grid than atmosphere.
!
! Land may have sun above the horizon (coszen > 0) but atmosphere may
! have sun below the horizon (forc_solad = 0 and forc_solai = 0). This is okay
! because all fluxes (absorbed, reflected, transmitted) are multiplied
! by the incoming flux and all will equal zero.
!
! Atmosphere may have sun above horizon (forc_solad > 0 and forc_solai > 0) but
! land may have sun below horizon. This is okay because fabd, fabi,
! ftdd, ftid, and ftii all equal zero so that sabv=sabg=fsa=0. Also,
! albd and albi equal one so that fsr=forc_solad+forc_solai. In other words, all
! the radiation is reflected. NDVI should equal zero in this case.
! However, the way the code is currently implemented this is only true
! if (forc_solad+forc_solai)|vis = (forc_solad+forc_solai)|nir.
!
! Output variables are parsun,parsha,sabv,sabg,fsa,fsr,ndvi
! 
! Author:
! Gordon Bonan
! 
!-----------------------------------------------------------------------
! $Id: SurfaceRadiation.F90,v 1.7 2004/11/24 22:56:46 jim Exp $
!-----------------------------------------------------------------------

  use clmtype
  use clm_varpar,  only : numrad
  implicit none

!----Arguments----------------------------------------------------------

  type (clm1d), intent(inout) :: clm    !CLM 1-D Module

!----Local Variables----------------------------------------------------

  integer ib              ! waveband number (1=vis, 2=nir)
  integer  nband          ! number of solar radiation waveband classes              
  real(r8) abs            ! absorbed solar radiation (W/m**2) 
!  real(r8) rnir           ! reflected solar radiation [nir] (W/m**2)
!  real(r8) rvis           ! reflected solar radiation [vis] (W/m**2)
  real(r8) laifra         ! leaf area fraction of canopy
  real(r8) trd            ! transmitted solar radiation: direct (W/m**2)
  real(r8) tri            ! transmitted solar radiation: diffuse (W/m**2)
  real(r8) cad(numrad)    ! direct beam absorbed by canopy (W/m**2)
  real(r8) cai(numrad)    ! diffuse radiation absorbed by canopy (W/m**2)
  real(r8) fsha           ! shaded fraction of canopy
  real(r8) vai            ! total leaf area index + stem area index, one sided
  real(r8) mpe            ! prevents overflow for division by zero                  

!----End Variable List--------------------------------------------------

  mpe   = 1.e-06
  nband = numrad

  fsha = 1.-clm%fsun
  clm%laisun = clm%elai*clm%fsun
  clm%laisha = clm%elai*fsha
  vai = clm%elai+ clm%esai

!
! Zero summed solar fluxes
!

  clm%sabg = 0.
  clm%sabv = 0.
  clm%fsa  = 0.

!
! Loop over nband wavebands
!

  do ib = 1, nband

!
! Absorbed by canopy
!
     cad(ib)  = clm%forc_solad(ib)*clm%fabd(ib)
     cai(ib)  = clm%forc_solai(ib)*clm%fabi(ib)
     clm%sabv = clm%sabv + cad(ib) + cai(ib)
     clm%fsa  = clm%fsa  + cad(ib) + cai(ib)

!
! Transmitted = solar fluxes incident on ground
!

     trd = clm%forc_solad(ib)*clm%ftdd(ib)
     tri = clm%forc_solad(ib)*clm%ftid(ib) + clm%forc_solai(ib)*clm%ftii(ib)

!
! Solar radiation absorbed by ground surface
!

     abs = trd*(1.-clm%albgrd(ib)) + tri*(1.-clm%albgri(ib)) 
     clm%sabg = clm%sabg + abs
     clm%fsa  = clm%fsa  + abs

  end do

!
! Partition visible canopy absorption to sunlit and shaded fractions
! to get average absorbed par for sunlit and shaded leaves
!

  laifra = clm%elai / max(vai,mpe)
  if (clm%fsun > 0.) then
     clm%parsun = (cad(1) + cai(1)) * laifra
!     clm%parsha = 0._r4 
  else
     clm%parsun = 0._r4 
!     clm%parsha = 0._r4 
  endif

!
! NDVI and reflected solar radiation
!

!  rvis = clm%albd(1)*clm%forc_solad(1) + clm%albi(1)*clm%forc_solai(1) 
!  rnir = clm%albd(2)*clm%forc_solad(2) + clm%albi(2)*clm%forc_solai(2)
!  clm%fsr = rvis + rnir
!  clm%ndvi = (rnir-rvis) / max(rnir+rvis,mpe)
  return
end subroutine SurfaceRadiation 

!=================================================
!=================================================

subroutine SurfaceAlbedo (clm)

!-----------------------------------------------------------------------
!
!  CLMCLMCLMCLMCLMCLMCLMCLMCL  A community developed and sponsored, freely
!  L                        M  available land surface process model.
!  M --COMMON LAND MODEL--  C
!  C                        L  CLM WEB INFO: http://clm.gsfc.nasa.gov
!  LMCLMCLMCLMCLMCLMCLMCLMCLM  CLM ListServ/Mailing List:
!
!-----------------------------------------------------------------------
! Purpose: 
! Surface albedo and two-stream fluxes
! 
! Method: 
! Surface albedos. Also fluxes (per unit incoming direct and diffuse
! radiation) reflected, transmitted, and absorbed by vegetation. 
! Also sunlit fraction of the canopy. 
!
! The calling sequence is:
!   -> SurfaceAlbedo:     albedos for next time step
!        -> SnowAge:      snow age
!        -> SnowAlbedo:   snow albedos: direct beam
!        -> SnowAlbedo:   snow albedos: diffuse
!        -> SoilAlbedo:   soil/lake/glacier/wetland albedos
!        -> TwoStream:    absorbed, reflected, transmitted solar fluxes (vis dir)
!        -> TwoStream:    absorbed, reflected, transmitted solar fluxes (vis dif)
!        -> TwoStream:    absorbed, reflected, transmitted solar fluxes (nir dir)
!        -> TwoStream:    absorbed, reflected, transmitted solar fluxes (nir dif)
! 
! Author: 
! Gordon Bonan
! 
!-----------------------------------------------------------------------
! $Id: SurfaceAlbedo.F90,v 1.7 2004/11/24 22:56:45 jim Exp $
!-----------------------------------------------------------------------
  use clm_varpar,  only : numrad
  use clmtype
  implicit none

!----Arguments----------------------------------------------------------

  type (clm1d), intent(inout) :: clm  !CLM 1-D Module

!----Local Variables----------------------------------------------------

  integer  :: ib              ! band index
  integer  :: ic              ! direct beam: ic=0; diffuse: ic=1
  integer  :: nband = numrad  ! number of solar radiation wave bands
  real(r8) :: wl              ! fraction of LAI+SAI that is LAI
  real(r8) :: ws              ! fraction of LAI+SAI that is SAI
  real(r8) :: mpe = 1.e-06    ! prevents overflow for division by zero 
  real(r8) :: vai             ! elai+esai
  real(r8) :: rho(numrad)     ! leaf/stem refl weighted by fraction LAI and SAI
  real(r8) :: tau(numrad)     ! leaf/stem tran weighted by fraction LAI and SAI
  real(r8) :: ftdi(numrad)    ! down direct flux below veg per unit dif flux = 0
  real(r8) :: albsnd(numrad)  ! snow albedo (direct)
  real(r8) :: albsni(numrad)  ! snow albedo (diffuse)
  real(r8) :: gdir            ! aver projected leaf/stem area in solar direction
  real(r8) :: ext             ! optical depth direct beam per unit LAI+SAI
  real(r8) :: delta           ! solar declination angle in radians
  real(r8) :: eccf            ! earth orbit eccentricity factor
  real(r8) :: coszen          ! cosine solar zenith angle for next time step
  real(r8) :: albd(2)
  real(r8) :: albi(2)

!----End Variable List--------------------------------------------------

!
! Solar declination  for next time step
!

!  call shr_orb_decl (caldayp1, eccen, mvelpp, lambm0, obliqr, &
!                     delta   , eccf )

!
! Cosine solar zenith angle for next time step
!

#if ( defined COUPLED )
!!  COSZ calculated in WRF and passed to LIS
    coszen = clm%forc_cosz 
#else
!   coszen = shr_orb_cosz(caldayp1, clm%lat, clm%lon, delta)
#endif

    coszen = clm%forc_cosz 

!!print*, 'coszen:', coszen

!!

!
! Initialize output because solar radiation only done if coszen > 0
!

  do ib = 1, nband
     albd(ib)   = 1.0
     albi(ib)   = 1.0
     clm%albgrd(ib) = 0._r8
     clm%albgri(ib) = 0._r8
     clm%fabd(ib)   = 0._r8
     clm%fabi(ib)   = 0._r8
     clm%ftdd(ib)   = 0._r8
     clm%ftid(ib)   = 0._r8
     clm%ftii(ib)   = 0._r8
     if (ib==1) clm%fsun = 0.
  end do

!===LDAS modification
!Original CLM2 code for solar zenith angle less than 0 is replaced with that in CLM1
!
! Return if coszen is not positive
!! Joe S. Quick Fix for Zenith Angle in LIS-WRF
!!  if (coszen <= 0._r8) then
!!  clm%surfalb = 0.0
!!   clm%snoalb =  0.0
!!  RETURN
!!  endif
!!
!===LDAS modification
! IF COSZEN IS NOT POSITIVE - NEW for CLM_OFFLINE
! NOTE: All NEW changes for CLM offline assumes that the incoming solar
! radiation is 70% direct, 30% diffuse, and 50% visible, 50% near-infrared
! (as is assumed in atmdrvMod.F90).

  do ib = 1, nband
    albsnd(ib)     = 0._r8
    albsni(ib)     = 0._r8
  enddo
  if (coszen <= 0._r8) then

   clm%surfalb = 35./100.*(albd(1)+albd(2)) &
                 +15./100.*(albi(1)+albi(2))
   clm%snoalb = 35./100.*(albsnd(1)+albsnd(2)) &
                 +15./100.*(albsni(1)+albsni(2))

#if ( defined COUPLED )
   clm%surfalb = 0.0
   clm%snoalb = 0.0
#else
!   clm%surfalb = LIS_rc%udef
!   clm%snoalb = LIS_rc%udef
#endif
   RETURN
  endif

!
! weight reflectance/transmittance by lai and sai
!

  do ib = 1, nband
     vai = clm%elai + clm%esai
     wl = clm%elai / max( vai,mpe )
     ws = clm%esai / max( vai,mpe )
     rho(ib) = max( clm%rhol(ib)*wl + clm%rhos(ib)*ws, mpe )
     tau(ib) = max( clm%taul(ib)*wl + clm%taus(ib)*ws, mpe )
  end do

!
! Snow albedos: only if h2osno > 0
!

  if ( clm%h2osno > 0._r8 ) then
     ic=0; call SnowAlbedo (clm, coszen, nband, ic, albsnd)
     ic=1; call SnowAlbedo (clm, coszen, nband, ic, albsni)  
  else
     albsnd(:) = 0._r8
     albsni(:) = 0._r8
  endif

!===LDAS modification
!===NEW for CLM offline

!  clm%snoalb = 35./100.*(albsnd(1)+albsnd(2)) + 15./100.*(albsni(1)+albsni(2))

!
! Ground surface albedos
!
  call SoilAlbedo (clm, coszen, nband, albsnd, albsni)      
  if (vai /= 0.) then  ! vegetated patch
!   if(clm%itypveg .ne.12) then 
!
! Loop over nband wavebands to calculate surface albedos and solar 
! fluxes for vegetated patch for unit incoming direct 
! (ic=0) and diffuse flux (ic=1)
!
     do ib = 1, nband
        ic = 0
        call TwoStream (clm,      ib,  ic,       coszen,   vai,      &
                        rho,      tau, clm%fabd, albd, clm%ftdd, &
                        clm%ftid, gdir )
        ic = 1
        call TwoStream (clm,      ib,  ic,       coszen,   vai,  &
                        rho,      tau, clm%fabi, albi, ftdi, &
                        clm%ftii, gdir )
     end do
!     print*, 'aft twostream.. ',clm%kpatch, clm%fabd,clm%fabi
!
! Sunlit fraction of canopy. Set fsun = 0 if fsun < 0.01.
!
     
     ext = gdir/coszen * sqrt(1.-rho(1)-tau(1))
     clm%fsun = (1.-exp(-ext*vai)) / max(ext*vai,mpe)
     ext = clm%fsun                                       !temporary fsun
     if (ext < 0.01) then 
        wl = 0._r8                                        !temporary fsun
     else
        wl = ext                                          !temporary fsun
     end if
     clm%fsun = wl

  else     ! non-vegetated patch

     do ib = 1,numrad
        clm%fabd(ib) = 0.
        clm%fabi(ib) = 0.
        clm%ftdd(ib) = 1.
        clm%ftid(ib) = 0.
        clm%ftii(ib) = 1.
        albd(ib) = clm%albgrd(ib)
        albi(ib) = clm%albgri(ib)
        clm%fsun     = 0.
     end do

  endif
!===LDAS modification
!===NEW for CLM offline:
  
    clm%surfalb = 35./100.*(albd(1)+albd(2)) +15./100.*(albi(1)+albi(2))
!!print*, 'SurfaceAlbedo:', clm%surfalb 
!    print*, 'end twostream.. ',clm%kpatch, clm%fabd,clm%fabi
  return
end subroutine SurfaceAlbedo
!=================================================
!=================================================

subroutine SnowAlbedo (clm, coszen, nband, ind, alb)

!-----------------------------------------------------------------------
!
!  CLMCLMCLMCLMCLMCLMCLMCLMCL  A community developed and sponsored, freely
!  L                        M  available land surface process model.
!  M --COMMON LAND MODEL--  C
!  C                        L  CLM WEB INFO: http://clm.gsfc.nasa.gov
!  LMCLMCLMCLMCLMCLMCLMCLMCLM  CLM ListServ/Mailing List:
!
!-----------------------------------------------------------------------
! Purpose: 
! Determine snow albedos
! 
! Method: 
! 
! Author:
! Gordon Bonan
! 
!-----------------------------------------------------------------------
! $Id: SnowAlbedo.F90,v 1.6 2004/11/24 22:56:36 jim Exp $
!-----------------------------------------------------------------------

  use clmtype
  use clm_varpar,   only : numrad
  implicit none

!----Arguments----------------------------------------------------------

  type (clm1d), intent(inout) :: clm       !CLM 1-D Module

  real(r8), intent(in) :: coszen ! cosine solar zenith angle for next time step
  integer, intent(in) :: nband   ! number of solar radiation waveband classes
  integer, intent(in) :: ind     ! 0=direct beam, 1=diffuse radiation

  real(r8), intent(out):: alb(numrad) !snow albedo by waveband

!----Local Variables----------------------------------------------------

  integer  :: ib          !waveband class

! CLM variables and constants for snow albedo calculation

  real(r8) :: snal0 = 0.95 !vis albedo of new snow for sza<60
  real(r8) :: snal1 = 0.65 !nir albedo of new snow for sza<60
  real(r8) :: conn  = 0.5  !constant for visible snow alb calculation [-]
  real(r8) :: cons  = 0.2  !constant (=0.2) for nir snow albedo calculation [-]
  real(r8) :: sl    = 2.0  !factor that helps control alb zenith dependence [-]
  real(r8) :: age          !factor to reduce visible snow alb due to snow age [-]
  real(r8) :: albs         !temporary vis snow albedo
  real(r8) :: albl         !temporary nir snow albedo
  real(r8) :: cff          !snow alb correction factor for zenith angle > 60 [-]
  real(r8) :: czf          !solar zenith correction for new snow albedo [-]

!----End Variable List--------------------------------------------------

!
! Zero albedos
!

  do ib = 1, nband
     alb(ib) = 0._r4
  end do

!
! CLM Albedo for snow cover.
! Snow albedo depends on snow-age, zenith angle, and thickness of snow,
! age gives reduction of visible radiation
!

!
! Correction for snow age

!  print*,'sn..',iam, clm%snowage
  age = 1.-1./(1.+clm%snowage)
  albs = snal0*(1.-cons*age)
  albl = snal1*(1.-conn*age)

  if (ind == 0) then

!
! czf corrects albedo of new snow for solar zenith
!

    cff    = ((1.+1./sl)/(1.+max(0.001_r4,coszen)*2.*sl )- 1./sl)
    cff    = max(cff,0._r4)
    czf    = 0.4*cff*(1.-albs)
    albs = albs+czf
    czf    = 0.4*cff*(1.-albl)
    albl = albl+czf

  endif

  alb(1) = albs
  alb(2) = albl
  return
end subroutine SnowAlbedo

!=================================================
!=================================================
subroutine SoilAlbedo (clm, coszen, nband, albsnd, albsni)

!-----------------------------------------------------------------------
!
!  CLMCLMCLMCLMCLMCLMCLMCLMCL  A community developed and sponsored, freely
!  L                        M  available land surface process model.
!  M --COMMON LAND MODEL--  C
!  C                        L  CLM WEB INFO: http://clm.gsfc.nasa.gov
!  LMCLMCLMCLMCLMCLMCLMCLMCLM  CLM ListServ/Mailing List:
!
!-----------------------------------------------------------------------
! Purpose: 
! Determine ground surface albedo, accounting for snow
! 
! Method: 
! 
! Author: 
! Gordon Bonan
! 
!-----------------------------------------------------------------------
! $Id: SoilAlbedo.F90,v 1.7 2004/11/24 22:56:39 jim Exp $
!-----------------------------------------------------------------------
  use clmtype
  use clm_varpar, only : numrad
  use clm_varcon, only : albsat, albdry, alblak, albice, tfrz, istice, istsoil
!  use spmdMod
  implicit none

!----Arguments----------------------------------------------------------

  type (clm1d), intent(inout) :: clm            !CLM 1-D Module

  integer, intent(in) :: nband           ! number of solar radiation waveband classes
  real(r8), intent(in) :: coszen         ! cosine solar zenith angle for next time step
  real(r8), intent(in) :: albsnd(numrad) ! snow albedo (direct)
  real(r8), intent(in) :: albsni(numrad) ! snow albedo (diffuse)

!----Local Variables----------------------------------------------------

  integer  ib      !waveband number (1=vis, 2=nir)
  real(r8) inc     !soil water correction factor for soil albedo
  real(r8) albsod  !soil albedo (direct)
  real(r8) albsoi  !soil albedo (diffuse)
!----End Variable List--------------------------------------------------

  do ib = 1, nband
     if (clm%itypwat == istsoil)  then               !soil
        inc    = max(0.11-0.40*clm%h2osoi_vol(1), 0._r4)
        albsod = min(albsat(clm%isoicol,ib)+inc, albdry(clm%isoicol,ib))
        albsoi = albsod
     else if (clm%itypwat == istice)  then           !land ice
        albsod = albice(ib)
        albsoi = albsod
     else if (clm%t_grnd > tfrz) then                !unfrozen lake, wetland
        albsod = 0.05/(max(0.001_r4,coszen) + 0.15)
        albsoi = albsod
     else                                            !frozen lake, wetland
        albsod = alblak(ib)
        albsoi = albsod
     end if
     clm%albgrd(ib) = albsod*(1.-clm%frac_sno) + albsnd(ib)*clm%frac_sno
     clm%albgri(ib) = albsoi*(1.-clm%frac_sno) + albsni(ib)*clm%frac_sno
  end do
  return
end subroutine SoilAlbedo

!=================================================
!=================================================

subroutine TwoStream     (clm, ib , ic , coszen, vai, &
                          rho, tau, fab, fre   , ftd, &
                          fti, gdir )
       
!-----------------------------------------------------------------------
!
!  CLMCLMCLMCLMCLMCLMCLMCLMCL  A community developed and sponsored, freely
!  L                        M  available land surface process model.
!  M --COMMON LAND MODEL--  C
!  C                        L  CLM WEB INFO: http://clm.gsfc.nasa.gov
!  LMCLMCLMCLMCLMCLMCLMCLMCLM  CLM ListServ/Mailing List:
!
!-----------------------------------------------------------------------
! Purpose: 
! Two-stream fluxes for canopy radiative transfer
! 
! Method: 
! Use two-stream approximation of Dickinson (1983) Adv Geophysics
! 25:305-353 and Sellers (1985) Int J Remote Sensing 6:1335-1372
! to calculate fluxes absorbed by vegetation, reflected by vegetation,
! and transmitted through vegetation for unit incoming direct or diffuse 
! flux given an underlying surface with known albedo.
! 
! Author: 
! Gordon Bonan
! 
!-----------------------------------------------------------------------
! $Id: TwoStream.F90,v 1.8 2004/11/24 22:56:49 jim Exp $
!-----------------------------------------------------------------------

  use clmtype
  use clm_varpar, only : numrad
  use clm_varcon, only : omegas, tfrz, betads, betais
  implicit none

!----Arguments----------------------------------------------------------

  type (clm1d), intent(inout) :: clm   !CLM 1-D Module

  integer , intent(in)  :: ib          ! waveband number 
  integer , intent(in)  :: ic          ! 0=unit incoming direct; 1=unit incoming diffuse
  real(r8), intent(in)  :: coszen      ! cosine solar zenith angle for next time step
  real(r8), intent(in)  :: vai         ! elai+esai
  real(r8), intent(in)  :: rho(numrad) ! leaf/stem refl weighted by fraction LAI and SAI
  real(r8), intent(in)  :: tau(numrad) ! leaf/stem tran weighted by fraction LAI and SAI

  real(r8), intent(out) :: fab(numrad) ! flux abs by veg layer (per unit incoming flux)   
  real(r8), intent(out) :: fre(numrad) ! flux refl above veg layer (per unit incoming flx)
  real(r8), intent(out) :: ftd(numrad) ! down dir flux below veg layer (per unit in flux) 
  real(r8), intent(out) :: fti(numrad) ! down dif flux below veg layer (per unit in flux) 
  real(r8), intent(out) :: gdir        ! aver projected leaf/stem area in solar direction

!----Local Variables----------------------------------------------------

!  integer i,j       ! array indices
  real(r8) cosz     ! 0.001 <= coszen <= 1.000
  real(r8) asu      ! single scattering albedo
  real(r8) chil     ! -0.4 <= xl <= 0.6
  real(r8) tmp0,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9
  real(r8) p1,p2,p3,p4,s1,s2,u1,u2,u3
  real(r8) b,c,d,d1,d2,f,h,h1,h2,h3,h4,h5,h6,h7,h8,h9,h10
  real(r8) phi1,phi2,sigma
  real(r8) ftds,ftis,fres
  real(r8) twostext ! optical depth of direct beam per unit leaf area 
  real(r8) avmu     ! average diffuse optical depth
  real(r8) omega    ! fraction of intercepted radiation that is scattered
  real(r8) omegal   ! omega for leaves
  real(r8) betai    ! upscatter parameter for diffuse radiation 
  real(r8) betail   ! betai for leaves
  real(r8) betad    ! upscatter parameter for direct beam radiation 
  real(r8) betadl   ! betad for leaves

!----End Variable List--------------------------------------------------

!
! Calculate two-stream parameters omega, betad, betai, avmu, gdir, twostext.
! Omega, betad, betai are adjusted for snow. Values for omega*betad 
! and omega*betai are calculated and then divided by the new omega
! because the product omega*betai, omega*betad is used in solution. 
! Also, the transmittances and reflectances (tau, rho) are linear 
! weights of leaf and stem values.
!
  cosz = max(0.001_r4, coszen)
  chil = min( max(clm%xl, -0.4_r4), 0.6_r4 )
  if (abs(chil) <= 0.01) chil = 0.01
  phi1 = 0.5 - 0.633*chil - 0.330*chil*chil
  phi2 = 0.877 * (1.-2.*phi1)
  gdir = phi1 + phi2*cosz
  twostext = gdir/cosz
  avmu = ( 1. - phi1/phi2 * log((phi1+phi2)/phi1) ) / phi2
  omegal = rho(ib) + tau(ib)
  tmp0 = gdir + phi2*cosz
  tmp1 = phi1*cosz
  asu = 0.5*omegal*gdir/tmp0 * ( 1. - tmp1/tmp0 * log((tmp1+tmp0)/tmp1) )
  betadl = (1.+avmu*twostext)/(omegal*avmu*twostext)*asu
  betail = 0.5 * ((rho(ib)+tau(ib)) + (rho(ib)-tau(ib)) * ((1.+chil)/2.)**2) / omegal

!
! Adjust omega, betad, and betai for intercepted snow
!

  if (clm%t_veg > tfrz) then                             !no snow
     tmp0 = omegal           
     tmp1 = betadl 
     tmp2 = betail 
  else
     tmp0 =   (1.-clm%fwet)*omegal        + clm%fwet*omegas(ib)           
     tmp1 = ( (1.-clm%fwet)*omegal*betadl + clm%fwet*omegas(ib)*betads ) / tmp0
     tmp2 = ( (1.-clm%fwet)*omegal*betail + clm%fwet*omegas(ib)*betais ) / tmp0
  end if
  omega = tmp0           
  betad = tmp1 
  betai = tmp2  

!
! Absorbed, reflected, transmitted fluxes per unit incoming radiation
!

  b = 1. - omega + omega*betai
  c = omega*betai
  tmp0 = avmu*twostext
  d = tmp0 * omega*betad
  f = tmp0 * omega*(1.-betad)
  tmp1 = b*b - c*c
  h = sqrt(tmp1) / avmu
  sigma = tmp0*tmp0 - tmp1
  if(sigma.lt.1e-10) then
     omega = omega*0.98
     b = 1. - omega + omega*betai
     c = omega*betai
     tmp0 = avmu*twostext
     d = tmp0 * omega*betad
     f = tmp0 * omega*(1.-betad)
     tmp1 = b*b - c*c
     h = sqrt(tmp1) / avmu
     sigma = tmp0*tmp0 - tmp1
  endif

  p1 = b + avmu*h
  p2 = b - avmu*h
  p3 = b + tmp0
  p4 = b - tmp0
  s1 = exp(-h*vai)
  s2 = exp(-twostext*vai)
  if (ic == 0) then
     u1 = b - c/clm%albgrd(ib)
     u2 = b - c*clm%albgrd(ib)
     u3 = f + c*clm%albgrd(ib)
  else
     u1 = b - c/clm%albgri(ib)
     u2 = b - c*clm%albgri(ib)
     u3 = f + c*clm%albgri(ib)
  end if
  tmp2 = u1 - avmu*h
  tmp3 = u1 + avmu*h
  d1 = p1*tmp2/s1 - p2*tmp3*s1
  tmp4 = u2 + avmu*h
  tmp5 = u2 - avmu*h
  d2 = tmp4/s1 - tmp5*s1
  h1 = -d*p4 - c*f
  tmp6 = d - h1*p3/sigma
  tmp7 = ( d - c - h1/sigma*(u1+tmp0) ) * s2
  h2 = ( tmp6*tmp2/s1 - p2*tmp7 ) / d1
   h3 = - ( tmp6*tmp3*s1 - p1*tmp7 ) / d1
  h4 = -f*p3 - c*d
  tmp8 = h4/sigma
  tmp9 = ( u3 - tmp8*(u2-tmp0) ) * s2
  h5 = - ( tmp8*tmp4/s1 + tmp9 ) / d2
  h6 = ( tmp8*tmp5*s1 + tmp9 ) / d2
  h7 = (c*tmp2) / (d1*s1)
  h8 = (-c*tmp3*s1) / d1
  h9 = tmp4 / (d2*s1)
  h10 = (-tmp5*s1) / d2

!
! Downward direct and diffuse fluxes below vegetation
!

  if (ic == 0) then
     ftds = s2
     ftis = h4*s2/sigma + h5*s1 + h6/s1
  else
     ftds = 0.
     ftis = h9*s1 + h10/s1
  end if
  ftd(ib) = ftds
  fti(ib) = ftis

!
! Flux reflected by vegetation
!

  if (ic == 0) then
     fres = h1/sigma + h2 + h3
  else
     fres = h7 + h8
  end if
  fre(ib) = fres

!
! Flux absorbed by vegetation

  fab(ib) = 1. - fre(ib) - (1.-clm%albgrd(ib))*ftd(ib) - (1.-clm%albgri(ib))*fti(ib)

  return
end subroutine TwoStream

!=================================================
!=================================================

subroutine Hydrology1 (clm)

!-----------------------------------------------------------------------
!
!  CLMCLMCLMCLMCLMCLMCLMCLMCL  A community developed and sponsored, freely   
!  L                        M  available land surface process model.  
!  M --COMMON LAND MODEL--  C  
!  C                        L  CLM WEB INFO: http://clm.gsfc.nasa.gov
!  LMCLMCLMCLMCLMCLMCLMCLMCLM  CLM ListServ/Mailing List: 
!
!-----------------------------------------------------------------------
! Purpose:
! Calculation of 
! (1) water storage of intercepted precipitation
! (2) direct throughfall and canopy drainage of precipitation
! (3) the fraction of foliage covered by water and the fraction
!     of foliage that is dry and transpiring. 
! (4) snow layer initialization if the snow accumulation exceeds 10 mm.
!
! Note:  The evaporation loss is taken off after the calculation of leaf 
! temperature in the subroutine clm_leaftem.f90, not in this subroutine.
!
! Method:
!
! Author:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision 
!-----------------------------------------------------------------------
! $Id: Hydrology1.F90,v 1.6 2004/11/24 22:56:26 jim Exp $
!-----------------------------------------------------------------------

  use clmtype
  use clm_varcon, only : tfrz, istice, istwet, istsoil
  implicit none

!----Arguments----------------------------------------------------------

  type (clm1d), intent(inout) :: clm        !CLM 1-D Module

!----Local Variables----------------------------------------------------

  integer newnode       ! flag when new snow node is set, (1=yes, 0=no)
  real(r8) prcp         ! precipitation rate [mm/s]
  real(r8) h2ocanmx     ! maximum allowed water on canopy [mm]
  real(r8) fpi          ! coefficient of interception
  real(r8) xrun         ! excess water that exceeds the leaf capacity [mm/s]
  real(r8) qflx_candrip ! rate of canopy runoff and snow falling off canopy [mm/s]
  real(r8) qflx_through ! direct throughfall [mm/s]
  real(r8) dz_snowf     ! layer thickness rate change due to precipitation [mm/s]
  real(r8) flfall       ! fraction of liquid water within falling precip.
  real(r8) bifall       ! bulk density of newly fallen dry snow [kg/m3]
!  real(r8) coefb        ! Slope of "Alta" expression for dependence of flfall on temp
!  real(r8) coefa        ! Offset of  of "Alta" expression for dependence of flfall on temp
  real(r8) :: dewmx
 
  real(r8) :: qflx_prec_intr
  real(r8) :: qflx_prec_grnd
  real(r8) :: qflx_snow_grnd
!----End Variable List--------------------------------------------------

!
! [1] Canopy interception and precipitation onto ground surface
!

!
! 1.1 Add precipitation to leaf water 
! 
  dewmx = 0.1
  if (clm%itypwat==istsoil .OR. clm%itypwat==istwet) then

     qflx_candrip = 0.                 ! rate of canopy runoff
     qflx_through = 0.                 ! precipitation direct through canopy
     qflx_prec_intr = 0.           ! intercepted precipitation  

     prcp = clm%forc_rain + clm%forc_snow  ! total precipitation

     if (clm%frac_veg_nosno == 1 .AND. prcp > 0.) then

!
! The leaf water capacities for solid and liquid are different, 
! generally double for snow, but these are of somewhat less significance
! for the water budget because of lower evap. rate at lower temperature.
! Hence, it is reasonable to assume that vegetation storage of solid water 
! is the same as liquid water.
!

        h2ocanmx = dewmx * (clm%elai + clm%esai)

!
! Direct throughfall
!

        fpi = 1. - exp(-0.5*(clm%elai + clm%esai))
        qflx_through  = prcp*(1.-fpi)

!
! Water storage of intercepted precipitation and dew
!

        qflx_prec_intr = prcp*fpi
        clm%h2ocan = max(0._r4, clm%h2ocan + clm%dtime*qflx_prec_intr)

!
! Initialize rate of canopy runoff and snow falling off canopy
!

        qflx_candrip = 0.0

!
! Excess water that exceeds the leaf capacity
!

        xrun = (clm%h2ocan - h2ocanmx)/clm%dtime

!
! Test on maximum dew on leaf
!

        if (xrun > 0.) then
           qflx_candrip = xrun
           clm%h2ocan = h2ocanmx
        endif

     endif

  else if (clm%itypwat == istice) then

     qflx_prec_intr = 0.
     clm%h2ocan = 0.
     qflx_candrip = 0.
     qflx_through = 0.  

  endif

!
! 1.2 Precipitation onto ground (kg/(m2 s))
!

  if (clm%frac_veg_nosno == 0) then
     qflx_prec_grnd = clm%forc_rain + clm%forc_snow
  else
     qflx_prec_grnd = qflx_through + qflx_candrip  
  endif
  
!
! 1.3 The percentage of liquid water by mass, which is arbitrarily set to 
!     vary linearly with air temp, from 0% at 273.16 to 40% max at 275.16.
!

  if (clm%itypprc <= 1) then
     flfall = 1.                        ! fraction of liquid water within falling precip.
     if (clm%do_capsnow) then
        clm%qflx_snowcap = qflx_prec_grnd
        qflx_snow_grnd = 0.
        clm%qflx_rain_grnd = 0.
     else
        clm%qflx_snowcap = 0.
        qflx_snow_grnd = 0.                  ! ice onto ground (mm/s)
        clm%qflx_rain_grnd = qflx_prec_grnd  ! liquid water onto ground (mm/s)
     endif
     dz_snowf = 0.                            ! rate of snowfall, snow depth/s (m/s)
  else
#if (defined PERGRO)
     coefb = 0.4/2.0
     coefa = -coefb*tfrz
     if (clm%forc_t <= tfrz) then
        flfall = 0.0
     else if (clm%forc_t <= tfrz+2.) then
        flfall = coefa + coefb*clm%forc_t
     else
        flfall = coefa + coefb*(tfrz+2.)
     endif
#else
     if (clm%forc_t <= tfrz) then
        flfall = 0.
     else if (clm%forc_t <= tfrz+2.) then
        flfall = -54.632 + 0.2*clm%forc_t
     else
        flfall = 0.4
     endif
#endif
     
!
! Use Alta relationship, Anderson(1976); LaChapelle(1961), 
! U.S.Department of Agriculture Forest Service, Project F, 
! Progress Rep. 1, Alta Avalanche Study Center:Snow Layer Densification.
!
#if (defined PERGRO)
     if (clm%forc_t > tfrz + 2.) then
        bifall=50. + 1.7*(17.0)**1.5
     else if (clm%forc_t > tfrz - 15.) then
        bifall=50. + 1.7*(clm%forc_t - tfrz + 15.)**1.5
     else
        bifall=50.
     endif
#else
     if (clm%forc_t > tfrz + 2.) then
        bifall =189.
     else if (clm%forc_t > tfrz - 15.) then
        bifall=50. + 1.7*(clm%forc_t - tfrz + 15.)**1.5
     else
        bifall=50.
     endif
#endif
     
     if (clm%do_capsnow) then
        clm%qflx_snowcap = qflx_prec_grnd
        qflx_snow_grnd = 0.
        clm%qflx_rain_grnd = 0.
        dz_snowf = 0.
     else
        clm%qflx_snowcap = 0.
        qflx_snow_grnd = qflx_prec_grnd*(1.-flfall)                 
        clm%qflx_rain_grnd = qflx_prec_grnd*flfall
        dz_snowf = qflx_snow_grnd/bifall                
        clm%snowdp = clm%snowdp + dz_snowf*clm%dtime         
        clm%h2osno = clm%h2osno + qflx_snow_grnd*clm%dtime  ! snow water equivalent (mm)
     endif

     if (clm%itypwat==istwet .AND. clm%t_grnd>tfrz) then
        clm%h2osno=0. 
        clm%snowdp=0. 
        clm%snowage=0.
     endif

  endif

!
! [2] Determine the fraction of foliage covered by water and the 
!     fraction of foliage that is dry and transpiring.
!

  call Fwet(clm)

!
! [3] When the snow accumulation exceeds 10 mm, initialize snow layer
! Currently, the water temperature for the precipitation is simply set 
! as the surface air temperature
!

  newnode = 0    ! flag for when snow node will be initialized
  if (clm%snl == 0 .AND. qflx_snow_grnd > 0.0 .AND. clm%snowdp >= 0.01) then  
     newnode = 1
     clm%snl = -1
     clm%dz(0) = clm%snowdp                       ! meter
     clm%z(0) = -0.5*clm%dz(0)
     clm%zi(-1) = -clm%dz(0)
     clm%snowage = 0.                             ! snow age
     clm%t_soisno (0) = min(tfrz, clm%forc_t)     ! K
     clm%h2osoi_ice(0) = clm%h2osno               ! kg/m2
     clm%h2osoi_liq(0) = 0.                       ! kg/m2
     clm%frac_iceold(0) = 1.
  endif

!
! The change of ice partial density of surface node due to precipitation.
! Only ice part of snowfall is added here, the liquid part will be added later
!

  if (clm%snl < 0 .AND. newnode == 0) then
     clm%h2osoi_ice(clm%snl+1) = clm%h2osoi_ice(clm%snl+1)+clm%dtime*qflx_snow_grnd
     clm%dz(clm%snl+1) = clm%dz(clm%snl+1)+dz_snowf*clm%dtime
  endif

end subroutine Hydrology1


!========================================================================

subroutine Fwet(clm)

  use clmtype
  implicit none

!=== Arguments ===========================================================

  type (clm1d), intent(inout) :: clm        !CLM 1-D Module

!=== Local Variables =====================================================

! fwet is the fraction of all vegetation surfaces which are wet 
! including stem area which contribute to evaporation.
! fdry is the fraction of elai which is dry because only leaves
! can transpire.  Adjusted for stem area which does not transpire.

  real(r8) vegt             ! frac_veg_nosno*lsai
  real(r8) dewmxi           ! inverse of maximum allowed dew [1/mm]
  real(r8) dewmx
!=== End Variable List ===================================================
  dewmx = 0.1
  if (clm%frac_veg_nosno == 1) then
     if (clm%h2ocan > 0.) then
        vegt     = clm%frac_veg_nosno*(clm%elai + clm%esai)
        dewmxi   = 1.0/dewmx
        clm%fwet = ((dewmxi/vegt)*clm%h2ocan)**.666666666666
        clm%fwet = min (clm%fwet,1.0_r4)     ! Check for maximum limit of fwet
     else
        clm%fwet = 0.
     endif
     clm%fdry = (1.-clm%fwet)*clm%elai/(clm%elai+clm%esai)
#if (defined PERGRO)
     clm%fwet = 0.
     clm%fdry = clm%elai/(clm%elai+clm%esai)
#endif
  else
     clm%fwet = 0.
     clm%fdry = 0.
  endif

end subroutine Fwet

subroutine updatelai( clm )
  use clmtype
  implicit none
  type (clm1d), intent(inout) :: clm        !CLM 1-D Module
  real(r8) :: ol      !thickness of canopy layer covered by snow (m)
  real(r8) :: fb      !fraction of canopy layer covered by snow

 
     ol = min( max(clm%snowdp-clm%hbot,0._r4), clm%htop-clm%hbot)
     fb = 1. - ol / max(1.e-06, clm%htop-clm%hbot)
     clm%elai = clm%tlai*fb
     clm%esai = clm%tsai*fb
     if (clm%elai < 0.05) clm%elai = 0._r4
     if (clm%esai < 0.05) clm%esai = 0._r4

! Fraction of vegetation free of snow

     if ((clm%elai + clm%esai) >= 0.05) then
        clm%frac_veg_nosno_alb = 1
     else
        clm%frac_veg_nosno_alb = 0
     endif

!     if (clm%itypveg .eq.LIS_rc%bareclass .and. (clm%elai + clm%esai) >= 0.05) then
!       clm%frac_veg_nosno_alb = 0
!     endif

     if (clm%itypveg.eq.0)  clm%frac_veg_nosno_alb = 0

end subroutine updatelai

!====================================
! canopy
!====================================

subroutine CanopyFluxes (z0mv,   z0hv,  z0qv,  thm,   th,     &
                         thv,    tg,    qg,    dqgdT, htvp,   &
                         emv,    emg,   dlrad, ulrad, cgrnds, &
                         cgrndl, cgrnd, clm    )

!-----------------------------------------------------------------------
!
!  CLMCLMCLMCLMCLMCLMCLMCLMCL  A community developed and sponsored, freely   
!  L                        M  available land surface process model.  
!  M --COMMON LAND MODEL--  C  
!  C                        L  CLM WEB INFO: http://clm.gsfc.nasa.gov
!  LMCLMCLMCLMCLMCLMCLMCLMCLM  CLM ListServ/Mailing List: 
!
!-----------------------------------------------------------------------
! Purpose:
! This subroutine:
! 1. Calculates the leaf temperature: 
! 2. Calculates the leaf fluxes, transpiration, photosynthesis and 
!    updates the dew accumulation due to evaporation.
!
! Method:
! Use the Newton-Raphson iteration to solve for the foliage 
! temperature that balances the surface energy budget:
!
! f(t_veg) = Net radiation - Sensible - Latent = 0
! f(t_veg) + d(f)/d(t_veg) * dt_veg = 0     (*)
!
! Note:
! (1) In solving for t_veg, t_grnd is given from the previous timestep.
! (2) The partial derivatives of aerodynamical resistances, which cannot 
!     be determined analytically, are ignored for d(H)/dT and d(LE)/dT
! (3) The weighted stomatal resistance of sunlit and shaded foliage is used 
! (4) Canopy air temperature and humidity are derived from => Hc + Hg = Ha
!                                                          => Ec + Eg = Ea
! (5) Energy loss is due to: numerical truncation of energy budget equation
!     (*); and "ecidif" (see the code) which is dropped into the sensible 
!     heat 
! (6) The convergence criteria: the difference, del = t_veg(n+1)-t_veg(n) and 
!     del2 = t_veg(n)-t_veg(n-1) less than 0.01 K, and the difference of 
!     water flux from the leaf between the iteration step (n+1) and (n) 
!     less than 0.1 W/m2; or the iterative steps over 40.
!
! Author:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision 
!
!-----------------------------------------------------------------------
! $Id: CanopyFluxes.F90,v 1.8 2004/11/24 22:56:20 jim Exp $
!-----------------------------------------------------------------------

  use clmtype
  use clm_varcon, only : sb, cpair, hvap, vkc, grav, denice, denh2o, tfrz
  use clm_varpar, only : nlevsoi

  implicit none

!----Arguments----------------------------------------------------------

  type (clm1d), intent(inout) :: clm !CLM 1-D Module

  real(r8), intent(in) :: z0mv  ! roughness length, momentum [m]
  real(r8), intent(in) :: z0hv  ! roughness length, sensible heat [m]
  real(r8), intent(in) :: z0qv  ! roughness length, latent heat [m]
  real(r8), intent(in) :: thm   ! intermediate variable (forc_t+0.0098*forc_hgt_t)
  real(r8), intent(in) :: th    ! potential temperature (kelvin)
  real(r8), intent(in) :: thv   ! virtual potential temperature (kelvin)
  real(r8), intent(in) :: tg    ! ground surface temperature [K]
  real(r8), intent(in) :: qg    ! specific humidity at ground surface [kg/kg]
  real(r8), intent(in) :: dqgdT ! temperature derivative of "qg"
  real(r8), intent(in) :: htvp  ! latent heat of evaporation (/sublimation) [J/kg]
  real(r8), intent(in) :: emv   ! ground emissivity
  real(r8), intent(in) :: emg   ! vegetation emissivity

  real(r8), intent(inout) :: cgrnd  ! deriv. of soil energy flux wrt to soil temp [w/m2/k]
  real(r8), intent(inout) :: cgrndl ! deriv, of soil sensible heat flux wrt soil temp 
                                    ! [w/m2/k]
  real(r8), intent(inout) :: cgrnds ! deriv of soil latent heat flux wrt soil temp 
                                    ! [w/m2/k]

  real(r8), intent(out) :: dlrad    ! downward longwave radiation below the canopy [W/m2]
  real(r8), intent(out) :: ulrad    ! upward longwave radiation above the canopy [W/m2]
  real(r8) :: parsha
  real(r8) :: btran
  real(r8) :: csoilc
!----Local Variables----------------------------------------------------

  real(r8) zldis       ! reference height "minus" zero displacement height [m]
  real(r8) zii         ! convective boundary layer height [m]
  real(r8) zeta        ! dimensionless height used in Monin-Obukhov theory
  real(r8) beta        ! coefficient of conective velocity [-]
  real(r8) wc          ! convective velocity [m/s]
  real(r8) dth         ! diff of virtual temp. between ref. height and surface 
  real(r8) dthv        ! diff of vir. poten. temp. between ref. height and surface
  real(r8) dqh         ! diff of humidity between ref. height and surface
  real(r8) obu         ! Monin-Obukhov length (m)
  real(r8) um          ! wind speed including the stablity effect [m/s]
  real(r8) ur          ! wind speed at reference height [m/s]
  real(r8) uaf         ! velocity of air within foliage [m/s]
  real(r8) temp1       ! relation for potential temperature profile
  real(r8) temp2       ! relation for specific humidity profile
  real(r8) temp3       ! relation for 2m potential temperature profile
  real(r8) temp4       ! relation for 2m specific humidity profile
  real(r8) ustar       ! friction velocity [m/s]
  real(r8) tstar       ! temperature scaling parameter
  real(r8) qstar       ! moisture scaling parameter
  real(r8) thvstar     ! virtual potential temperature scaling parameter
  real(r8) taf         ! air temperature within canopy space [K]
  real(r8) qaf         ! humidity of canopy air [kg/kg]
  real(r8) rpp         ! fraction of potential evaporation from leaf [-]
  real(r8) rppdry      ! fraction of potential evaporation through transp [-]
  real(r8) cf          ! heat transfer coefficient from leaves [-]
  real(r8) rb          ! leaf boundary layer resistance [s/m]
  real(r8) ram(2)      ! aerodynamical resistance [s/m]
  real(r8) rah(2)      ! thermal resistance [s/m]
  real(r8) raw(2)      ! moisture resistance [s/m]
  real(r8) wta         ! heat conductance for air [m/s]
  real(r8) wtg         ! heat conductance for ground [m/s]
  real(r8) wtl         ! heat conductance for leaf [m/s]
  real(r8) wta0        ! normalized heat conductance for air [-]
  real(r8) wtl0        ! normalized heat conductance for leaf [-]
  real(r8) wtg0        ! normalized heat conductance for ground [-]
  real(r8) wtal        ! normalized heat conductance for air and leaf [-]
  real(r8) wtgl        ! normalized heat conductance for leaf and ground [-]
  real(r8) wtga        ! normalized heat cond. for air and ground  [-]
  real(r8) wtaq        ! latent heat conductance for air [m/s]
  real(r8) wtlq        ! latent heat conductance for leaf [m/s]
  real(r8) wtgq        ! latent heat conductance for ground [m/s]
  real(r8) wtaq0       ! normalized latent heat conductance for air [-]
  real(r8) wtlq0       ! normalized latent heat conductance for leaf [-]
  real(r8) wtgq0       ! normalized heat conductance for ground [-]
  real(r8) wtalq       ! normalized latent heat cond. for air and leaf [-]
  real(r8) wtglq       ! normalized latent heat cond. for leaf and ground [-]
  real(r8) wtgaq       ! normalized latent heat cond. for air and ground [-]
  real(r8) el          ! vapor pressure on leaf surface [pa]
  real(r8) deldT       ! derivative of "el" on "t_veg" [pa/K]
  real(r8) qsatl       ! leaf specific humidity [kg/kg]
  real(r8) qsatldT     ! derivative of "qsatl" on "t_veg"
  real(r8) air,bir,cir ! atmos. radiation temporay set
  real(r8) dc1,dc2     ! derivative of energy flux [W/m2/K]
  real(r8) delt        ! temporary
  real(r8) delq        ! temporary
  real(r8) del         ! absolute change in leaf temp in current iteration [K]
  real(r8) del2        ! change in leaf temperature in previous iteration [K]
  real(r8) dele        ! change in latent heat flux from leaf [K]
  real(r8) delmax      ! maximum change in  leaf temperature [K]
  real(r8) dels        ! change in leaf temperature in current iteration [K]
  real(r8) det         ! maximum leaf temp. change in two consecutive iter [K]
  real(r8) dlemin      ! max limit for energy flux convergence [w/m2]
  real(r8) dtmin       ! max limit for temperature convergence [K]
  real(r8) efeb        ! latent heat flux from leaf (previous iter) [mm/s]
  real(r8) efeold      ! latent heat flux from leaf (previous iter) [mm/s]
  real(r8) efpot       ! potential latent energy flux [kg/m2/s]
  real(r8) efe         ! water flux from leaf [mm/s]
  real(r8) efsh        ! sensible heat from leaf [mm/s]
  real(r8) obuold      ! monin-obukhov length from previous iteration
  real(r8) tlbef       ! leaf temperature from previous iteration [K]
  real(r8) ecidif      ! excess energies [W/m2]
  real(r8) err         ! balance error
  real(r8) erre        ! balance error
  real(r8) co2         ! atmospheric co2 concentration (pa)
  real(r8) o2          ! atmospheric o2 concentration (pa)
  real(r8) svpts       ! saturation vapor pressure at t_veg (pa)
  real(r8) eah         ! canopy air vapor pressure (pa)
  real(r8) s_node      ! vol_liq/eff_porosity
  real(r8) smp_node    ! matrix potential
  real(r8) vol_ice(1:nlevsoi)      ! partial volume of ice lens in layer
  real(r8) eff_porosity(1:nlevsoi) ! effective porosity in layer
  real(r8) vol_liq(1:nlevsoi)      ! partial volume of liquid water in layer
  real(r8) rresis(1:nlevsoi)       ! soil water contribution to root resistance

! Constant atmospheric co2 and o2
  real(r8) po2                   ! partial pressure  o2 (mol/mol)
  real(r8) pco2                  ! partial pressure co2 (mol/mol)
  data po2,pco2 /0.209,355.e-06/

  real(r8) :: mpe = 1.e-6        ! prevents overflow error if division by zero

  integer i       ! loop index
  integer itlef   ! counter for leaf temperature iteration [-]
  integer itmax   ! maximum number of iteration [-]
  integer itmin   ! minimum number of iteration [-]
  integer nmozsgn ! number of times stability changes sign
  real(r8):: dt_veg
  real(r8) :: smpmax
!----End Variable List--------------------------------------------------

!
! Initialization
!
  smpmax = -1.5e5
  csoilc = 0.004

  del   = 0.0  ! change in leaf temperature from previous iteration
  itlef = 0    ! counter for leaf temperature iteration
  efeb  = 0.0  ! latent head flux from leaf for previous iteration

  wtlq = 0.0
  wtlq0 = 0.0
  wtgq0 = 0.0
  wtalq = 0.0
  wtgaq = 0.0
  wtglq = 0.0
  wtaq = 0.0
  wtgq = 0.0
  wtaq0 = 0.0
  wtlq0 = 0.0
  wtgq0 = 0.0
  wtalq = 0.0
  wtgaq = 0.0
  wtglq = 0.0

!
! Assign iteration parameters
!

  delmax = 1.0  ! maximum change in  leaf temperature
  itmax  = 40   ! maximum number of iteration
  itmin  = 2    ! minimum number of iteration
  dtmin  = 0.01 ! max limit for temperature convergence
  dlemin = 0.1  ! max limit for energy flux convergence

!
! Effective porosity of soil, partial volume of ice and liquid (needed
! for btran)
!

  do i = 1,nlevsoi
     vol_ice(i) = min(clm%watsat(i), clm%h2osoi_ice(i)/(clm%dz(i)*denice))
     eff_porosity(i) = clm%watsat(i)-vol_ice(i)
     vol_liq(i) = min(eff_porosity(i), clm%h2osoi_liq(i)/(clm%dz(i)*denh2o))
  enddo

!
! Root resistance factors
!

  btran = 1.e-10
  do i = 1,nlevsoi
     if (clm%t_soisno(i) > tfrz) then
        s_node = max(vol_liq(i)/eff_porosity(i),0.01_r4)
        smp_node = max(smpmax, -clm%sucsat(i)*s_node**(-clm%bsw(i)))
        rresis(i) = (1.-smp_node/smpmax)/(1.+clm%sucsat(i)/smpmax)
        clm%rootr(i) = clm%rootfr(i)*rresis(i)
        btran = btran + clm%rootr(i)
     else
        clm%rootr(i) = 0.
     endif
  enddo

!
! Normalize root resistances to get layer contribution to ET
!

  do i = 1,nlevsoi
     clm%rootr(i)  = clm%rootr(i)/btran
  enddo

!
! Net absorbed longwave radiation by canopy and ground
! =air+bir*t_veg**4+cir*t_grnd**4
!

  air =   emv * (1.+(1.-emv)*(1.-emg)) * clm%forc_lwrad
  bir = - (2.-emv*(1.-emg)) * emv * sb
  cir =   emv*emg*sb

!
! Saturated vapor pressure, specific humidity, and their derivatives
! at the leaf surface
!

  call QSatclm (clm%t_veg, clm%forc_pbot, el, deldT, qsatl, &
             qsatldT)

!
! Determine atmospheric co2 and o2
!

  co2 = pco2*clm%forc_pbot
  o2  = po2*clm%forc_pbot

!
! Initialize flux profile
!

  nmozsgn = 0
  obuold = 0.
  zii=1000.         ! m  (pbl height)
  beta=1.           ! -  (in computing W_*)

  taf = (tg + thm)/2.
  qaf = (clm%forc_q+qg)/2.

  ur = max(1.0_r8,sqrt(clm%forc_u*clm%forc_u+clm%forc_v*clm%forc_v))
  dth = thm-taf
  dqh = clm%forc_q-qaf
  dthv = dth*(1.+0.61*clm%forc_q)+0.61*th*dqh

  zldis = clm%forc_hgt_u - clm%displa

!
! Initialize Monin-Obukhov length and wind speed
!
  call MoninObukIni(ur, thv, dthv, zldis, z0mv, &
                    um, obu  )

!
! Begin stability iteration
!

  ITERATION : do while (itlef <= itmax) 

     tlbef = clm%t_veg
     del2 = del

!
! Determine friction velocity, and potential temperature and humidity
! profiles of the surface boundary layer
!


     call FrictionVelocity (clm%displa, z0mv,  z0hv,  z0qv,  obu, &
                            itlef+1, ur, um, ustar, temp1, temp2, &
                            temp3,temp4,clm)
    
!
! Determine aerodynamic resistances
!

!=== LDAS modification

      clm%ch = ustar*ustar/um ! Add-in for ALMA output
      clm%chs2 = temp3*ustar ! For calc of 2M T
      clm%cqs2 = temp4*ustar ! For calc of 2M Q

     ram(1)=1./(ustar*ustar/um)
     rah(1)=1./(temp1*ustar) 
     raw(1)=1./(temp2*ustar) 
     
!
! Bulk boundary layer resistance of leaves
!

     uaf = um*sqrt( 1./(ram(1)*um) )
     cf = 0.01/(sqrt(uaf)*sqrt(clm%dleaf))
     rb = 1./(cf*uaf)

!
! Aerodynamic resistances raw and rah between heights zpd+z0h and z0hg.
! if no vegetation, rah(2)=0 because zpd+z0h = z0hg.
! (Dickinson et al., 1993, pp.54)
!

     ram(2) = 0.               ! not used
     rah(2) = 1./(csoilc*uaf)
     raw(2) = rah(2) 

!
! Stomatal resistances for sunlit and shaded fractions of canopy.
! Done each iteration to account for differences in eah, tv.
!
     parsha = 0._r4 
     svpts = el                        ! pa
     eah = clm%forc_pbot * qaf / 0.622 ! pa

     call Stomata(mpe,        clm%parsun, svpts,     eah,      thm,        &
                  o2,         co2,        btran, rb,       clm%rssun,  &
                  clm%psnsun, clm%qe25,   clm%vcmx25,clm%mp,   clm%c3psn,  &
                  clm       )
                                                      
!write(*,*) rb,  clm%rssun,  clm%psnsun
     call Stomata(mpe,        parsha, svpts,     eah,      thm,        &
                  o2,         co2,        btran, rb,       clm%rssha,  &
                  clm%psnsha, clm%qe25,   clm%vcmx25,clm%mp,   clm%c3psn,  &
                  clm       )

!write(*,*) rb,  clm%rssha,  clm%psnsun
!
! Heat conductance for air, leaf and ground  
!

     call SensibleHCond(rah(1), rb,   rah(2), wta,  wtl,  &
                        wtg,    wta0, wtl0,   wtg0, wtal, &
                        wtga,   wtgl, clm     )

!
! Fraction of potential evaporation from leaf
!

     if (clm%fdry .gt. 0.0) then
        rppdry  = clm%fdry*rb*(clm%laisun/(rb+clm%rssun) + clm%laisha/(rb+clm%rssha))/ &
                  clm%elai
     else
        rppdry = 0.0
     endif
     efpot = clm%forc_rho*wtl*(qsatl-qaf)

     if (efpot > 0.) then
       if (btran > 1.e-10) then

        clm%qflx_tran_veg = efpot*rppdry
        rpp = rppdry + clm%fwet

!
! No transpiration if btran below 1.e-10
!

       else
        rpp = clm%fwet
        clm%qflx_tran_veg = 0.
       endif

!
! Check total evapotranspiration from leaves
!

       rpp = min(rpp, (clm%qflx_tran_veg+clm%h2ocan/clm%dtime)/efpot)

     else

!
! No transpiration if potential evaporation less than zero
!

       rpp = 1.
       clm%qflx_tran_veg = 0.

     endif

!
! Update conductances for changes in rpp 
! Latent heat conductances for ground and leaf.
! Air has same conductance for both sensible and latent heat.
!

     call LatentHCond(raw(1), rb,    raw(2), rpp,   wtaq,  &
                      wtlq,   wtgq,  wtaq0,  wtlq0, wtgq0, &
                      wtalq,  wtgaq, wtglq,  clm    ) 

     dc1 = clm%forc_rho*cpair*wtl
     dc2 = hvap*clm%forc_rho*wtlq

     efsh = dc1*(wtga*clm%t_veg-wtg0*tg-wta0*thm)
     efe = dc2*(wtgaq*qsatl-wtgq0*qg-wtaq0*clm%forc_q)

!
! Evaporation flux from foliage
!

     erre = 0.
     if (efe*efeb < 0.0) then
        efeold = efe
        efe  = 0.1*efeold
        erre = efe - efeold
     endif
     dt_veg = (clm%sabv + air + bir*clm%t_veg**4 + cir*tg**4 - efsh - efe) &
          / (- 4.*bir*clm%t_veg**3 +dc1*wtga +dc2*wtgaq*qsatldT)
     clm%t_veg = tlbef + dt_veg

     dels = clm%t_veg-tlbef
     del  = abs(dels)
     err = 0.
     if (del > delmax) then
        dt_veg = delmax*dels/del
        clm%t_veg = tlbef + dt_veg
        err = clm%sabv + air + bir*tlbef**3*(tlbef + 4.*dt_veg) &
             + cir*tg**4 - (efsh + dc1*wtga*dt_veg)          &
             - (efe + dc2*wtgaq*qsatldT*dt_veg)
     endif

!
! Fluxes from leaves to canopy space
! "efe" was limited as its sign changes frequently.  This limit may
! result in an imbalance in "hvap*qflx_evap_veg" and "efe + dc2*wtgaq*qsatldT*dt_veg" 
!

     efpot = clm%forc_rho*wtl*(wtgaq*(qsatl+qsatldT*dt_veg) &
          -wtgq0*qg-wtaq0*clm%forc_q)
     clm%qflx_evap_veg = rpp*efpot

!
! Calculation of evaporative potentials (efpot) and
! interception losses; flux in kg m**-2 s-1.  ecidif 
! holds the excess energy if all intercepted water is evaporated
! during the timestep.  This energy is later added to the
! sensible heat flux.
!

     ecidif = 0.
     if (efpot > 0. .AND. btran > 1.e-10) then
        clm%qflx_tran_veg = efpot*rppdry
     else
        clm%qflx_tran_veg = 0.
     endif
     ecidif = max(0._r4, clm%qflx_evap_veg-clm%qflx_tran_veg-clm%h2ocan/clm%dtime)
     clm%qflx_evap_veg = min(clm%qflx_evap_veg,clm%qflx_tran_veg+clm%h2ocan/clm%dtime)

!
! The energy loss due to above two limits is added to 
! the sensible heat flux.
!

     clm%eflx_sh_veg = efsh + dc1*wtga*dt_veg + err + erre +hvap*ecidif

!
! Re-calculate saturated vapor pressure, specific humidity, and their
! derivatives at the leaf surface
!

     call QSatclm(clm%t_veg, clm%forc_pbot, el, deldT, qsatl, &
               qsatldT    )

!
! Update vegetation/ground surface temperature, canopy air temperature, 
! canopy vapor pressure, aerodynamic temperature, and
! Monin-Obukhov stability parameter for next iteration. 
!

     taf = wtg0*tg + wta0*thm + wtl0*clm%t_veg
     qaf = wtlq0*qsatl+wtgq0*qg+clm%forc_q*wtaq0

!
! Update Monin-Obukhov length and wind speed including the stability effect
!

     dth = thm-taf       
     dqh = clm%forc_q-qaf

     tstar=temp1*dth
     qstar=temp2*dqh

     dthv=dth*(1.+0.61*clm%forc_q)+0.61*th*dqh

     thvstar=tstar*(1.+0.61*clm%forc_q) + 0.61*th*qstar
     zeta=zldis*vkc*grav*thvstar/(ustar**2*thv)

     if (zeta >= 0.) then     !stable
        zeta = min(2._r4,max(zeta,0.01_r4))
        um = max(ur,0.1_r4)
     else                     !unstable
        zeta = max(-100._r4,min(zeta,-0.01_r4))
        wc = beta*(-grav*ustar*thvstar*zii/thv)**0.333
        um = sqrt(ur*ur+wc*wc)
     endif
     obu = zldis/zeta

     if (obuold*obu < 0.) nmozsgn = nmozsgn+1
     if (nmozsgn >= 4) then 
        obu = zldis/(-0.01)
     endif

     obuold = obu

!
! Test for convergence
!

     itlef = itlef+1
     if (itlef > itmin) then
        dele = abs(efe-efeb)
        efeb = efe
        det  = max(del,del2)
        if (det < dtmin .AND. dele < dlemin) exit 
     endif

  enddo ITERATION     ! End stability iteration

!
! Energy balance check in canopy
!

  err = clm%sabv + air + bir*tlbef**3*(tlbef + 4.*dt_veg) &
       + cir*tg**4 - clm%eflx_sh_veg - hvap*clm%qflx_evap_veg
  if (abs(err) > 0.1) then
     write(6,*) 'energy balance in canopy X',err
  endif

!
! Fluxes from ground to canopy space 
!

  delt  = wtal*tg-wtl0*clm%t_veg-wta0*thm
  delq  = wtalq*qg-wtlq0*qsatl-wtaq0*clm%forc_q
  clm%taux  = -clm%forc_rho*clm%forc_u/ram(1)
  clm%tauy  = -clm%forc_rho*clm%forc_v/ram(1)
  clm%eflx_sh_grnd = cpair*clm%forc_rho*wtg*delt
  clm%qflx_evap_soi = clm%forc_rho*wtgq*delq

!
! 2 m height air temperature
!

  clm%t_ref2m   = clm%t_ref2m + (taf + temp1*dth * &
       1./vkc *log((2.+z0hv)/z0hv))

!
! Downward longwave radiation below the canopy    
!

  dlrad = (1.-emv)*emg*clm%forc_lwrad &
       + emv*emg * sb * &
       tlbef**3*(tlbef + 4.*dt_veg)

!
! Upward longwave radiation above the canopy    
!

  ulrad = ( (1.-emg)*(1.-emv)*(1.-emv)*clm%forc_lwrad &
       + emv*(1.+(1.-emg)*(1.-emv))*sb * tlbef**3 &
       *(tlbef + 4.*dt_veg) + emg *(1.-emv) *sb * tg**4)

!
! Derivative of soil energy flux with respect to soil temperature (cgrnd) 
!

  cgrnds = cgrnds + cpair*clm%forc_rho*wtg*wtal
  cgrndl = cgrndl + clm%forc_rho*wtgq*wtalq*dqgdT
  cgrnd  = cgrnds  + cgrndl*htvp

!
! Update dew accumulation (kg/m2) 
!

  clm%h2ocan = max(0._r4,clm%h2ocan + (clm%qflx_tran_veg-clm%qflx_evap_veg)*clm%dtime)

end subroutine CanopyFluxes

!======================================

subroutine MoninObukIni(ur, thv, dthv, zldis, z0m, &
                        um, obu  )

!-----------------------------------------------------------------------
!
!  CLMCLMCLMCLMCLMCLMCLMCLMCL  A community developed and sponsored, freely   
!  L                        M  available land surface process model.  
!  M --COMMON LAND MODEL--  C  
!  C                        L  CLM WEB INFO: http://clm.gsfc.nasa.gov
!  LMCLMCLMCLMCLMCLMCLMCLMCLM  CLM ListServ/Mailing List: 
!
!-----------------------------------------------------------------------
! Purpose:
! Initialization of the Monin-Obukhov length.
!
! Method:
! The scheme is based on the work of Zeng et al. (1998): 
! Intercomparison of bulk aerodynamic algorithms for the computation 
! of sea surface fluxes using TOGA CORE and TAO data. J. Climate, 
! Vol. 11, 2628-2644.
!
! Author:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision 
!
!-----------------------------------------------------------------------
! $Id: MoninObukIni.F90,v 1.6 2004/11/24 22:56:31 jim Exp $
!-----------------------------------------------------------------------

  use clm_varcon, only : grav
  implicit none

!----Arguments----------------------------------------------------------

  real(r8), intent(in) :: ur    ! wind speed at reference height [m/s]
  real(r8), intent(in) :: thv   ! virtual potential temperature (kelvin)
  real(r8), intent(in) :: dthv  ! diff of vir. poten. temp. between ref. height and surface
  real(r8), intent(in) :: zldis ! reference height "minus" zero displacement heght [m]
  real(r8), intent(in) :: z0m   ! roughness length, momentum [m]

  real(r8), intent(out) :: um   ! wind speed including the stability effect [m/s]
  real(r8), intent(out) :: obu  ! monin-obukhov length (m)

!----Local Variables----------------------------------------------------

  real(r8)  wc    ! convective velocity [m/s]
  real(r8)  rib   ! bulk Richardson number
  real(r8)  zeta  ! dimensionless height used in Monin-Obukhov theory
!  real(r8)  ustar ! friction velocity [m/s]     

!----End Variable List--------------------------------------------------

!
! Initial values of u* and convective velocity
!

!  ustar=0.06
  wc=0.5
  if (dthv >= 0.) then
     um=max(ur,0.1_r4)
  else
     um=sqrt(ur*ur+wc*wc)
  endif

  rib=grav*zldis*dthv/(thv*um*um)
#if (defined PERGRO)
  rib = 0.
#endif
  
  if (rib >= 0.) then      ! neutral or stable
     zeta = rib*log(zldis/z0m)/(1.-5.*min(rib,0.19_r4))
     zeta = min(2._r4,max(zeta,0.01_r4 ))
  else                    !unstable
     zeta=rib*log(zldis/z0m)
     zeta = max(-100._r4,min(zeta,-0.01_r4 ))
  endif

  obu=zldis/zeta

end subroutine MoninObukIni

!===============================

subroutine FrictionVelocity (displa, z0m,   z0h,   z0q,   obu, &
                             iter, ur, um, ustar, temp1, temp2, &
                             temp3, temp4, clm)

!-----------------------------------------------------------------------
!
!  CLMCLMCLMCLMCLMCLMCLMCLMCL  A community developed and sponsored, freely   
!  L                        M  available land surface process model.  
!  M --COMMON LAND MODEL--  C  
!  C                        L  CLM WEB INFO: http://clm.gsfc.nasa.gov
!  LMCLMCLMCLMCLMCLMCLMCLMCLM  CLM ListServ/Mailing List: 
!
!-----------------------------------------------------------------------
! Purpose:
! Calculation of the friction velocity, relation for potential 
! temperature and humidity profiles of surface boundary layer. 
!
! Method:
! The scheme is based on the work of Zeng et al. (1998): 
! Intercomparison of bulk aerodynamic algorithms for the computation 
! of sea surface fluxes using TOGA CORE and TAO data. J. Climate, 
! Vol. 11, 2628-2644.
!
! Author:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision 
!
!-----------------------------------------------------------------------
! $Id: FrictionVelocity.F90,v 1.6 2004/11/24 22:56:25 jim Exp $
!-----------------------------------------------------------------------

  use clmtype
  use clm_varcon, only : vkc
  implicit none

!----Arguments----------------------------------------------------------

  type (clm1d), intent(inout) :: clm !CLM 1-D Module

  real(r8), intent(in) :: displa ! displacement height [m]
  real(r8), intent(in) :: z0m    ! roughness length, momentum [m]
  real(r8), intent(in) :: z0h    ! roughness length, sensible heat [m]
  real(r8), intent(in) :: z0q    ! roughness length, latent heat [m]
  real(r8), intent(in) :: obu    ! monin-obukhov length (m)
  real(r8), intent(in) :: um     ! wind speed including the stablity effect [m/s]
  real(r8), intent(in) :: ur
  integer , intent(in) :: iter

  real(r8), intent(out) :: ustar ! friction velocity [m/s]
  real(r8), intent(out) :: temp1 ! relation for potential temperature profile
  real(r8), intent(out) :: temp2 ! relation for specific humidity profile
  real(r8), intent(out) :: temp3 ! relation for 2m potential temperature profile
  real(r8), intent(out) :: temp4 ! relation for 2m specific humidity profile
 
!----Local Variables----------------------------------------------------

  real(r8) zldis   ! reference height "minus" zero displacement heght [m]
!  real(r8) StabilityFunc ! stability function for unstable case
  real(r8) zetam   ! transition point of flux-gradient relation (wind profile)
  real(r8) zetat   ! transition point of flux-gradient relation (temp. profile)
  real(r8) zeta    ! dimensionless height used in Monin-Obukhov theory
  real(r8) tmp_stb1, tmp_stb2
#if (defined BGC)
! Variables used to diagnose the 10 meter wind
  real(r8) :: tmp1,tmp2,tmp3,tmp4
  real(r8) :: fmnew
  real(r8) :: fm10
  real(r8) :: zeta10
  real(r8), save :: fm
#endif

!----End Variable List--------------------------------------------------

!
! Adjustment factors for unstable (moz < 0) or stable (moz > 0) conditions.
! Wind profile
!

  zldis=clm%forc_hgt_u-displa
  zeta=zldis/obu
  zetam=1.574

  if (zeta < -zetam) then           ! zeta < -1
     call StabilityFunc(1,-zetam, tmp_stb1) 
     call StabilityFunc(1,z0m/obu, tmp_stb2) 
     ustar=vkc*um/(log(-zetam*obu/z0m)- &
          tmp_stb1 + tmp_stb2 &
          +1.14*((-zeta)**0.333-(zetam)**0.333))
  else if (zeta < 0.) then         ! -1 <= zeta < 0
     call StabilityFunc(1,zeta, tmp_stb1) 
     call StabilityFunc(1,z0m/obu, tmp_stb2) 
     ustar=vkc*um/(log(zldis/z0m)- &
           tmp_stb1+tmp_stb2)
  else if (zeta <= 1.) then        !  0 <= ztea <= 1
     ustar=vkc*um/(log(zldis/z0m) + &
          5.*zeta -5.*z0m/obu)
  else                             !  1 < zeta, phi=5+zeta
     ustar=vkc*um/(log(obu/z0m)+5.-5.*z0m/obu &
          +(5.*log(zeta)+zeta-1.))
  endif

#if (defined PERGRO)
  if (zeta < -zetam) then           ! zeta < -1
     ustar=vkc*um/log(-zetam*obu/z0m)
  else if (zeta < 0.) then         ! -1 <= zeta < 0
     ustar=vkc*um/log(zldis/z0m)
  else if (zeta <= 1.) then        !  0 <= ztea <= 1
     ustar=vkc*um/log(zldis/z0m)
  else                             !  1 < zeta, phi=5+zeta
     ustar=vkc*um/log(obu/z0m)
  endif
#endif


#if (defined BGC)
!
! diagnose 10-m wind for dust model (dstmbl.F)
! Notes from C. Zender's dst.F:
! According to Bon96 p. 62, the displacement height d (here displa) is
! 0.0 <= d <= 0.34 m in dust source regions (i.e., regions w/o trees).
! Therefore d <= 0.034*z1 and may safely be neglected.
! Code from LSM routine SurfaceTemperature was used to obtain u10
!
  if (min(zeta,1.) < 0.) then
     tmp1 = (1. - 16.*min(zeta,1.))**0.25
     tmp2 = log((1.+tmp1*tmp1)/2.)
     tmp3 = log((1.+tmp1)/2.)
     fmnew = 2.*tmp3 + tmp2 - 2.*atan(tmp1) + 1.5707963
  else
     fmnew = -5.*min(zeta,1.)
  endif
  if (iter == 1) then
     fm = fmnew
  else
     fm = 0.5 * (fm+fmnew)
  endif

  zeta10 = min(10./obu, 1.)
  if (zeta == 0.) zeta10 = 0.

  if (zeta10 < 0.) then
     tmp1 = (1.0 - 16.0 * zeta10)**0.25
     tmp2 = log((1.0 + tmp1*tmp1)/2.0)
     tmp3 = log((1.0 + tmp1)/2.0)
     fm10 = 2.0*tmp3 + tmp2 - 2.0*atan(tmp1) + 1.5707963
  else  ! not stable
     fm10 = -5.0 * zeta10
  endif ! not stable
  tmp4 = log(clm%forc_hgt / 10.)

!  clm%u10 = ur - ustar/vkc * (tmp4 - fm + fm10)
!  clm%fv  = ustar
#endif

!
! Temperature profile
!

  zldis=clm%forc_hgt_t-displa
  zeta=zldis/obu
  zetat=0.465
  if (zeta < -zetat) then           ! zeta < -1
     call StabilityFunc(2,-zetat, tmp_stb1) 
     call StabilityFunc(2,z0h/obu, tmp_stb2) 
     temp1=vkc/(log(-zetat*obu/z0h)- tmp_stb1 &
          + tmp_stb2 &
          + 0.8*((zetat)**(-0.333)-(-zeta)**(-0.333)))
  else if (zeta < 0.) then         ! -1 <= zeta < 0
     call StabilityFunc(2,zeta, tmp_stb1) 
     call StabilityFunc(2,z0h/obu, tmp_stb2) 
     temp1=vkc/(log(zldis/z0h) - tmp_stb1 + &
          tmp_stb2)
  else if (zeta <= 1.) then        !  0 <= ztea <= 1
     temp1=vkc/(log(zldis/z0h) + 5.*zeta - 5.*z0h/obu)
  else                             !  1 < zeta, phi=5+zeta
     temp1=vkc/(log(obu/z0h) + 5. - 5.*z0h/obu &
          + (5.*log(zeta)+zeta-1.))
  endif

  ! 2M Calculation
  zldis=2.0 + z0h
  zeta=zldis/obu
  zetat=0.465
  if (zeta < -zetat) then           ! zeta < -1
     call StabilityFunc(2,-zetat, tmp_stb1) 
     call StabilityFunc(2,z0h/obu, tmp_stb2) 
     temp3=vkc/(log(-zetat*obu/z0h)-tmp_stb1 &
          + tmp_stb2 &
          + 0.8*((zetat)**(-0.333)-(-zeta)**(-0.333)))
  else if (zeta < 0.) then         ! -1 <= zeta < 0
     call StabilityFunc(2,zeta, tmp_stb1) 
     call StabilityFunc(2,z0h/obu, tmp_stb2) 
     temp3=vkc/(log(zldis/z0h) - tmp_stb1 + &
          tmp_stb2)
  else if (zeta <= 1.) then        !  0 <= ztea <= 1
     temp3=vkc/(log(zldis/z0h) + 5.*zeta - 5.*z0h/obu)
  else                             !  1 < zeta, phi=5+zeta
     temp3=vkc/(log(obu/z0h) + 5. - 5.*z0h/obu &
          + (5.*log(zeta)+zeta-1.))
  endif

#if (defined PERGRO)
  if (zeta < -zetat) then           ! zeta < -1
     temp1=vkc/log(-zetat*obu/z0h)
  else if (zeta < 0.) then         ! -1 <= zeta < 0
     temp1=vkc/log(zldis/z0h)
  else if (zeta <= 1.) then        !  0 <= ztea <= 1
     temp1=vkc/log(zldis/z0h)
  else                             !  1 < zeta, phi=5+zeta
     temp1=vkc/log(obu/z0h)
  endif
#endif

!
! Humidity profile
!

  zldis=clm%forc_hgt_q-displa
  zeta=zldis/obu
  zetat=0.465
  if (zeta < -zetat) then          ! zeta < -1
     call StabilityFunc(2,-zetat, tmp_stb1) 
     call StabilityFunc(2,z0q/obu, tmp_stb2) 
     temp2=vkc/(log(-zetat*obu/z0q) - &
          tmp_stb1 + tmp_stb2 &
          + 0.8*((zetat)**(-0.333)-(-zeta)**(-0.333)))
  else if (zeta < 0.) then         ! -1 <= zeta < 0
     call StabilityFunc(2,zetat, tmp_stb1) 
     call StabilityFunc(2,z0q/obu, tmp_stb2) 
     temp2=vkc/(log(zldis/z0q) - &
          tmp_stb1+tmp_stb2)
  else if (zeta <= 1.) then        !  0 <= ztea <= 1
     temp2=vkc/(log(zldis/z0q)+5.*zeta-5.*z0q/obu)
  else                             !  1 < zeta, phi=5+zeta
     temp2=vkc/(log(obu/z0q) + 5. - 5.*z0q/obu &
          + (5.*log(zeta)+zeta-1.))
  endif

 !2m Calculation
  zldis=2.0 + z0q
  zeta=zldis/obu
  zetat=0.465
  if (zeta < -zetat) then          ! zeta < -1
     call StabilityFunc(2,-zetat, tmp_stb1) 
     call StabilityFunc(2,z0q/obu, tmp_stb2) 
     temp4=vkc/(log(-zetat*obu/z0q) - &
          tmp_stb1 + tmp_stb2 &
          + 0.8*((zetat)**(-0.333)-(-zeta)**(-0.333)))
  else if (zeta < 0.) then         ! -1 <= zeta < 0
     call StabilityFunc(2,zetat, tmp_stb1) 
     call StabilityFunc(2,z0q/obu, tmp_stb2) 
     temp4=vkc/(log(zldis/z0q) - &
          tmp_stb1+tmp_stb2)
  else if (zeta <= 1.) then        !  0 <= ztea <= 1
     temp4=vkc/(log(zldis/z0q)+5.*zeta-5.*z0q/obu)
  else                             !  1 < zeta, phi=5+zeta
     temp4=vkc/(log(obu/z0q) + 5. - 5.*z0q/obu &
          + (5.*log(zeta)+zeta-1.))
  endif

#if (defined PERGRO)
  if (zeta < -zetat) then          ! zeta < -1
     temp2=vkc/log(-zetat*obu/z0q)
  else if (zeta < 0.) then         ! -1 <= zeta < 0
     temp2=vkc/log(zldis/z0q)
  else if (zeta <= 1.) then        !  0 <= ztea <= 1
     temp2=vkc/log(zldis/z0q)
  else                             !  1 < zeta, phi=5+zeta
     temp2=vkc/log(obu/z0q)
  endif
#endif

end subroutine FrictionVelocity

!==================================

subroutine StabilityFunc(k, zeta, tmp_out)

!-----------------------------------------------------------------------
! Purpose:
! Stability function for rib < 0.
!
! Method:
!
! Author:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision 
!
!-----------------------------------------------------------------------
! $Id: StabilityFunc.F90,v 1.6 2004/11/24 22:56:43 jim Exp $
!-----------------------------------------------------------------------

  use shr_const_mod, only: SHR_CONST_PI
  implicit none

!----Local Variables----------------------------------------------------

  integer k         !
  real(r8) zeta     ! dimensionless height used in Monin-Obukhov theory
!  real(r8) StabilityFunc  ! stability function for unstable case
  real(r8) chik, tmp_out     ! 

!=== End Variable List ===================================================

  chik = (1.-16.*zeta)**0.25
  if (k == 1) then
     tmp_out = 2.*log((1.+chik)*0.5) &
          + log((1.+chik*chik)*0.5)-2.*atan(chik)+SHR_CONST_PI*0.5
  else
     tmp_out = 2.*log((1.+chik*chik)*0.5)
  endif

end subroutine StabilityFunc

!====================================

subroutine Stomata(mpe,  apar,   ei,    ea,   tgcm,   &
                   o2,   co2,    btran, rb,   rs,     &
                   psn,  qe25,   vcmx25,mp,   c3psn,  &
                   clm   ) 

!-----------------------------------------------------------------------
!
!  CLMCLMCLMCLMCLMCLMCLMCLMCL  A community developed and sponsored, freely   
!  L                        M  available land surface process model.  
!  M --COMMON LAND MODEL--  C  
!  C                        L  CLM WEB INFO: http://clm.gsfc.nasa.gov
!  LMCLMCLMCLMCLMCLMCLMCLMCLM  CLM ListServ/Mailing List: 
!
!-----------------------------------------------------------------------
! Purpose:
! Leaf stomatal resistance and leaf photosynthesis.
!
! Method:
!
! Author:
! author:            Gordon Bonan
! standardized:      J. Truesdale, Feb. 1996
! reviewed:          G. Bonan, Feb. 1996
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision 
!
!-----------------------------------------------------------------------
! $Id: Stomata.F90,v 1.6 2004/11/24 22:56:44 jim Exp $
!-----------------------------------------------------------------------

  use clmtype
  use clm_varcon   , only : tfrz
  use shr_const_mod, only : SHR_CONST_TKFRZ,SHR_CONST_RGAS
  implicit none

!----Arguments----------------------------------------------------------

  type (clm1d), intent(inout) :: clm    !CLM 1-D Module

  real(r8), intent(in) :: mpe    ! prevents division by zero errors
  real(r8), intent(in) :: ei     ! vapor pressure inside leaf (sat vapor press at t_veg) (pa)
  real(r8), intent(in) :: ea     ! vapor pressure of canopy air (pa)
  real(r8), intent(in) :: apar   ! par absorbed per unit lai (w/m**2)
  real(r8), intent(in) :: o2     ! atmospheric o2 concentration (pa)
  real(r8), intent(in) :: co2    ! atmospheric co2 concentration (pa)
  real(r8), intent(in) :: tgcm   ! air temperature at agcm reference height (kelvin)
  real(r8), intent(in) :: btran  ! soil water transpiration factor (0 to 1)
  real(r8), intent(in) :: qe25   ! quantum efficiency at 25c (umol co2 / umol photon)
  real(r8), intent(in) :: vcmx25 ! maximum rate of carboxylation at 25c (umol co2/m**2/s)
  real(r8), intent(in) :: mp     ! slope for conductance-to-photosynthesis relationship 
  real(r8), intent(in) :: c3psn  ! photosynthetic pathway: 0. = c4, 1. = c3

  real(r8), intent(inout) :: rb  ! boundary layer resistance (s/m)

  real(r8), intent(out)   :: rs  ! leaf stomatal resistance (s/m)
  real(r8), intent(out)   :: psn ! foliage photosynthesis (umol co2 /m**2/ s) [always +]

!----Local Variables----------------------------------------------------

  integer, parameter :: niter = 3  ! number of iterations
  integer  iter                    ! iteration index

  real(r8) ab      ! used in statement functions
  real(r8) bc      ! used in statement functions
  real(r8) f1      ! generic temperature response (statement function)
  real(r8) f2      ! generic temperature inhibition (statement function)
  real(r8) tc      ! foliage temperature (degree celsius)
  real(r8) cs      ! co2 concentration at leaf surface (pa)
  real(r8) kc      ! co2 michaelis-menten constant (pa)
  real(r8) ko      ! o2 michaelis-menten constant (pa)
  real(r8) a,b,c,q ! intermediate calculations for rs
  real(r8) r1,r2   ! roots for rs
  real(r8) ppf     ! absorb photosynthetic photon flux (umol photons/m**2/s)
  real(r8) wc      ! rubisco limited photosynthesis (umol co2/m**2/s)
  real(r8) wj      ! light limited photosynthesis (umol co2/m**2/s)
  real(r8) we      ! export limited photosynthesis (umol co2/m**2/s)
  real(r8) cp      ! co2 compensation point (pa)
  real(r8) ci      ! internal co2 (pa)
  real(r8) awc     ! intermediate calcuation for wc
  real(r8) vcmx    ! maximum rate of carboxylation (umol co2/m**2/s)
#if (defined DGVM)
  real(r8) vcmxpot ! potential maximum rate of carboxylation (umol co2/m**2/s)
#endif
  real(r8) j       ! electron transport (umol co2/m**2/s)
  real(r8) cea     ! constrain ea or else model blows up
  real(r8) cf      ! s m**2/umol -> s/m
  real(r8) rsmax0  ! maximum stomatal resistance [s/m]

  real(r8) kc25    ! co2 michaelis-menten constant at 25c (pa)
  real(r8) akc     ! q10 for kc25
  real(r8) ko25    ! o2 michaelis-menten constant at 25c (pa)
  real(r8) ako     ! q10 for ko25
  real(r8) avcmx   ! q10 for vcmx25
  real(r8) bp      ! minimum leaf conductance (umol/m**2/s)

!----End Variable List--------------------------------------------------

  f1(ab,bc) = ab**((bc-25.)/10.)
  f2(ab) = 1. + exp((-2.2e05+710.*(ab+SHR_CONST_TKFRZ))/(SHR_CONST_RGAS*0.001*(ab+SHR_CONST_TKFRZ)))

  kc25 = 30.
  akc = 2.1
  ko25 = 30000.
  ako = 1.2
  avcmx = 2.4
  bp = 2000.

!
! Initialize rs=rsmax and psn=0 because calculations are performed only
! when apar > 0, in which case rs <= rsmax and psn >= 0
! Set constants
!

  rsmax0 = 2.e4
  cf = clm%forc_pbot/(SHR_CONST_RGAS*0.001*tgcm)*1.e06 

  if (apar <= 0.) then          ! night time
     rs = min(rsmax0, 1./bp * cf)
     psn = 0.
     return
  else                          ! day time
     tc = clm%t_veg-tfrz                            
     ppf = 4.6*apar                  
     j = ppf*qe25
     kc = kc25 * f1(akc,tc)       
     ko = ko25 * f1(ako,tc)
     awc = kc * (1.+o2/ko)
     cp = 0.5*kc/ko*o2*0.21
     vcmx = vcmx25 * f1(avcmx,tc) / f2(tc) * btran
#if (defined DGVM)
     vcmxpot = vcmx25 * f1(avcmx,tc) / f2(tc) 
#endif

!
! First guess ci
!

     ci = 0.7*co2*c3psn + 0.4*co2*(1.-c3psn)  

!
! rb: s/m -> s m**2 / umol
!

     rb = rb/cf 

!
! Constrain ea
!

     cea = max(0.25*ei*c3psn+0.40*ei*(1.-c3psn), min(ea,ei) ) 

#if (defined DGVM)
!
! ci iteration for 'potential' photosynthesis
!

     do iter = 1, niter
        wj = max(ci-cp,0._r4)*j/(ci+2.*cp)*c3psn + j*(1.-c3psn)
        wc = max(ci-cp,0._r4)*vcmxpot/(ci+awc)*c3psn + vcmxpot*(1.-c3psn)
        we = 0.5*vcmxpot*c3psn + 4000.*vcmxpot*ci/clm%forc_pbot*(1.-c3psn) 
        psn = min(wj,wc,we) 
        cs = max( co2-1.37*rb*clm%forc_pbot*psn, mpe )
        a = mp*psn*clm%forc_pbot*cea / (cs*ei) + bp
        b = ( mp*psn*clm%forc_pbot/cs + bp ) * rb - 1.
        c = -rb
        if (b >= 0.) then
           q = -0.5*( b + sqrt(b*b-4.*a*c) )
        else
           q = -0.5*( b - sqrt(b*b-4.*a*c) )
        endif
        r1 = q/a
        r2 = c/q
        rs = max(r1,r2)
        ci = max( cs-psn*clm%forc_pbot*1.65*rs, 0._r4 )
     enddo

     clm%annpsnpot = clm%annpsnpot + psn

! First guess ci

     ci = 0.7*co2*c3psn + 0.4*co2*(1.-c3psn)  
#endif

!
! ci iteration for 'actual' photosynthesis
!

     do iter = 1, niter
        wj = max(ci-cp,0._r4)*j/(ci+2.*cp)*c3psn + j*(1.-c3psn)
        wc = max(ci-cp,0._r4)*vcmx/(ci+awc)*c3psn + vcmx*(1.-c3psn)
        we = 0.5*vcmx*c3psn + 4000.*vcmx*ci/clm%forc_pbot*(1.-c3psn) 
        psn = min(wj,wc,we) 
!write(*,*) psn
        cs = max( co2-1.37*rb*clm%forc_pbot*psn, mpe )
        a = mp*psn*clm%forc_pbot*cea / (cs*ei) + bp
        b = ( mp*psn*clm%forc_pbot/cs + bp ) * rb - 1.
        c = -rb
        if (b >= 0.) then
           q = -0.5*( b + sqrt(b*b-4.*a*c) )
        else
           q = -0.5*( b - sqrt(b*b-4.*a*c) )
        endif
        r1 = q/a
        r2 = c/q
        rs = max(r1,r2)
        ci = max( cs-psn*clm%forc_pbot*1.65*rs, 0._r4 )
     enddo

#if (defined DGVM)
     clm%annpsn = clm%annpsn + psn
#endif

     ! rs, rb:  s m**2 / umol -> s/m 
     rs = min(rsmax0, rs*cf)
     rb = rb*cf 

  endif

end subroutine Stomata

!================================

subroutine SensibleHCond (ra,   rb,   rd,   wta,  wtl,  &
                          wtg,  wta0, wtl0, wtg0, wtal, &
                          wtga, wtgl, clm   ) 

!-----------------------------------------------------------------------
!
!  CLMCLMCLMCLMCLMCLMCLMCLMCL  A community developed and sponsored, freely   
!  L                        M  available land surface process model.  
!  M --COMMON LAND MODEL--  C  
!  C                        L  CLM WEB INFO: http://clm.gsfc.nasa.gov
!  LMCLMCLMCLMCLMCLMCLMCLMCLM  CLM ListServ/Mailing List: 
!
!-----------------------------------------------------------------------
! Purpose:
! Provides dimensional and non-dimensional sensible heat
! conductances for canopy and soil flux calculations.
!
! Method:
!
! Author:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision 
!
!-----------------------------------------------------------------------
! $Id: SensibleHCond.F90,v 1.6 2004/11/24 22:56:34 jim Exp $
!-----------------------------------------------------------------------

  use clmtype
  implicit none

!----Arguments----------------------------------------------------------

  type (clm1d), intent(inout)  :: clm !CLM 1-D Module

  real(r8), intent(in) :: ra ! aerodynamical resistance [s/m]
  real(r8), intent(in) :: rb ! leaf boundary layer resistance [s/m]
  real(r8), intent(in) :: rd ! thermal resistance between ground and bottom of canopy

  real(r8), intent(out) :: wta  ! heat conduactance for air [m/s]
  real(r8), intent(out) :: wtg  ! heat conduactance for ground [m/s]
  real(r8), intent(out) :: wtl  ! heat conduactance for leaf [m/s]
  real(r8), intent(out) :: wta0 ! normalized heat conductance for air [-]
  real(r8), intent(out) :: wtl0 ! normalized heat conductance for air [-]
  real(r8), intent(out) :: wtg0 ! normalized heat conductance for ground [-]
  real(r8), intent(out) :: wtal ! normalized heat conductance for air and leaf [-]
  real(r8), intent(out) :: wtgl ! normalized heat conductance for leaf and ground [-]
  real(r8), intent(out) :: wtga ! normalized heat conductance for air and ground  [-]

!----Local Variables----------------------------------------------------

  real(r8) wtshi                ! heat resistance for air, ground and leaf [s/m]

!----End Variable List--------------------------------------------------

  wta   = 1./ra                     ! air
  wtl   = (clm%elai+clm%esai)/rb    ! leaf
  wtg   = 1./rd                     ! ground
  wtshi = 1./(wta+wtl+wtg)

  wtl0  = wtl*wtshi         ! leaf
  wtg0  = wtg*wtshi         ! ground
  wta0  = wta*wtshi         ! air

  wtgl  = wtl0+wtg0         ! ground + leaf
  wtga  = wta0+wtg0         ! ground + air
  wtal  = wta0+wtl0         ! air + leaf

end subroutine SensibleHCond

!==================================

subroutine LatentHCond (raw,   rbw,   rdw,   rpp,   wtaq,  &
                        wtlq,  wtgq,  wtaq0, wtlq0, wtgq0, &
                        wtalq, wtgaq, wtglq, clm    )

!-----------------------------------------------------------------------
!
!  CLMCLMCLMCLMCLMCLMCLMCLMCL  A community developed and sponsored, freely   
!  L                        M  available land surface process model.  
!  M --COMMON LAND MODEL--  C  
!  C                        L  CLM WEB INFO: http://clm.gsfc.nasa.gov
!  LMCLMCLMCLMCLMCLMCLMCLMCLM  CLM ListServ/Mailing List: 
!
!-----------------------------------------------------------------------
! Purpose:
! Provides dimensional and non-dimensional latent heat 
! conductances for canopy and soil flux calculations.  Latent fluxes 
! differs from the sensible heat flux due to stomatal resistance.
!
! Method:
!
! Author:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision 
!
!-----------------------------------------------------------------------
! $Id: LatentHCond.F90,v 1.6 2004/11/24 22:56:30 jim Exp $
!-----------------------------------------------------------------------

  use clmtype
  implicit none

!----Arguments----------------------------------------------------------

  type (clm1d), intent(inout) :: clm !CLM 1-D Module

  real(r8), intent(in) :: raw    ! aerodynamical resistance [s/m]
  real(r8), intent(in) :: rbw    ! leaf boundary layer resistance [s/m]
  real(r8), intent(in) :: rdw    ! latent heat resistance between ground and bottom 
                                 ! of canopy
  real(r8), intent(in) :: rpp    ! fraction of potential evaporation from leaf [-]

  real(r8), intent(out) :: wtaq  ! latent heat conduactance for air [m/s]
  real(r8), intent(out) :: wtlq  ! latent heat conduactance for leaf [m/s]
  real(r8), intent(out) :: wtgq  ! latent heat conduactance for ground [m/s]
  real(r8), intent(out) :: wtaq0 ! normalized latent heat conduactance for air [-]
  real(r8), intent(out) :: wtlq0 ! normalized latent heat conduactance for leaf [-]
  real(r8), intent(out) :: wtgq0 ! normalized heat conduactance for ground [-]
  real(r8), intent(out) :: wtalq ! normalized latent heat cond. for air and leaf [-]
  real(r8), intent(out) :: wtglq ! normalized latent heat cond. for leaf and ground [-]
  real(r8), intent(out) :: wtgaq ! normalized latent heat cond. for air and ground [-]

!----Local Variables----------------------------------------------------

  real(r8) wtsqi                 ! latent heat resistance for air, grd and leaf [-]

!----End Variable List--------------------------------------------------

  wtaq  = clm%frac_veg_nosno/raw                                ! air
  wtlq  = clm%frac_veg_nosno*(clm%elai+clm%esai)/rbw * rpp      ! leaf
  wtgq  = clm%frac_veg_nosno/rdw                                ! ground
  wtsqi = 1./(wtaq+wtlq+wtgq)

  wtgq0 = wtgq*wtsqi                    ! ground
  wtlq0 = wtlq*wtsqi                    ! leaf
  wtaq0 = wtaq*wtsqi                    ! air

  wtglq = wtgq0+wtlq0                   ! ground + leaf
  wtgaq = wtaq0+wtgq0                   ! air + ground
  wtalq = wtaq0+wtlq0                   ! air + leaf

end subroutine LatentHCond

end module module_biogeophysics
