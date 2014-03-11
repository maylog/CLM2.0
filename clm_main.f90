use module_biogeophysics
use clmtype
use initialize
use clm_varpar, only : nlevsoi
use shr_const_mod, only : r8 => shr_kind_r8

implicit none
  real(r8) :: T, Pa, es, esdT, qs, qsdT 
  type (clm1d) :: clm

  real(r8) eg      ! water vapor pressure at temperature T [pa]
  real(r8) qsatg   ! saturated humidity [kg/kg]
  real(r8) degdT   ! d(eg)/dT
  real(r8) qsatgdT ! d(qsatg)/dT

  real(r8) :: cgrnd
  real(r8) :: cgrndl
  real(r8) :: cgrnds
  real(r8) :: tg
  real(r8) :: emg
  real(r8) :: htvp
  real(r8) :: dlrad
  real(r8) :: ulrad
  real(r8) :: tssbef(-5:10)

  integer :: nstep, j
!-----------------------
  
  clm%forc_rain              = 10.0
  clm%forc_snow              = 10.0
  clm%forc_pbot              = 101325.
  clm%forc_t                 = 280.
  clm%forc_solad(1)          = 0.7*600
  clm%forc_solad(2)          = 0.3*600
  clm%forc_solai(1)          = 0.7*300
  clm%forc_solai(2)          = 0.3*300
  clm%forc_cosz              = 0.5
  clm%dtime                  = 1.
  clm%forc_hgt               = 10.
  clm%forc_t                 = 296.
  clm%forc_th                = 298.
  clm%forc_u                 = 6.
  clm%forc_v                 = 2.
  clm%forc_rho               = 1.29

  call initialize_clm( clm )
  call iniTimeVar(clm)
!write(*,*) clm%fsun

!========================================
!========================================
   clm%nstep = nstep
   clm%h2osno_old = clm%h2osno  ! snow mass at previous time step
   clm%frac_veg_nosno = clm%frac_veg_nosno_alb
   if (clm%h2osno > 1000.) then
      clm%do_capsnow = .true.
   else
      clm%do_capsnow = .false.
   endif


   if (.not. clm%lakpoi) then

      do j = clm%snl+1, 0       ! ice fraction of snow at previous time step
         clm%frac_iceold(j) = clm%h2osoi_ice(j)/(clm%h2osoi_liq(j)+clm%h2osoi_ice(j))
      enddo

! Determine beginning water balance (water balance at previous time step)
     clm%begwb = clm%h2ocan + clm%h2osno
     do j = 1, nlevsoi
        clm%begwb = clm%begwb +  clm%h2osoi_ice(j) + clm%h2osoi_liq(j)
     enddo

! Determine canopy interception and precipitation onto ground surface.
! Determine the fraction of foliage covered by water and the fraction
! of foliage that is dry and transpiring. Initialize snow layer if the
! snow accumulation exceeds 10 mm.
        call Hydrology1(clm)
!        call updatelai(clm)
!        write(*,*) 'ini',clm%elai, clm%esai

! Determine leaf temperature and surface fluxes based on ground
! temperature from previous time step.
        call Biogeophysics1(clm,cgrnd,cgrndl,cgrnds,tg,emg,htvp, dlrad,ulrad,tssbef)


   endif






!-----------------------


!  tg = 290.
!  clm%forc_pbot = 101325 
!  call QSatclm( tg, clm%forc_pbot, eg, degdT, qsatg, qsatgdT )
!  write(*,*) 'qsatclm', eg, degdT, qsatg, qsatgdT

!  call SurfaceRadiation (clm)
!  write(*,*) 'sur_rad',clm%sabv, clm%sabg

!  call CanopyFluxes (z0mv,   z0hv,  z0qv,  thm,   clm%forc_th, &
!                        thv,    tg,    qg,    dqgdT, htvp,        &
!                        emv,    emg,   dlrad, ulrad, cgrnds,      &
!                        cgrndl, cgrnd, clm    )



!write(*,*) 'in clm'

END
