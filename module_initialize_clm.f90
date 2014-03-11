module initialize

use clmtype
use shr_const_mod, only : r8 => shr_kind_r8
contains

!-------------------------
subroutine initialize_clm( clm )

use clm_varpar   , only : nlevsoi, nlevlak, numrad
use clm_varsur   , only : zsoi, dzsoi, zisoi, sand, clay
use pft_varcon   , only : ncorn, nwheat, roota_par, rootb_par,  &
                          z0mr, displar, dleaf, rhol, rhos, taul, taus, xl, &
                          qe25, vcmx25, mp, c3psn
use shr_const_mod, only : SHR_CONST_PI
implicit none
type (clm1d) :: clm
integer :: j, ivt, ib
integer :: isp  ! soil type, 19 kinds
real(r8) :: scalez = 0.025  !Soil layer thickness discretization (m) 
real(r8) :: hkdepth = 0.5   !Length scale for Ksat decrease (m)
real(r8) :: bd              !bulk density of dry soil material [kg/m^3]
real(r8) :: tkm             !mineral conductivity
real(r8) :: xksat           !maximum hydraulic conductivity of soil [mm/s]

!--------------------
  clm%fsun                   = 0.7
  clm%elai                   = 2.0
  clm%esai                   = 1.0
  clm%tlai                   = 2.0
  clm%tsai                   = 1.0

  clm%htop                   = 20.  
  clm%hbot                   = 1.

  clm%itypwat                = 1
  clm%itypveg                = 10
  clm%isoicol                = 1
!..........................................

  clm%fabd(1)                = 0.3 
  clm%fabd(2)                = 0.1 
  clm%fabi(1)                = 0.3  
  clm%fabi(2)                = 0.1  

  clm%ftdd(1)                = 0.6
  clm%ftdd(2)                = 0.4
  clm%ftid(1)                = 0.1
  clm%ftid(2)                = 0.1
  clm%ftii(1)                = 0.1
  clm%ftii(2)                = 0.1

  clm%albgrd(1)              = 0.4
  clm%albgrd(2)              = 0.4
  clm%albgri(1)              = 0.4
  clm%albgri(2)              = 0.4

!..........................................
! Define layer structure for soil
!..........................................

  do j = 1, nlevsoi
     zsoi(j) = scalez*(exp(0.5*(j-0.5))-1.)    !node depths
  enddo

  dzsoi(1) = 0.5*(zsoi(1)+zsoi(2))             !thickness b/n two interfaces
  do j = 2,nlevsoi-1
     dzsoi(j)= 0.5*(zsoi(j+1)-zsoi(j-1))
  enddo
  dzsoi(nlevsoi) = zsoi(nlevsoi)-zsoi(nlevsoi-1)

  zisoi(0) = 0.
  do j = 1, nlevsoi-1
     zisoi(j) = 0.5*(zsoi(j)+zsoi(j+1))         !interface depths
  enddo
  zisoi(nlevsoi) = zsoi(nlevsoi) + 0.5*dzsoi(nlevsoi)
 
  clm%z(1:nlevsoi)  = zsoi(1:nlevsoi)
  clm%dz(1:nlevsoi) = dzsoi(1:nlevsoi)
  clm%zi(0:nlevsoi) = zisoi(0:nlevsoi)


! --------------------------------------------------------------------
! Initialize root fraction (computing from surface, d is depth in meter):
! Y = 1 -1/2 (exp(-ad)+exp(-bd) under the constraint that
! Y(d =0.1m) = 1-beta^(10 cm) and Y(d=d_obs)=0.99 with beta & d_obs
! given in Zeng et al. (1998).
! --------------------------------------------------------------------
   ivt = clm%itypveg
   do j = 1, nlevsoi-1
       clm%rootfr(j) = .5*( exp(-roota_par(ivt)*clm%zi(j-1))  &
                       + exp(-rootb_par(ivt)*clm%zi(j-1))  &
                       - exp(-roota_par(ivt)*clm%zi(j  ))  &
                       - exp(-rootb_par(ivt)*clm%zi(j  )) )
   end do
   clm%rootfr(nlevsoi) = .5*( exp(-roota_par(ivt)*clm%zi(nlevsoi-1))  &
                         + exp(-rootb_par(ivt)*clm%zi(nlevsoi-1)) )

! --------------------------------------------------------------------
! Initialize soil thermal and hydraulic properties 
! --------------------------------------------------------------------
  isp = 8
  do j = 1, nlevsoi
      clm%bsw(j)    = 2.91 + 0.159*clay(isp)
      clm%watsat(j) = 0.489 - 0.00126*sand(isp)

      xksat         = 0.0070556 *( 10.**(-0.884+0.0153*sand(isp)) ) ! mm/s

      clm%hksat(j)  = xksat * exp(-clm%zi(j)/hkdepth)
      clm%sucsat(j) = 10. * ( 10.**(1.88-0.0131*sand(isp)) )
      if(sand(isp).eq.0.and.clay(isp).eq.0) then
         sand(isp) = 0.001
         clay(isp) = 0.001
      endif
      tkm              = (8.80*sand(isp)+2.92*clay(isp))/(sand(isp)+clay(isp)) ! W/(m K)
      bd               = (1.-clm%watsat(j))*2.7e3
      clm%tkmg(j)   = tkm ** (1.- clm%watsat(j))
      clm%tksatu(j) = clm%tkmg(j)*0.57**clm%watsat(j)
      clm%tkdry(j)  = (0.135*bd + 64.7) / (2.7e3 - 0.947*bd)
      clm%csol(j)   = (2.128*sand(isp)+2.385*clay(isp))/ (sand(isp)+clay(isp))*1.e6  ! J/(m3 K)
  end do

! --------------------------------------------------------------------
! Initialize clm derived type components from pft_varcon 
! --------------------------------------------------------------------
     ivt = clm%itypveg
     clm%z0mr    = z0mr(ivt)
     clm%displar = displar(ivt)
     clm%dleaf  = dleaf(ivt)
     clm%xl     = xl(ivt)
     do ib = 1,numrad
        clm%rhol(ib) = rhol(ivt,ib)
        clm%rhos(ib) = rhos(ivt,ib)
        clm%taul(ib) = taul(ivt,ib)
        clm%taus(ib) = taus(ivt,ib)
     end do
     clm%qe25   = qe25(ivt)      ! quantum efficiency at 25c (umol co2 / umol photon)
     clm%vcmx25 = vcmx25(ivt)    ! maximum rate of carboxylation at 25c (umol co2/m**2/s)
     clm%mp     = mp(ivt)        ! slope for conductance-to-photosynthesis relationship
     clm%c3psn  = c3psn(ivt)     ! photosynthetic pathway: 0. = c4, 1. = c3



end subroutine initialize_clm

!...................................
!
!...................................

subroutine iniTimeVar(clm)
  use clm_varpar,   only : nlevsno, nlevsoi
  use clm_varcon  , only : bdsno, istice, istwet, istsoil, denice, denh2o, tfrz, spval, doalb
  use  module_biogeophysics
implicit none
type (clm1d) :: clm
integer  :: i,k 
real :: clm_it, clm_ism

! ========================================================================
! Set snow water 
! ========================================================================
  clm%h2ocan = 0.
  if (clm%itypwat == istice) then
     clm%h2osno = 1000.
  else
     clm%h2osno = 200.
  endif
  clm%snowdp  = clm%h2osno/bdsno
  clm%snowage = 0.

! ========================================================================
! Set snow layer number, depth and thickiness 
! ========================================================================
  call snowdp2lev (clm)

! ========================================================================
! Set snow/soil temperature
! ========================================================================
     clm_it = 288.
     clm%t_soisno(-nlevsno+1:nlevsoi) = 0
     clm%t_veg = clm_it
     if (.not. clm%lakpoi) then  !not lake
        clm%t_soisno(-nlevsno+1:0) = spval
        if (clm%snl < 0) then    !snow layer temperatures
           do i = clm%snl+1, 0
              if (clm_it  < 273.16) then
                 clm%t_soisno(i) = clm_it
              else
                 clm%t_soisno(i) = 273.16 - 1.
              endif
           enddo
        endif
        do i = 1, nlevsoi
           if (clm%itypwat == istice) then
              clm%t_soisno(i) = clm_it
           else if (clm%itypwat == istwet) then
              clm%t_soisno(i) = clm_it
           else
              clm%t_soisno(i) = clm_it
           endif
        enddo
        clm%t_grnd = clm%t_soisno(clm%snl+1)
     else                           !lake
        clm%t_grnd = clm_it
        clm%t_lake = clm_it
     endif

! ========================================================================
! Set snow/soil ice and liquid mass
! ========================================================================
     clm%h2osoi_vol(         1:nlevsoi) = spval
     clm%h2osoi_liq(-nlevsno+1:nlevsoi) = spval
     clm%h2osoi_ice(-nlevsno+1:nlevsoi) = spval
     clm_ism = 5.
     
     if (.not. clm%lakpoi) then  !not lake
        if (clm%snl < 0) then    !snow 
           do i = clm%snl+1, 0
              clm%h2osoi_ice(i) = clm%dz(i)*250.
              clm%h2osoi_liq(i) = 0.
           enddo
        endif
        do i = 1, nlevsoi           !soil layers
           if (clm%t_soisno(i) <= tfrz) then       
              clm%h2osoi_ice(i) = clm%dz(i)*clm_ism*clm%watsat(i)*denice
              clm%h2osoi_liq(i) = 0.
              if (clm%itypwat==istwet .or. clm%itypwat==istice) & 
                   clm%h2osoi_ice(i)=clm%dz(i)*denice
           else
              if (clm%itypwat == istsoil) then
                 clm%h2osoi_liq(i) = clm%dz(i)*clm_ism*clm%watsat(i)*denh2o
                 clm%h2osoi_ice(i) = 0.
              elseif (clm%itypwat==istwet .or. clm%itypwat==istice) then 
                 clm%h2osoi_liq(i)=clm%dz(i)*denh2o
                 clm%h2osoi_ice(i) = 0.                
              endif
           endif
        enddo
        
        do i = 1,nlevsoi
           if (clm%itypwat == istsoil) then
              clm%h2osoi_vol(i) = 0.5_r8
              clm%h2osoi_vol(i) = clm%h2osoi_liq(i)/(clm%dz(i)*denh2o) + clm%h2osoi_ice(i)/(clm%dz(i)*denice)
           else
              clm%h2osoi_vol(i) = 1.0_r8
           endif
           clm%h2osoi_vol(i) = min(clm%h2osoi_vol(i),clm%watsat(i))
        end do
     endif

! ========================================================================
! Remaining variables are initialized by calls to ecosystem dynamics and
! albedo subroutines. 
! Note: elai, esai, frac_veg_nosno are computed in Ecosysdyn and needed
! by Fwet and SurfaceAlbedo
! Note: fwet is needed in routine clm_twostream (called by clm_surfalb)
! ========================================================================

     clm%frac_sno = clm%snowdp/(0.1 + clm%snowdp)  
     clm%frac_veg_nosno = clm%frac_veg_nosno_alb
     call Fwet(clm)
     call updatelai(clm)
     call SurfaceAlbedo (clm)


end subroutine iniTimeVar

!===============================

subroutine snowdp2lev(clm)

  use clm_varpar, only : nlevsoi, nlevsno, nlevlak
  implicit none

  type (clm1d) :: clm
! ------------------- local variables -----------------------------
  integer i,k    !indices
! -----------------------------------------------------------------

     clm%dz(-nlevsno+1:0) = 1.e36 
     clm%z (-nlevsno+1:0) = 1.e36 
     clm%zi(-nlevsno:-1)  = 1.e36 
     if (.not. clm%lakpoi) then  !not lake
        if (clm%snowdp < 0.01) then
           clm%snl = 0
           clm%dz(-nlevsno+1:0) = 0.
           clm%z (-nlevsno+1:0) = 0.
           clm%zi(-nlevsno+0:0) = 0.
        else
           if ((clm%snowdp >= 0.01) .AND. (clm%snowdp <= 0.03)) then
              clm%snl = -1
              clm%dz(0)  = clm%snowdp
           else if ((clm%snowdp > 0.03) .AND. (clm%snowdp <= 0.04)) then
              clm%snl = -2
              clm%dz(-1) = clm%snowdp/2.
              clm%dz( 0) = clm%dz(-1)
           else if ((clm%snowdp > 0.04) .AND. (clm%snowdp <= 0.07)) then
              clm%snl = -2
              clm%dz(-1) = 0.02
              clm%dz( 0) = clm%snowdp - clm%dz(-1)
           else if ((clm%snowdp > 0.07) .AND. (clm%snowdp <= 0.12)) then
              clm%snl = -3
              clm%dz(-2) = 0.02
              clm%dz(-1) = (clm%snowdp - 0.02)/2.
              clm%dz( 0) = clm%dz(-1)
           else if ((clm%snowdp > 0.12) .AND. (clm%snowdp <= 0.18)) then
              clm%snl = -3
              clm%dz(-2) = 0.02
              clm%dz(-1) = 0.05
              clm%dz( 0) = clm%snowdp - clm%dz(-2) - clm%dz(-1)
           else if ((clm%snowdp > 0.18) .AND. (clm%snowdp <= 0.29)) then
              clm%snl = -4
              clm%dz(-3) = 0.02
              clm%dz(-2) = 0.05
              clm%dz(-1) = (clm%snowdp - &
                                            clm%dz(-3) - &
                                            clm%dz(-2))/2.
              clm%dz( 0) = clm%dz(-1)
           else if ((clm%snowdp > 0.29) .and. &
                    (clm%snowdp <= 0.41)) then
              clm%snl = -4
              clm%dz(-3) = 0.02
              clm%dz(-2) = 0.05
              clm%dz(-1) = 0.11
              clm%dz( 0) = clm%snowdp - &
                                           clm%dz(-3) - &
                                           clm%dz(-2) - &
                                           clm%dz(-1)
           else if ((clm%snowdp > 0.41) .and. &
                    (clm%snowdp <= 0.64)) then
              clm%snl = -5
              clm%dz(-4) = 0.02
              clm%dz(-3) = 0.05
              clm%dz(-2) = 0.11
              clm%dz(-1) = (clm%snowdp - &
                                            clm%dz(-4) - &
                                            clm%dz(-3) - &
                                            clm%dz(-2))/2.
              clm%dz( 0) = clm%dz(-1)
           else if (clm%snowdp > 0.64) then 
              clm%snl = -5
              clm%dz(-4) = 0.02
              clm%dz(-3) = 0.05
              clm%dz(-2) = 0.11
              clm%dz(-1) = 0.23
              clm%dz( 0) = clm%snowdp - &
                                           clm%dz(-4) - &
                                           clm%dz(-3) - &
                                           clm%dz(-2) - &
                                           clm%dz(-1)
           endif
           do i = 0, clm%snl+1, -1
              clm%z(i)    = clm%zi(i) - 0.5*clm%dz(i)
              clm%zi(i-1) = clm%zi(i) - clm%dz(i)
           enddo
        endif 
     else   !lake points
        clm%snl = 0
        clm%dz(-nlevsno+1:0) = 0.
        clm%z (-nlevsno+1:0) = 0.
        clm%zi(-nlevsno+0:0) = 0.
     endif

end subroutine snowdp2lev

end module initialize
