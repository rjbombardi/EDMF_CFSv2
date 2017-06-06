      !**************************************************************
      !***               :::    DESCRIPTION    :::                ***
      !***    Defines a new flag variable TRIGGERFLG used to tell ***
      !*** the model to use the convection test below.            ***
      !***    Call a new test for convection based on the buoyant ***
      !*** condensation level (BCL).                              ***
      !**************************************************************
!---------------------------------------------------------------------------------
! Purpose:
!
! Calls a new test for convection based on the buoyant condensation level (BCL)
! This method assumes moist deep convection is triggered when saturation occurs
! at the top of the mixed layer.  Mixed layer is defined by the neutral buoyancy 
! level of a surface parcel (usually 2-meter parcels).
! Reference: Tawfik and Dirmeyer 2014 GRL: A processed based framework for 
!            quantifying the atmospheric background state of surface triggered
!            convection
!
! Output variables:
! TDEF = potential temperature deficit (K)
! TBM  = buoyant mixing potential temperature (K)
! BCLH = buoyant condensation level height (m)
! BCLP = buoyant condensation pressure level (Pa)
!
! Author and Revision History: 
! A.B. Tawfik on Aug 2014
!
!---------------------------------------------------------------------------------
      module BCLTRIGGER
      !    machine restores: machine dependent constants
      USE MACHINE , ONLY : kind_phys
      !    funcphys FPVS elementally computes saturation vapor pressure
      USE FUNCPHYS , ONLY : fpvs
      !    physcons restores: some of the most used math and phycis
      !    constants
      USE PHYSCONS, ONLY : grav => con_g, CP => con_CP, RD => con_rd
     &,             RV => con_rv, EPS => con_eps, EPSM1 => con_epsm1

      implicit none

      ! Private
      private
      ! Public Variables
      logical, public :: TRIGGERFLG = .false.
      logical, public :: hybedmf    = .false.
      logical, public :: dspheat    = .false.

      public bcl_mixed_layer

      contains
!------------------------------------------------------------------------------
!      subroutine bcl_mixed_layer(nlev,tmp_in,press_in,qhum_in 
!     &,hgt_in, t2m, psfc, q2m, h2m, TBM, BCLH, BCLP, TDEF)
      subroutine bcl_mixed_layer(nlev,tmp_in,press_in,qhum_in
     &,                           TDEF, BCLP)

!------------------------------------------------------------------------------
      !**************************************************************
      !***               :::    DESCRIPTION    :::                ***
      !***    Call a new test for convection based on the buoyant ***
      !*** condensation level (BCL). This method assumes moist    ***
      !*** deep convection is triggered when saturation occurs at ***
      !*** the top of the mixed layer.  Mixed layer is defined by ***
      !*** the neutral buoyancy level of a surface parcel         ***
      !*** (usually 2-meter parcels).                             ***
      !***                                                        ***   
      !***               :::       NOTES       :::                ***
      !***    Vertical profiles are assumed to be bottom-up       ***
      !*** (surface to tropopause). Assumes first model level is  ***
      !*** either 2-m or surface layer.                           ***
      !***                                                        ***
      !***               :::     REFERENCE     :::                ***
      !***    #If you use this routine please reference the       ***
      !*** publication: Tawfik and Dirmeyer 2013 GRL: A processed ***
      !*** based framework for quantifying the atmospheric        ***
      !*** background state of surface triggered convection       ***
      !***                                                        ***
      !***               :::  INPUT VARIABLES  :::                ***
      !*** press  =  pressure levels in [Pascals]                 ***
      !*** tmp_in  =  temperature profile [Kelvin]                 ***
      !*** pot    =  potential temperature profile Kelvin]        ***
      !*** qhum_in   =  specific humidity profile in [kg H2O/kg air] ***
      !*** hgt    =  geopotential height levels  in [m]           ***
      !*** n      =  number of vertical levels                    ***
      !***                                                        ***
      !***               :::  OUTPUT VARIABLES :::                ***
      !*** qdef =  saturation deficit at height of neutral        ***
      !***         buoyancy given some pot. temp. [kg/kg]         ***
      !*** bclh =  height of the buoyant condensatin level [m]    ***
      !*** bclp =  pressure at the buoyant condensatin level [Pa] ***
      !***                                                        ***
      !**************************************************************
      !**************************************************************
      !***                                                        ***
      !***    This routine computes q (mixing ratio or specific   ***
      !*** humidity) from rh and temp.                            ***
      !***                                                        ***
      !**************************************************************
      !**************************************************************
      !***               Input/Output arguments:                  ***
      !**************************************************************
      implicit none

      integer              , intent(in   ) :: nlev
      real(kind=kind_phys) , intent(in   ) :: press_in(nlev)
      real(kind=kind_phys) , intent(in   ) :: tmp_in(nlev)
      real(kind=kind_phys) , intent(in   ) :: qhum_in(nlev)
      real(kind=kind_phys) , intent(inout) :: TDEF
      real(kind=kind_phys) , intent(inout) :: BCLP
      real(kind=kind_phys)                 :: TBM

      !*************************************************************
      !***              LOCAL VARIABLES AND ARRAYS               ***
      !*************************************************************

      real(kind=kind_phys) :: pbar(nlev)
      real(kind=kind_phys) :: tbar(nlev)
      real(kind=kind_phys) :: qbar(nlev)
      real(kind=kind_phys) :: qdef(nlev)
      real(kind=kind_phys) :: qmix(nlev)
      real(kind=kind_phys) :: qsat(nlev)
      real(kind=kind_phys) :: logp(nlev)
      real(kind=kind_phys) :: rhoh(nlev)
      real(kind=kind_phys) :: logp0
      real(kind=kind_phys) :: logp1
      real(kind=kind_phys) :: logp2
      real(kind=kind_phys) :: dpress
      real(kind=kind_phys) :: pot2m
      real(kind=kind_phys) :: qmix1
      real(kind=kind_phys) :: qmix2
      real(kind=kind_phys) :: rhoh1
      real(kind=kind_phys) :: p_up
      real(kind=kind_phys) :: t_up
      real(kind=kind_phys) :: q_up
      real(kind=kind_phys) :: p_lo
      real(kind=kind_phys) :: t_lo
      real(kind=kind_phys) :: q_lo

      integer                :: zz
      real(kind=kind_phys)   :: es, prime

      logical   :: first_negative
      real(kind=kind_phys), parameter :: RDRV = 1./(RD/RV)
     &, RDCP=RD/CP, grav1=1./grav, p_ref=1e5 

!--- Allocating variables

!-----------------------------------------------------------------------------
      first_negative=.false.
      do zz = 1,nlev-1

        logp0 = log(press_in(zz))
        logp1 = log(press_in(zz+1))
        logp2 = log(press_in(zz+1)*press_in(zz))

      !-- Calculate middle layer temperature
        tbar(zz) = ( (tmp_in(zz+1)*logp1     ) +
     &                tmp_in(zz)*logp0    ) /
     &                logp2

      !-- Calculate pressure difference between layers
        dpress  =  press_in(zz) - press_in(zz+1)

      !-- Calculate middle layer specific humidity (kg/kg)
        qbar(zz) = ( (qhum_in(zz+1)*logp1     ) +
     &                qhum_in(zz)*logp0    ) /
     &                logp2

      !-- Calculate middle layer pressure (Pa)
        pbar(zz) = ( (press_in(zz+1)*logp1  ) +
     &                press_in(zz)*logp0 ) /
     &                logp2

      !-- Calculate log of pressure
        logp(zz) =  log(pbar(zz))

      !-- Calculate calculate mixed moisture (kg/m2) and column density (kg/m2)
        qmix1 = qbar(zz)*dpress/grav
        rhoh1 = dpress/grav

      !-- Calculate calculate mixed humidity (kg/kg)
        if(zz .eq. 1)then
          qmix(zz) = qmix1
          rhoh(zz) = rhoh1
        else
          qmix(zz) = qmix(zz-1) + qmix1
          rhoh(zz) = rhoh(zz-1) + rhoh1
        endif

      !-- Calculate calculate mixed humidity (kg/kg)
        ES = 0.01 * fpvs(tbar(zz))      ! fpvs is in Pa
        PRIME = pbar(zz)/100. + EPSM1 * ES
        qsat(zz) = EPS * ES / PRIME

      !-- Calucalte the saturation specific humidity (kg/kg)
        qmix2 = qmix(zz) / rhoh(zz)
        qdef(zz) = qsat(zz) - qmix2

        if(pbar(zz) .le. 20000.) qdef(zz) = -1.

      !***********************************************************
      !***   Calculate slope of each variable to find the      ***
      !***   y-intercept for each variable;                    ***
      !***   Meaning locate the two data points surrounding    ***
      !***   the sign change in qdef and linearly interpolate  ***
      !***   to find the "zero point"                          ***
      !***********************************************************

      !----- Get the higher (from ground) level point in the qdef sign transition
      !----- Find the point where the sign first turns negative from the ground up


        if(qdef(zz) .le. 0.0)then
            !--- find where index equals the first level of negative
          p_up = logp(zz)        !--- Pressure (Pa)
          t_up = tbar(zz)        !--- Temperature (K)
          q_up = qdef(zz)        !--- Specific Humidity Deficit (kg/kg)

          if(zz .eq. 1)then
             p_lo = press_in(1)  !--- Pressure (Pa)
             t_lo = tmp_in(1)    !--- Temperature (K)
             q_lo = qhum_in(1)   !--- Specific Humidity Deficit (kg/kg)
          else
             p_lo = logp(zz-1)   !--- Pressure (Pa)
             t_lo = tbar(zz-1)   !--- Temperature (K)
             q_lo = qdef(zz-1)   !--- Specific Humidity Deficit (kg/kg)
          endif
      !----- Get the lower (closest to ground) level before qdef sign transition
          first_negative=.true.
        endif


      !--- Calculate output variables; BCL height, BCL pressure, 
      !--- Buoyant Mixing Potential Temp, and Potential Temperature Deficit
      !*** Note ***
      !--- Check to see if first level is satured; if so then set BCL and other
      !--- output variables to near-surface quantities and return
      !--- !!!!  This is the Virtual Potential Temperature (K) !!!
        pot2m  =  (tbar(1)*(1.+0.61*qbar(1))) * 
     &            ((p_ref/pbar(1))**(RDCP))
        if( qdef(1).le.0. ) then
          BCLP  =  pbar(1)
          TBM   =  pot2m
          TDEF  =  0.
        else
          BCLP  =  exp( p_up - ((p_up-p_lo) /
     &                    (q_up-q_lo))*q_up )
          TBM   = (t_up - ((t_up-t_lo) /
     &                (q_up-q_lo))*q_up ) * 
     &                ((p_ref/BCLP)**(RDCP))
          TDEF  =     TBM  - pot2m
        endif

        if(first_negative) return

      enddo

!--- Deallocating variables

      return

      end subroutine bcl_mixed_layer


      end module bcltrigger
