! ===================================================================== !
!  description:                                                         !
!                                                                       !
!     gbphys is main program to invoke model physics (except radiation) !
!     at model time steps                                               !
!                                                                       !
!  usage:                                                               !
!                                                                       !
!    call gbphys                                                        !
!       inputs:                                                         !
!         ( im, ix, levs, lsoil, lsm, ntrac, ncld, ntoz, ntcw,          !
!           nmtvr, nrcm, ko3, lonf, latg, jcap, num_p3d, num_p2d,       !
!           kdt, lat, me, pl_coeff, nlons, ncw, flgmin, crtrh,          !
!           ccwf, ctei_rm, clstp, dtp, dtf, fhour, solhr,               !
!           slag, sdec, cdec, sinlat, coslat, pgr, ugrs, vgrs,          !
!           tgrs, qgrs, vvel, prsi, prsl, prslk, prsik, phii, phil,     !
!           xkt2, prdout, poz, rcs2, dpshc, hprime, xlon, xlat,         !
!           slope, shdmin, shdmax, snoalb, tg3, slmsk, vfrac,           !
!           vtype, stype, uustar, oro, coszen, sfcdsw, sfcnsw,          !
!           sfcdlw, tsflw, sfcemis, sfalb, swh, hlw, ras, pre_rad,      !
!--- rjb            
!!          ldiag3d, trans_trac, lssav, lssav_cc, xkzm, flipv,          !
!           ldiag3d, trans_trac, lssav, lssav_cc, xkzm, xkzm_m,         !
!           xkzm_h, xkzm_s, flipv,                                      !
!--- rjb            
!           old_monin, cnvgwd, sashal, newsas, mom4ice, mstrat,         !
!       input/outputs:                                                  !
!           hice, fice, tisfc, tsea, tprcp, cv, cvb, cvt,               !
!           acv, acvb, acvt, cldwrk, phy_f3d, phy_f2d, ep,              !
!           slc, smc, stc, upd_mf, dwn_mf, det_mf,                      !
!           srflag, snwdph, sheleg, sncovr, zorl, canopy,               !
!           ffmm, ffhh, f10m, srunoff, evbsa, evcwa, snohfa,            !
!           transa, sbsnoa, snowca, soilm, tmpmin, tmpmax,              !
!           dugwd, dvgwd, psmean, bengsh, spfhmin, spfhmax,             !
!           dusfc, dvsfc, dtsfc, dqsfc, geshem, gflux,                  !
!           dt3dt, dq3dt, du3dt, dv3dt, runoff,                         !
!           dlwsfc, ulwsfc, dlwsfc_cc, ulwsfc_cc, swsfc_cc,             !
!           dusfc_cc, dvsfc_cc, dtsfc_cc, dqsfc_cc, precr_cc,           !
!       outputs:                                                        !
!           gt0, gq0, gu0, gv0, t2m, q2m, u10m, v10m,                   !
!           zlvl, psurf, hpbl, pwat, t1, q1, u1, v1,                    !
!           chh, cmm, dlwsfci, ulwsfci, dswsfci, uswsfci,               !
!           dtsfci, dqsfci, gfluxi, epi,                                !
!           xmu_cc, dlw_cc, dsw_cc, snw_cc, lprec_cc )                  !
!                                                                       !
!  subprograms called:                                                  !
!                                                                       !
!     get_prs,  dcyc2t2_pre_rad (testing),    dcyc2t3,  sfc_diff,       !
!     sfc_ocean,sfc_drv,  sfc_land, sfc_sice, sfc_diag, moninp1,        !
!     moninp,   moninq1,  moninq,   gwdps,    ozphys,   get_phi,        !
!     sascnv,   sascnvn,  rascnv,   gwdc,     shalcvt3, shalcv,         !
!     shalcnv,  cnvc90,   lrgscl,   gsmdrive, gscond,   precpd,         !
!     progt2                                                            !
!                                                                       !
!                                                                       !
!  program history log:                                                 !
!           19xx  - ncep mrf/gfs                                        !
!      nov  2002  - s. moorthi  modify and restructure                  !
!           200x  - h. juang    modify (need description)               !
!      nov  2004  - x. wu       modify sea-ice model                    !
!      may  2005  - s. moorthi  modify and restructure                  !
!           2005  - s. lu       modify (need description)               !
!      oct  2006  - h. wei      modify land sfc mdl including options   !
!                     for noah land sfc mdl and osu land sfc mdl.       !
!           2008  - jun wang    added spfhmax/spfhmin as input/output.  !
!      apr  2008  - y.-t. hou   added lw sfc emissivity var (sfcemis),  !
!                     define the lw sfc dn/up fluxes in two forms: atmos!
!                     and ground. also changed sw sfc net flux direction!
!                     (positive) from ground -> atmos to the direction  !
!                     of atmos -> ground. recode the program and add    !
!                     program documentation block.                      !
!                                                                       !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                       size   !
!     ix, im   - integer, horiz dimention and num of used pts      1    !
!     levs     - integer, vertical layer dimension                 1    !
!     lsoil    - integer, number of soil layers                    1    !
!     lsm      - integer, flag for land surface model scheme       1    !
!                =0: use osu scheme; =1: use noah scheme                !
!     ntrac    - integer, number of tracers                        1    !
!     ncld     - integer, indicate cloud types                     1    !
!     ntoz     - integer, flag for ozone generation schemes        1    !
!     ntcw     - integer, flag for cloud condensate computations   1    !
!     nmtvr    - integer, number of montain variables              1    !
!     nrcm     - integer, number of random clouds                  1    !
!     ko3      - integer, number of layers for ozone data          1    !
!     lonf,latg- integer, number of lon/lat points                 1    !
!     jcap     - integer, number of spectral wave trancation       1    !
!     num_p3d  - integer,                                          1    !
!     num_p2d  - integer,                                          1    !
!     kdt,lat,me-integer, for check print debug control purpose    1    !
!     pl_coeff - integer,                                          1    !
!     nlons    - integer,                                          im   !
!     ncw      - integer,                                          2    !
!     flgmin   - real, min large ice fraction                      2    !
!     crtrh    - real,                                             3    !
!     ccwf     - real,                                             1    !
!     ctei_rm  - real,                                             1    !
!     clstp    - real,                                             1    !
!     dtp,dtf  - real, time interval (second)                      1    !
!     fhour    - real, forecast hour                               1    !
!     solhr    - real, fcst hour at the end of prev time step      1    !
!     slag     - real, equation of time ( radian )                 1    !
!     sdec,cdec- real, sin and cos of the solar declination angle  1    !
!     sinlat   - real, sin of latitude                             im   !
!     coslat   - real, cos of latitude                             im   !
!     pgr      - real, surface pressure                            im   !
!     ugrs,vgrs- real, u/v component of layer wind              ix,levs !
!     tgrs     - real, layer mean temperature ( k )             ix,levs !
!     qgrs     - real, layer mean tracer concentration     ix,levs,ntrac!
!     vvel     - real, layer mean vertical velocity             ix,levs !
!     prsi     - real, pressure at layer interfaces             ix,levs+1
!     prsl     - real, mean layer presure                       ix,levs !
!     prsik    - real,                                          ix,levs+1
!     prslk    - real,                                          ix,levs !
!     phii     - real,                                          ix,levs+1
!     phil     - real,                                          ix,levs !
!     xkt2     - real,                                          ix,nrcm !
!     prdout   - real,                                   ix,ko3,pl_coeff!
!     poz      - real,                                             ko3  !
!     rcs2     - real,                                             im   !
!     dpshc    - real,                                             im   !
!     hprime   - real, orographic std dev                       ix,nmtvr!
!     xlon,xlat- real, longitude and latitude ( radian )           im   !
!     slope    - real, sfc slope type for lsm                      im   !
!     shdmin   - real, min fractional coverage of green veg        im   !
!     shdmax   - real, max fractnl cover of green veg (not used)   im   !
!     snoalb   - real, max snow albedo over land (for deep snow)   im   !
!     tg3      - real, deep soil temperature                       im   !
!     slmsk    - real, sea/land/ice mask (=0/1/2)                  im   !
!     vfrac    - real,                                             im   !
!     vtype    - real, vegetation type                             im   !
!     stype    - real, soil type                                   im   !
!     uustar   - real,                                             im   !
!     oro      - real,                                             im   !
!     coszen   - real, avg cosz over daytime sw radiation interval im   !
!     sfcdsw   - real, total sky sfc downward sw flux ( w/m**2 )   im   !
!     sfcnsw   - real, total sky sfc netsw flx into ground(w/m**2) im   !
!     sfcdlw   - real, total sky sfc downward lw flux ( w/m**2 )   im   !
!     tsflw    - real, sfc air (layer 1) temp over lw interval (k) im   !
!     sfcemis  - real, sfc lw emissivity ( fraction )              im   !
!     sfalb    - real, mean sfc diffused sw albedo                 im   !
!     swh      - real, total sky sw heating rates ( k/s )       ix,levs !
!     hlw      - real, total sky lw heating rates ( k/s )       ix,levs !
!     ras      - logical, flag for ras convection scheme           1    !
!     pre_rad  - logical, flag for testing purpose                 1    !
!     ldiag3d  - logical, flag for 3d diagnostic fields            1    !
!     trans_trac-logical,                                          1    !
!     lssav    - logical, flag controls data store and output      1    !
!     lssav_cc - logical, flag for save data for ocean coupling    1    !
!     xkzm     - logical,                                          1    !
!     flipv    - logical, flag for vertical direction              1    !
!     old_monin- logical, flag for diff monin schemes              1    !
!     cnvgwd   - logical, flag for conv gravity wave drag          1    !
!     sashal   - logical, flag for shallow conv schemes            1    !
!     newsas   - logical, flag for sas conv schemes                1    !
!     mom4ice  - logical, flag controls mom4 sea-ice               1    !
!     mstrat   - logical,                                          1    !
!                                                                       !
!  input/outputs:                                                       !
!     hice     - real, sea-ice thickness                           im   !
!     fice     - real, sea-ice concentration                       im   !
!     tisfc    - real, sea-ice temperature                         im   !
!     tsea     - real, ground surface temperature ( k )            im   !
!     tprcp    - real, total precipitation                         im   !
!     cv,cvb,cvt-real, convective clouds amt, base and top         im   !
!     acv,acvb,acvt                                                     !
!              - real,                                             im   !
!     cldwrk   - real,                                             im   !
!     phy_f3d  - real,                                   ix,levs,num_p3d!
!     phy_f2d  - real,                                        ix,num_p2d!
!     ep       - real,                                             im   !
!     slc      - real, liquid soil moisture                     ix,lsoil!
!     smc      - real, total soil moisture content (fractional) ix,lsoil!
!     stc      - real,                                          ix,lsoil!
!     upd_mf   - real,                                          ix,levs !
!     dwn_mf   - real,                                          ix,levs !
!     det_mf   - real,                                          ix,levs !
!     srflag   - real, snow/rain flag for precipitation            im   !
!     snwdph   - real, snow depth (water equiv) over land          im   !
!     sheleg   - real, snow depth (water equiv)                    im   !
!     sncovr   - real, snow cover over land                        im   !
!     zorl     - real, surface roughness                           im   !
!     canopy   - real, canopy water                                im   !
!     ffmm     - real,                                             im   !
!     ffhh     - real,                                             im   !
!     f10m     - real,                                             im   !
!     srunoff  - real,                                             im   !
!     evbsa    - real, accumulated direct soil evaporation         im   !
!     evcwa    - real, accumulated canopy water evaporation        im   !
!     snohfa   - real, accumulated snow/freezing_rain latent flux  im   !
!     transa   - real, accumulated total plant transpiration (m/s) im   !
!     sbsnoa   - real, accumulated sublimation from snowpack       im   !
!     snowca   - real, accumulated fractional snow cover           im   !
!     soilm    - real, total soil column moisture content (m)      im   !
!     tmpmin   - real, min temperature at 2m height (k)            im   !
!     tmpmax   - real, max temperature at 2m height (k)            im   !
!     dugwd    - real,                                             im   !
!     dvgwd    - real,                                             im   !
!     psmean   - real,                                             im   !
!     bengsh   - real,                                             im   !
!     spfhmin  - real,                                             im   !
!     spfhmax  - real,                                             im   !
!     dusfc    - real,                                             im   !
!     dvsfc    - real,                                             im   !
!     dtsfc    - real,                                             im   !
!     dqsfc    - real,                                             im   !
!     geshem   - real,                                             im   !
!     gflux    - real, groudn conductive heat flux                 im   !
!     runoff   - real, total water runoff                          im   !
!     dt3dt    - real,                                        ix,levs,6 !
!     dq3dt    - real,                                ix,levs,5+pl_coeff!
!     du3dt    - real,                                        ix,levs,4 !
!     dv3dt    - real,                                        ix,levs,4 !
!     dlwsfc   - real, time accumulated sfc dn lw flux ( w/m**2 )  im   !
!     ulwsfc   - real, time accumulated sfc up lw flux ( w/m**2 )  im   !
!     dlwsfc_cc- real, sfc dnwd lw flux (w/m**2) for ocn coupling  im   !
!     ulwsfc_cc- real, sfc upwd lw flux (w/m**2) for ocn coupling  im   !
!     swsfc_cc - real, sfc net sw  flux (w/m**2) for ocn coupling  im   !
!     dusfc_cc - real, sfc u-wind                for ocn coupling  im   !
!     dvsfc_cc - real, sfc v-wind                for ocn coupling  im   !
!     dtsfc_cc - real, sfc temperature  (k)      for ocn coupling  im   !
!     dqsfc_cc - real, sfc pressure              for ocn coupling  im   !
!     precr_cc - real, total precipitation       for ocn coupling  im   !
!                                                                       !
!  outputs:                                                             !
!     gt0      - real,                                          ix,levs !
!     gq0      - real,                                    ix,levs,ntrac !
!     gu0,gv0  - real,                                          ix,levs !
!     t2m,q2m  - real, 2 meter temperature and humidity            im   !
!     u10m,v10m- real, 10 meater u/v wind speed                    im   !
!     zlvl     - real,                                             im   !
!     psurf    - real, surface pressure                            im   !
!     hpbl     - real,                                             im   !
!     pwat     - real, precipitable water                          im   !
!     t1, q1   - real,                                             im   !
!     u1, v1   - real,                                             im   !
!     chh,cmm  - real,                                             im   !
!     dlwsfci  - real, instantaneous sfc dnwd lw flux ( w/m**2 )   im   !
!     ulwsfci  - real, instantaneous sfc upwd lw flux ( w/m**2 )   im   !
!     dswsfci  - real, instantaneous sfc dnwd sw flux ( w/m**2 )   im   !
!     uswsfci  - real, instantaneous sfc upwd sw flux ( w/m**2 )   im   !
!     dtsfci   - real,                                             im   !
!     dqsfci   - real,                                             im   !
!     gfluxi   - real,                                             im   !
!     epi      - real,                                             im   !
!     xmu_cc   - real, cosine of zenith angle at time step         im   !
!     dlw_cc   - real, sfc dnwd lw flux at time step for ocn cpl   im   !
!     dsw_cc   - real, sfc dnwd sw flux at time step for ocn cpl   im   !
!     snw_cc   - real, lower atms snow fall rate for ocn cpl       im   !
!     lprec_cc - real, lower atms rain fall rate for ocn cpl       im   !
!                                                                       !
!  ====================    end of description    =====================  !

      subroutine gbphys                                                 &
!  ---  inputs:
     &    ( im, ix, levs, lsoil, lsm, ntrac, ncld, ntoz, ntcw,          &
     &      nmtvr, nrcm, ko3, lonf, latg, jcap, num_p3d, num_p2d,       &
     &      kdt, lat, me, pl_coeff, nlons, ncw, flgmin, crtrh,          &
     &      ccwf, ctei_rm, clstp, dtp, dtf, fhour, solhr,               &
     &      slag, sdec, cdec, sinlat, coslat, pgr, ugrs, vgrs,          &
     &      tgrs, qgrs, vvel, prsi, prsl, prslk, prsik, phii, phil,     &
     &      xkt2, prdout, poz, rcs2, dpshc, hprime, xlon, xlat,         &
     &      slope, shdmin, shdmax, snoalb, tg3, slmsk, vfrac,           &
     &      vtype, stype, uustar, oro, coszen, sfcdsw, sfcnsw,          &
     &      sfcdlw, tsflw, sfcemis, sfalb, swh, hlw, ras, pre_rad,      &
!--- rjb            
!    &      ldiag3d, trans_trac, lssav, lssav_cc, xkzm, flipv,          &
     &      ldiag3d, trans_trac, lssav, lssav_cc, xkzm, xkzm_m,         &
     &      xkzm_h, xkzm_s, flipv,                                      &
!--- rjb            

     &      old_monin, cnvgwd, sashal, newsas, mom4ice, mstrat,         &
!  ---  input/outputs:
     &      hice, fice, tisfc, tsea, tprcp, cv, cvb, cvt,               &
     &      acv, acvb, acvt, cldwrk, phy_f3d, phy_f2d, ep,              &
     &      slc, smc, stc, upd_mf, dwn_mf, det_mf,                      &
     &      srflag, snwdph, sheleg, sncovr, zorl, canopy,               &
     &      ffmm, ffhh, f10m, srunoff, evbsa, evcwa, snohfa,            &
     &      transa, sbsnoa, snowca, soilm, tmpmin, tmpmax,              &
     &      dugwd, dvgwd, psmean, bengsh, spfhmin, spfhmax,             &
     &      dusfc, dvsfc, dtsfc, dqsfc, geshem, gflux,                  &
     &      runoff, dt3dt, dq3dt, du3dt, dv3dt,                         &
     &      dlwsfc, ulwsfc, dlwsfc_cc, ulwsfc_cc, swsfc_cc,             &
     &      dusfc_cc, dvsfc_cc, dtsfc_cc, dqsfc_cc, precr_cc,           &
!  ---  outputs:
     &      gt0, gq0, gu0, gv0, t2m, q2m, u10m, v10m,                   &
     &      zlvl, psurf, hpbl, pwat, t1, q1, u1, v1,                    &
     &      chh, cmm, dlwsfci, ulwsfci, dswsfci, uswsfci,               &
     &      dtsfci, dqsfci, gfluxi, epi,                                &
     &      xmu_cc, dlw_cc, dsw_cc, snw_cc, lprec_cc                    &
     &    )

!
      use machine ,   only : kind_phys
      use physcons,   only : con_cp, con_fvirt, con_g, con_rd, con_rv,  &
     &                       con_hvap, con_hfus, con_rerth, con_pi
      USE BCLTRIGGER, ONLY : hybedmf, dspheat

!
      implicit none
!

!  ---  constant parameters:
      real (kind=kind_phys), parameter :: hocp    = con_hvap/con_cp
      real (kind=kind_phys), parameter :: fhourpr = 0.0
      real (kind=kind_phys), parameter :: rhc_max = 0.9999      ! 20060512
!     real (kind=kind_phys), parameter :: rhc_max = 0.999       ! for pry
      real (kind=kind_phys), parameter :: qmin    = 1.0e-10
      real (kind=kind_phys), parameter :: p850    = 85.0
      real (kind=kind_phys), parameter :: epsq    = 1.e-20
      real (kind=kind_phys), parameter :: hsub    = con_hvap+con_hfus
      real (kind=kind_phys), parameter :: cb2mb   = 10.0
      real (kind=kind_phys), parameter :: czmin   = 0.0001      ! cos(89.994)

!     real (kind=kind_phys), parameter ::                               &
!    & dxmax=log(1.0/7200.0),          dxmin=log(1.0/192.0)
!    & dxmax=ln(1.0/14000.0),          dxmin=ln(1.0/192.0)
!    & dxmax=ln(1.0/(14000.0*7000.0)), dxmin=ln(1.0/(192.0*94.0)
!    & dxmax=ln(1.0/(4000.0*2000.0)),  dxmin=ln(1.0/(192.0*94.0)
!    & dxmax=ln(1.0/(3000.0*1500.0)),  dxmin=ln(1.0/(192.0*94.0)
!    & dxmax=ln(1.0/(2500.0*1250.0)),  dxmin=ln(1.0/(192.0*94.0)
!    & dxmax=ln(1.0/(2000.0*1000.0)),  dxmin=ln(1.0/(192.0*94.0)
!    & dxmax=ln(1.0/(5000.0*2500.0)),  dxmin=ln(1.0/(192.0*94.0)

      real (kind=kind_phys), parameter ::                               &
!    & dxmax=-8.8818363,   dxmin=-5.2574954,   dxinv=1.0/(dxmax-dxmin)
!    & dxmax=-9.5468126,   dxmin=-5.2574954,   dxinv=1.0/(dxmax-dxmin)
!    & dxmax=-18.40047804, dxmin=-9.800790154, dxinv=1.0/(dxmax-dxmin)
!    & dxmax=-15.8949521,  dxmin=-9.800790154, dxinv=1.0/(dxmax-dxmin)
!    & dxmax=-15.31958795, dxmin=-9.800790154, dxinv=1.0/(dxmax-dxmin)
!    & dxmax=-14.95494484, dxmin=-9.800790154, dxinv=1.0/(dxmax-dxmin)
!    & dxmax=-14.50865774, dxmin=-9.800790154, dxinv=1.0/(dxmax-dxmin)
     & dxmax=-16.118095651,dxmin=-9.800790154, dxinv=1.0/(dxmax-dxmin)

!  ---  inputs:
      integer, intent(in) :: ix, im, levs, lsoil, lsm, ntrac,           &
     &      ncld, ntoz, ntcw, nmtvr, nrcm, ko3, lonf, latg, jcap,       &
     &      num_p3d, num_p2d, kdt, lat, me, pl_coeff


      integer, intent(in) :: nlons(im), ncw(2)

      logical, intent(in) :: ras, pre_rad, ldiag3d, flipv, old_monin,   &
     &                       cnvgwd, sashal, newsas, lssav, lssav_cc,   &
     &                       mom4ice, mstrat, trans_trac

      real (kind=kind_phys), dimension(im),           intent(in) ::     &
     &       sinlat, coslat, pgr, rcs2, dpshc, xlon, xlat,              &
     &       slope, shdmin, shdmax, snoalb, tg3, slmsk, vfrac,          &
     &       vtype, stype, uustar, oro, coszen, sfcnsw, sfcdsw,         &
     &       sfcdlw, tsflw, sfalb, sfcemis

      real (kind=kind_phys), dimension(ix,levs),      intent(in) ::     &
     &       ugrs, vgrs, tgrs, vvel, prsl, prslk, phil, swh, hlw

      real (kind=kind_phys), dimension(ix,levs+1),    intent(in) ::     &
     &       prsi, prsik, phii

      real (kind=kind_phys), intent(in) ::  hprime(ix,nmtvr),           &
     &       qgrs(ix,levs,ntrac), prdout(ix,ko3,pl_coeff),              &
     &       crtrh(3), flgmin(2), xkt2(ix,nrcm), poz(ko3)

      real (kind=kind_phys), intent(in) ::  ccwf, ctei_rm, clstp,       &
     &       dtp, dtf, fhour, solhr, slag, sdec, cdec, xkzm,            &
     &       xkzm_m, xkzm_h, xkzm_s

!  ---  input/output:
      real (kind=kind_phys), dimension(im),           intent(inout) ::  &
     &       hice, fice, tisfc, tsea, tprcp, cv, cvb, cvt, srflag,      &
     &       snwdph, sheleg, sncovr, zorl, canopy, ffmm, ffhh, f10m,    &
     &       srunoff, evbsa, evcwa, snohfa, transa, sbsnoa, snowca,     &
     &       soilm, tmpmin, tmpmax, dusfc, dvsfc, dtsfc, dqsfc,         &
     &       geshem, gflux, dlwsfc, ulwsfc, runoff, ep, cldwrk,         &
     &       dugwd, dvgwd, psmean, bengsh, spfhmin, spfhmax,            &
     &       acv, acvb, acvt, dlwsfc_cc, ulwsfc_cc, dtsfc_cc,           &
     &       swsfc_cc, dusfc_cc, dvsfc_cc, dqsfc_cc, precr_cc

      real (kind=kind_phys), dimension(ix,lsoil),     intent(inout) ::  &
     &       smc, stc, slc

      real (kind=kind_phys), dimension(ix,levs),      intent(inout) ::  &
     &       upd_mf, dwn_mf, det_mf

      real (kind=kind_phys),                          intent(inout) ::  &
     &       phy_f3d(ix,levs,num_p3d), phy_f2d(ix,num_p2d),             &
     &       dt3dt(ix,levs,6), du3dt(ix,levs,4), dv3dt(ix,levs,4),      &
     &       dq3dt(ix,levs,5+pl_coeff)

!  ---  output:
      real (kind=kind_phys), dimension(im),           intent(out) ::    &
     &       t2m, q2m, u10m, v10m, zlvl, psurf, hpbl, pwat, t1, q1,     &
     &       u1, v1, chh, cmm, dlwsfci, ulwsfci, dswsfci, uswsfci,      &
     &       dtsfci, dqsfci, gfluxi, epi, xmu_cc, dlw_cc, dsw_cc,       &
     &       snw_cc, lprec_cc

      real (kind=kind_phys), dimension(ix,levs),      intent(out) ::    &
     &       gt0, gu0, gv0

      real (kind=kind_phys), dimension(ix,levs,ntrac),intent(out) ::    &
     &       gq0

!  ---  local:
      real (kind=kind_phys), dimension(im)         :: ccwfac, garea,    &
     &       dlength, xncw, cumabs, qmax, cice, zice, tice, gflx,       &
     &       rain, rainc, rainl, rain1, raincs, snowmt, cd, cdq,        &
     &       qss, dusfcg, dvsfcg, dusfc1, dvsfc1, dtsfc1, dqsfc1, rb,   &
     &       rhscnpy, drain, cld1d, evap, hflx, stress, t850, ep1d,     &
     &       gamt, gamq, sigmaf, rcl, rcs, oc, theta, gamma, sigma,     &
     &       elvmax, wind, work1, work2, runof, xmu, oro_land, fm10,    &
     &       fh2, tsurf, tx1,      ctei_r, flgmin_l, evbs, evcw, trans, &
     &       sbsno, snowc, adjsfcdsw, adjsfcnsw, adjsfcdlw, adjsfculw,  &
     &       asfcdlw, asfculw, gsfcdlw, gsfculw, xcosz, snohf

      real (kind=kind_phys), dimension(ix,levs)    :: ud_mf, dd_mf,     &
     &       dt_mf, del
      real(kind=kind_phys), dimension(im,levs-1)   :: dkt            !rjb

      real (kind=kind_phys), dimension(im,levs)    :: rhc, sr, dtdt,    &
     &       dudt, dvdt, gwdcu, gwdcv, diagn1, diagn2, cuhr, cumchr,    &
     &       qr_col, fc_ice

      real (kind=kind_phys), dimension(im,lsoil)   :: smsoil, stsoil,   &
     &       ai, bi, cci, rhsmc, zsoil, slsoil

      real (kind=kind_phys) :: rhbbot, rhbtop, rhpbl, frain, f_rain,    &
     &       f_ice, qi, qw, qr, wc, tem, tem1, tem2,                    &
     &       dqdt(im,levs,ntrac), oa4(im,4), clx(im,4)

!  --- ...  in clw, the first two varaibles are cloud water and ice.
!           from third to ntrac are convective transportable tracers,
!           third being the ozone, when ntrac=3 (valid only with ras)

      real (kind=kind_phys), allocatable :: clw(:,:,:)

      integer, dimension(im) :: kbot, ktop, kcnv, soiltyp, vegtype,     &
     &          kpbl, slopetyp, kinver, lmh, levshc

      integer :: i, nvdiff, kk, ic, k, n, ipr, lv, k1, iter, levshcm,   &
     &          tracers, trc_shft, tottracer

      logical, dimension(im) :: flag_iter, flag_guess, invrsn

      logical :: lprnt
!
!===> ...  begin here
!

!  --- ...  set up check print point
!
!     lprnt = .true.
      lprnt = .false.

!     ipr = 1
!     lprnt = ( kdt>0 .and. ilon==1 )
!     do i = 1, im
!       work1(1) = xlon(i) * 57.29578
!       if ( work1(1) >= 180.0 ) work1(1) = work1(1) - 360.0
!       work2(1) = xlat(i) * 57.29578
!       print *,' me=',me,' work1=',work1(1),' work2=',work2(1),' i=',i

!       lprnt = kdt > 4320
!       lprnt = kdt > 0 .and. abs(work1(1)-103.0) < 1.0                 &
!    &                  .and. abs(work2(1)-20.0)  < 0.5
!       lprnt = kdt >= 14 .and. lat == 43
!       lprnt = kdt >= 1 .and. abs(xlon(i)*57.29578-114.325) < 0.101    &
!    &                   .and. abs(xlat(i)*57.29578+6.66)  < 0.101
!       print *,' i=',i,' xlon=',xlon(i)*57.29578,                      &
!    &                  ' xlat=',xlat(i)*57.29578,' i=',i,' me=',me
!       if ( lprnt ) then
!         ipr = i
!         exit
!       endif
!     enddo

!     lprnt = .false.
!     if ( lprnt ) then
!       print *,' im=',im,' ix=',ix,' levs=',levs,' lsoil=',lsoil,      &
!    &   ' ntrac=',ntrac,' ntoz=',ntoz,' ntcw=',ntcw,' me=',me,         &
!    &   ' xlat=',xlat(ipr),' kdt=',kdt,' slmsk=',slmsk(ipr),           &
!    &   ' nrcm=',nrcm,                                                 &
!    &   ' xlon=',xlon(ipr),' sfalb=',sfalb(ipr),' kdt=',kdt
!       print *,' pgr=',pgr(ipr),' kdt=',kdt,' ipr=',ipr
!       print *,' ipr=',ipr,' phy_f2d=',phy_f2d(ipr,1:num_p2d)
!       print *,' ugrs=',ugrs(ipr,:)
!       print *,' vgrs=',vgrs(ipr,:)
!       print *,' tgrs=',tgrs(ipr,:),' kdt=',kdt,' ipr=',ipr
!    &,         ' xlon=',xlon(ipr),' xlat=',xlat(ipr)
!       print *,' qgrs=',qgrs(ipr,:,1)
!       print *,' ozg=',qgrs(ipr,:,2)
!       print *,' clw=',qgrs(ipr,:,3)
!    &,         ' xlon=',xlon(ipr),' xlat=',xlat(ipr)
!     endif
!
      lmh(:) = levs
      nvdiff = ntrac    ! vertical diffusion of all tracers!
!
!  --- ...  figure out how many extra tracers are there
!
      if ( trans_trac ) then
        if ( ntcw > 0 ) then
          if ( ntoz < ntcw ) then
            trc_shft = ntcw + ncld
          else
            trc_shft = ntoz
          endif
        elseif ( ntoz > 0 ) then
          trc_shft = ntoz
        else
          trc_shft = 1
        endif

        tracers   = ntrac - trc_shft
        if ( ntoz > 0 ) tottracer = tracers + 1  ! ozone is added separately
      else
        tottracer = 0                            ! no convective transport of tracers
      endif

!     if ( lprnt ) print *,' trans_trac=',trans_trac,' tottracer=',     &
!    &                     tottracer,' trc_shft=',trc_shft,' kdt=',kdt

      allocate ( clw(ix,levs,tottracer+2) )

      if ( ras ) then
        if ( ccwf >= 0.0 ) then
          ccwfac = ccwf
        else
          ccwfac = -999.0
        endif
      endif
!
      call get_prs( im,ix,levs,ntrac,tgrs,qgrs                          &
     &,             prsi,prsik,prsl,prslk,phii,phil,del )
!
      rhbbot = crtrh(1)
      rhpbl  = crtrh(2)
      rhbtop = crtrh(3)

!  --- ...  for pry version

      do i = 1, im
        levshc(i) = 0
      enddo

      do k = 2, levs
        do i = 1, im
          if ( prsi(i,1)-prsi(i,k) <= dpshc(i) ) levshc(i) = k
        enddo
      enddo

      levshcm = 1
      do i = 1, im
        levshcm = max( levshcm, levshc(i) )
      enddo

      if ( mstrat ) then
        levshcm = max( levshcm, levs/2 )
      else
        levshcm = min( levshcm, levs/2 )
      endif

!  --- ...  frain=factor for centered difference scheme correction of rain amount.

      frain = dtf / dtp

      do i = 1, im
        soiltyp(i) = int( stype(i) + 0.5 )
        sigmaf(i)  = max( vfrac(i), 0.01 )
!       sigmaf(i)  = max( vfrac(i), 0.3 )
        if ( lsm == 0 ) sigmaf(i) = 0.5 + vfrac(i)*0.5

        vegtype(i)  = int( vtype(i) + 0.5 )
        slopetyp(i) = int( slope(i) + 0.5 ) !! clu: slope -> slopetyp

        if ( slmsk(i) == 2.0 ) then
          soiltyp(i) = 9
          vegtype(i) = 13
          slopetyp(i) = 9                 !! clu: qa(slopetyp)
        endif
      enddo

!  --- ...  transfer soil moisture and temperature from global to local variables

      do k = 1, lsoil
        do i = 1, im
          smsoil(i,k) = smc(i,k)
          stsoil(i,k) = stc(i,k)
          slsoil(i,k) = slc(i,k)          !! clu: slc -> slsoil
        enddo
      enddo

!  --- ...  xw: transfer ice thickness & concentration from global to local variables

      do i = 1, im
        zice(i) = hice(i)
        cice(i) = fice(i)
        tice(i) = tisfc(i)
      enddo

      do k = 1, levs
        do i = 1, im
          dudt(i,k) = 0.0
          dvdt(i,k) = 0.0
          dtdt(i,k) = 0.0
        enddo
      enddo

      do n = 1, ntrac
        do k = 1, levs
          do i = 1, im
            dqdt(i,k,n) = 0.0      ! dqdt may be dimensioned (levs,ntrac)
          enddo
        enddo
      enddo

      do i = 1, im
        rcl(i)   = rcs2(i)
        rcs(i)   = sqrt(rcl(i))
        psurf(i) = pgr(i)
        work1(i) = prsik(i,1) / prslk(i,1)
!       garea(i) = con_rerth * con_rerth                                &
!    &           * (con_pi+con_pi)*con_pi*coslat(i)/(nlons(i)*latg)
        tem1       = con_rerth * (con_pi + con_pi)*coslat(i) / nlons(i)
        tem2       = con_rerth * con_pi/latg
        garea(i)   = tem1 * tem2
        dlength(i) = sqrt( tem1*tem1 + tem2*tem2 )
      enddo

      if ( lssav ) then
        do i = 1, im
          psmean(i) = psmean(i) + pgr(i)*dtf
        enddo
      endif

!  --- ...  initialize dtdt with heating rate from dcyc2

!     if ( lprnt ) then
!       do ipr=1,im
!         print *,' before dcyc2: im=',im,' lsoil=',lsoil,' levs=',levs &
!    &,    ' sde=',sdec,' cdec=',cdec,' tsea=',tsea(ipr),' ipr=',ipr    &
!    &,    ' lat=',lat,' me=',me,' kdt=',kdt                            &
!    &,    ' sfcdlw=',sfcdlw(ipr),' sfcnsw=',sfcnsw(ipr)
!         print *,' hlw=',hlw(ipr,:),' me=',me,' lat=',lat,xlon(ipr)
!         print *,' swh=',swh(ipr,:),' me=',me,' lat=',lat,xlon(ipr)
!       enddo
!     endif

!  --- ...  adjust radiation fluxes and heating rates (longer calling period
!           of time to fit for shorter model time steps.
!      sw:  using cos of zenith angle as scaling factor
!      lw:  using surface air temperature as scaling factor

      if ( pre_rad ) then

        call dcyc2t3_pre_rad                                            &
!  ---  inputs:
     &     ( solhr,slag,sdec,cdec,sinlat,coslat,                        &
     &       xlon,coszen,tsea,tgrs(1,1),tgrs(1,1),                      &
     &       sfcdsw,sfcnsw,sfcdlw,swh,hlw,                              &
     &       ix, im, levs,                                              &
!  ---  input/output:
     &       dtdt,                                                      &
!  ---  outputs:
     &       adjsfcdsw,adjsfcnsw,adjsfcdlw,adjsfculw,xmu,xcosz          &
!old vars   ( dswsfc,    -radsl,  dlwsf1,   ulwsf1,  xmu,xcosz )
     &     )

      else

        call dcyc2t3                                                    &
!  ---  inputs:
     &     ( solhr,slag,sdec,cdec,sinlat,coslat,                        &
     &       xlon,coszen,tsea,tgrs(1,1),tsflw,                          &
     &       sfcdsw,sfcnsw,sfcdlw,swh,hlw,                              &
     &       ix, im, levs,                                              &
!  ---  input/output:
     &       dtdt,                                                      &
!  ---  outputs:
     &       adjsfcdsw,adjsfcnsw,adjsfcdlw,adjsfculw,xmu,xcosz          &
!old vars   ( dswsfc,    -radsl,  dlwsf1,   ulwsf1,  xmu,xcosz )
     &     )

      endif

!  ---  convert lw fluxes for land/ocean/sea-ice models
!  note: adjsfcdlw and adjsfculw are time-step adjusted lw sfc dn/up fluxes.
!        those fluxes from dcyc2t3 subprogram are not applied with emissivity
!        effect.  one needs to be cautious that, when emis<1, sfc dn/up lw
!        fluxes will expressed differently when looking from atmospheric
!        (as the output from radiation) to those applied for land/ocean/sea-ice
!        model processes. however, the net effects are the same.
!
!   - flux to/from lnd/oc/ice:  gsfcdlw=emis*adjsfcdlw; gsfculw=emis*adjsfculw
!   - flux to/from atmos mdl: asfcdlw=adjsfcdlw; asfculw=emis*adjsfculw+(1-emis)*adjsfcdlw

        do i = 1, im
          gsfcdlw(i) = sfcemis(i) * adjsfcdlw(i)                  ! for lnd/ocn/sice
          gsfculw(i) = sfcemis(i) * adjsfculw(i)                  ! for lnd/ocn/sice
          asfcdlw(i) = adjsfcdlw(i)                               ! for atmos output
          asfculw(i) = gsfculw(i) + adjsfcdlw(i) - gsfcdlw(i)     ! for atmos output
!org      asfculw(i) = gsfculw(i) + (1.0-sfcemis(i))*adjsfcdlw(i) ! for atmos output
        enddo

!  --- ...  coupling insertion
!  note: all radiation fluxes are possitive values.
!        sfc net sw is defined as sfcnsw = sfcdsw - sfcusw (*** positive downward
!        see subr grrad for def)

      if ( lssav_cc ) then
        do i = 1, im
          xmu_cc(i)    = max( 0.0, min( 1.0, xcosz(i) ) )
          dsw_cc(i)    = adjsfcdsw(i)
          dlw_cc(i)    = adjsfcdlw(i)
          dlwsfc_cc(i) = dlwsfc_cc(i) + gsfcdlw(i)
          ulwsfc_cc(i) = ulwsfc_cc(i) + gsfculw(i)
          swsfc_cc(i)  = swsfc_cc(i)  - adjsfcnsw(i) ! swsfc is positive upward
!         sfcems_cc(i) = sfcemis(i)                  ! for ocn mdl if needed
        enddo
      endif

!  --- ...  accumulate/save output variables

      if ( lssav ) then

!  --- ...  sfc lw fluxes used by atmospheric model are saved for output

        do i = 1, im
          dlwsfc(i) = dlwsfc(i) + asfcdlw(i)*dtf
          ulwsfc(i) = ulwsfc(i) + asfculw(i)*dtf
        enddo

        if ( ldiag3d ) then
          do k = 1, levs
            do i = 1, im
              dt3dt(i,k,1) = dt3dt(i,k,1) + hlw(i,k)*dtf
              dt3dt(i,k,2) = dt3dt(i,k,2) + swh(i,k)*dtf*xmu(i)
            enddo
          enddo
        endif

      endif    ! end if_lssav_block

      do i = 1, im
        kcnv(i)   = 0
        kinver(i) = levs
        invrsn(i) = .false.
        tx1(i)    = 0.0
        ctei_r(i) = 10.0
      enddo

!     ctei_rm = 0.20
!     ctei_rm = 0.25
!     ctei_rm = 0.50
!     ctei_rm = 0.70

      do k = 1, levs/2
        do i = 1, im

          if ( prsi(i,1)-prsi(i,k+1) < 0.35*prsi(i,1)                   &
     &        .and. (.not. invrsn(i)) ) then
            tem = (tgrs(i,k+1) - tgrs(i,k)) / (prsl(i,k) - prsl(i,k+1))

!           if ( tem > 0.02 .and. tx1(i) < 0.0 ) then
!           if ( tem > 0.10 .and. tx1(i) < 0.0 ) then
!           if ( tem > 0.05 .and. tx1(i) < 0.0 ) then
            if ( tem > 0.025 .and. tx1(i) < 0.0 ) then
              invrsn(i) = .true.

              if ( qgrs(i,k,1)/=qgrs(i,k+1,1) .and. .not.sashal ) then
                tem1 = (1.0 + hocp*max(qgrs(i,k+1,1),qmin)/tgrs(i,k+1))
                tem2 = (1.0 + hocp*max(qgrs(i,k,1),qmin)/tgrs(i,k))

                tem1 = tem1 * tgrs(i,k+1) / prslk(i,k+1)                &
     &               - tem2 * tgrs(i,k)   / prslk(i,k)

!  --- ...  (cp/l)(delthetae)/(deltotwater) < 0.7
                ctei_r(i) = (1.0/hocp)*tem1/(qgrs(i,k+1,1)-qgrs(i,k,1)  &
     &                    + qgrs(i,k+1,ntcw)-qgrs(i,k,ntcw))
              else
                ctei_r(i) = 10
              endif

              kinver(i) = k
!             kinver(i) = k + 1
            endif

            tx1(i) = tem
          endif

        enddo
      enddo

!  --- ...  check print

!     ipr = 1
!     if ( lprnt ) then
!       print *,' before progtm: im=',im,' lsoil=',lsoil                &
!    &,         ' nvdiff=',nvdiff,' adjsfcnsw=',adjsfcnsw(ipr)          &
!    &,         ' asfcdlw=',asfcdlw(ipr),' asfculw=',asfculw(ipr)       &
!    &,         ' gsfcdlw=',gsfcdlw(ipr),' gsfculw=',gsfculw(ipr)       &
!    &,         ' sfcemis=',sfcemis(ipr),' tsea2=',tsea(ipr)            &
!    &,         ' ipr=',ipr,' me=',me,' lat=',lat,' xlon=',xlon(ipr)    &
!    &,         ' kdt=',kdt
!       print *,' dtdth=',dtdt(ipr,:),' kdt=',kdt
!     endif

!  --- ...  lu: initialize flag_guess, flag_iter, tsurf

      do i = 1, im
        tsurf(i)      = tsea(i)
        flag_guess(i) = .false.
        flag_iter(i)  = .true.
        drain(i)      = 0.0
        ep1d(i)       = 0.0
        runof(i)      = 0.0
        hflx(i)       = 0.0
        gflx(i)       = 0.0
        evap(i)       = 0.0

        evbs(i)       = 0.0
        evcw(i)       = 0.0
        trans(i)      = 0.0
        sbsno(i)      = 0.0
        snowc(i)      = 0.0
        snohf(i)      = 0.0
      enddo

!  --- ...  lu: iter-loop over (sfc_diff,sfc_drv,sfc_ocean,sfc_sice)

      do iter = 1, 2

!  --- ...  surface exchange coefficients

        call sfc_diff(im,pgr,ugrs,vgrs,tgrs,qgrs,                       &
     &                tsea,zorl,cd,cdq,rb,                              &
     &                rcl,prsl(1,1),work1,slmsk,                        &
     &                stress,ffmm,ffhh,                                 &
     &                uustar,wind,phy_f2d(1,num_p2d),fm10,fh2,          &
     &                tsurf, flag_iter)

!       if ( lprnt ) print *,' cdq=',cdq(ipr),' iter=',iter             &
!    &,   ' wind=',wind(ipr),'phy_f2d=',phy_f2d(ipr,num_p2d),' ugrs='   &
!    &,   ugrs(ipr,1),' vgrs=',vgrs(ipr,1)

!  --- ...  lu: update flag_guess

        do i = 1, im
          if ( iter==1 .and. wind(i)<2.0 ) then
            flag_guess(i) = .true.
          endif
        enddo

!  --- ...  surface energy balance over ocean

        call sfc_ocean                                                  &
!  ---  inputs:
     &     ( im,pgr,ugrs,vgrs,tgrs,qgrs,tsea,cd,cdq,rcl,                &
     &       prsl(1,1),work1,slmsk,phy_f2d(1,num_p2d),flag_iter,        &
!  ---  outputs:
     &       qss,cmm,chh,evap,hflx                                      &
     &     )

!       if ( lprnt ) print *,' sfalb=',sfalb(ipr),' ipr=',ipr           &
!    &,   ' sheleg=',sheleg(ipr),' snwdph=',snwdph(ipr)                 &
!    &,   ' tprcp=',tprcp(ipr),' kdt=',kdt,' iter=',iter

!  --- ...  surface energy balance over land

        if ( lsm == 1 ) then                        ! noah lsm call

          call sfc_drv                                                  &
!  ---  inputs:
     &     ( im,lsoil,pgr,ugrs,vgrs,tgrs,qgrs,soiltyp,vegtype,sigmaf,   &
     &       sfcemis,gsfcdlw,adjsfcdsw,adjsfcnsw,dtf,tg3,cd,cdq,        &
     &       rcl,prsl(1,1),work1,slmsk,phy_f2d(1,num_p2d),slopetyp,     &
     &       shdmin,shdmax,snoalb,sfalb,flag_iter,flag_guess,           &
!  ---  in/outs:
     &       sheleg,snwdph,tsea,tprcp,srflag,smsoil,stsoil,slsoil,      &
     &       canopy,trans,tsurf,                                        &
!  ---  outputs:
     &       sncovr,qss,gflx,drain,evap,hflx,ep1d,runof,                &
     &       cmm,chh,zlvl,evbs,evcw,sbsno,snowc,soilm,snohf             &
     &     )

        else                                       ! osu lsm call

          call sfc_land                                                 &
!  ---  inputs:
     &     ( im,lsoil,pgr,ugrs,vgrs,tgrs,qgrs,smsoil,soiltyp,           &
     &       sigmaf,vegtype,sfcemis,gsfcdlw,adjsfcnsw,dtf,              &
     &       zorl,tg3,cd,cdq,rcl,prsl(1,1),work1,slmsk,                 &
     &       phy_f2d(1,num_p2d),flag_iter,flag_guess,                   &
!  ---  input/outputs:
     &       sheleg,tsea,tprcp,srflag,stsoil,canopy,tsurf,              &
!  ---  outputs:
     &       qss,snowmt,gflx,zsoil,rhscnpy,rhsmc,                       &
     &       ai,bi,cci,drain,evap,hflx,ep1d,cmm,chh,                    &
     &       zlvl,evbs,evcw,trans,sbsno,snowc,soilm,snohf               &
     &     )

        endif

!       if ( lprnt ) print *,' tseabeficemodel =',tsea(ipr),' me=',me   &
!    &,   ' kdt=',kdt

!  --- ...  surface energy balance over seaice

        call sfc_sice                                                   &
!  ---  inputs:
     &     ( im,lsoil,pgr,ugrs,vgrs,tgrs,qgrs,dtf,                      &
     &       sfcemis,gsfcdlw,adjsfcnsw,adjsfcdsw,srflag,                &
     &       cd,cdq,rcl,prsl(1,1),work1,slmsk,phy_f2d(1,num_p2d),       &
     &       flag_iter,mom4ice,lsm,                                     &
!  ---  input/outputs:
     &       zice,cice,tice,sheleg,tsea,tprcp,stsoil,ep1d,              &
!  ---  outputs:
     &       snwdph,qss,snowmt,gflx,cmm,chh,                            &
     &       zlvl,evap,hflx                                             &
     &     )

!  --- ...  lu: update flag_iter and flag_guess

        do i = 1, im
          flag_iter(i)  = .false.
          flag_guess(i) = .false.

          if ( slmsk(i)==1.0 .and. iter==1 ) then
            if ( wind(i) < 2.0 ) flag_iter(i) = .true.
          endif
        enddo

      enddo   ! end iter_loop

      do i = 1, im
        epi(i)     = ep1d(i)
        dlwsfci(i) = gsfcdlw(i)
        ulwsfci(i) = gsfculw(i)
        uswsfci(i) = adjsfcdsw(i) - adjsfcnsw(i)
        dswsfci(i) = adjsfcdsw(i)
        gfluxi(i)  = gflx(i)
        t1(i)      = tgrs(i,1)
        q1(i)      = qgrs(i,1,1)
        u1(i)      = ugrs(i,1)
        v1(i)      = vgrs(i,1)
      enddo

      if ( lsm == 0 ) then                        ! osu lsm call
        do i = 1, im
         sncovr(i) = 0.0
         if ( sheleg(i) > 0.0 ) sncovr(i) = 1.0
        enddo
      endif

!  --- ...  update near surface fields

      call sfc_diag(im,lsoil,pgr,ugrs,vgrs,tgrs,qgrs,                   &
     &              tsea,qss,f10m,u10m,v10m,t2m,q2m,rcl,work1,slmsk,    &
     &              evap,ffmm,ffhh,fm10,fh2)

      do i = 1, im
        phy_f2d(i,num_p2d) = 0.0
      enddo

!     if ( lprnt ) print *,' tseaim=',tsea(ipr),' me=',me,' kdt=',kdt

      if ( lssav ) then
        do i = 1, im
          gflux (i) = gflux (i) + gflx (i)*dtf
          evbsa (i) = evbsa (i) + evbs (i)*dtf
          evcwa (i) = evcwa (i) + evcw (i)*dtf
          transa(i) = transa(i) + trans(i)*dtf
          sbsnoa(i) = sbsnoa(i) + sbsno(i)*dtf
          snowca(i) = snowca(i) + snowc(i)*dtf
          snohfa(i) = snohfa(i) + snohf(i)*dtf
          ep    (i) = ep    (i) + ep1d (i)*dtf

          tmpmax(i) = max( tmpmax(i), t2m(i) )
          tmpmin(i) = min( tmpmin(i), t2m(i) )

          spfhmax(i) = max( spfhmax(i), q2m(i) )
          spfhmin(i) = min( spfhmin(i), q2m(i) )
        enddo
      endif

!  --- ...  vertical diffusion

!     if ( lprnt ) print *,' tsea3=',tsea(ipr),' slmsk=',slmsk(ipr)     &
!    &, ' kdt=',kdt,' evap=',evap(ipr)
!     if ( lprnt )  print *,' dtdtb=',dtdt(ipr,:)

      do i = 1, im
        if ( slmsk(i) == 0.0 ) then
          oro_land(i) = 0.0
        else
          oro_land(i) = oro(i)
        endif
      enddo



      if (hybedmf) then

        call moninedmf(ix,im,levs,nvdiff,ntcw,dvdt,dudt,dtdt,dqdt,      &
     &     ugrs,vgrs,tgrs,qgrs,swh,hlw,xmu,                             &
     &     prsik(1,1),rb,zorl,u10m,v10m,ffmm,ffhh,                      &
     &     tsea,qss,hflx,evap,stress,wind,kpbl,                         &
     &     prsi,del,prsl,prslk,phii,phil,rcs,dtp,dspheat,               &
     &     dusfc1,dvsfc1,dtsfc1,dqsfc1,hpbl,gamt,gamq,dkt,              &
     &     kinver, xkzm_m, xkzm_h, xkzm_s, lprnt, ipr)

      elseif (.not. old_monin) then

          call moninq(ix,im,levs,nvdiff,dvdt,dudt,dtdt,dqdt,            &
     &     ugrs,vgrs,tgrs,qgrs,swh,hlw,xmu,slmsk,                       &
     &     prsik(1,1),rb,ffmm,ffhh,tsea,qss,hflx,evap,stress,wind,kpbl, &
     &     prsi,del,prsl,prslk,phii,phil,rcs,dtp,                       &
     &     dusfc1,dvsfc1,dtsfc1,dqsfc1,hpbl,gamt,gamq)

      else

        if ( mstrat ) then
          call moninp1(ix,im,levs,nvdiff,dvdt,dudt,dtdt,dqdt,           &
     &     ugrs,vgrs,tgrs,qgrs,                                         &
     &     prsik(1,1),rb,ffmm,ffhh,tsea,qss,hflx,evap,stress,wind,kpbl, &
     &     prsi,del,prsl,prslk,phii,phil,rcl,dtp,                       &
     &     dusfc1,dvsfc1,dtsfc1,dqsfc1,hpbl,gamt,gamq,                  &
     &     kinver, oro_land, xkzm)
        else
          call moninp(ix,im,levs,nvdiff,dvdt,dudt,dtdt,dqdt,            &
     &     ugrs,vgrs,tgrs,qgrs,                                         &
     &     prsik(1,1),rb,ffmm,ffhh,tsea,qss,hflx,evap,stress,wind,kpbl, &
     &     prsi,del,prsl,prslk,phii,phil,rcl,dtp,                       &
     &     dusfc1,dvsfc1,dtsfc1,dqsfc1,hpbl,gamt,gamq,xkzm)
        endif

      endif   ! end if_hydedmf


!     if ( ntoz > 0 ) dqdt(:,:,ntoz) = 0.0
!
!     if ( lprnt ) then
!       print *,' dusfc1=',dusfc1(ipr)
!       print *,' dtsfc1=',dtsfc1(ipr)
!       print *,' dqsfc1=',dqsfc1(ipr)
!       print *,' dtdt=',dtdt(ipr,:)
!       print *,' dudtm=',dudt(ipr,:)
!     endif

!  --- ...  coupling insertion

      if ( lssav_cc ) then
        do i = 1, im
          dusfc_cc(i) = dusfc_cc(i) + dusfc1(i) !*dth <-na h.(?)
          dvsfc_cc(i) = dvsfc_cc(i) + dvsfc1(i) !*dth <-na h.(?)
          dtsfc_cc(i) = dtsfc_cc(i) + dtsfc1(i) !*dth <-na h.(?)
          dqsfc_cc(i) = dqsfc_cc(i) + dqsfc1(i) !*dth <-na h.(?)
        enddo
      endif

      if ( lssav ) then
        do i = 1, im
          dusfc(i)  = dusfc(i) + dusfc1(i)*dtf
          dvsfc(i)  = dvsfc(i) + dvsfc1(i)*dtf
          dtsfc(i)  = dtsfc(i) + dtsfc1(i)*dtf
          dqsfc(i)  = dqsfc(i) + dqsfc1(i)*dtf
          dtsfci(i) = dtsfc1(i)
          dqsfci(i) = dqsfc1(i)
        enddo

        if ( ldiag3d ) then
          do k = 1, levs
            do i = 1, im
              tem1 = rcs(i) * dtf
              tem  = dtdt(i,k) - (hlw(i,k)+swh(i,k)*xmu(i))
              dt3dt(i,k,3) = dt3dt(i,k,3) + tem*dtf
              dq3dt(i,k,1) = dq3dt(i,k,1) + dqdt(i,k,1)*dtf
              du3dt(i,k,1) = du3dt(i,k,1) + dudt(i,k)*tem1
              du3dt(i,k,2) = du3dt(i,k,2) - dudt(i,k)*tem1
              dv3dt(i,k,1) = dv3dt(i,k,1) + dvdt(i,k)*tem1
              dv3dt(i,k,2) = dv3dt(i,k,2) - dvdt(i,k)*tem1
            enddo
          enddo
          if ( ntoz > 0 ) then
            do k = 1, levs
              do i = 1, im
                dq3dt(i,k,5) = dq3dt(i,k,5) + dqdt(i,k,ntoz)*dtf
              enddo
            enddo
          endif
        endif   ! end if_ldiag3d

      endif   ! end if_lssav

      if ( nmtvr == 6 ) then

        do i = 1, im
          oc(i) = hprime(i,2)
        enddo

        do k = 1, 4
          do i = 1, im
            oa4(i,k) = hprime(i,k+2)
            clx(i,k) = 0.0
          enddo
        enddo

      elseif ( nmtvr == 10 ) then

        do i = 1, im
          oc(i) = hprime(i,2)
        enddo

        do k = 1, 4
          do i = 1, im
            oa4(i,k) = hprime(i,k+2)
            clx(i,k) = hprime(i,k+6)
          enddo
        enddo

      elseif ( nmtvr == 14 ) then

!  --- ...  get the kim fields (until this is changed)

        do i = 1, im
          oc(i) = hprime(i,2)
        enddo

        do k = 1, 4
          do i = 1, im
            oa4(i,k) = hprime(i,k+2)
            clx(i,k) = hprime(i,k+6)
          enddo
        enddo

        do i = 1, im
          theta(i)  = hprime(i,11)
          gamma(i)  = hprime(i,12)
          sigma(i)  = hprime(i,13)
          elvmax(i) = hprime(i,14)
        enddo

      else

        oc    = 0.0
        oa4   = 0.0
        clx   = 0.0
        theta = 0.0
        gamma = 0.0
        sigma = 0.0
        elvmax= 0.0

      endif   ! end if_nmtvr

      call gwdps(im, ix, im,   levs,  dvdt, dudt,                       &
     &           ugrs,   vgrs, tgrs,  qgrs,                             &
     &           kpbl,   prsi, del,   prsl, prslk,                      &
     &           phii,   phil, rcl, dtp,                                &
     &           kdt,    hprime(1,1), oc, oa4, clx,                     &
     &           theta,sigma,gamma,elvmax,dusfcg, dvsfcg,               &
     &           con_g,con_cp,con_rd,con_rv, lonf, nmtvr, me, lprnt,ipr)

      if ( lprnt )  print *,' dudtg=',dudt(ipr,:)

      if ( lssav ) then
        do i = 1, im
          dugwd(i) = dugwd(i) + dusfcg(i)*dtf
          dvgwd(i) = dvgwd(i) + dvsfcg(i)*dtf
        enddo

!       if ( lprnt ) print *,' dugwd=',dugwd(ipr),' dusfcg=',dusfcg(ipr)
!       if ( lprnt ) print *,' dvgwd=',dvgwd(ipr),' dvsfcg=',dvsfcg(ipr)

        if ( ldiag3d ) then
          do k = 1, levs
            do i = 1, im
              tem = rcs(i) * dtf
              du3dt(i,k,2) = du3dt(i,k,2) + dudt(i,k)*tem
              dv3dt(i,k,2) = dv3dt(i,k,2) + dvdt(i,k)*tem
            enddo
          enddo
        endif
      endif

      do k = 1, levs
        do i = 1, im
          gt0(i,k) = tgrs(i,k) + dtdt(i,k)*dtp
          gu0(i,k) = ugrs(i,k) + dudt(i,k)*dtp
          gv0(i,k) = vgrs(i,k) + dvdt(i,k)*dtp
        enddo
      enddo

      do n = 1, ntrac
        do k = 1, levs
          do i = 1, im
            gq0(i,k,n) = qgrs(i,k,n) + dqdt(i,k,n)*dtp
          enddo
        enddo
      enddo

!  --- ...  check print

!     if ( me == 0 ) then
!       sumq = 0.0
!       do k = 1, levs
!         do i = 1, im
!           sumq = sumq + (dqdt(i,k,1) + dqdt(i,k,ntcw))*del(i,k)
!         enddo
!       enddo

!       sume = 0.0
!       do i = 1, im
!         sume = sume + dqsfc1(i)
!       enddo

!       sumq = sumq * 1000.0 / con_g
!       sume = sume / con_hvap
!       print *,' after mon: sumq=',sumq,' sume=',sume, ' kdt=',kdt
!     endif

!  --- ...  ozone physics

      if ( ntoz>0 .and. ntrac>=ntoz ) then

        call ozphys(ix,im,levs,ko3,dtp,gq0(1,1,ntoz),gq0(1,1,ntoz)      &
     &,             gt0, poz, prsl, prdout, pl_coeff, del, ldiag3d      &
     &,             dq3dt(1,1,6), me)

      endif

!  --- ...  to side-step the ozone physics

!      if ( ntrac >= 2 ) then
!        do k = 1, levs
!          gq0(k,ntoz) = qgrs(k,ntoz)
!        enddo
!      endif

!     if ( lprnt ) then
!       print *,' levs=',levs,' jcap=',jcap,' dtp',dtp                  &
!    *,  ' slmsk=',slmsk(ilon,ilat),' rcs=',rcs,' kdt=',kdt
!       print *,' xkt2=',xkt2,' ncld=',ncld,' iq=',iq,' lat=',lat
!       print *,' pgr=',pgr
!       print *,' del=',del(1,:)
!       print *,' prsl=',prsl(1,:)
!       print *,' prslk=',prslk(1,:)
!       print *,' xkt2=',xkt2(ipr,1)
!       print *,' gt0=',gt0(ipr,:)                                      &
!    &,         ' kdt=',kdt,' xlon=',xlon(ipr),' xlat=',xlat(ipr)
!       print *,' dtdt=',dtdt(ipr,:)
!       print *,' gu0=',gu0(ipr,:)
!       print *,' gv0=',gv0(ipr,:)
!       print *,' gq0=',(gq0(ipr,k,1),k=1,levs)
!       print *,' gq1=',(gq0(ipr,k,ntcw),k=1,levs)
!       print *,' vvel=',vvel
!     endif

      if ( ldiag3d ) then

        do k = 1, levs
          do i = 1, im
            dtdt(i,k)   = gt0(i,k)
            dqdt(i,k,1) = gq0(i,k,1)
            dudt(i,k)   = gu0(i,k)
            dvdt(i,k)   = gv0(i,k)
          enddo
        enddo

      elseif ( cnvgwd ) then

        do k = 1, levs
          do i = 1, im
            dtdt(i,k) = gt0(i,k)
          enddo
        enddo

      endif   ! end if_ldiag3d/cnvgwd

      call get_phi(im,ix,levs,ntrac,gt0,gq0,                            &
     &             prsi,prsik,prsl,prslk,phii,phil)

!     if ( lprnt ) then
!       print *,' phii2=',phii(ipr,:)
!       print *,' phil2=',phil(ipr,:)
!     endif

      do k = 1, levs
        do i = 1, im
          clw(i,k,1) = 0.0
          clw(i,k,2) = -999.9
        enddo
      enddo

!  --- ...  for convective tracer transport (while using ras)

      if ( ras ) then
        if ( tottracer > 0 ) then

          if ( ntoz > 0 ) then
            clw(:,:,3) = gq0(:,:,ntoz)

            if ( tracers > 0 ) then
              do n = 1, tracers
                clw(:,:,3+n) = gq0(:,:,n+trc_shft)
              enddo
            endif
          else
            do n = 1, tracers
              clw(:,:,2+n) = gq0(:,:,n+trc_shft)
            enddo
          endif

        endif
      endif   ! end if_ras

      do i = 1, im
        ktop(i)  = 1
        kbot(i)  = levs
      enddo

!  --- ...  calling precipitation processes

      do i = 1, im
!pry    work1(i) = (log(1.0 / (rcs(i)*nlons(i))) - dxmin) * dxinv
!       work1(i) = (log(1.0 / (rcs(i)*nlons(i)*latg)) - dxmin) * dxinv
        work1(i) = (log(coslat(i) / (nlons(i)*latg)) - dxmin) * dxinv
        work1(i) = max( 0.0, min( 1.0, work1(i) ) )
        work2(i) = 1.0 - work1(i)
      enddo

!  --- ...  calling convective parameterization

      if ( ntcw > 0 ) then

        do k = 1, levs
          do i = 1, im
            rhc(i,k) = rhbbot - (rhbbot - rhbtop) * (1.0 - prslk(i,k))
            rhc(i,k) = rhc_max*work1(i) + rhc(i,k)*work2(i)
            rhc(i,k) = max( 0.0, min( 1.0, rhc(i,k) ) )
          enddo
        enddo

        if ( num_p3d == 3 ) then    ! call brad ferrier's microphysics
          do i = 1, im
            flgmin_l(i) = flgmin(1)*work1(i) + flgmin(2)*work2(i)
          enddo

!  --- ...  algorithm to separate different hydrometeor species

          do k = 1, levs
            do i = 1, im
              wc     = gq0(i,k,ntcw)
              qi     = 0.0
              qr     = 0.0
              qw     = 0.0
              f_ice  = max( 0.0, min( 1.0, phy_f3d(i,k,1) ) )
              f_rain = max( 0.0, min( 1.0, phy_f3d(i,k,2) ) )

              qi = f_ice * wc
              qw = wc - qi
              if ( qw > 0.0 ) then
                qr = f_rain * qw
                qw = qw - qr
              endif

!             if ( f_ice >= 1.0 ) then
!               qi = wc
!             elseif ( f_ice <= 0.0 ) then
!               qw = wc
!             else
!               qi = f_ice * wc
!               qw = wc - qi
!             endif

!             if ( qw>0.0 .and. f_rain>0.0 ) then
!               if ( f_rain >= 1.0 ) then
!                 qr = qw
!                 qw = 0.0
!               else
!                 qr = f_rain * qw
!                 qw = qw - qr
!               endif
!             endif

              qr_col(i,k) = qr
!             clw(i,k)    = qi + qw
              clw(i,k,1)  = qi
              clw(i,k,2)  = qw

!  --- ...  array to track fraction of "cloud" in the form of ice

!             if ( qi+qw > epsq ) then
!               fc_ice(i,k) = qi / (qi + qw)
!             else
!               fc_ice(i,k) = 0.0
!             endif
            enddo
          enddo

        else   ! if_num_p3d

          do k = 1, levs
            do i = 1, im
              clw(i,k,1) = gq0(i,k,ntcw)
            enddo
          enddo

        endif  ! end if_num_p3d

      else    ! if_ntcw

        rhc(:,:) = 1.0

      endif   ! end if_ntcw
!
      if ( .not. ras ) then

        if ( .not. newsas ) then

          call sascnv(im,ix,levs,jcap,dtp,del,prsl,pgr,phil,            &
     &                clw,gq0,gt0,gu0,gv0,rcs,cld1d,                    &
     &                rain1,kbot,ktop,kcnv,slmsk,                       &
     &                vvel,xkt2,ncld)

        else

          call sascnvn(im,ix,levs,jcap,dtp,del,prsl,pgr,phil,           &
     &                clw,gq0,gt0,gu0,gv0,rcs,cld1d,                    &
     &                rain1,kbot,ktop,kcnv,slmsk,                       &
     &                vvel,ncld,ud_mf,dd_mf,dt_mf)

        endif

!       if ( lprnt ) print *,' rain1=',rain1(ipr),' xkt2=',xkt2(ipr,1)

      else    ! if_not_newsas

!       if ( lprnt ) print *,' calling ras for kdt=',kdt,' me=',me      &
!    &,                    ' lprnt=',lprnt

        call rascnv(im,     ix,    levs,   dtp, dtf, xkt2               &
     &,             gt0,    gq0,   gu0,    gv0, clw, tottracer          &
     &,             prsi,   prsl,   prsik,  prslk, phil,  phii          &
     &,             kpbl,   cd,     rain1,  kbot,  ktop,  kcnv          &
     &,             phy_f2d(1,num_p2d), flipv, cb2mb                    &
     &,             me, garea, lmh, ccwfac, nrcm, rhc                   &
     &,             ud_mf, dd_mf, dt_mf, lprnt, ipr)

!  --- ...  check print

!       do i = 1, im
!         if ( tsea(i)>380.0 .or. tsea(i)<10 ) then
!           print *,' tsea=', tsea(i),' i=',i,' lat=',lat,              &
!    &        ' kdt=',kdt,' xlon=',xlon(i),' xlat=',xlat(i),' slmsk=',  &
!    &        slmsk(i),' me=',me
!           stop
!         endif
!       enddo

!       do k = 1, levs
!         do i = 1, im
!           if ( gt0(i,k)>330.0 .or. gt0(i,k)<80.0 ) then               &
!             print *,' gt0=', gt0(i,k),' i=',i,' k=',k,' lat=',lat,    &
!    &          ' kdt=',kdt,' xlon=',xlon(i),' xlat=',xlat(i)
!             stop
!           endif

!           if ( gq0(i,k,1) > 1.0 ) then
!             print *,' gq0=', gq0(i,k,1),' i=',i,' k=',k,' lat=',lat,  &
!    &          ' kdt=',kdt
!             stop
!           endif
!         enddo
!       enddo

!       if ( lprnt ) print *,' returning from ras for kdt=', kdt,       &
!    &                     ' me=',me,' lat=',lat

!  --- ...  end check print

        cld1d = 0

!  --- ...  update the tracers due to convective transport

        if ( tottracer > 0 ) then
          if ( ntoz > 0 ) then                       ! for ozone
            gq0(:,:,ntoz) = clw(:,:,3)

            if ( tracers > 0 ) then                  ! for other tracers
              do n = 1, tracers
                gq0(:,:,n+trc_shft) = clw(:,:,3+n)
              enddo
            endif
          else
            do n = 1, tracers
              gq0(:,:,n+trc_shft) = clw(:,:,2+n)
            enddo
          endif
        endif
      endif   ! end if_not_ras

!     if ( me == 0 ) then
!       do k = 1, levs
!         do i = 1, im
!           sumq = sumq                                                 &
!    &           + (gq0(i,k,1) + clw(i,k,1) + clw(i,k,2)) * del(i,k)
!         enddo
!       enddo

!       sumr = 0.0
!       do i = 1, im
!         sumr = sumr + rain1(i)
!       enddo
!       sumq = sumq * 1000.0 / con_g
!       sumr = sumr *1000

!       print *,' after ras: sumq=',sumq,' sumr=',sumr, ' kdt=',kdt
!     endif

      if ( lssav ) then
        do i = 1, im
          cldwrk(i) = cldwrk(i) + cld1d(i)*dtf
        enddo

        if ( ldiag3d ) then
          do k = 1, levs
            do i = 1, im
              tem = rcs(i) * frain
              dt3dt(i,k,4) = dt3dt(i,k,4) + (gt0(i,k) - dtdt(i,k))*frain
              dq3dt(i,k,2) = dq3dt(i,k,2) + (gq0(i,k,1)-dqdt(i,k,1))    &
     &                                                            *frain
              du3dt(i,k,3) = du3dt(i,k,3) + (gu0(i,k) - dudt(i,k))*tem
              dv3dt(i,k,3) = dv3dt(i,k,3) + (gv0(i,k) - dvdt(i,k))*tem

              upd_mf(i,k)  = upd_mf(i,k)  + ud_mf(i,k) * frain
              dwn_mf(i,k)  = dwn_mf(i,k)  + dd_mf(i,k) * frain
              det_mf(i,k)  = det_mf(i,k)  + dt_mf(i,k) * frain
            enddo
          enddo
        endif   ! end if_ldiag3d
      endif   ! end if_lssav

      do i = 1, im
        rainc(i) = frain * rain1(i)
      enddo

      if ( lssav ) then
        do i = 1, im
          bengsh(i) = bengsh(i) + rainc(i)
        enddo
      endif
!
      if ( cnvgwd ) then       !        call convective gravity wave drag

!  --- ...  calculate maximum convective heating rate            qmax [k/s]
!           cuhr = temperature change due to deep convection

        do i = 1, im
          qmax(i)   = 0.0
          cumabs(i) = 0.0
        enddo

        do k = 1, levs
          do i = 1, im
!           cuhr(i,k) = (gt0(i,k) - dtdt(i,k)) / dtf
            cuhr(i,k) = (gt0(i,k) - dtdt(i,k)) / dtp    ! moorthi

            cumchr(i,k) = 0.0
            gwdcu(i,k)  = 0.0
            gwdcv(i,k)  = 0.0
            diagn1(i,k) = 0.0
            diagn2(i,k) = 0.0

            if ( k>=kbot(i) .and. k<=ktop(i) ) then
              qmax(i)     = max( qmax(i), cuhr(i,k) )
              cumabs(i)   = cuhr(i,k) + cumabs(i)
            endif
          enddo
        enddo

        do i = 1, im
          do k = kbot(i), ktop(i)
            do k1 = kbot(i), k
              cumchr(i,k) = cuhr(i,k1) + cumchr(i,k)
            enddo

            cumchr(i,k) = cumchr(i,k) / cumabs(i)
          enddo
        enddo

!  --- ...  check print

        if ( lprnt ) then
          if ( kbot(ipr) <= ktop(ipr) ) then
            write(*,*) 'kbot <= ktop     for (lat,lon) = ',             &
     &            xlon(ipr)*57.29578,xlat(ipr)*57.29578
            write(*,*) 'kcnv kbot ktop qmax dlength  ',kcnv(ipr),       &
     &            kbot(ipr),ktop(ipr),(86400.*qmax(ipr)),dlength(ipr)
            write(*,9000) kdt
 9000       format(/,3x,'k',5x,'cuhr(k)',4x,'cumchr(k)',5x,             &
     &            'at kdt = ',i4,/)

            do k = ktop(ipr), kbot(ipr),-1
              write(*,9010) k,(86400.*cuhr(ipr,k)),(100.*cumchr(ipr,k))
 9010         format(2x,i2,2x,f8.2,5x,f6.0)
            enddo
          endif

          print *,' Before gwdc in gbphys fhour ',fhour

          if ( fhour >= fhourpr ) then
            print *,' before gwdc in gbphys start print'
            write(*,*) 'fhour ix im levs = ',fhour,ix,im,levs
            print *,'dtp  dtf  rcs = ',dtp,dtf,rcs(ipr)

            write(*,9100)
 9100       format(//,14x,'pressure levels',//                          &
     &             ' ilev',7x,'prsi',8x,'prsl',8x,'delp',/)

            k = levs + 1
            write(*,9110) k,(10.*prsi(ipr,k))
 9110       format(i4,2x,f10.3)

            do k = levs, 1, -1
              write(*,9120) k,(10.*prsl(ipr,k)),(10.*del(ipr,k))
              write(*,9110) k,(10.*prsi(ipr,k))
            enddo
 9120       format(i4,12x,2(2x,f10.3))

            write(*,9130)
 9130       format(//,10x,'before gwdc in gbphys',//,' ilev',6x,        &
     &             'ugrs',9x,'gu0',8x,'vrgs',9x,'gv0',8x,               &
     &             'tgrs',9x,'gt0',8x,'gt0b',8x,'dudt',8x,'dvdt',/)

            do k = levs, 1, -1
              write(*,9140) k,ugrs(ipr,k),gu0(ipr,k),                   &
     &                        vgrs(ipr,k),gv0(ipr,k),                   &
     &                        tgrs(ipr,k),gt0(ipr,k),dtdt(ipr,k),       &
     &                        dudt(ipr,k),dvdt(ipr,k)
            enddo
 9140       format(i4,9(2x,f10.3))

            print *,' before gwdc in gbphys end print'
          endif
        endif   ! end if_lprnt

!  --- ...  end check print

        call gwdc(im, ix, im, levs, lat, ugrs, vgrs, tgrs, qgrs,        &
     &          rcs, prsl, prsi, del, qmax, cumchr, ktop, kbot, kcnv,   &
     &          gwdcu, gwdcv,con_g,con_cp,con_rd,con_fvirt, dlength,    &
     &          lprnt, ipr, fhour,                                      &
     &          dusfcg,dvsfcg,diagn1,diagn2)

        if ( lprnt ) then
          if ( fhour >= fhourpr ) then
            print *,' after gwdc in gbphys start print'

            write(*,9131)
 9131       format(//,10x,'after gwdc in gbphys',//,' ilev',6x,         &
     &             'ugrs',9x,'gu0',8x,'vgrs',9x,'gv0',8x,               &
     &             'tgrs',9x,'gt0',8x,'gt0b',7x,'gwdcu',7x,'gwdcv',/)

            do k = levs, 1, -1
              write(*,9141) k,ugrs(ipr,k),gu0(ipr,k),                   &
     &                        vgrs(ipr,k),gv0(ipr,k),                   &
     &                        tgrs(ipr,k),gt0(ipr,k),dtdt(ipr,k),       &
     &                        gwdcu(ipr,k),gwdcv(ipr,k)
            enddo
 9141       format(i4,9(2x,f10.3))

!           print *,' after gwdc in gbphys end print'
          endif
        endif

!  --- ...  write out cloud top stress and wind tendencies

        if ( lssav ) then
          do i = 1, im
            dugwd(i) = dugwd(i) + dusfcg(i)*dtf
            dvgwd(i) = dvgwd(i) + dvsfcg(i)*dtf
          enddo

          if ( ldiag3d ) then
            do k = 1, levs
              do i = 1, im
                tem = rcs(i) * dtf
                du3dt(i,k,4) = du3dt(i,k,4) + gwdcu(i,k)*tem
                dv3dt(i,k,4) = dv3dt(i,k,4) + gwdcv(i,k)*tem
!               du3dt(i,k,2) = du3dt(i,k,2) + diagn1(i,k)*tem
!               dv3dt(i,k,2) = dv3dt(i,k,2) + diagn2(i,k)*tem
              enddo
            enddo
          endif
        endif   ! end if_lssav

!  --- ...  update the wind components with  gwdc tendencies

        do k = 1, levs
          do i = 1, im
            gu0(i,k) = gu0(i,k) + gwdcu(i,k)*dtp
            gv0(i,k) = gv0(i,k) + gwdcv(i,k)*dtp
          enddo
        enddo

        if ( lprnt ) then
          if ( fhour >= fhourpr ) then
            write(*,9132)
 9132       format(//,10x,'after tendency gwdc in gbphys',//,' ilev',6x,&
     &             'ugrs',9x,'gu0',8x,'vgrs',9x,'gv0',8x,               &
     &             'tgrs',9x,'gt0',8x,'gt0b',7x,'gwdcu',7x,'gwdcv',/)

            do k = levs, 1, -1
              write(*,9142) k,ugrs(ipr,k),gu0(ipr,k),vgrs(ipr,k),       &
     &              gv0(ipr,k),tgrs(ipr,k),gt0(ipr,k),dtdt(ipr,k),      &
     &              gwdcu(ipr,k),gwdcv(ipr,k)
            enddo
 9142       format(i4,9(2x,f10.3))
          endif
        endif

      endif   ! end if_cnvgwd (convective gravity wave drag)

      if ( ldiag3d ) then
        do k = 1, levs
          do i = 1, im
            dtdt(i,k)   = gt0(i,k)
            dqdt(i,k,1) = gq0(i,k,1)
          enddo
        enddo
      endif

      if ( .not. sashal ) then

!       if ( lprnt ) print *,' levshcm=',levshcm,' gt0sh=',gt0(ipr,:)

        if ( .not. mstrat ) then
          call shalcvt3(im,ix,levshcm,dtp,del,prsi,prsl,prslk,          &
     &                  kcnv,gq0,gt0)                             ! for pry
        else
          call shalcv(im,ix,levshcm,dtp,del,prsi,prsl,prslk,kcnv,       &
     &                gq0,gt0,levshc,phil,kinver,ctei_r,ctei_rm)

!         if ( lprnt ) print *,' levshcm=',levshcm,' gt0sha=',gt0(ipr,:)
        endif

      else

        call shalcnv(im,ix,levs,jcap,dtp,del,prsl,pgr,phil,             &
     &               clw,gq0,gt0,gu0,gv0,rcs,                           &
     &               rain1,kbot,ktop,kcnv,slmsk,                        &
     &               vvel,ncld,hpbl,hflx,evap,ud_mf,dt_mf)

        do i = 1, im
          raincs(i) = frain * rain1(i)
        enddo

        if ( lssav ) then
          do i = 1, im
            bengsh(i) = bengsh(i) + raincs(i)
          enddo
        endif

        do i = 1, im
          rainc(i) = rainc(i) + raincs(i)
        enddo

        if ( lssav ) then
          do i = 1, im
            cldwrk(i) = cldwrk(i) + cld1d(i)*dtf
          enddo
        endif

      endif   ! end if_not_sashal

      if ( lssav ) then
        if ( ldiag3d ) then
          do k = 1, levs
            do i = 1, im
              dt3dt(i,k,5) = dt3dt(i,k,5) + (gt0(i,k) - dtdt(i,k))*frain
              dq3dt(i,k,3) = dq3dt(i,k,3)                               &
     &                     + (gq0(i,k,1) - dqdt(i,k,1)) * frain
              upd_mf(i,k)  = upd_mf(i,k)  + ud_mf(i,k) * FRAIN
              det_mf(i,k)  = det_mf(i,k)  + dt_mf(i,k) * FRAIN
            enddo
          enddo
          do k=1,levs
            do i=1,im
              dtdt(i,k)   = GT0(i,k)
              dqdt(i,k,1) = gq0(i,k,1)
            enddo
          enddo
        endif
      endif   ! end if_lssav

      do k = 1, levs
        do i = 1, im
          if ( clw(i,k,2) <= -999.0 ) clw(i,k,2) = 0.0
        enddo
      enddo

      if ( ntcw > 0 ) then

        if ( num_p3d == 3 ) then    ! call brad ferrier's microphysics

!  --- ...  extract cloud water & ice from fc_ice

          do k = 1, levs
            do i = 1, im

!             qi = clw(i,k) * fc_ice(i,k)
!             qw = clw(i,k) - qi
              qi = clw(i,k,1)
              qw = clw(i,k,2)

!  --- ...  algorithm to combine different hydrometeor species

!             gq0(i,k,ntcw) = max( epsq, qi+qw+qr_col(i,k) )
              gq0(i,k,ntcw) = qi + qw + qr_col(i,k)

              if ( qi <= epsq ) then
                phy_f3d(i,k,1) = 0.0
              else
                phy_f3d(i,k,1) = qi / gq0(i,k,ntcw)
              endif

              if ( qr_col(i,k) <= epsq ) then
                phy_f3d(i,k,2) = 0.0
              else
                phy_f3d(i,k,2) = qr_col(i,k) / (qw+qr_col(i,k))
              endif

            enddo
          enddo

        else    ! if_num_p3d

          do k = 1, levs
            do i = 1, im
!             gq0(i,k,ntcw) = clw(i,k) + gq0(i,k,ntcw)
              gq0(i,k,ntcw) = clw(i,k,1) + clw(i,k,2)
!             gq0(i,k,ntcw) = clw(i,k,1)               ! for pry
            enddo
          enddo

        endif   ! end if_num_p3d

      else    ! if_ntcw

        do k = 1, levs
          do i = 1, im
            clw(i,k,1) = clw(i,k,1) + clw(i,k,2)
          enddo
        enddo

      endif   ! end if_ntcw

      call cnvc90(clstp, im,   ix,   rainc, kbot, ktop, levs, prsi,     &
     &            acv,   acvb, acvt, cv,    cvb,  cvt)

      if ( ncld == 0 ) then

        call lrgscl(ix,im,levs,dtp,gt0,gq0,prsl,del,prslk,rain1,clw)

      elseif ( ncld == 1 ) then

        if ( num_p3d == 3 ) then    ! call brad ferrier's microphysics

          do i = 1, im
            xncw(i) = ncw(2) * work1(i) + ncw(1) * work2(i)
          enddo

          if ( kdt==1 .and. abs(xlon(1))<0.0001 ) then
            print *,' xncw=',xncw(1),' rhc=',rhc(1,1),' work1=',work1(1)&
     &,         ' work2=',work2(1),' flgmin=',flgmin_l(1)               &
     &,         ' lon=',xlon(1) * 57.29578,' lat=',lat,' me=',me
!    &,         ' lon=',xlon(1) * 57.29578,' lat=',xlat(1) * 57.29578
!    &,         ' kinver=',kinver(1)
          endif

!         if ( lprnt ) print *,' ipr=',ipr,' gt0_gsmb=',gt0(ipr,:)      &
!    &,         ' xlon=',xlon(ipr),' xlat=',xlat(ipr)

          call gsmdrive(im, ix, levs, dtp, prsl, del                    &
     &,                 gt0, gq0(1,1,1), gq0(1,1,ntcw), slmsk           &
     &,                 phy_f3d(1,1,1),  phy_f3d(1,1,2)                 &
     &,                 phy_f3d(1,1,3), rain1, sr, con_g                &
     &,                 con_hvap, hsub,con_cp, rhc, xncw, flgmin_l      &
     &,                 me, lprnt, ipr)

!         if ( lprnt ) print *,' ipr=',ipr,' gt0_gsma=',gt0(ipr,:)

        elseif ( num_p3d == 4 ) then ! call zhao/carr/sundqvist microphysics

          call gscond(im, ix, levs, dtp, prsl, pgr,                     &
     &                gq0(1,1,1), gq0(1,1,ntcw), gt0,                   &
     &                phy_f3d(1,1,1), phy_f3d(1,1,2), phy_f2d(1,1),     &
     &                phy_f3d(1,1,3), phy_f3d(1,1,4), phy_f2d(1,2),     &
     &                rhc,lprnt, ipr)

          call precpd(im, ix, levs, dtp, del, prsl, pgr,                &
     &                gq0(1,1,1), gq0(1,1,ntcw), gt0, rain1,            &
     &                 rhc, lprnt, ipr)

        endif   ! end if_num_p3d

      endif   ! end if_ncld

!     if ( lprnt ) print *,' rain1=',rain1(ipr),' rainc=',rainc(ipr)

      do i = 1, im
        rainl(i) = frain * rain1(i)
        rain(i) = rainc(i) + rainl(i)
      enddo

!  --- ...  coupling insertion

      if ( lssav_cc ) then
        do i = 1, im
          precr_cc(i) = precr_cc(i) + rain(i)
        enddo
      endif

      if ( lssav ) then
        do i = 1, im
          geshem(i) = geshem(i) + rain(i)
        enddo

        if ( ldiag3d ) then
          do k = 1, levs
            do i = 1, im
              dt3dt(i,k,6) = dt3dt(i,k,6) + (gt0(i,k) - dtdt(i,k))*frain
              dq3dt(i,k,4) = dq3dt(i,k,4)                               &
     &                     + (gq0(i,k,1) - dqdt(i,k,1)) * frain
            enddo
          enddo
        endif
      endif   ! end if_lssav

!  --- ...  estimate t850 for rain-snow decision

      do i = 1, im
        t850(i) = gt0(i,1)
      enddo

      do k = 1, levs - 1
        do i = 1, im
          if ( prsl(i,k)>p850 .and. prsl(i,k+1)<=p850 ) then
            t850(i) = gt0(i,k) - (prsl(i,k) - p850)                     &
     &              / (prsl(i,k)-prsl(i,k+1)) * (gt0(i,k)-gt0(i,k+1))
          endif
        enddo
      enddo

!  --- ...  lu: snow-rain detection is performed in land/sice module
!           factor=weighted mean tep.

      do i = 1, im
        tprcp(i) = rain(i)               ! clu: rain -> tprcp
        srflag(i) = 0.0                  ! clu: default srflag as 'rain'
        if ( t850(i) <= 273.16 ) then
          srflag(i) = 1.0                ! clu: set srflag to 'snow'

!  --- ...  wei: when call osu lsm, neutral impact, for code consistency

          if ( lsm==0 .and. slmsk(i)/=0.0 ) then
            sheleg(i) = sheleg(i) + 1.0e3*rain(i)  ! neutral impact, for code consistency
            tprcp(i)  = 0.0
          endif
        endif
      enddo

!     if ( lprnt ) print *,' tprcp=',tprcp(ipr),' rain=',rain(ipr)

!  --- ...  cpl insertion

!     snw_cc = 0.
      do i = 1, im
!       if ( t850(i)<=273.16 .and. slmsk(i)/=0.0 ) then
        if ( t850(i) <= 273.16 ) then
          lprec_cc(i) = 0.0
          snw_cc(i) = rain(i)
        else
          lprec_cc(i) = rain(i)
          snw_cc(i) = 0.0
        endif
      enddo

!  --- ...  wei: when call osu lsm
!           update soil moisture and canopy water after precipitation computaion

      if ( lsm == 0 ) then

        call progt2(im,lsoil,rhscnpy,rhsmc,ai,bi,cci,smsoil,            &
     &              slmsk,canopy,tprcp,runof,snowmt,                    &
     &              zsoil,soiltyp,sigmaf,dtf,me)

!  --- ...  wei: let soil liquid water equal to soil total water

        do k = 1, lsoil
          do i = 1, im
            if ( slmsk(i) == 1.0 ) then
              slsoil(i,k) = smsoil(i,k)
            endif
          enddo
        enddo

      endif   ! end if_lsm

!  --- ...  total runoff is composed of drainage into water table and
!           runoff at the surface and is accumulated in unit of meters

      if ( lssav ) then
        do i = 1, im
          runoff(i) = runoff(i) + (drain(i) + runof(i)) * dtf / 1000.0
          srunoff(i)= srunoff(i) + runof(i)*dtf / 1000.0
        enddo
      endif

!  --- ...  xw: return updated ice thickness & concentration to global array

      do i = 1, im
        if ( slmsk(i) == 2.0 ) then
          hice(i)  = zice(i)
          fice(i)  = cice(i)
          tisfc(i) = tice(i)
        else
          hice(i)  = 0.0
          fice(i)  = 0.0
          tisfc(i) = tsea(i)
        endif
      enddo

!  --- ...  return updated smsoil and stsoil to global arrays

      do k = 1, lsoil
        do i = 1, im
          smc(i,k) = smsoil(i,k)
          stc(i,k) = stsoil(i,k)
          slc(i,k) = slsoil(i,k)
        enddo
      enddo

!  --- ...  calc. integral of moistue in pwat

      do i = 1, im
        pwat(i) = 0.0
      enddo

      do k = 1, levs
        do i = 1, im
          work1(i) = 0.0
        enddo

        if ( ncld > 0 ) then
          do ic = ntcw, ntcw+ncld-1
            do i = 1, im
!             work1(i) = work1(i) + max( gq0(i,k,ic), qmin )
              work1(i) = work1(i) + gq0(i,k,ic)
            enddo
          enddo
        endif

        do i = 1, im
!  --- ...  inside this routine, let t as dry temperature only   ! hmhj
!         work2(i) = 1.0 + con_fvirt * max( gq0(i,k,1),qmin )    ! hmhj
!         gt0(i,k) = gt0(i,k) * work2(i)                         ! hmhj
          pwat(i)  = pwat(i)  + del(i,k)*(gq0(i,k,1) + work1(i))
        enddo
      enddo

      do i = 1, im
        pwat(i)  = pwat(i) * (1.0e3/con_g)
      enddo
!
      deallocate (clw)

!     if ( lprnt ) call mpi_quit(7)
!     if ( kdt > 1 )  call mpi_quit(7)

      return
!...................................
      end subroutine gbphys
!-----------------------------------

