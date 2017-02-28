      !
      ! Hardwired sizes.
      !
      integer, parameter :: maxreac  = 472
      integer, parameter :: maxspec  = 72
      integer, parameter :: maxelts  = 6
      integer, parameter :: maxthrdb = 20
      integer, parameter :: maxtp    = 3
      integer, parameter :: maxsp    = 12
      integer, parameter :: maxspnml = 16
      !
      ! The only ChemKin stuff we need.
      !
      integer, parameter :: ckbi = 1
      integer, parameter :: ckbr = 1
      !
      ! EGLib package for transport properties.
      !
      ! For normal LMC usage:
      !
      integer, parameter :: eg_nodes=1, eg_IFLAG=3, eg_ITLS=1
      !
      ! If you're using the DD* routines:
      !
      ! integer, parameter :: eg_nodes=1, eg_IFLAG=7, eg_ITLS=3
      !
      ! DVODE & DVODE driver workspace requirements.
      !
      integer, parameter :: dvbr = 1
      integer, parameter :: dvbi = 1

      integer dvr, dvi, dvdr, dvdbr, dvder
      common / wrkidx / dvr, dvi, dvdr, dvdbr, dvder
      save   / wrkidx /

      integer,          pointer :: IWRK(:)
      double precision, pointer :: RWRK(:)
      common / cdwork / RWRK, IWRK
      save   / cdwork /

!$omp threadprivate(/CDWORK/)
      !
      ! EGlib workspace -- allocated once in FORT_INITCHEM()
      !
      integer,          pointer :: EGIWRK(:)
      double precision, pointer :: EGRWRK(:)
      common / EGWORK / EGRWRK, EGIWRK
      save   / EGWORK /
      !
      ! Tranlib workspace -- allocated once in FORT_INITCHEM()
      !
      integer,          pointer :: MCIWRK(:)
      double precision, pointer :: MCRWRK(:)
      common / MCWORK / MCRWRK, MCIWRK
      save   / MCWORK /

!$omp threadprivate(/EGWORK/,/MCWORK/)
      !
      ! Actual sizes
      !
      integer            Nelt, Nspec, Nreac, Nfit
      common / ckdsize / Nelt, Nspec, Nreac, Nfit
      save   / ckdsize /
      !
      ! I/O
      !
      integer          LLINKMC
      common / ckdio / LLINKMC
      save   / ckdio /

      integer verbose_vode, max_vode_subcycles, vode_itol
      common / ckdio1 / verbose_vode, max_vode_subcycles, vode_itol
      save   / ckdio1 /
      !
      ! Integrator
      !
      double precision  spec_scalY(maxspec), vode_rtol, vode_atol
      common /vode1/ spec_scalY, vode_rtol, vode_atol
      save /vode1/

      double precision  thickFacCH
      common / vode2 / thickFacCH
      save   / vode2 /
      !            
      ! Transport library
      !
      double precision  TMIN_TRANS
      common / tran2 / TMIN_TRANS
      save   / tran2 /

      integer Ncoefs,Doffset,PTCoffset,DMIXoffset,TDoffset,use_mc,use_eg
      common / tran3 / Ncoefs,Doffset,PTCoffset,DMIXoffset,TDoffset,
     &     use_mc,use_eg
      save   / tran3 /

#ifdef LMC_SDC
      ! SDC additions
      integer negative_Y_test
      double precision c_0(maxspec+1), c_1(maxspec+1), rhoH_INIT, T_cell, rhoY_INIT(maxspec)
      common /sdc_function/ c_0, c_1, rhoH_INIT, T_cell, rhoY_INIT, negative_Y_test
      save /sdc_function/
!$omp threadprivate(/sdc_function/)
#endif
