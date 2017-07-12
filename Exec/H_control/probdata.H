!
!     This means that FORT_ACTIVECONTROL is callable.
!     Otherwise FORT_ACTIVECONTROL will abort if it is called.
!
#define ACTIVE_CONTROL_IS_USABLE
      
      REAL_T vorterr, temperr, adverr, tempgrad, flametracval, twall
      REAL_T splitx, splity, traceSpecVal, xfrontw, yfrontw
      REAL_T standoff
      REAL_T domnlo(3), domnhi(3)
      
      integer max_vort_lev, max_temp_lev, max_trac_lev, max_nozzle_lev
      integer fuelID, oxidID, prodID, nspecies
      character*50 probtype

      common /prob1/ vorterr, temperr, adverr, tempgrad, flametracval, twall,
     &               splitx, splity, traceSpecVal, xfrontw, yfrontw, standoff
      common /prob2/ domnlo, domnhi
      common /prob3/ max_vort_lev, max_temp_lev, max_trac_lev,
     &               max_nozzle_lev
      common /prob4/ fuelID, oxidID, prodID, nspecies
      common /prob5/ probtype
      save /prob1/, /prob2/, /prob3/, /prob4/, /prob5/

      integer refine_nozzle
      REAL_T refine_nozzle_x, refine_nozzle_y, refine_nozzle_z, blobx, bloby, blobz, blobr, blobT, Tfrontw, xcen
      common /myprob1/ refine_nozzle
      common /myprob2/ refine_nozzle_x, refine_nozzle_y, refine_nozzle_z, blobx, bloby, blobz, 
     &     blobr, blobT, Tfrontw, xcen

      REAL_T v_strength,v_width,v_xcen,v_ycen, v_cl_x
      common /myprob3/ v_strength,v_width,v_xcen,v_ycen,v_cl_x

      REAL_T v_blob_r, v_blob_T, v_blob_airfrac, turb_scale,phi_in
      common /myprob4/ v_blob_r, v_blob_T, v_blob_airfrac, turb_scale,phi_in

      REAL_T stTh, Rf, V_in, V_co, T_in, T_co, R_hot, R_hotBL, stBL, pipeTh, pipeBL
      REAL_T tV_in_l, tV_in_r, V_in_new, tV_co_l, tV_co_r, V_co_new, H2_frac
      REAL_T T_stick, rho_stick, h_stick, fuel_N2_vol_percent,T_switch,V_jet
      common /myprob5/ stTh, Rf, V_in, V_co, T_in, T_co, R_hot, R_hotBL, stBL,
     &     pipeTh, pipeBL, tV_in_l, tV_in_r, V_in_new, tV_co_l, tV_co_r,
     &     V_co_new, T_stick, rho_stick, h_stick, fuel_N2_vol_percent,
     &     H2_frac,T_switch,V_jet
      save /myprob1/, /myprob2/, /myprob3/, /myprob4/, /myprob5/ 

      REAL_T pertmag,  pert_scale
      common /pert1/ pertmag, pert_scale
      save /pert1/

      integer dimFile
      logical forceInflow
      REAL_T convVelInit, probSizeFile(3), dxFile(3)
      common /inflow1/ convVelInit, forceInflow
      common /inflow2/ probSizeFile, dxFile
      common /inflow3/ dimFile
      save /inflow1/, /inflow2/, /inflow3/
