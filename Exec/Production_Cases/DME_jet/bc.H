c     Note: Requires ckdwrk.H to obtain the value of maxspec
      
      logical bcinit
      common / stats/ bcinit

      integer Nzones
      parameter (Nzones=5)      

      REAL_T u_bc(Nzones), v_bc(Nzones), w_bc(Nzones), rho_bc(Nzones)
      REAL_T Y_bc(0:maxspec-1, Nzones), T_bc(Nzones)
      REAL_T h_bc(Nzones)

      integer MAXPNTS
      parameter(MAXPNTS = 50)
      REAL_T time_points(0:MAXPNTS),vel_points(0:MAXPNTS),cntl_points(0:MAXPNTS)

      common / bcs / u_bc, v_bc, w_bc, rho_bc, Y_bc, T_bc, h_bc

      character(50) ac_hist_file
      REAL_T tau_control, cfix, coft_old, sest, V_in_old, corr,
     &     changeMax_control, tbase_control, dV_control, scale_control,
     &     zbase_control, h_control, controlVelMax
      integer navg_pnts

      common / control1 / tau_control, cfix, coft_old, sest, V_in_old,
     &     corr, changeMax_control, tbase_control, dV_control, scale_control,
     &     zbase_control, h_control, controlVelMax,
     &     time_points,vel_points,cntl_points,navg_pnts,ac_hist_file

      integer pseudo_gravity
      common / control2 / pseudo_gravity

      save /bcs/, /control1/, /control2/
