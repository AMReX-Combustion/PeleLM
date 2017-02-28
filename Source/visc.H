      logical use_constant_mu, use_constant_lambda, use_constant_rhoD
      REAL_T  constant_mu_val, constant_lambda_val, constant_rhoD_val,
     &     Pr, Sc, thickFacTR
      logical LeEQ1
      common /xxvisc1/ use_constant_mu, use_constant_lambda,
     &     use_constant_rhoD
      common /xxvisc2/ constant_mu_val, constant_lambda_val,
     &     constant_rhoD_val, Pr, Sc, thickFacTR
      common /xxvisc3/ LeEQ1
      save /xxvisc1/,/xxvisc2/,/xxvisc3/
      

