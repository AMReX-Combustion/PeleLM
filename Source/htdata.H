      REAL_T pamb, dpdt_factor
      integer closed_chamber
      common /htdata1/ pamb, dpdt_factor, closed_chamber
      
      integer Density, Temp, RhoH, Trac, FirstSpec, LastSpec
      common /scalID/ Density, Temp, RhoH, Trac, FirstSpec, LastSpec

      save /htdata1/,/scalID/
