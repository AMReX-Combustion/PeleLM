% epsilon=0.6
Th = 1200
Tc = 300

mu = 1.68e-5
Pr = 0.71
T0 = 600
g = 9.81

P0 = 101325
R = 287
rho = P0/(R*T0)  %0.585981 %

Ra=1e6


Lcube = (Ra*T0*mu*mu)/(Pr*g*rho*rho*(Th-Tc))

L=nthroot(Lcube,3)






epsilon = 0.6
Tc = 300
Th = Tc*((1+epsilon)/(1-epsilon))