Convected Vortex (CoVo) - README
================================

Description
-----------

The CoVo case is a classical test case for computational fluid dynamics solvers. It allows to directly evaluate the numerical characteristics of the solver by comparison with analytical solution. It consists in convecting an isentropic vortex across a periodic 2D box. The analytic solution of the isentropic vortex is given by the stream function:

<a href="https://www.codecogs.com/eqnedit.php?latex=\large&space;\psi(x,y)&space;=&space;\psi_0&space;e^{&space;-&space;{(x-x_c)^2&plus;(y-y_c)^2&space;\over&space;2&space;r_c^2}&space;}=\psi_0&space;e^{-r^2/2r_c^2}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\large&space;\psi(x,y)&space;=&space;\psi_0&space;e^{&space;-&space;{(x-x_c)^2&plus;(y-y_c)^2&space;\over&space;2&space;r_c^2}&space;}=\psi_0&space;e^{-r^2/2r_c^2}" title="\large \psi(x,y) = \psi_0 e^{ - {(x-x_c)^2+(y-y_c)^2 \over 2 r_c^2} }=\psi_0 e^{-r^2/2r_c^2}" /></a>


Current implementation
----------------------


