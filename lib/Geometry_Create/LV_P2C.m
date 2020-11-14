function [x1,x2,x3]=LV_P2C(z1,z2,z3,focus)
% Adapted by Thien-Khoi N. Phung from SLH_CMI_P2C.m (Feb 16, 2017)
% Transformation from prolate spheroid to 
% rectangular Cartesian coordinates.  
% z1 = lambda
% z2 = mu
% z3 = theta
% ALTERED by Thien-Khoi N. Phung to allow for vectors or matrices of points
% to be passed through. Feb. 23, 2017 (I just added the dots)

x1=focus.*cosh(z1).*cos(z2);
x2=focus.*sinh(z1).*sin(z2).*cos(z3);
x3=focus.*sinh(z1).*sin(z2).*sin(z3);
