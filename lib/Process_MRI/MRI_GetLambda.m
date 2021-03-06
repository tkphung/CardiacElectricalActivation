function L=MRI_GetLambda(C,E)
% Adapted by Thien-Khoi N. Phung from SLH_CMI_GetLambda.m Feb. 23, 2017
% This function actually computes the lambda values and derivatives 
% according to the cubic Hermite function.
%
H00(1)=1-3*(E(1)^2)+2*(E(1)^3);
H00(2)=1-3*(E(2)^2)+2*(E(2)^3);
H10(1)=E(1)*((E(1)-1)^2);
H10(2)=E(2)*((E(2)-1)^2);
H01(1)=(E(1)^2)*(3-2*E(1));
H01(2)=(E(2)^2)*(3-2*E(2));
H11(1)=(E(1)^2)*(E(1)-1);
H11(2)=(E(2)^2)*(E(2)-1);

dH00(1)=6*((E(1)^2)-E(1));
dH00(2)=6*((E(2)^2)-E(2));
dH10(1)=3*(E(1)^2)-4*E(1)+1;
dH10(2)=3*(E(2)^2)-4*E(2)+1;
dH01(1)=6*(E(1)-E(1)^2);
dH01(2)=6*(E(2)-E(2)^2);
dH11(1)=3*(E(1)^2)-2*E(1);
dH11(2)=3*(E(2)^2)-2*E(2);

L1=H00(1)*H00(2)*C(1,1)+H01(1)*H00(2)*C(2,1)+H00(1)*H01(2)*C(3,1)+H01(1)*H01(2)*C(4,1);
L2=H10(1)*H00(2)*C(1,4)+H11(1)*H00(2)*C(2,4)+H10(1)*H01(2)*C(3,4)+H11(1)*H01(2)*C(4,4);
L3=H00(1)*H10(2)*C(1,5)+H01(1)*H10(2)*C(2,5)+H00(1)*H11(2)*C(3,5)+H01(1)*H11(2)*C(4,5);
L4=H10(1)*H10(2)*C(1,6)+H11(1)*H10(2)*C(2,6)+H10(1)*H11(2)*C(3,6)+H11(1)*H11(2)*C(4,6);
L(:,:,1)=L1+L2+L3+L4;

dL1wrt1=dH00(1)*H00(2)*C(1,1)+dH01(1)*H00(2)*C(2,1)+dH00(1)*H01(2)*C(3,1)+dH01(1)*H01(2)*C(4,1);
dL2wrt1=dH10(1)*H00(2)*C(1,4)+dH11(1)*H00(2)*C(2,4)+dH10(1)*H01(2)*C(3,4)+dH11(1)*H01(2)*C(4,4);
dL3wrt1=dH00(1)*H10(2)*C(1,5)+dH01(1)*H10(2)*C(2,5)+dH00(1)*H11(2)*C(3,5)+dH01(1)*H11(2)*C(4,5);
dL4wrt1=dH10(1)*H10(2)*C(1,6)+dH11(1)*H10(2)*C(2,6)+dH10(1)*H11(2)*C(3,6)+dH11(1)*H11(2)*C(4,6);
L(:,:,2)=dL1wrt1+dL2wrt1+dL3wrt1+dL4wrt1;

dL1wrt2=H00(1)*dH00(2)*C(1,1)+H01(1)*dH00(2)*C(2,1)+H00(1)*dH01(2)*C(3,1)+H01(1)*dH01(2)*C(4,1);
dL2wrt2=H10(1)*dH00(2)*C(1,4)+H11(1)*dH00(2)*C(2,4)+H10(1)*dH01(2)*C(3,4)+H11(1)*dH01(2)*C(4,4);
dL3wrt2=H00(1)*dH10(2)*C(1,5)+H01(1)*dH10(2)*C(2,5)+H00(1)*dH11(2)*C(3,5)+H01(1)*dH11(2)*C(4,5);
dL4wrt2=H10(1)*dH10(2)*C(1,6)+H11(1)*dH10(2)*C(2,6)+H10(1)*dH11(2)*C(3,6)+H11(1)*dH11(2)*C(4,6);
L(:,:,3)=dL1wrt2+dL2wrt2+dL3wrt2+dL4wrt2;

dL1wrt12=dH00(1)*dH00(2)*C(1,1)+dH01(1)*dH00(2)*C(2,1)+dH00(1)*dH01(2)*C(3,1)+dH01(1)*dH01(2)*C(4,1);
dL2wrt12=dH10(1)*dH00(2)*C(1,4)+dH11(1)*dH00(2)*C(2,4)+dH10(1)*dH01(2)*C(3,4)+dH11(1)*dH01(2)*C(4,4);
dL3wrt12=dH00(1)*dH10(2)*C(1,5)+dH01(1)*dH10(2)*C(2,5)+dH00(1)*dH11(2)*C(3,5)+dH01(1)*dH11(2)*C(4,5);
dL4wrt12=dH10(1)*dH10(2)*C(1,6)+dH11(1)*dH10(2)*C(2,6)+dH10(1)*dH11(2)*C(3,6)+dH11(1)*dH11(2)*C(4,6);
L(:,:,4)=dL1wrt12+dL2wrt12+dL3wrt12+dL4wrt12;
