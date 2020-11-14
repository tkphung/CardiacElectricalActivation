function [Q] = LV_ElementtoGlobal(n)
%LV_ElementtoGlobal(n) takes nodes for an element (each row is a node and
%each column is X Y Z)
%   From SAMANTHA CLARK, adapted by TNPhung to add extra comments
%   INPUT: n (nodes matrix for single element)
%   OUTPUT: Q is the transformation matrix to rotate a vector v(3,1) into 
%   element coordinate system
%       v_rot = Q*v

n = n(:,1:3);
% radial unit vector originating at node 4
a = (n(3,:)-n(4,:))'; b = (n(8,:)-n(4,:))';
c = cross(a,b);
er = c/(c(1)^2+c(2)^2+c(3)^2)^0.5;

% local longitudinal vector
qp = [-1 0 0]'; % CHANGED by TNP. This should be the Apex to Base vector
nnorm = [n(4,:)/((n(4,1)^2+n(4,2)^2+n(4,3)^2)^0.5)]';
m = cross(qp,nnorm); % unit vector perpendicular to bisecting plane
mnorm = m/((m(1)^2+m(2)^2+m(3)^2)^0.5);
d = cross(er,mnorm);
el = d/((d(1)^2+d(2)^2+d(3)^2)^0.5);

% circumferential unit vector originating at node 4
% should be perpendicular to both er and el
f = cross(er,el);
ec = f/((f(1)^2+f(2)^2+f(3)^2)^0.5);

Q = -1.*[el ec er];
end

