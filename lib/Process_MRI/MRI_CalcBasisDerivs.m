function H=MRI_CalcBasisDerivs(E,deriv_num,order)
% Adapted from SLH_CMI_CalcBasisDerivs.m by TK Phung (Feb 14, 2017)
% Author: CM Ingrassia
%
% This function computes the derivatives of the bicubic Hermite and Lagrange
% basis function coefficients at E.  Order is the order of the time polynomial.
% The coefficient derivatives are returned in H.  H potentially has 3 rows:
% row 1 for either space-only coefficients (if order=0) or space-time
% coefficients with derivatives in space only (if order>0); row 2 for
% first-order time derivatives (order>0); row 3 for second-order time
% derivatives (order>1).
%
% deriv_num is between 0 and 5:  0 for no derivative, 1 (wrt 1), 2 (wrt 11), 
% 3 (wrt 2), 4 (wrt 22), 5 (wrt 12).  This specifies the SPATIAL
% derivatives.  If order is 0, just these coefficients are specified:  TIME IS
% NOT USED AT ALL.  However, if order>0/(1), first/(and second) time derivative 
% coefficients are combined with the spatial derivative coefficients.  To obtain 
% derivatives wrt time only, use deriv_num=0 and order>0.  To obtain
% derivatives wrt space only, but combined with 0-order time coefficients,
% use deriv_num>0 and order>0.  0-order time derivatives are output in
% addition to velocity and acceleration.
%

% 1D shape functions.
H00(1)=1-3*(E(1)^2)+2*(E(1)^3);
H00(2)=1-3*(E(2)^2)+2*(E(2)^3);
H10(1)=E(1)*((E(1)-1)^2);
H10(2)=E(2)*((E(2)-1)^2);
H01(1)=(E(1)^2)*(3-2*E(1));
H01(2)=(E(2)^2)*(3-2*E(2));
H11(1)=(E(1)^2)*(E(1)-1);
H11(2)=(E(2)^2)*(E(2)-1);

% First derivatives.
dH00(1)=-6*E(1)+6*(E(1)^2);
dH00(2)=-6*E(2)+6*(E(2)^2);
dH10(1)=2*E(1)*(E(1)-1)+(E(1)-1)^2;
dH10(2)=2*E(2)*(E(2)-1)+(E(2)-1)^2;
dH01(1)=-2*(E(1)^2)+2*E(1)*(3-2*E(1));
dH01(2)=-2*(E(2)^2)+2*E(2)*(3-2*E(2));
dH11(1)=(E(1)^2)+2*E(1)*(E(1)-1);
dH11(2)=(E(2)^2)+2*E(2)*(E(2)-1);

% Second derivatives.
d2H00(1)=12*E(1)-6;
d2H00(2)=12*E(2)-6;
d2H10(1)=6*E(1)-4;
d2H10(2)=6*E(2)-4;
d2H01(1)=-12*E(1)+6;
d2H01(2)=-12*E(2)+6;
d2H11(1)=6*E(1)-2;
d2H11(2)=6*E(2)-2;

% Assemble spatial derivative coefficients.
if deriv_num==0 % No derivative
    H_init=[H00(1)*H00(2) H01(1)*H00(2) H00(1)*H01(2) H01(1)*H01(2) H10(1)*H00(2) H11(1)*H00(2) H10(1)*H01(2) H11(1)*H01(2) H00(1)*H10(2) H01(1)*H10(2) H00(1)*H11(2) H01(1)*H11(2) H10(1)*H10(2) H11(1)*H10(2) H10(1)*H11(2) H11(1)*H11(2)];    
elseif deriv_num==1 % wrt 1
    H_init=[dH00(1)*H00(2) dH01(1)*H00(2) dH00(1)*H01(2) dH01(1)*H01(2) dH10(1)*H00(2) dH11(1)*H00(2) dH10(1)*H01(2) dH11(1)*H01(2) dH00(1)*H10(2) dH01(1)*H10(2) dH00(1)*H11(2) dH01(1)*H11(2) dH10(1)*H10(2) dH11(1)*H10(2) dH10(1)*H11(2) dH11(1)*H11(2)];    
elseif deriv_num==2 % wrt 11
    H_init=[d2H00(1)*H00(2) d2H01(1)*H00(2) d2H00(1)*H01(2) d2H01(1)*H01(2) d2H10(1)*H00(2) d2H11(1)*H00(2) d2H10(1)*H01(2) d2H11(1)*H01(2) d2H00(1)*H10(2) d2H01(1)*H10(2) d2H00(1)*H11(2) d2H01(1)*H11(2) d2H10(1)*H10(2) d2H11(1)*H10(2) d2H10(1)*H11(2) d2H11(1)*H11(2)];
elseif deriv_num==3 % wrt 2
    H_init=[H00(1)*dH00(2) H01(1)*dH00(2) H00(1)*dH01(2) H01(1)*dH01(2) H10(1)*dH00(2) H11(1)*dH00(2) H10(1)*dH01(2) H11(1)*dH01(2) H00(1)*dH10(2) H01(1)*dH10(2) H00(1)*dH11(2) H01(1)*dH11(2) H10(1)*dH10(2) H11(1)*dH10(2) H10(1)*dH11(2) H11(1)*dH11(2)];
elseif deriv_num==4 % wrt 22
    H_init=[H00(1)*d2H00(2) H01(1)*d2H00(2) H00(1)*d2H01(2) H01(1)*d2H01(2) H10(1)*d2H00(2) H11(1)*d2H00(2) H10(1)*d2H01(2) H11(1)*d2H01(2) H00(1)*d2H10(2) H01(1)*d2H10(2) H00(1)*d2H11(2) H01(1)*d2H11(2) H10(1)*d2H10(2) H11(1)*d2H10(2) H10(1)*d2H11(2) H11(1)*d2H11(2)];
elseif deriv_num==5 % wrt 12
    H_init=[dH00(1)*dH00(2) dH01(1)*dH00(2) dH00(1)*dH01(2) dH01(1)*dH01(2) dH10(1)*dH00(2) dH11(1)*dH00(2) dH10(1)*dH01(2) dH11(1)*dH01(2) dH00(1)*dH10(2) dH01(1)*dH10(2) dH00(1)*dH11(2) dH01(1)*dH11(2) dH10(1)*dH10(2) dH11(1)*dH10(2) dH10(1)*dH11(2) dH11(1)*dH11(2)];
end

% Take time derivatives, if any.
if order==0 % SPACE ONLY, static problem
    H=H_init;
else % At least linear Lagrange, output 0-order and first-order time derivatives.
    %
    % 0-order (no) time derivative.
    %
    
    % Normalized time coordinates for an element of this order.
    t=[0:1/order:1];
    num_t=length(t);
    
    % Compute Lagrange polynomial, 0-order time derivative.
    for i=1:num_t
        LG(i)=1;
        for j=1:num_t
            if j~=i
                LG(i)=LG(i)*(E(3)-t(j))/(t(i)-t(j));
            end
        end
    end
    
    % IMPORTANT: Rearrange to put ELEMENT nodes in the first 2 positions.
    if order>1
        temp=LG(2:num_t-1);
        LG(2)=LG(num_t);
        LG(3:num_t)=temp;
    end
    
    % Assemble first row of H (0-order time derivative).
    H_temp=[];
    for i=1:num_t
        H_temp=[H_temp LG(i).*H_init];
    end
    H(1,:)=H_temp;
    
    % 
    % First-order time derivative (velocity).
    %
    
    % Compute first order time derivatives (velocity).
    for i=1:num_t
        vel(i)=0;
        for j=1:num_t
            if j~=i
                vel_temp=1;
                for k=1:num_t
                    if k~=j & k~=i
                        vel_temp=vel_temp*(E(3)-t(k));
                    end
                end
                vel(i)=vel(i)+vel_temp;
            end
        end
        for m=1:num_t
            if m~=i
                vel(i)=vel(i)/(t(i)-t(m));
            end
        end
    end
    
    % IMPORTANT: Rearrange to put ELEMENT nodes in the first 2 positions.
    if order>1
        temp=vel(2:num_t-1);
        vel(2)=vel(num_t);
        vel(3:num_t)=temp;
    end
    
    % Assemble second row of H (first-order time derivative).
    H_temp=[];
    for i=1:num_t
        H_temp=[H_temp vel(i).*H_init];
    end
    H(2,:)=H_temp;
    
    %
    % Second-order time derivative (acceleration), if applicable.
    %
        
    if order>1 % At least quadratic Lagrange, output 0-order, first-order, and second-order time derivatives.
        % Compute first order time derivatives (velocity).
        for i=1:num_t
            acc(i)=0;
            for j=1:num_t
                if j~=i
                    for k=1:num_t
                        if k~=j & k~=i
                            acc_temp=1;
                            for m=1:num_t
                                if m~=k & m~=j & m~=i
                                    acc_temp=acc_temp*(E(3)-t(m));
                                end
                            end
                            acc(i)=acc(i)+acc_temp;
                        end
                    end
                end
            end
            for n=1:num_t
                if n~=i
                    acc(i)=acc(i)/(t(i)-t(n));
                end
            end
        end
        
        % IMPORTANT: Rearrange to put ELEMENT nodes in the first 2 positions.
        if order>1
            temp=acc(2:num_t-1);
            acc(2)=acc(num_t);
            acc(3:num_t)=temp;
        end
        
        % Assemble third row of H (second-order time derivative).
        H_temp=[];
        for i=1:num_t
            H_temp=[H_temp acc(i).*H_init];
        end
        H(3,:)=H_temp;            
    end
end
