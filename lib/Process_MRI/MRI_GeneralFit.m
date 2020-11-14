function [L,H] = MRI_GeneralFit(D,E,order)
%Adapted from SLH_CMI_GeneralFit.m by Thien-Khoi N. Phung (Feb. 14, 2017)
% This function actually computes the lambda values according to the bicubic 
% Hermite-Lagrange basis functions.
% The second output parameter returns the values of the shape functions
% at the local coordinates.  Local coordinates are stored in E.
% Order is the temporal order.  For space-only (no time) fits, use order=0.
%

% 1D shape functions
H00(1)=1-3*(E(1)^2)+2*(E(1)^3);
H00(2)=1-3*(E(2)^2)+2*(E(2)^3);
H10(1)=E(1)*((E(1)-1)^2);
H10(2)=E(2)*((E(2)-1)^2);
H01(1)=(E(1)^2)*(3-2*E(1));
H01(2)=(E(2)^2)*(3-2*E(2));
H11(1)=(E(1)^2)*(E(1)-1);
H11(2)=(E(2)^2)*(E(2)-1);

% Assemble 2D shape functions.
H_init=[H00(1)*H00(2) H01(1)*H00(2) H00(1)*H01(2) H01(1)*H01(2) H10(1)*H00(2) H11(1)*H00(2) H10(1)*H01(2) H11(1)*H01(2) H00(1)*H10(2) H01(1)*H10(2) H00(1)*H11(2) H01(1)*H11(2) H10(1)*H10(2) H11(1)*H10(2) H10(1)*H11(2) H11(1)*H11(2)];

% For temporal basis functions.
if order>0
    % Normalized time coordinates for an element of this order.
    t=[0:1/order:1];
    num_t=length(t); % =order+1

    % Compute Lagrange polynomial.
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

    % Assemble H and get result, L.
    L=0;
	H=[];
	for i=1:num_t
        L=L+LG(i)*dot(H_init,D(:,i));
        H=[H LG(i).*H_init];
	end
else % Space-only fit.
    H=H_init;
    L=dot(H,D);
end
