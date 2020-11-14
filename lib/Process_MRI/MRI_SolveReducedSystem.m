function displacement=MRI_SolveReducedSystem(LHS_Matrix,Global_RHS,C,nd,nn,num_dof_mesh,num_time_nodes)
% Adapted from SLH_CMI_SolveReducedSystem.m by TK Phung (Feb 14, 2017)
% This function will use the constraints encoded in C to reduce the linear
% system of equations and solve.  Constraints can be free dofs (default if not
% listed), coupled dofs, or fixed parameters.  C is set up as follows:
% column 1: node number (i)
% column 2: node number (j) of lambda coupling, -1 if lambda(i) is free, 0
%           if lambda(i) is fixed (no displacement)
% column 3: coupling factor (x) for lambda such that lambda(i)=x*lambda(j)
% column 4: node number (dj1) of derivative (wrt 1) coupling, -1 if free, 0
%           if fixed
% column 5: coupling factor (dx1) for derivative (wrt 1) 
% column 6: node number (dj2) of derivative (wrt 2) coupling, -1 if free, 0
%           if fixed
% column 7: coupling factor (dx2) for derivative (wrt 2) 
% column 8: node number (dj12) of derivative (wrt 12) coupling, -1 if free,
%           0 if fixed
% column 9: coupling factor (dx12) for derivative (wrt 12)
%
% Each row specifies a new constraint equation.  Order of rows should not
% matter and indirect coupling is allowed.  
%
% IMPORTANT: Node numbers refer to array position AFTER sorting by theta
% and mu.
%
% IMPORTANT: Every coupled parameter must be linked, directly or indirectly, 
% to a free dof.  Cannot couple to a fixed parameter.  The free dof must appear 
% in C and must contain -1 in column 2 (for lambda coupling), 4 (for derivative 
% wrt 1 coupling), 6 (for derivative wrt 2 coupling), or 8 (for derivative
% wrt 12 coupling).
%
% All derivatives (remember lambda is a 0-order derivative) must be either linked
% to derivatives of the same order or fixed.  
%
% At the end of the code, the fixed and coupled parameters are merged with
% the optimized dofs to produce displacement vector which is of length nd.
%

% Initialize V, the mapping vector.
% V(i,1)=j with i=j implies a free dof (no coupling or fixing).
% V(i,1)=0 implies dof i is fixed.
% Otherwise, V(i,1)=j with i~=j implies dof i=V(i,2)*dof j (coupling).
% Weights are stored in V(:,2) and are 1 by default.
V=zeros(nd,2);
V(:,1)=[1:nd]';
V(:,2)=1;

% Arrange dofs.
[m,n]=size(C);
k=nd;

for i=1:m % Cycle through each node listed as constraint.
    for j=2:2:8 % Cycle through lambda and 3 derivative values.
        if C(i,j)>0 % Coupling.
            k=k-num_time_nodes;
            mult=C(i,j+1);
            search_flag=1;
            target=C(i,j);
            while search_flag % Cycle through until you find the end of the link.
                p=find(target==C(:,1));
                if C(p,j)==-1 % End of the trail.
                    for tnn=1:num_time_nodes
                        V(num_dof_mesh*(tnn-1)+C(i,1)+((j/2)-1)*nn,1)=num_dof_mesh*(tnn-1)+C(p,1)+((j/2)-1)*nn;
                        V(num_dof_mesh*(tnn-1)+C(i,1)+((j/2)-1)*nn,2)=mult;
                    end
                    search_flag=0;
                else
                    target=C(p,j);
                    mult=mult*C(p,j+1);
                end
            end
        elseif C(i,j)==0 % Fixing.
            k=k-num_time_nodes;
            for tnn=1:num_time_nodes
                V(num_dof_mesh*(tnn-1)+C(i,1)+((j/2)-1)*nn,1)=0;
            end
        end
    end
end

[m,n,NoFixDOF]=find(V(:,1)); % Eliminate zeros (fixed).
TrueDOF=unique(NoFixDOF); % Eliminate redundancies (coupled) and order.

% Reduction and solve.
Red_RHS(1:k,1)=0;
Red_LHS(1:k,1:k)=0;
ind1=1;

for i=1:nd
    p1=V(i,1);
    mult1=V(i,2);
    
    if p1~=0 % Make sure p1 is not fixed.
        if p1==i % p1 is free dof.
            Red_RHS(ind1,1)=Red_RHS(ind1,1)+Global_RHS(i,1);
            ind2=1;
            for j=1:nd
                p2=V(j,1);
                mult2=V(j,2);
                if p2~=0 % Make sure p2 is not fixed.
                    if p2==j % p2 is free dof.
                        Red_LHS(ind1,ind2)=Red_LHS(ind1,ind2)+LHS_Matrix(i,j);
                        ind2=ind2+1;
                    else % p2 is coupled.
                        m=find(p2==TrueDOF);
                        Red_LHS(ind1,m)=Red_LHS(ind1,m)+LHS_Matrix(i,j)*mult2;
                    end
                end
            end
            ind1=ind1+1;
        else % p1 is coupled.
            m=find(p1==TrueDOF);
            Red_RHS(m,1)=Red_RHS(m,1)+Global_RHS(i,1)*mult1;
            ind2=1;
            for j=1:nd
                p2=V(j,1);
                mult2=V(j,2);
                if p2~=0 % Make sure p2 is not fixed.
                    if p2==j % p2 is free dof.
                        Red_LHS(m,ind2)=Red_LHS(m,ind2)+LHS_Matrix(i,j)*mult1;
                        ind2=ind2+1;
                    else % p2 is coupled.
                        n=find(p2==TrueDOF);
                        Red_LHS(m,n)=Red_LHS(m,n)+LHS_Matrix(i,j)*mult1*mult2;
                    end
                end
            end
        end
    end
end

Red_disp=Red_LHS\Red_RHS;

% Reassemble.
for i=1:nd
    p=V(i,1);
    mult=V(i,2);
    if p~=0 % p is not fixed.
        ind=find(p==TrueDOF);
        displacement(i,1)=Red_disp(ind)*mult;
    else % p is fixed
        displacement(i,1)=0;
    end
end
