function Global_Damp=MRI_CalcDamping(smooth_weights,num_gp,num_dof_total,num_dof_elem,size_nodal_mu,size_nodal_theta,num_nodes,num_dof_mesh,num_time_nodes)
% SLH_CMI_CalcDamping Adapted by TK Phung (Feb 14, 2017)
% Calculate the combined global damping and mass matrix for smoothing.
% This is achieved by integrating the weighted derivatives as in Hashima
% et al.  ASSUMES THERE ARE 5 DERIVATIVES: 1/11/2/22/12.
%

% Initialize output matrix.
Global_Damp=zeros(num_dof_total,num_dof_total);

% Set up matricies for Gauss point position and weight for integration.
% For now, each is 7x7 since 7 is the maximum number of Gauss points in
% any direction.  See gauss1.f in cmiss code.
P=zeros(7,7);
W=zeros(7,7);

P(2,1)=-0.2886751345948130;
P(2,2)=-P(2,1);
P(3,1)=-0.3872983346207410;
P(3,3)=-P(3,1);
P(4,1)=-0.4305681557970260;
P(4,2)=-0.1699905217924280;
P(4,3)=-P(4,2);
P(4,4)=-P(4,1);
P(5,1)=-0.4530899229693320;
P(5,2)=-0.2692346550528410;
P(5,4)=-P(5,2);
P(5,5)=-P(5,1);
P(6,1)=-0.4662347571015760;
P(6,2)=-0.3306046932331330;
P(6,3)=-0.1193095930415990;
P(6,4)=-P(6,3);
P(6,5)=-P(6,2);
P(6,6)=-P(6,1);
P(7,1)=-0.4745539561713800;
P(7,2)=-0.3707655927996970;
P(7,3)=-0.2029225756886990;
P(7,5)=-P(7,3);
P(7,6)=-P(7,2);
P(7,7)=-P(7,1);
P=P+0.5000000000000000;

W(1,1)=1.0000000000000000;
W(2,1)=0.5000000000000000;
W(2,2)=W(2,1);
W(3,1)=0.2777777777777780;
W(3,2)=0.4444444444444440;
W(3,3)=W(3,1);
W(4,1)=0.1739274225687270;
W(4,2)=0.3260725774312730;
W(4,3)=W(4,2);
W(4,4)=W(4,1);
W(5,1)=0.1184634425280940;
W(5,2)=0.2393143352496830;
W(5,3)=0.2844444444444440;
W(5,4)=W(5,2);
W(5,5)=W(5,1);
W(6,1)=0.0856622461895850;
W(6,2)=0.1803807865240700;
W(6,3)=0.2339569672863460;
W(6,4)=W(6,3);
W(6,5)=W(6,2);
W(6,6)=W(6,1);
W(7,1)=0.0647424830844350;
W(7,2)=0.1398526957446390;
W(7,3)=0.1909150252525600;
W(7,4)=0.2089795918367350;
W(7,5)=W(7,3);
W(7,6)=W(7,2);
W(7,7)=W(7,1);

% Step through each element and build Global_Damp by integrating at Gauss
% points.  The mu and theta coordinates do not matter.  The same number of
% Gauss points are used for all elements.
for m=1:size_nodal_mu-1
    for t=1:size_nodal_theta
        % Calculate element number.
        element_number=(m-1)*size_nodal_theta+t;
        
        % Get theta indicies.
        corner24theta=t;
        if t==size_nodal_theta
            corner13theta=1;
        else
            corner13theta=corner24theta+1;
        end
        
        % Get mu indicies.
        corner12mu=m;
        corner34mu=corner12mu+1;
        
        % Build vector of the relevant dofs.  Use corner13theta,
        % corner24theta, corner12mu, and corner34mu to get indicies.  Ind
        % is index vector.  It is of length 16.  Ind has the order:
        % [lam_corner1 lam_corner2 lam_corner3 lam_corner4
        % d(lam)/d(e1)_corner1 d(lam)/d(e1)_corner2 etc]
        Ind(1)=size_nodal_mu*(corner13theta-1)+corner12mu; % 1:4 are lambda vals.
        Ind(2)=size_nodal_mu*(corner24theta-1)+corner12mu;
        Ind(3)=size_nodal_mu*(corner13theta-1)+corner34mu;
        Ind(4)=size_nodal_mu*(corner24theta-1)+corner34mu;
        Ind(5)=Ind(1)+num_nodes; % 5:8 are first derivatives wrt e1.
        Ind(6)=Ind(2)+num_nodes;
        Ind(7)=Ind(3)+num_nodes;
        Ind(8)=Ind(4)+num_nodes;
        Ind(9)=Ind(5)+num_nodes; % 9:12 are first derivatives wrt e2.
        Ind(10)=Ind(6)+num_nodes;
        Ind(11)=Ind(7)+num_nodes;
        Ind(12)=Ind(8)+num_nodes;
        Ind(13)=Ind(9)+num_nodes; % 13:16 are second derivatives wrt e1,e2.
        Ind(14)=Ind(10)+num_nodes;
        Ind(15)=Ind(11)+num_nodes;
        Ind(16)=Ind(12)+num_nodes;
        
        % Go through each Gauss point.
        for e1=1:num_gp(1)
            for e2=1:num_gp(2)
                % Calculate Gauss point weight.
                wgp=W(num_gp(1),e1)*W(num_gp(2),e2);
                
                % Get Gauss point position.
                E=[P(num_gp(1),e1) P(num_gp(2),e2)];
                
                % ASSUME 5 DERIVATIVES (1/11/2/22/12).
                for der=1:5 
                    % Obtain smoothing weighting factor for this derivative.
                    wder=smooth_weights(element_number,der);
                    
                    % H has length 16 and is in same order as Ind.
                    H=MRI_CalcBasisDerivs(E,der,0);
                    
                    for h1=1:num_dof_elem
                        for h2=1:num_dof_elem
                            for i=1:num_time_nodes
                                Global_Damp(num_dof_mesh*(i-1)+Ind(h1),num_dof_mesh*(i-1)+Ind(h2))=Global_Damp(num_dof_mesh*(i-1)+Ind(h1),num_dof_mesh*(i-1)+Ind(h2))+H(h1)*H(h2)*wgp*wder;
                            end
                        end
                    end
                end
            end
        end
    end
end
