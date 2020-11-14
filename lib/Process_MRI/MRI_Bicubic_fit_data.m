function [node_matrix,rms_err] = MRI_Bicubic_fit_data(Data,focus,fileout,filename,meshdensity)
% from SLH_CMI_Bicubic_fit_data.m 02/14/17
% from SLH_CMI_Bicubic.m version 10/5/04
% 
% Bicubic fit of (x,y,z) data to a prolate mesh.
%
% 
% Modified:
% uses focal length input by user rather than focal length from starting_mesh_file
% does not read in data from input .IPDATA file,
% instead uses data [x y z] matrix supplied as input argument
% does not output a node file,
% instead returns a matrix 'node_matrix': a 72 x 6 matrix with nodal values 
% and columns listed in folowing order: lambda, mu, theta, dlds1, dlds2, d2lds1ds2
%
% Adapted by Thien-Khoi Phung (February 14, 2016 <3)

format short;

%% Flags
write_file     = fileout;
smooth         = true; % Kept from SLH_CMI_Bicubic_fit_data.m
constraints    = true; % Kept from SLH_CMI_Bicubic_fit_data.m
compute_errors = true; % Kept from SLH_CMI_Bicubic_fit_data.m

%% Set Mesh Desnity-dependent Start Files and Options
starting_mesh_file = ['start_mesh_' meshdensity 'elem.ipnode'];

if smooth
    switch meshdensity 
        case '4x2'
            elem_file='surf_8.ipelem';
            ipfit_file='surf_8_smooth2.ipfit';
        case '4x4'
            elem_file='surf_16.ipelem';
            ipfit_file='surf_16_smooth1.ipfit';
        case '4x8'
            elem_file='surf_32.ipelem';
            ipfit_file='surf_32_smooth1.ipfit';
        otherwise
            error('ERROR: Check bicubic mesh density input')
    end
end

% If requested, set up C, the nodal constraint matrix as follows:
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
if constraints
    % Link lambda for apex nodes to #1, derivatives wrt 1 and 12 are
    % fixed, derivatives wrt 2 are coupled by -1 to the node 180 degrees
    % apart in theta.
    switch meshdensity 
        case '4x2'
            % for a 4x2 mesh:
            C=[ 1 -1 1 0 1 -1  1 0 1;
                4  1 1 0 1 -1  1 0 1;
                7  1 1 0 1  1 -1 0 1;
                10 1 1 0 1  4 -1 0 1];
        case '4x4'
            % for a 4x4 mesh:
            C = [ 1 -1 1 0 1 -1  1 0 1;
                  6  1 1 0 1 -1  1 0 1;
                 11  1 1 0 1  1 -1 0 1;
                 16  1 1 0 1  6 -1 0 1];
        case '4x8'
            % for a 4x8 mesh (more elements in circ direction):
            C=[ 1 -1 1 0 1 -1  1 0 1;
                6  1 1 0 1 -1  1 0 1;
               11  1 1 0 1 -1  1 0 1;
               16  1 1 0 1 -1  1 0 1;
               21  1 1 0 1  1 -1 0 1;
               26  1 1 0 1  6 -1 0 1;
               31  1 1 0 1 11 -1 0 1;
               36  1 1 0 1 16 -1 0 1];
        otherwise
            error('ERROR: Check bicubic mesh density')
    end
end

%% Read nodal mesh data and sort
[nodal_mesh,~] = MRI_IPNodeRead(starting_mesh_file);
unsorted_nodal_mesh_deg = nodal_mesh(:,2:3);
nodal_mesh = sortrows(nodal_mesh,[3 2 1]);

% Get sizes
[m,~] = size(nodal_mesh);
num_0_theta = length(find(nodal_mesh(:,3)==0));
size_nodal_mu = num_0_theta;
size_nodal_theta =(m/size_nodal_mu);
num_elem = size_nodal_theta*(size_nodal_mu-1);

% Put into radians
rads = pi/180;
nodal_mu_deg = nodal_mesh(1:size_nodal_mu,2);
nodal_theta_deg = nodal_mesh(1:size_nodal_mu:end,3);
nodal_mu = nodal_mu_deg.*rads;
nodal_theta = nodal_theta_deg.*rads;

% Build initial guess vector (collection of nodal lambda and derivatives)
init_guess = [nodal_mesh(:,1); nodal_mesh(:,4); nodal_mesh(:,5); nodal_mesh(:,6)];
num_dof_total = length(init_guess);
num_dof_elem = 16;

if smooth>0
    % Set number of Gauss points in psi1 and psi2 directions.
    gp_limit = 7;
    num_gp(1) = 3;
    num_gp(2) = 3;
    if num_gp(1)>gp_limit || num_gp(2)>gp_limit
        error('Gauss point limit exceeded');
    end
end

% Initialize
pts_elem(1:num_elem,1) = 0;
Global_Stiff(1:num_dof_total,1:num_dof_total) = 0;
Global_RHS(1:num_dof_total,1) = 0;
LHS_Matrix(1:num_dof_total,1:num_dof_total) = 0;

% Read data points
data_x = Data(:,1);
data_y = Data(:,2);
data_z = Data(:,3);
count  = size(Data,1);
data_w = ones(count,3);

% Convert to prolate, ignoring any values outside nodal mu range, and find
% the 4 nearest nodal points using convention in Hashima et al.  We assume
% that min_nodal_theta=0.
j = 0;
min_nodal_mu = min(nodal_mu);
max_nodal_mu = max(nodal_mu);
% min_nodal_theta = min(nodal_theta); % Not Used
max_nodal_theta = max(nodal_theta);

for i = 1:count
    [z1,z2,z3] = LV_C2P(data_x(i),data_y(i),data_z(i),focus);
    % Replaced SLH_CMI_C2P with LV_C2P
    if z2 >= min_nodal_mu && z2 <= max_nodal_mu
        j = j+1;
        
        % Special theta cases 
        if z3>=2*pi % If SLH_CMI_C2P happens to give us a theta>=2*pi.
            z3=z3-2*pi;
            t13=nodal_theta(2);
            t24=nodal_theta(1);
            E(1)=(t13-z3)/(t13-t24);    
            % indicies for theta vals
            corner13theta=2;
            corner24theta=1;
        elseif z3>=max_nodal_theta % If theta=0 is E(1)=0.           
            t13=2*pi;
            t24=max_nodal_theta;
            E(1)=(t13-z3)/(t13-t24);
            % indicies for theta vals
            corner13theta=1;
            corner24theta=size_nodal_theta;
        else % General theta case.
            min_t=length(find(nodal_theta<=z3));
            t13=nodal_theta(min_t+1);
            t24=nodal_theta(min_t);
            E(1)=(t13-z3)/(t13-t24);
            % indicies for theta vals
            corner13theta=min_t+1;
            corner24theta=min_t;
        end
        % Special mu cases
        if z2==min_nodal_mu % Smallest mu value.
            m12=nodal_mu(1);
            m34=nodal_mu(2);
            E(2)=0;
            % indicies for mu vals
            corner12mu=1;
            corner34mu=2;
        elseif z2==max_nodal_mu % Maximum mu value.
            m12=nodal_mu(size_nodal_mu-1);
            m34=nodal_mu(size_nodal_mu);
            E(2)=1;
            % indicies for mu vals
            corner12mu=size_nodal_mu-1;
            corner34mu=size_nodal_mu;
        else % General mu case.
            min_m=length(find(nodal_mu<z2));
            m12=nodal_mu(min_m);
            m34=nodal_mu(min_m+1);
            E(2)=(z2-m12)/(m34-m12);
            % indicies for mu vals
            corner12mu=min_m;
            corner34mu=min_m+1;
        end
        
        % Determine element number.
        element_number=(corner12mu-1)*size_nodal_theta+corner24theta;
        pts_elem(element_number)=pts_elem(element_number)+1;
        
        % Build vector of the relevant dofs.  Use corner13theta,
        % corner24theta, corner12mu, and corner34mu to get indicies.  Ind
        % is index vector.  It is of length 16.  Ind has the order:
        % [lam_corner1 lam_corner2 lam_corner3 lam_corner4
        % d(lam)/d(e1)_corner1 d(lam)/d(e1)_corner2 etc]
        % It will be used to access init_guess, so it follows the same
        % order.
        Ind(1)=size_nodal_mu*(corner13theta-1)+corner12mu; % 1:4 are lambda vals.
        Ind(2)=size_nodal_mu*(corner24theta-1)+corner12mu;
        Ind(3)=size_nodal_mu*(corner13theta-1)+corner34mu;
        Ind(4)=size_nodal_mu*(corner24theta-1)+corner34mu;
        Ind(5)=Ind(1)+m; % 5:8 are first derivatives wrt e1.
        Ind(6)=Ind(2)+m;
        Ind(7)=Ind(3)+m;
        Ind(8)=Ind(4)+m;
        Ind(9)=Ind(5)+m; % 9:12 are first derivatives wrt e2.
        Ind(10)=Ind(6)+m;
        Ind(11)=Ind(7)+m;
        Ind(12)=Ind(8)+m;
        Ind(13)=Ind(9)+m; % 13:16 are second derivatives wrt e1,e2.
        Ind(14)=Ind(10)+m;
        Ind(15)=Ind(11)+m;
        Ind(16)=Ind(12)+m;
        
        % Get vector from init_guess containing data from 4 corner nodes.
        dof_model = init_guess(Ind);
        
        % Get lambda value.  H returns the Hermite products (interpolation coefficients).
        [lam_model,H] = MRI_GeneralFit(dof_model,E,0);

        % Compute Global_Stiff and Global_RHS entries.
        lam_diff=z1-lam_model;
        for h1=1:num_dof_elem
            Global_RHS(Ind(h1),1)=Global_RHS(Ind(h1),1)+H(h1)*lam_diff*data_w(i);
            for h2=1:num_dof_elem
                Global_Stiff(Ind(h1),Ind(h2))=Global_Stiff(Ind(h1),Ind(h2))+H(h1)*H(h2)*data_w(i);
            end
        end
    end
end
size_data=j;

%% SMOOTHING
% Add smoothing functional if requested.  Global_Damp assumes dof vector is
% ordered.  Global_Damp is arranged the same way as Global_Stiff.
if smooth
    % Get derivative coefficients from file.
    smooth_weights=MRI_IPFitParse(ipfit_file,elem_file,num_elem,unsorted_nodal_mesh_deg,nodal_mu_deg,nodal_theta_deg,size_nodal_theta);
    % This is actually a combination of the damping and mass matricies.
    Global_Damp=MRI_CalcDamping(smooth_weights,num_gp,num_dof_total,num_dof_elem,size_nodal_mu,size_nodal_theta,m,num_dof_total,1);
    LHS_Matrix=Global_Stiff+Global_Damp;
else
    LHS_Matrix=Global_Stiff;    
end

%% Solve for displacement using \.
if constraints
    % This function will return displacement vector of the proper size.
    displacement=MRI_SolveReducedSystem(LHS_Matrix,Global_RHS,C,num_dof_total,m,num_dof_total,1);
else
    displacement=LHS_Matrix\Global_RHS;
end

% These are the fitted dof values sorted by mu and theta.
optimized_dof=init_guess+displacement;

%% Write out a new node file in same order as original mesh.
if write_file
    unsorted_optimized_dof=MRI_DOF_NodeWrite(filename,starting_mesh_file,optimized_dof,unsorted_nodal_mesh_deg,nodal_mu_deg,nodal_theta_deg,m,size_nodal_mu,focus);
end

%% Instead of writing out node file, calculate values for function return

for i=1:m
    m1=length(find(nodal_mu_deg<=unsorted_nodal_mesh_deg(i,1)));
    t=length(find(nodal_theta_deg<=unsorted_nodal_mesh_deg(i,2)));
    lam_ind=size_nodal_mu*(t-1)+m1;
    d1_ind=lam_ind+m;
    d2_ind=d1_ind+m;
    d3_ind=d2_ind+m;
    dof_data(i,1)=optimized_dof(lam_ind);
    dof_data(i,2)=optimized_dof(d1_ind);
    dof_data(i,3)=optimized_dof(d2_ind);
    dof_data(i,4)=optimized_dof(d3_ind);
end

% 72 x 6 matrix with nodal values - columns in folowing order: lambda, mu
% theta, dlds1, dlds2, d2lds1ds2
node_matrix = [dof_data(:,1) unsorted_nodal_mesh_deg dof_data(:,2:end)];

%
% The following code performs error estimate.  It is the 2-norm of the
% difference between lambda and the corresponding lambda on the model.
% Until we come up with something better, we must cycle through each
% data point again to calculate this.
%

if compute_errors>0
	% Initialize global error.
	err=0;
	
	% Go through points again.
	for i=1:count
        [z1,z2,z3]=LV_C2P(data_x(i),data_y(i),data_z(i),focus);
        if z2>=min_nodal_mu && z2<=max_nodal_mu
            j=j+1;
            
            % Special theta cases 
            if z3>=2*pi % If SLH_CMI_C2P happens to give us a theta>=2*pi.
                z3=z3-2*pi;
                t13=nodal_theta(2);
                t24=nodal_theta(1);
                E(1)=(t13-z3)/(t13-t24);    
                % indicies for theta vals
                corner13theta=2;
                corner24theta=1;
            elseif z3>=max_nodal_theta % If theta=0 is E(1)=0.           
                t13=2*pi;
                t24=max_nodal_theta;
                E(1)=(t13-z3)/(t13-t24);
                % indicies for theta vals
                corner13theta=1;
                corner24theta=size_nodal_theta;
            else % General theta case.
                min_t=length(find(nodal_theta<=z3));
                t13=nodal_theta(min_t+1);
                t24=nodal_theta(min_t);
                E(1)=(t13-z3)/(t13-t24);
                % indicies for theta vals
                corner13theta=min_t+1;
                corner24theta=min_t;
            end
            % Special mu cases
            if z2==min_nodal_mu % Smallest mu value.
                m12=nodal_mu(1);
                m34=nodal_mu(2);
                E(2)=0;
                % indicies for mu vals
                corner12mu=1;
                corner34mu=2;
            elseif z2==max_nodal_mu % Maximum mu value.
                m12=nodal_mu(size_nodal_mu-1);
                m34=nodal_mu(size_nodal_mu);
                E(2)=1;
                % indicies for mu vals
                corner12mu=size_nodal_mu-1;
                corner34mu=size_nodal_mu;
            else % General mu case.
                min_m=length(find(nodal_mu<z2));
                m12=nodal_mu(min_m);
                m34=nodal_mu(min_m+1);
                E(2)=(z2-m12)/(m34-m12);
                % indicies for mu vals
                corner12mu=min_m;
                corner34mu=min_m+1;
            end
            
            % Build vector of the relevant dofs.  Use corner13theta,
            % corner24theta, corner12mu, and corner34mu to get indicies.  Ind
            % is index vector.  It is of length 16.  Ind has the order:
            % [lam_corner1 lam_corner2 lam_corner3 lam_corner4
            % d(lam)/d(e1)_corner1 d(lam)/d(e1)_corner2 etc]
            % It will be used to access init_guess, so it follows the same
            % order.
            Ind(1)=size_nodal_mu*(corner13theta-1)+corner12mu; % 1:4 are lambda vals.
            Ind(2)=size_nodal_mu*(corner24theta-1)+corner12mu;
            Ind(3)=size_nodal_mu*(corner13theta-1)+corner34mu;
            Ind(4)=size_nodal_mu*(corner24theta-1)+corner34mu;
            Ind(5)=Ind(1)+m; % 5:8 are first derivatives wrt e1.
            Ind(6)=Ind(2)+m;
            Ind(7)=Ind(3)+m;
            Ind(8)=Ind(4)+m;
            Ind(9)=Ind(5)+m; % 9:12 are first derivatives wrt e2.
            Ind(10)=Ind(6)+m;
            Ind(11)=Ind(7)+m;
            Ind(12)=Ind(8)+m;
            Ind(13)=Ind(9)+m; % 13:16 are second derivatives wrt e1,e2.
            Ind(14)=Ind(10)+m;
            Ind(15)=Ind(11)+m;
            Ind(16)=Ind(12)+m;
            
            % Get vector from optimized_dof containing data from 4 corner nodes.
            dof_model=optimized_dof(Ind);
            
            % Get lambda value.  H returns the Hermite products (interpolation coefficients).
            [lam_model,H] = MRI_GeneralFit(dof_model,E,0);
            
            % For calculating RMS error in lambda.
            err=err+((lam_model-z1)*data_w(i))^2;
        end
    end
	
	% Finish error computations.
    % RMS error in lambda.
	rms_err=sqrt(err/size_data);
    %avg_err=err/size_data
    %avg_err_norm=avg_err/focus
end
