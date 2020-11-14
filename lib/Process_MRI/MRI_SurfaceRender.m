function H = MRI_SurfaceRender(ipnode_file,varargin)
% Adapted from SLH_CMI_Render.m by Thien-Khoi N. Phung (Feb 16, 2017)
% added varargin to pass figure handle (to plot on same graph)
% 
% function H=SLH_CMI_Render(ipnode_file,ipdata_file)
% SLH_CMI_Render.m
%
% Author: CM Ingrassia
%
% Render elements and data points in 3D.  Must specify an '.ipnode' file in
% prolate spheroidal coordinates with unscaled derivatives.  The mesh must
% represent a bicubic surface (no wall thickness).  Must be in the same format
% as the Continuity 5.5 node files.  The output of SLH_CMI_DOF_NodeWrite.m is
% acceptable.  A rectangular mu-theta grid is assumed.  Also assume that 
% min(theta)=0.  Second parameter is optional.  It is an '.IPDATA' file.  Use
% SLH_CMI_DataRead.m to read in the Cartesian coordinates.  The return is a 
% handle to the figure window.
%
% added color option (Thien-Khoi N. Phung October 25, 2018)
% If one varargin- that is the figure handle
% If two- start with a name-pair combination
% ie MRI_SurfaceRender(ipnode_file)
%    MRI_SurfaceRender(ipnode_file,H)
%    MRI_SurfaceRender(ipnode_file,'color',[253,141,60]./255,'handle',H)

handle_flag = false;
color_flag = false;
if ~isempty(varargin) && numel(varargin)>1
    for jz = 1:2:numel(varargin)
        switch varargin{jz}
            case 'handle'
                figure(varargin{jz+1});
                hold on;
                handle_flag = true;
            case 'color'
                color_flag = true;
                col = varargin{jz+1};
            otherwise
                error('ERROR: Check your varargins.')
        end
    end
elseif numel(varargin)==1
    figure(varargin{1});
    hold on;
    handle_flag = true;
end

if ~handle_flag
    H = figure('WindowStyle','docked');
    hold on
    handle_flag = true;
end

% Set number of divisions per arc or element boundary.
num_div=20;

E=0:1/num_div:1;

rads=pi/180;

% Read node file.
[nodal_mesh,focus]=MRI_IPNodeRead(ipnode_file);
nodal_mesh=sortrows(nodal_mesh,[3 2 1]);
nodal_mesh(:,2:3)=nodal_mesh(:,2:3).*rads;

% Get sizes.
[m,n]=size(nodal_mesh);
num_0_theta=length(find(nodal_mesh(:,3)==0));
size_nodal_mu=num_0_theta;
size_nodal_theta=(m/size_nodal_mu);
num_elem=size_nodal_theta*(size_nodal_mu-1);

% Get vectors.
nodal_mu=nodal_mesh(1:size_nodal_mu,2);
nodal_theta=nodal_mesh(1:size_nodal_mu:end,3);

% Plot circumferential element boundaries.  Use 1D interpolation to
% find additional points along the arc.  If min(mu)=0, just plot the
% point.  No need to interpolate around the circumference there.  All
% interpolation is done here, rather than calling a function, for speed.
% The interpolation is very simple since we only need 1 dimension.  Plot
% node points here.  Use blue for nodes and lines.  Use diamonds for nodes.
if min(nodal_mu)==0
    % Find and plot apex node.
    [x,y,z]=LV_P2C(nodal_mesh(1,1),nodal_mesh(1,2),nodal_mesh(1,3),focus);
    plot3(x,y,z,'bd');
    start_mu=2;
else
    start_mu=1;
end

for i=start_mu:size_nodal_mu
    for j=1:size_nodal_theta
        % Define nodal values needed for the interpolation.
        if j==size_nodal_theta
            ind0=i;
            t0=2*pi;
        else
            ind0=j*size_nodal_mu+i;
            t0=nodal_theta(j+1);
        end
        ind1=(j-1)*size_nodal_mu+i;
        t1=nodal_theta(j);
        
        % Get lambda and dL/de1 at nodes.
        L0=nodal_mesh(ind0,1); 
        dL0=nodal_mesh(ind0,4);
        L1=nodal_mesh(ind1,1);
        dL1=nodal_mesh(ind1,4);
        
        % Plot the node at E=0.
        [N0x,N0y,N0z]=LV_P2C(nodal_mesh(ind0,1),nodal_mesh(ind0,2),nodal_mesh(ind0,3),focus);
        plot3(N0x,N0y,N0z,'bd');
        
        % Build arc (like Noah did).
        for k=2:length(E)
            % We need a second point to plot a line.  Use the node (at E=0) or
            % the point from the previous iteration.
            if k==2
                pt_x=N0x;
                pt_y=N0y;
                pt_z=N0z;
            else
                pt_x=x_here;
                pt_y=y_here;
                pt_z=z_here;
            end
            
            % Get lambda values.
            HL0=1-3*(E(k)^2)+2*(E(k)^3);
            HdL0=E(k)*((E(k)-1)^2);
            HL1=(E(k)^2)*(3-2*E(k));
            HdL1=(E(k)^2)*(E(k)-1);
            L=HL0*L0+HdL0*dL0+HL1*L1+HdL1*dL1;
            
            % Get theta and calculate Cartesian point.
            t_here=t0-E(k)*(t0-t1);            
            [x_here,y_here,z_here]=LV_P2C(L,nodal_mu(i),t_here,focus);
            
            % Make vectors for plotting.
            x=[pt_x; x_here];
            y=[pt_y; y_here];
            z=[pt_z; z_here];
            
            if ~color_flag
                plot3(x,y,z,'k.-','LineWidth',3);
            else
                plot3(x,y,z,'.-','color',col,'LineWidth',3)
            end
        end
    end
end

% Plot longitudinal element boundaries.
for i=1:size_nodal_theta
    for j=1:size_nodal_mu-1
        % Define nodal values needed for the interpolation.
        ind0=(i-1)*size_nodal_mu+j;
        ind1=ind0+1;
        m0=nodal_mu(j);
        m1=nodal_mu(j+1);
        
        % Get lambda and dL/de2 at nodes.
        L0=nodal_mesh(ind0,1);
        dL0=nodal_mesh(ind0,5);
        L1=nodal_mesh(ind1,1);
        dL1=nodal_mesh(ind1,5);
        
        % Find the node at E=0.
        [N0x,N0y,N0z]=LV_P2C(nodal_mesh(ind0,1),nodal_mesh(ind0,2),nodal_mesh(ind0,3),focus);
        
        % Build arc (like Noah did).
        for k=2:length(E)
            % We need a second point to plot a line.  Use the node (at E=0) or
            % the point from the previous iteration.
            if k==2
                pt_x=N0x;
                pt_y=N0y;
                pt_z=N0z;
            else
                pt_x=x_here;
                pt_y=y_here;
                pt_z=z_here;
            end
            
            % Get lambda values.
            HL0=1-3*(E(k)^2)+2*(E(k)^3);
            HdL0=E(k)*((E(k)-1)^2);
            HL1=(E(k)^2)*(3-2*E(k));
            HdL1=(E(k)^2)*(E(k)-1);
            L=HL0*L0+HdL0*dL0+HL1*L1+HdL1*dL1;
            
            % Get mu and calculate Cartesian point.
            m_here=m0+E(k)*(m1-m0);            
            [x_here,y_here,z_here]=LV_P2C(L,m_here,nodal_theta(i),focus);
            
            % Make vectors for plotting.
            x=[pt_x; x_here];
            y=[pt_y; y_here];
            z=[pt_z; z_here];
            
            if ~color_flag
                plot3(x,y,z,'k.-','LineWidth',3);
            else
                plot3(x,y,z,'.-','color',col,'LineWidth',3)
            end
        end
    end
end

% % Plot data points, if present.
% if nargin==2
%     [x,y,z,count]=SLH_CMI_DataRead(ipdata_file);
%     for i=1:count
%         [z1,z2,z3]=SLH_CMI_C2P(x(i),y(i),z(i),focus);
%         if z2<=max(nodal_mu)
%             plot3(x(i),y(i),z(i),'bo','MarkerSize',1);
%         else
%             plot3(x(i),y(i),z(i),'g.');
%         end
%     end
% end

% Clean up figure and send back handle.
axis equal;
xlabel('X'); 
ylabel('Y'); 
zlabel('Z');

set(gcf,'color','w');
view([-90 -60])
