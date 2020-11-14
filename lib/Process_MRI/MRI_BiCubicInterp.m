function A=MRI_BiCubicInterp(xData,yData,NodeData,scale_der)
% Adapted by Thien-Khoi N. Phung Feb. 23, 2017 from SLH_CMI_Interp.m
% Interpolate new lambda values and derivatives
%See Hashima et al. 1992
%
% xData and yData are the new grid points
%
% NodeData has following format:
% (:,:,1) is lambda
% (:,:,2) is mu
% (:,:,3) is theta
% (:,:,4) is D lambda/D theta
% (:,:,5) is D lambda/D mu
% (:,:,6) is D^2 lambda/ D theta D mu
%

% Get all mu and theta values
NewMuVec=xData(1,:);
NewThetaVec=yData(:,1);
OldMuVec=NodeData(:,1,2);
OldThetaVec=NodeData(1,:,3);

SizeNewMu=length(NewMuVec);
SizeNewTheta=length(NewThetaVec);
[SizeOldMu,SizeOldTheta,n]=size(NodeData);

rads = pi/180;

% Set up A (output matrix of lambda, mu, theta, and derivatives in same
% order as NodeData).  For now, interchange lambda and theta to make call
% to GetLambda easier.  See below.
A(1:SizeNewMu,1:SizeNewTheta,1:n)=0;
A(:,:,2)=xData';
A(:,:,1)=yData';

% First, search to see if we're at a nodal point, else do interpolation 
for i=1:SizeNewMu
    for j=1:SizeNewTheta
        mu=NewMuVec(i);
        theta=NewThetaVec(j);
        
        % Search
        found_nodal_point=0;
        for k=1:SizeOldMu
            for m=1:SizeOldTheta
                if(mu==OldMuVec(k) & theta==OldThetaVec(m))
                    A(i,j,3)=NodeData(k,m,1);
                    A(i,j,4:n)=NodeData(k,m,4:n);
                    found_nodal_point=1;
                    break;
                end
            end
        end
        
        % If we're not at a nodal position, do interpolation
        if(found_nodal_point==0)
            % Find the appropriate lambda values via mu and theta
            % Handle special cases at mu=min and theta=max
            if(i==1) % mu=min
                min_t=length(find(OldThetaVec<=theta));
                corner(1,:)=NodeData(1,min_t+1,:);
                corner(2,:)=NodeData(1,min_t,:);
                corner(3,:)=NodeData(2,min_t+1,:);
                corner(4,:)=NodeData(2,min_t,:);
                E(1)=(corner(1,3)-theta)/(corner(1,3)-corner(2,3));
                E(2)=0;
            elseif(j==SizeNewTheta) % theta=max
                min_m=length(find(OldMuVec<mu));
                corner(1,:)=NodeData(min_m,2,:);
                corner(2,:)=NodeData(min_m,1,:);
                corner(3,:)=NodeData(min_m+1,2,:);
                corner(4,:)=NodeData(min_m+1,1,:);
                E(1)=1;
                E(2)=(mu-corner(2,2))/(corner(4,2)-corner(2,2));
            else
                min_t=length(find(OldThetaVec<=theta));
                min_m=length(find(OldMuVec<mu));
                corner(1,:)=NodeData(min_m,min_t+1,:);
                corner(2,:)=NodeData(min_m,min_t,:);
                corner(3,:)=NodeData(min_m+1,min_t+1,:);
                corner(4,:)=NodeData(min_m+1,min_t,:);
                E(1)=(corner(1,3)-theta)/(corner(1,3)-corner(2,3));
                E(2)=(mu-corner(2,2))/(corner(4,2)-corner(2,2));
            end
            
            if scale_der>0
                % Rescale derivatives by angle change to put deriv in terms of
                % psi 1 and psi 2 (local coordinates)
                corner(:,4)=corner(:,4)*((corner(1,3)-corner(2,3))*rads);
                corner(:,5)=corner(:,5)*((corner(4,2)-corner(2,2))*rads);
                corner(:,6)=corner(:,6)*((pi/180)*(corner(1,3)-corner(2,3))*(corner(4,2)-corner(2,2))*rads);
            end
            
            % Bicubic interpolation
            A(i,j,3:n)=MRI_GetLambda(corner,E);
            
            if scale_der>0
                % rescale interpolated derivatives from psi1 and psi2 back to
                % degrees
                A(i,j,4)=A(i,j,4)/((corner(1,3)-corner(2,3))*rads);
                A(i,j,5)=A(i,j,5)/((corner(4,2)-corner(2,2))*rads);
                A(i,j,6)=A(i,j,6)/((corner(1,3)-corner(2,3))*(corner(4,2)-corner(2,2))*(pi/180)*rads);
            end    
        end
    end
end

% Interchange lambda and theta matrices
temp=A(:,:,1);
A(:,:,1)=A(:,:,3);
A(:,:,3)=temp;