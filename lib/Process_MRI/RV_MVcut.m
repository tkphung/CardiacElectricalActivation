function [GEOM] = RV_MVcut(GEOM,mvX)
%RV_MVcut: cuts RV data above MV
% [GEOM] = RV_MVcut(GEOM,mvX)
%   INPUT VARIABLE: 
%       GEOM: structure including fields x, y, z
%             x, y, and z have dimensions:
%                   (rings base-apex +1, elems per ring, layers(Endo2Epi))
%       mvX: the x value of the Mitral Valve plane (negative X is above the
%            mitral valve)
%   METHOD:
%         for each longitudinal column of nodes (in RVgeom this is each row of values)
%         if the node(@base) > mvX
%             Interpolate (x,mu) for mvMU given the mvX
%             In Prolate Spheroid Coordinates
%             Create vector for mu = 0 to mvMU
%             Interpolate new Lambda and Theta values
%             Convert them to Cartesian Coordinates
%           Store both of those in the RVgeom structure
%   OUTPUT VARIABLE:
%       GEOM
%
% Thien-Khoi N. Phung (February 20, 2018)

[c, l, r] = size(GEOM.x); % number of points [circ. long. rad.]

for ci = 1:c        % Along Circumferential Direction
    for ri = 1:r    % Along Radial Direction
        % PULL X VALUES (longitudinal direction)
        x  = squeeze(GEOM.x(ci,:,ri));
            
        if x(end) < mvX % if base of RV is above MV plane
            
            % PULL MU VALUES
            mu = squeeze(GEOM.mu(ci,:,ri));
            
            mvMU = interp1(x,mu,mvX,'pchip');
            
            % Create new points in PROLATE SPHEROID COORDINATES
            newmu = linspace(mu(1),mvMU,l);
            
            % Interpoalte lambdas in PSC
            theta  = squeeze(GEOM.theta(ci,:,ri));
            lambda = squeeze(GEOM.lambda(ci,:,ri));
            
            newlambda = interp1(mu,lambda,newmu,'pchip');
            % theta is the same for each longitudinal set of nodes
            
            % PSC to CARTESIAN COORDINATES
            [newx,newy,newz] = LV_P2C(newlambda,newmu,theta,GEOM.d);
                        
            % Store data
            GEOM.x(ci,:,ri) = newx;
            GEOM.y(ci,:,ri) = newy;
            GEOM.z(ci,:,ri) = newz;
            GEOM.mu(ci,:,ri) = newmu;
            GEOM.lambda(ci,:,ri) = newlambda;
        end
    end
end
            
        