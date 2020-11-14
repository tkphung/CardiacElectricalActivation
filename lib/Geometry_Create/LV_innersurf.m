function [QUAD, TRI] = LV_innersurf(Le, Ce, Re)
%LV_innersurf: creates connectivity matrices for the QUAD and TRI surfaces
%of the inner LV cavity
% [QUAD, TRI] = LV_innersurf(Le, Ce, Re)
%   INPUT VARIABLES:
%       Le: number of long. elements excluding apex cap elem.
%       Ce: number of circ. elements around each ring
%       Re: number of rad. elements through wall
%   OUTPUT VARIABLE:
%       QUAD: matrix with each row being one elements connectivity with
%       respect to nodes defined by LV_nodenumb QUAD surfaces are for the
%       HEX elements
%		TRI: connectivity for nodes for the PENT elements
% 
% Thien-Khoi N. Phung (April 26, 2016)

NPL = (Le+1)*Ce + 1; % number of nodes per layer
EPL = Le*Ce;         % number of elements per layer (excludes apex cap)

% QUAD node connectivity
QUAD = zeros(EPL, 4);
for q = 1:EPL
    if mod(q,Ce)~= 0
        QUAD(q,:) = [q, q+1, (q+1)+Ce, q+Ce];
    else
        QUAD(q,:) = [q, (q+1)-Ce, q+1, q+Ce];
    end
end

% TRI node connectivity
TRI = zeros(Ce,3);
for q = 1:Ce
    if mod(q,Ce)~= 0
        TRI(q,:) = [NPL, (NPL-Ce)+q-1, (NPL-Ce)+q];
    else
        TRI(q,:) = [NPL, (NPL-Ce)+q-1, (NPL-Ce)];
    end
end