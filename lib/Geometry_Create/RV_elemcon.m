function [HEX] = RV_elemcon(Le, Ce, Re)
%RV_elemcon: creates element connnectivity matrix for given number of
%elements in the longitudinal (excluding cap), circumferential, and radial
%directions
% [HEX] = RV_elemcon(Le, Ce, Re)
%   INPUT VARIABLES:
%       Le: number of long. elements excluding apex cap elem.
%       Ce: number of circ. elements around each ring
%       Re: number of rad. elements through wall
%   OUTPUT VARIABLE:
%       HEX: matrix with each row being one elements connectivity with
%       respect to nodes defined by RV_nodenumb
% 
% Thien-Khoi N. Phung (September 20, 2016)

% Using Sliding Window similar to LV_elemnbr
% Create indices of nodes (associated with RV_nodenumb)
nodeindex = reshape(1:((Ce+1)*(Le+1)*(Re+1)),(Ce+1),(Le+1),(Re+1));

% Nodes
nONE = nodeindex(1:end-1, 1:end-1, 1:end-1);
nTWO = nodeindex(2:end,   1:end-1, 1:end-1);
nTHR = nodeindex(2:end,   1:end-1, 2:end);
nFOU = nodeindex(1:end-1, 1:end-1, 2:end);
nFIV = nodeindex(1:end-1, 2:end,   1:end-1);
nSIX = nodeindex(2:end,   2:end,   1:end-1);
nSEV = nodeindex(2:end,   2:end,   2:end);
nEIG = nodeindex(1:end-1, 2:end,   2:end);

% HEX elements
HEX = [nONE(:) nTWO(:) nTHR(:) nFOU(:) nFIV(:) nSIX(:) nSEV(:) nEIG(:)];