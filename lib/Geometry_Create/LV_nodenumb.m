function [NODES] = LV_nodenumb(GEOM)
%LV_nodenumb: Numbers nodes from LV_geom function
% NODES = LV_nodenumbering(GEOM)
%   INPUT VARIABLE: (from LV_geom)
%       GEOM: structure including fields x, y, z
%             x, y, and z have dimensions:
%                   (rings base-apex +1, elems per ring, layers(Endo2Epi))
%   METHOD:
%       1. Numbering from theta increasing to 2pi 
%       2. Starting at the endocardial surface
%       3. From base to apical ring with apex node last
%   OUTPUT VARIABLE:
%       NODES: list of X Y Z coordinates with row number corresponding to
%       node number
%
% Thien-Khoi N. Phung (April 19, 2016)

% STEP 1: Identify number of element (rings, slices, peels)
[c, l, r] = size(GEOM.x); % number of points [circ. long. rad.]

% STEP 2: Parsing together X Y Z for nodes
NODES = [];
for lyr = 1:r % for loop for each peel (starting endo to epi)
    % Pull X Y Z and flip so that first node in peel is index 1
    % Reshape to make each matrix a vector going from theta 0 to 360 and 
    % base to apex
    x = reshape(fliplr(GEOM.x(:,:,lyr)), c*l, 1);
    y = reshape(fliplr(GEOM.y(:,:,lyr)), c*l, 1);
    z = reshape(fliplr(GEOM.z(:,:,lyr)), c*l, 1);
    
    % Store X Y Z
    % Keep only one APEX node (s*(r-1)+1 nodes per peel)
    NODES = [NODES; x(1:(end-(c-1))) y(1:(end-(c-1))) z(1:(end-(c-1)))];
end