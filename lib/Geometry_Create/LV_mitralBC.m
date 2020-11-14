function [BCx, BCxyz] = LV_mitralBC(Le, Ce, Re)
%LV_mitralBC: returns the node indices for Boundary Conditions in the FE
%model
% [BCx, BCxyz] = LV_mitralBC(Le, Ce, Re)
%   INPUT VARIABLES:
%       Le: number of long. elements excluding apex cap elem.
%       Ce: number of circ. elements around each ring
%       Re: number of rad. elements through wall
%   OUTPUT VARIABLE:
%       BCx: fixed in X nodes
%       BCxyz: fixed in X, Y, Z nodes
% 
% Thien-Khoi N. Phung (April 26, 2016)

NPL = (Le+1)*Ce + 1; % number of nodes per layer

% FIX XYZ
% IDENTIFY OUTER RING OF NODES IN BASE
BCxyz = (1:Ce)' + NPL*Re;

% FIX X
% IDENTIFY OTHER NODES IN BASE
BCx = [];
for ly = 1:(Re+1) % Fixed to Re+1 credit to DJB
    BCx = [BCx; (1:Ce)' + NPL*(ly-1)];
end