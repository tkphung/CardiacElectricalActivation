function HEXv = LV_hexvol(nodes,HEX)
%LV_hexvol: Calculates the volumes of all HEX elements
%   HEXv = LV_hexvol(nodes,HEX)
%   INPUTS:
%       NODES: all nodes from model (rows referring to node indices)
%       HEX: HEX element connectivity matrix (rows are each element)
%   OUTPUTS:
%       HEXv: volumes of each element (rows are elements)
% 
% Thien-Khoi Phung (Aug 08, 2016)

HEXv = zeros(length(HEX),1);
for jz = 1:length(HEX)
    [~,v] = convhull(nodes(HEX(jz,:),:));
    HEXv(jz) = v;
end