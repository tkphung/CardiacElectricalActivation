function HEXc = LV_hexcent(nodes,HEX)
%LV_hexcent: Calculates the centroids of all HEX elements
%   HEXc = LV_hexcent(nodes,HEX)
%   INPUTS:
%       NODES: all nodes from model (rows referring to node indices)
%       HEX: HEX element connectivity matrix (rows are each element)
%   OUTPUTS:
%       HEXc: X Y Z for centroids of each element (rows are elements)
% 
% Thien-Khoi Phung (May 3, 2016)

HEXc = zeros(length(HEX),3);
for jz = 1:length(HEX)
    HEXc(jz,:) = mean(nodes(HEX(jz,:),:));
end