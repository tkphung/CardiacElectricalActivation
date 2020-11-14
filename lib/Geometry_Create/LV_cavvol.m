function CAVv = LV_cavvol(nodes,QUAD)
%LV_cavvol: Calculates the volumes of all pie wedge elements filling the
%cavity. This ASSUMES the base-apex axis is the X-axis.
%   CAVv = LV_cavvol(nodes,QUADS)
%   INPUTS:
%       nodes: all nodes from model (rows referring to node indices) (mm)
%       QUAD: QUAD element connectivity matrix (rows are each element) made
%             from LV_innersurf
%   OUTPUTS:
%       CAVv: volume of the cavity (mL)
% 
% Thien-Khoi Phung (February 7, 2018)
% Corrected coplanar quad error (June 17, 2019)

CAVv = zeros(length(QUAD),1);
for jz = 1:size(QUAD,1)
    
    % pull the quad points
    piepts = nodes(QUAD(jz,:),:);
    
    % define X axis points (pie center points)
    piepts = [piepts;
              max(piepts(:,1)) 0 0;
              min(piepts(:,1)) 0 0];
              
    if max(piepts(:,1)) == min(piepts(:,1))
        continue % catch if quad is coplanar in Y-Z
    end
    
    [~,v] = convhull(piepts);
    CAVv(jz) = v;
end

CAVv = sum(CAVv)/1000; % mL