function [eFIB, Qfcr] = LV_fibers(HEX,NODES,Le,Ce,Re,enA,epA)
%LV_fibers: generates fiber vectors for each HEX element according to the
%60 to -60 deg endo to epi distribution through the wall
%   [eFIB, Qfcr] = LV_fibers(HEX,nodes,Le,Ce,Re)
%   INPUTS:
%       HEX: Hexahedral elements (rows are elements, columns are node indices)
%       NODES: all nodes from model (rows referring to node indices)
%       Le: number of longitudinal elements
%       Ce: number of circumferential elements
%       Re: number of radial elements
%       enA: endo Angle (normal: 60 degrees)
%       epA: epi Angle (normal: -60 degrees)
%           NOTE: 0 degrees is circ, 90 degrees is long
%   OUTPUTS:
%       eFIB: fiber vectors in cartesian for each HEX element
%       Qfcr: fiber rotation matrix [F C R] for each element
% 
% ADAPTED from portion of A.C.Estrada NodeNumbering_General.m
% Thien-Khoi Phung (April 26, 2016)

EPL = Le*Ce;         % number of elements per layer (excludes apex cap)

FIBa = linspace(enA,epA,Re).*pi/180; % Fiber angles for each layer (Endo --> Epi)

eFIB = zeros(length(HEX),3);
Qfcr = zeros(3,3,length(HEX));
for ly = 1:Re
    % Set fiber vector in (long, circ, rad) for each layer
    layfib = [sin(FIBa(ly)); cos(FIBa(ly)); 0];
    
    % Cycle through elements in the layer
    for ei = (1:EPL) + EPL*(ly-1)
        NI = HEX(ei,:);   % index for nodes
        no = NODES(NI,:); % matrix of node coordinates for element
        Q  = LV_ElementtoGlobal(no); % rotation matrix for XYZ -> LCR
                
        eFIB(ei,:) = [Q*layfib]'; % assign fiber vector in cartesian
        
        Rv = Q(:,3); % radial vector
        Cv = cross(Rv, eFIB(ei,:)'); % cross-fiber vector
            Cv = Cv./norm(Cv); % unit vector

        % Store rotation matrix (Columns are vectors)
        Qfcr(:,:,ei) = [eFIB(ei,:)' Cv Rv];
%         Qfcr(:,:,ei) = Q; % TROUBLESHOOT LINE TO EXPORT LOCAL COORD (LCR)
    end
end
        