function regs = AHASegment(LOCS,focus)
%AHASegment: Takes a LV model and maps it onto a 16 segment AHA
% AHASeg = AHASegmentData(LOCS,DATA,focus,varargin)
%   INPUTS:
%       LOCS- cartesian points in cardiac coordinates
%       focus- prolate spheroidal coord focus 
%   OUTPUT:
%       regs- 16 seg bin for each element
% 
% Created by Thien-Khoi N. Phung (June 27, 2019)

%% Convert Data from Cartesian to PSC
[~,M,T] = LV_C2P(LOCS(:,1),LOCS(:,2),LOCS(:,3),focus);

% NOTE: theta (T) starts at mid-septum 0 rad and follows SPLAS to 2pi rads
theta = T;
% Reverse to SALPS
theta = -theta + 2*pi;

rad = M.*(180/pi);

%% Sort data in segments - data is in SALPS with theta = 0 at septum

% Construct topology, detect which of the points are contained in each
% segment (see Pim's NB Growth>UVA dog study>Infarct Sizes
AHAseg(1).topo  = find((theta >= 1/3*pi) & (theta < 2/3*pi) & (rad >= 86));
AHAseg(2).topo  = find((theta >= 0)      & (theta < 1/3*pi) & (rad >= 86));
AHAseg(3).topo  = find((theta >= 5/3*pi) & (theta < 6/3*pi) & (rad >= 86));
AHAseg(4).topo  = find((theta >= 4/3*pi) & (theta < 5/3*pi) & (rad >= 86));
AHAseg(5).topo  = find((theta >= 3/3*pi) & (theta < 4/3*pi) & (rad >= 86));
AHAseg(6).topo  = find((theta >= 2/3*pi) & (theta < 3/3*pi) & (rad >= 86));

AHAseg(7).topo  = find((theta >= 1/3*pi) & (theta < 2/3*pi) & (rad >= 45) & (rad < 86));
AHAseg(8).topo  = find((theta >= 0)      & (theta < 1/3*pi) & (rad >= 45) & (rad < 86));
AHAseg(9).topo  = find((theta >= 5/3*pi) & (theta < 6/3*pi) & (rad >= 45) & (rad < 86));
AHAseg(10).topo  = find((theta >= 4/3*pi) & (theta < 5/3*pi) & (rad >= 45) & (rad < 86));
AHAseg(11).topo  = find((theta >= 3/3*pi) & (theta < 4/3*pi) & (rad >= 45) & (rad < 86));
AHAseg(12).topo  = find((theta >= 2/3*pi) & (theta < 3/3*pi) & (rad >= 45) & (rad < 86));

AHAseg(13).topo  = find((theta >= 1/4*pi) & (theta < 3/4*pi) & (rad < 45));
AHAseg(15).topo  = find((theta >= 5/4*pi) & (theta < 7/4*pi) & (rad < 45));
AHAseg(16).topo  = find((theta >= 3/4*pi) & (theta < 5/4*pi) & (rad < 45));
AHAseg(14).topo  = find((  ((theta >= 7/4*pi) & (theta < 8/4*pi)) | ...
                            ((theta >= 0)      & (theta < 1/4*pi)) )...
                         & (rad < 45));
                     
regs = zeros(size(LOCS,1),1);
for jz = 1:16
    regs(AHAseg(jz).topo) = jz;
end