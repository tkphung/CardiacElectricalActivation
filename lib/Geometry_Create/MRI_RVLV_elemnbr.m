function [HN,varargout] = MRI_RVLV_elemnbr(LVnodes,LVHEX,LVLCR,RVnodes,RVHEX,RVLCR)
%MRI_RVLV_elemnbr: determines H-N pairs for the MRI-derived LV+RV, also
%outputs the NODE and HEX information for the combined model
% [HN,VARARGOUT] = MRI_RVLV_elemnbr(LVLCR, RVLCR)
%   INPUT VARIABLES:
%       LVnodes- from LV_nodenumb
%       LVHEX- from LV_elemcon
%       LVLCR- direct input
%       RVnodes- from LV_nodenumb
%       RVHEX- from LV_elemcon
%       RVLCR- direct input
%   CALCULATIONS:
%       MRI_RVinLVmask function
% 		STEP 1: LV Epicardial Elements
% 		STEP 2: RV Insertion elements (Based on LV_elemnbr.m)
%             STEP 2a: Make 3D grid of elements for RV
%             STEP 2b: Pad with NaN on all borders
%             STEP 2c: Pad circumferentially
%             STEP 2d: Logical for existing elements (outside of LV)
%             STEP 2e: Move logical Circumferentially +1 & -1 and Radial towards Endo
%             STEP 2f: Do circplus and circminus neighbors exist?
%             STEP 2g: RV insertion elements
% 		STEP 3: Use nearestneighbour.m to find associated neighbors
% 		STEP 4: Organize Home Neighbor pairing
%       STEP 5: VARAGOUT formatting
%   OUTPUT VARIABLES:
%       HN- Home Neighbor pairs (unique) for the LV+RV element and node num
%       VARAGOUT: NODES- combo LV+RV nodes (note some RV nodes will not be
%                        in the HEX because of LV-RV overlap)
%                 HEX-   combo LV+RV HEX elements
%                 LVEpiInsertion elements
%
% Developed in script RVLV_MRIFit/MRI_LVRV_EP.m
% Thien-Khoi N. Phung (May 02, 2017)

[RVinLV, RVnbr] = MRI_RVinLVmask(LVnodes,LVLCR,RVnodes,RVHEX,RVLCR);

% Define LV neighbors (unique HN pairs)
[~, LVnbr] = LV_elemnbr(LVLCR(1),LVLCR(2),LVLCR(3));

% Define LV+RV connection
% STEP 1: LV Epicardial Elements
LVEpi = (LVLCR(1)*LVLCR(2)*(LVLCR(3)-1) + 1):(LVLCR(1)*LVLCR(2)*LVLCR(3));

% STEP 2:RV Insertion elements (Based on LV_elemnbr.m)
% STEP 2a: Make 3D grid of elements for RV
egrid = reshape(1:(RVLCR(2)*RVLCR(1)*RVLCR(3)),RVLCR(2),RVLCR(1),RVLCR(3));
% STEP 2b: Pad with NaN on all borders
eg = nan(RVLCR(2)+2, RVLCR(1)+2, RVLCR(3)+ 2);
eg(2:end-1, 2:end-1, 2:end-1) = egrid;
% STEP 2c: Pad circumferentially
eg(1,:,:)   = eg(end-1,:,:);
eg(end,:,:) = eg(2,:,:);
% STEP 2d: Logical for existing elements (outside of LV)
L = ismember(eg,find(~RVinLV==1));
% STEP 2e: Move logical Circumferentially +1 & -1 and Radial towards Endo
% Shift logical matrix                    C  L  R
circplus  = eg(L(    [2:end 1],:,:)); % [+1, 0, 0]
circminus = eg(L([end 1:end-1],:,:)); % [-1, 0, 0]
radplus   = eg(L(:,:,    [2:end 1])); % [ 0, 0,+1]
% Home element
H = eg(L(:));
% STEP 2g: Do circplus and circminus neighbors exist?
% Are any of the circumferential neighbors in the LV?
LL = ismember(circplus,find(RVinLV)) | ismember(circminus,find(RVinLV))...
     | ismember(radplus,find(RVinLV));
% STEP 2h: RV insertion elements
RVins = unique(H(LL));

% STEP 3: Use nearestneighbour.m to find associated neighbors
% Define HEX centers for the LVEpi elements
LVEpiHC = LV_hexcent(LVnodes,LVHEX(LVEpi,:));

% Sampling Points from RV Insertion
RVinsHC = LV_hexcent(RVnodes,RVHEX(RVins,:));

% Nearest neighbour function
LV2RVins = nearestneighbour(RVinsHC',LVEpiHC');
% LVEpi elements that are nearest neighbor to RVins
LV2RVins = LVEpi(LV2RVins)';

% STEP 4: Organize Home Neighbor pairings
% Shift RV element numbering
RVnbr = RVnbr + (LVLCR(1)*LVLCR(2)*LVLCR(3));
RVins = RVins + (LVLCR(1)*LVLCR(2)*LVLCR(3));

% Home Neighbor Pairs
HN = [LVnbr; RVnbr; RVins LV2RVins];

% VARARGOUT
% Combine Nodes & Hex definitions
varargout{1} = [LVnodes;RVnodes];
varargout{2} = [LVHEX; (RVHEX + size(LVnodes,1))];
varargout{3} = LV2RVins; % idx of LV elements on RV insertion