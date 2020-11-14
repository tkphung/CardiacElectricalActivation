function [RVINLVMASK, RVnbr] = MRI_RVinLVmask(LVNODES,LVLCR,RVNODES,RVHEX,RVLCR)
%MRI_RVinLVmask: creates a logical mask that refers to the elements listed
%by RVHEX (RV element connectivity) for the RV elements within the LV
%myocardium due to segmentation. This information is used to separate the
%RV myocardium for electrical modeling
% RVINLVMASK = MRI_RVinLVmask(LVNODES,LVHEX,LVLCR,RVNODES,RVHEX,RVLCR)
%   INPUT VARIABLES:
%       LVNODES- NODES of LV (LV_nodenumb)
%       LVLCR- LV mesh density
%       RVNODES- NODES of RV (LV_nodenumb)
%       RVHEX- elem connectivity for HEX elements (LV_elemcon)
%       RVLCR- RV mesh density
%   CALCULATE:
% 		STEP 1: Calculate HEX centers of RV
% 		STEP 2: Limit search space of RV to the LV bounds
% 		STEP 2a: Find bounds of LV lateral to septal walls (y)
% 		STEP 2b: Find bounds of LV anterior to posterior walls (x)
% 		STEP 2c: Window RV elements that are within LV bounds
% 		STEP 3: Determine which RV HEX centers are inside the LV myocardium
% 		STEP 3a: Select LV epi NODES
% 		STEP 3b: Create a convexhull
% 		STEP 3c: Calculate which RV NODES are inhull
% 		STEP 3d: Connected Components- eliminate any single elements not surrounded by its nbrs
%   OUTPUT VARIABLE:
%       RVINLVMASK- logical same length as RVHEX referring to which
%       elements are in the LV myocardium
% 		RVnbr- Unique Home Neighbor pairs
% Thien-Khoi N. Phung (March 29, 2017)
% TNP (May 1, 2017)- added conncomp function

% STEP 1: Calculate HEX centers of RV
RVHEXc = LV_hexcent(RVNODES,RVHEX);

% STEP 2: Limit search space of RV to the LV bounds
% STEP 2a: Find bounds of LV lateral to septal walls (y)
LVylim = [min(LVNODES(:,2)) max(LVNODES(:,2))];

% STEP 2b: Find bounds of LV anterior to posterior walls (x)
LVzlim = [min(LVNODES(:,3)) max(LVNODES(:,3))];

% STEP 2c: Window RV elements that are within LV bounds
RVinLV = RVHEXc(:,2)>LVylim(1) & RVHEXc(:,2)<LVylim(2)...
       & RVHEXc(:,3)>LVzlim(1) & RVHEXc(:,3)<LVzlim(2);
RVinLVind = find(RVinLV==1);

% STEP 3: Determine which RV HEX centers are inside the LV myocardium
% STEP 3a: Select LV epi NODES
%          Number of Layer NODES + Apex NODES
idLVepi = ((LVLCR(1)+1)*LVLCR(2)*LVLCR(3) + 1*LVLCR(3)):size(LVNODES,1);

% STEP 3b: Create a convexhull
K = convhulln(LVNODES(idLVepi,:));

% Visualize convex hull and windowed RV hex centers
% figure, hold on
% trisurf(K,LVNODES(idLVepi,1),LVNODES(idLVepi,2),LVNODES(idLVepi,3))
% alpha(0.25)
% plot3(RVHEXc(RVinLV,1),RVHEXc(RVinLV,2),RVHEXc(RVinLV,3),'r.','MarkerSize',15)
% view(90,0), axis equal

% STEP 3c: Calculate which RV NODES are inhull
RVinLV = inhull(RVHEXc(RVinLV,:),LVNODES(idLVepi,:),K);
RVinLVind = RVinLVind(RVinLV);
% % Convert back to a logical
% RVinLV = ismember(1:size(RVHEXc,1),RVinLVind);

% STEP 3d: Connected Components- eliminate any single elements not surrounded by its
% neighbors
% Define all of the neighbor pairs
[~,RVnbr] = LV_elemnbr(RVLCR(1),RVLCR(2),RVLCR(3));
% Filter out any non-existent elements from the RVinLV mask
RVnbrTRUE = ~ismember(RVnbr(:,1),RVinLVind) & ~ismember(RVnbr(:,2),RVinLVind);
RVnbr = RVnbr(RVnbrTRUE,:);

% Connected Components Function
G = graph(RVnbr(:,1),RVnbr(:,2));
BINS = conncomp(G);
connected = find(BINS == 1); % Assume CC 1 is the RV myocardium
RVnbrTRUE = ismember(RVnbr(:,1),connected) & ismember(RVnbr(:,2),connected);
RVnbr = RVnbr(RVnbrTRUE,:);

% Sort into structure- for each element- what are its neighbors?
H = [RVnbr(:,1); RVnbr(:,2)];
N = [RVnbr(:,2); RVnbr(:,1)];
for el = 1:(RVLCR(1)*RVLCR(2)*RVLCR(3))
    E{el} = N(H==el);
end
% If an element has no neighbors- eliminate is from the RVinLVind
RVINLVMASK = cellfun(@isempty,E);


% Visualize which RV centers are inside(red) and outside(blue) the LV
% figure, hold on
% trisurf(K,LVNODES(idLVepi,1),LVNODES(idLVepi,2),LVNODES(idLVepi,3))
% alpha(0.25)
% plot3(RVHEXc(RVinLV,1),RVHEXc(RVinLV,2),RVHEXc(RVinLV,3),'r.','MarkerSize',15)
% plot3(RVHEXc(~RVinLV,1),RVHEXc(~RVinLV,2),RVHEXc(~RVinLV,3),'b.','MarkerSize',15)
% view(90,0), axis equal