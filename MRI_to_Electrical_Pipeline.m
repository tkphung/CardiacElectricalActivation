% MRI_to_Electrical_Pipeline.m
% This scipt takes cardiac MRI segmented using Medviso Segment and
% generates a biventricular finite element model suitable for electrical
% model.
% The script then simulates electrical activation using a shortest path
% tree algorithm
%
% Thien-Khoi N. Phung (December 6, 2020)

% Steps for MRI model creation
%     Step 1: Load in MRI segmentation
%     Step 2: Fit segmentations to bicubic hermite surfaces
%     Step 3: Create Finite Element mesh
% Steps for Electrical Activation model
%     Step 4: Define fiber orientations
%     Step 5: Set element electrical groupings and mechanical surfaces
%     Step 6: Map scar on LV (if scar_flag is set to TRUE)
%             This is based on work from (Phung+ 2020 JBME)
%     Step 7: Set the velocities for the electrical model
%     Step 8: Set which element(s) initiate the electrical propagation
%     Step 9: Electrical Simulation


%% Step 1: Load in MRI segmentation
% The segmentation was all processed using Medviso's Segment. There should
% be 4 files:
%   In the short axis stack, at matched frame throughout the cardiac cycle:
% 1. LVepifile: the LV endocardium and epicardium segmented with RV
%               insertion points indicates in one slice
% 2. RVendofile: a copy of LVepifile except the LV epicardial segmentation
%                is replaced by segmentation of the RV endocardium
% 3. RVepifile: a copy of LVepifile except the LV epicardial segmentation is
%               replaced by segmentation of the RV epicardium
%   In one long axis view, at matched frame as short axis segmentation:
% 4. LAfile: points placed at endocardial Apex location and at base of LV
%            at the midpoint of the mitral valve

% Does this data set include LGE scar?
scar_flag = false; % CRT006 did not have LGE scar data

% The data was segmented at End Diastole- indicate the frame number because
% we segmented CINE MRI
frameED = 25; % HARDCODED for data subject CRT006

% Path to the segmentation files
SApath     = 'data/CRT006/CINE'; % Short axis segmentations
LApath     = LApath;             % Long axis segmentation

% Filenames
LVepifile  = 'SA.mat';
RVendofile = 'SA_RVendo.mat';
RVepifile  = 'SA_RVepi.mat';
LAfile     = 'LA.mat';

% Combine the file names
SAname = {LVepifile RVendofile RVepifile};


%% Step 2: Fit segmentations to bicubic hermite surfaces
% Labels for naming processed files
segs = {'LVepi', 'RVEndo', 'RVEpi'};
lvendo = 'LVendo';

% Shift and Rotate MRI segmentations to orient them in LV (Base-apex,
% RV insertions) coordinates; do this for each of the short axis files
for jz = 1:3
    ProcessedStack = MRI_ProcessStack_withRV([SApath SAname{jz}],[LApath LAfile],'time',frameED);
    
    save([SApath 'ProcessedStack_' segs{jz} '.mat'],'ProcessedStack');

    % FE Fitting- Bicubic Fit for Processed Contours
    % INPUT Take ProcessedStack and fit surface to Endo & Epi
    %       SApath- saves the Endo & Epi IPNODE files in the same place as MRI
    %       Meshdensity- '4x2' '4x4' '4x8' options
    % OUTPUT IPNODE files for Endo & Epi
    MRI_FitContours(ProcessedStack,SApath,'4x8','epiname',segs{jz},'endoname',lvendo);
end

% Visualize the bicubic hermite surfaces for the segmentation
H = figure('Name','Bicubic Hermite Fit','NumberTitle','off',...
    'WindowStyle','docked');
col = [214,96,77; 103,169,207; 33,102,172]./255;
MRI_SurfaceRender([SApath lvendo '.IPNODE'],H);
for jz = 1:3
    MRI_SurfaceRender([SApath segs{jz} '.IPNODE'],'color',col(jz,:),'handle',H);
end


%% Step 3: Create Finite Element mesh
% LV SURFACES
% IPNODE Files
ENDOfile = [SApath lvendo '.IPNODE'];
EPIfile  = [SApath segs{1} '.IPNODE'];

% Mesh Density for Left Ventricle
LVLe = 28;    % rings base to apex (EXCLUDES APEX CAP element)
LVCe = 30;    % elements per ring
LVRe = 5;     % number of elements thru wall
LVLCR = [LVLe LVCe LVRe];

% Create Nodes
LVgeom = MRI_FEMeshRender(ENDOfile,EPIfile,LVLCR);

% Meshing Geometry- Reorder geom into NODES
[LVnodes] = LV_nodenumb(LVgeom);

% Element connectivity matrix with regards to element numbers
[LVHEX,LVPENT] = LV_elemcon(LVLe,LVCe,LVRe);

% RV SURFACES
% IPNODE Files
ENDOfile = [SApath segs{2} '.IPNODE'];
EPIfile  = [SApath segs{3} '.IPNODE'];

% Mesh Density for Right Ventricle
RVLe = LVLe;
RVCe = 96;
RVRe = LVRe;
RVLCR = [RVLe RVCe RVRe];

% Create Nodes
RVgeom = MRI_FEMeshRender(ENDOfile,EPIfile,RVLCR);
RVgeom = RV_MVcut(RVgeom,min(LVnodes(:,1)));

% Meshing Geometry- Reorder geom into NODES
[RVnodes] = LV_nodenumb(RVgeom);

% Element connectivity matrix with regards to element numbers
[RVHEX,RVPENT] = LV_elemcon(RVLe,RVCe,RVRe);

% Connect the LV+RV
[HN,NODES,HEX,LV2RVins] = MRI_RVLV_elemnbr(LVnodes,LVHEX,LVLCR,RVnodes,RVHEX,RVLCR);

% Index for elements used in the Electrical Model
EPnodeIDX = unique(HN(:));

% Calculate center point of each elements
HEXc = LV_hexcent(NODES,HEX);
HEXv = LV_hexvol(NODES,HEX);

% Store model information
MODEL.LVLCR = LVLCR; % Mesh dimensions
MODEL.RVLCR = RVLCR;
MODEL.HEX   = HEX;   % Element connectivity matrix
MODEL.NODES = NODES; % Node coordinates
MODEL.HEXc  = HEXc;  % Element Centers
MODEL.HEXv  = HEXv;  % Element volumes
MODEL.HN    = HN;    % Home-Neighbor pairs
MODEL.EPnodeIDX = EPnodeIDX; % Index of BiV elements (for electrical model)
MODEL.focus = LVgeom.d; % Focus
MODEL.LVPENT = LVPENT;


%% Step 4: Define fiber orientations
% Fiber directions (linear interpolation between epi-to-endo walls)
enfib = 60;  % Fiber angle in the plane of the wall at endocardium
epfib = -60; % Fiber angle in the plane of the wall at epicardium
[LVfib,LVfcr] = LV_fibers(LVHEX,LVnodes,LVLe,LVCe,LVRe,enfib,epfib);
[RVfib,RVfcr] = LV_fibers(RVHEX,RVnodes,RVLe,RVCe,RVRe,enfib,epfib);
fibvec = [LVfib;RVfib];
Qfcr   = cat(3,LVfcr,RVfcr);

% Store model information
MODEL.fiber = LVfib; % Store only LV fibers (for mechanical model)
MODEL.Qfcr  = Qfcr;  % Store Rotation Matrix (for electrical model)


%% Step 5: Set element electrical groupings and mechanical surfaces
% Element groupings (for electrical model)
% Purkinje Element Index
% LV Endocardial Layer Purkinje
LVPurk = 1:LVLe*LVCe; % LV Endocardial Elements

% RV Freewall Endocardial Purkinje
RVPurk = (LVLe*LVCe*LVRe) + (1:RVLe*RVCe); % RV Endocardial Elements

% RV Septal Purkinje
% Restructure logical of LV epi elements for the RV insertion mask
LV2RVinsLogical = zeros(LVLe*LVCe,1);
LV2RVinsLogical(LV2RVins-(LVLe*LVCe*(LVRe-1))) = 1;
LVEpiLogicalMatrix = logical(reshape(LV2RVinsLogical,LVCe,LVLe));
% Left column is base- Right Column is Apex
% Fill in the Base-Apex Purkinje mask on the RV Septum
LVEpiLogicalMatrix = imfill(LVEpiLogicalMatrix,[1 1]);
LVEpiLogicalMatrix = imfill(LVEpiLogicalMatrix,[LVCe 1]);
      
LVEpiMatrix = reshape((LVLe*LVCe*(LVRe-1)+1):(LVLe*LVCe*LVRe),LVCe,LVLe);
RVSeptalPurk = LVEpiMatrix(logical(LVEpiLogicalMatrix));
RVPurk = RVPurk(:);
RVPurk = [RVPurk(ismember(RVPurk,EPnodeIDX)); RVSeptalPurk(:)];

% Create Indices for Epicardial Elements
% Find Epicardial LV Wall
LVEpiWall = LVEpiMatrix(~logical(LVEpiLogicalMatrix));
% Find Epicardial RV wall
RVEpiWall = (LVLe*LVCe*LVRe) + (RVLe*RVCe*(RVRe-1)) + (1:RVLe*RVCe);
iRVEpiWall = ismember(RVEpiWall,EPnodeIDX);
RVEpiWall = RVEpiWall(iRVEpiWall);
Epicardium = [LVEpiWall(:); RVEpiWall(:)];

LV_WALL = 1:(LVLe*LVCe*LVRe);
removePurk = ~ismember(LV_WALL,[LVPurk(:) ;RVSeptalPurk(:)]);
LV_WALL = LV_WALL(removePurk);

% Groupings: LVPurk
%            RVPurk
%            -- The rest is bulk myocardium
%           
%            Epicardium (for plotting)

% Store model information (for electrical model)
MODEL.elem.LVP = LVPurk;
MODEL.elem.RVP = RVPurk;
MODEL.elem.LVM = LV_WALL;

% Surfaces for loading (for mechanical model)
RVSe = RVPurk(RVPurk<prod(LVLCR));

% Define Quad Surface Connectivity for RV septum
RVquad = LVHEX(RVSe,[3 4 8 7]);

% Define surfaces connectivity for LV endocardium
[LVQUAD,LVTRI] = LV_innersurf(LVLe,LVCe,LVRe);

% Store model information (for mechanical model)
MODEL.surf.RVQUAD = RVquad;
MODEL.surf.LVQUAD = LVQUAD;
MODEL.surf.LVTRI  = LVTRI;


%% Step 6: Map scar on LV (if scar_flag is set to TRUE)
if scar_flag
% EPI node numbers
nEPI = HEX(LVEpiMatrix(:),[3 4 7 8]);

% Load Short Axis LGE
% [LGEfile, SApath2] = uigetfile('*.mat','Pick short axis LGE');
LGEfile = 'SA.mat';
load([SApath2 LGEfile])

% Pre-process struct of data
% Add KeptSlices and parse Scar information 
setstruct = MRI_ProcessSegment(setstruct);

% Process RV insertion points from Short Axis
setstruct = MRI_RVInsertion(setstruct);

% Rotate all data into MRI imaging coordinates
SAxyz = MRI_Rotate2MRI(setstruct);

% slices with scar
slices = unique(SAxyz.scar(:,end))';
% keep only slices with scar
for jz = 1:numel(slices)
    slice = slices(jz);
    SAxis{jz}.Endo = SAxyz.endo(SAxyz.endo(:,end)==slice,1:3);
    SAxis{jz}.Epi = SAxyz.epi(SAxyz.epi(:,end)==slice,1:3);
    SAxis{jz}.Scar = SAxyz.scar(SAxyz.scar(:,end)==slice,1:3);
end
SAxis{jz+1}.pts  = SAxyz.pnpt;

% DATA long axis segmentation with scar
LAfiles = {[SApath2 'LA_2CH_scar.mat'];
           [SApath2 'LA_3CH_scar.mat'];
           [SApath2 'LA_4CH_scar.mat']};
for jz = 1:numel(LAfiles)

    load(LAfiles{jz})

    % Pre-process struct of data
    % Add KeptSlices and parse Scar information 
    setstruct = MRI_ProcessSegment(setstruct);
    
    % Rotate all data into MRI imaging coordinates
    LAxyz = MRI_Rotate2MRI(setstruct);
    
    LAxis{jz}.Endo = LAxyz.endo;
    LAxis{jz}.Epi  = LAxyz.epi;
    LAxis{jz}.Scar = LAxyz.scar;
end

% Long Axis pin points into Imaging Coordinates
LAfilepnpt = [SApath2 'LA_4CH_pnpt.mat'];
load(LAfilepnpt)
LAxis{jz+1}.pts = MRI_Rotate2MRI(setstruct);

% ROTATION of all data to cardiac coordinate system
% Pinpoints
A = LAxis{end}.pts(1,:);
B = LAxis{end}.pts(2,:);
Si = SAxis{end}.pts(3,:);

    % First basis vector, subtract Base from Apex and divide by magnitude
    C = A-B;
    e1 = C ./ norm(C);

    % Calculate Origin location
    origin = B + (1/3)*C;
 
    % Second basis vector using NEW METHOD- plane intersects septal point &
    % e1
    D = Si(1,:) - origin;
    D2 = D - dot(D,e1)*e1;
    e2 = D2 ./ norm(D2);

    % Third basis vector
    E = cross(e1,e2);
    e3 = E ./ norm(E);

    % Transformation matrix
    Transform = [e1; e2; e3];
    
G = figure('WindowStyle','docked'); hold on
axis equal tight,xlabel('X'),ylabel('Y'),zlabel('Z')
title('Data in Cardiac Coordinates')

for jz = 1:numel(SAxis)-1
    SAxis{jz}.REndo(:,1:3) = (SAxis{jz}.Endo(:,1:3) - origin)*Transform';
    SAxis{jz}.REpi(:,1:3)  = (SAxis{jz}.Epi(:,1:3) - origin)*Transform';
    SAxis{jz}.RScar(:,1:3) = (SAxis{jz}.Scar(:,1:3) - origin)*Transform';
    
    plot3(SAxis{jz}.REndo(:,1),SAxis{jz}.REndo(:,2),SAxis{jz}.REndo(:,3),'r-')
    plot3(SAxis{jz}.REpi(:,1),SAxis{jz}.REpi(:,2),SAxis{jz}.REpi(:,3),'k-')
    plot3(SAxis{jz}.RScar(:,1),SAxis{jz}.RScar(:,2),SAxis{jz}.RScar(:,3),'m.','MarkerSize',20)
end
SAxis{jz+1}.Rpts(:,1:3)  = (SAxis{jz+1}.pts(:,1:3) - origin)*Transform';
plot3(SAxis{jz+1}.Rpts(:,1),SAxis{jz+1}.Rpts(:,2),SAxis{jz+1}.Rpts(:,3),'k.','MarkerSize',20)

for jz = 1:numel(LAxis)-1
    LAxis{jz}.REndo(:,1:3) = (LAxis{jz}.Endo(:,1:3) - origin)*Transform';
    LAxis{jz}.REpi(:,1:3)  = (LAxis{jz}.Epi(:,1:3) - origin)*Transform';
    LAxis{jz}.RScar(:,1:3) = (LAxis{jz}.Scar(:,1:3) - origin)*Transform';
    
    plot3(LAxis{jz}.REndo(:,1),LAxis{jz}.REndo(:,2),LAxis{jz}.REndo(:,3),'r-')
    plot3(LAxis{jz}.REpi(:,1),LAxis{jz}.REpi(:,2),LAxis{jz}.REpi(:,3),'k-')
    plot3(LAxis{jz}.RScar(:,1),LAxis{jz}.RScar(:,2),LAxis{jz}.RScar(:,3),'m.','MarkerSize',20)
end
LAxis{jz+1}.Rpts(:,1:3)  = (LAxis{jz+1}.pts(:,1:3) - origin)*Transform';
plot3(LAxis{jz+1}.Rpts(:,1),LAxis{jz+1}.Rpts(:,2),LAxis{jz+1}.Rpts(:,3),'k.-','MarkerSize',20)

% CONVERT Cardiac Data to Prolate coordinates
% Short Axis
for jz = 1:numel(SAxis)-1
    [L,M,T] = LV_C2P(SAxis{jz}.REndo(:,1),SAxis{jz}.REndo(:,2),SAxis{jz}.REndo(:,3),LVgeom.d);
    SAxis{jz}.LMTEndo = [L M T];
    [L,M,T] = LV_C2P(SAxis{jz}.REpi(:,1),SAxis{jz}.REpi(:,2),SAxis{jz}.REpi(:,3),LVgeom.d);
    SAxis{jz}.LMTEpi = [L M T];
    [L,M,T] = LV_C2P(SAxis{jz}.RScar(:,1),SAxis{jz}.RScar(:,2),SAxis{jz}.RScar(:,3),LVgeom.d);
    SAxis{jz}.LMTScar = [L M T];
end

% Long Axis
for jz = 1:numel(LAxis)-1
    [L,M,T] = LV_C2P(LAxis{jz}.REndo(:,1),LAxis{jz}.REndo(:,2),LAxis{jz}.REndo(:,3),LVgeom.d);
    LAxis{jz}.LMTEndo = [L M T];
    [L,M,T] = LV_C2P(LAxis{jz}.REpi(:,1),LAxis{jz}.REpi(:,2),LAxis{jz}.REpi(:,3),LVgeom.d);
    LAxis{jz}.LMTEpi = [L M T];
    [L,M,T] = LV_C2P(LAxis{jz}.RScar(:,1),LAxis{jz}.RScar(:,2),LAxis{jz}.RScar(:,3),LVgeom.d);
    LAxis{jz}.LMTScar = [L M T];
end

% SAMPLE Geometry
[qL,qM,qT] = LV_C2P(LVnodes(:,1),LVnodes(:,2),LVnodes(:,3),LVgeom.d);

T = figure('WindowStyle','docked'); hold on
xlabel('\theta'),ylabel('\mu')
title('Data Locations on Epi')
colormap(flipud(colormap('parula')))
colorbar
% Plot SAMPLE points
plot(qT(nEPI(:)),qM(nEPI(:)),'ko','MarkerSize',6)
cardiacdata = [];
for jz = 1:numel(SAxis)-1
    % Short Axis use THETA bins
    numbins  = 50;
    angedges = linspace(0,2*pi,numbins+1);
    angcents = (angedges(1:end-1) + angedges(2:end))./2;
    
    % Interpolate LAMBDA and MU for EPI
    LM = interp1(SAxis{jz}.LMTEpi(:,3),SAxis{jz}.LMTEpi(:,1:2),angcents,'linear');
    
    % Group data into THETA bins
    [enN,~,enBIN] = histcounts(SAxis{jz}.LMTEndo(:,3),angedges);
    [epN,~,epBIN] = histcounts(SAxis{jz}.LMTEpi(:,3),angedges);
    [scN,~,scBIN] = histcounts(SAxis{jz}.LMTScar(:,3),angedges);
    
    % Wall Thickness | Scar Transmurality | Scar Depth from EPI
    WTHSCAR = nan(numbins,3);
    for wh = 1:numbins
        % Average Endo and Epi points in Cartesian Coordinates
        endopoint = mean(SAxis{jz}.REndo(enBIN==wh,:),1);
        epipoint  = mean(SAxis{jz}.REpi(epBIN==wh,:),1);
        
        % Store Wall Thickness
        WTHSCAR(wh,1) = sqrt(sum((epipoint-endopoint).^2));
        
        % Scar Points (only if more than 2 pixels scar)
        if sum(scBIN==wh)>2
            % Scar Points Cartesian and Prolate Sph Values
            scarpoints = [SAxis{jz}.RScar(scBIN==wh,:) SAxis{jz}.LMTScar(scBIN==wh,:)];
            
            % endo and epi most scar points determined by LAMDA
            scarpoints = sortrows(scarpoints,4);
            endosp = scarpoints(1,1:3);
            episp  = scarpoints(end,1:3);
            
            % Scar transmurality
            WTHSCAR(wh,2) = sqrt(sum((episp-endosp).^2))/WTHSCAR(wh,1);
            
            % Scar depth as percent from epi depth
            WTHSCAR(wh,3) = sqrt(sum((epipoint-episp).^2))/WTHSCAR(wh,1);
        end
    end
    
    % Store Information
    % [L M T WTH SCARTRANS SCARDEPTH]
    SAxis{jz}.wthscar = [LM angcents(:) WTHSCAR];
    
    data = [LM angcents(:) WTHSCAR];
    data(isnan(data(:,5)),5) = 0;
    cardiacdata = [cardiacdata; data];
    
    % Plot SA data points
    scatter(SAxis{jz}.wthscar(:,3),SAxis{jz}.wthscar(:,2),100,SAxis{jz}.wthscar(:,5),'filled')
end

% Assign data to epicardium using Prolate coordinates (LONG AXES)
for jz = 1:numel(LAxis)-1
    % Long Axis use MU bins
    numbins  = 20;
    angedges = linspace(0,120*(pi/180),numbins+1);
    angcents = (angedges(1:end-1) + angedges(2:end))./2;
    
    % Group LA data in two theta bins (since there are two sets of lambdas)
    mdpt = mean(LAxis{jz}.LMTEpi(:,3));
    edges = [mdpt-pi mdpt mdpt+pi];
    [~,~,Hepi]  = histcounts(LAxis{jz}.LMTEpi(:,3),edges);
    [~,~,Hendo] = histcounts(LAxis{jz}.LMTEndo(:,3),edges);
    [~,~,Hscar] = histcounts(LAxis{jz}.LMTScar(:,3),edges);
    
    % For each Hemisphere
    data = [];
    for by = 1:2
        LMTendo = LAxis{jz}.LMTEndo(Hendo==by,:);
        LMTepi  = LAxis{jz}.LMTEpi(Hepi==by,:);
        LMTscar = LAxis{jz}.LMTScar(Hscar==by,:);
        Rendo = LAxis{jz}.REndo(Hendo==by,:);
        Repi  = LAxis{jz}.REpi(Hepi==by,:);
        Rscar = LAxis{jz}.RScar(Hscar==by,:);
        
        % Interpolate LAMBDA and THETA for EPI
        LT = interp1(LMTepi(:,2),LMTepi(:,[1 3]),angcents,'linear','extrap');

        % Group data into MU bins
        [enN,~,enBIN] = histcounts(LMTendo(:,2),angedges);
        [epN,~,epBIN] = histcounts(LMTepi(:,2),angedges);
        [scN,~,scBIN] = histcounts(LMTscar(:,2),angedges);
    
        % Wall Thickness | Scar Transmurality | Scar Depth from EPI
        WTHSCAR = nan(numbins,3);
        for wh = 1:numbins
            % Average Endo and Epi points in Cartesian Coordinates
            endopoint = mean(Rendo(enBIN==wh,:),1);
            epipoint  = mean(Repi(epBIN==wh,:),1);

            % Store Wall Thickness
            WTHSCAR(wh,1) = sqrt(sum((epipoint-endopoint).^2));

            % Scar Points (only if more than 2 pixels scar)
            if sum(scBIN==wh)>2
                % Scar Points Cartesian and Prolate Sph Values
                scarpoints = [Rscar(scBIN==wh,:) LMTscar(scBIN==wh,:)];

                % endo and epi most scar points determined by LAMDA
                scarpoints = sortrows(scarpoints,4);
                endosp = scarpoints(1,1:3);
                episp  = scarpoints(end,1:3);

                % Scar transmurality
                WTHSCAR(wh,2) = sqrt(sum((episp-endosp).^2))/WTHSCAR(wh,1);

                % Scar depth as percent from epi depth
                WTHSCAR(wh,3) = sqrt(sum((epipoint-episp).^2))/WTHSCAR(wh,1);
            end
        end
        % Store Information
        % [L M T WTH SCARTRANS SCARDEPTH]
        data = [data; LT(:,1) angcents(:) LT(:,2) WTHSCAR];
    end
    LAxis{jz}.wthscar = data;
    data(isnan(data(:,5)),5) = 0;
    cardiacdata = [cardiacdata; data(~isnan(data(:,5)),:)];
    
    % Plot LA data points
    figure(T); hold on
    scatter(LAxis{jz}.wthscar(:,3),LAxis{jz}.wthscar(:,2),100,LAxis{jz}.wthscar(:,5),'filled')
end

% INTERPOLATE Scar Transmurality
% Pad Cardiac Data Circumferentially
interpdata = cardiacdata(~isnan(cardiacdata(:,2)),[3 2 5 6]); % Theta Mu Transmurality
interpdata = [interpdata;
              interpdata(:,1)+2*pi interpdata(:,2:end);
              interpdata(:,1)-2*pi interpdata(:,2:end)];

% Interpolant Function
transmurality = scatteredInterpolant(interpdata(:,1),interpdata(:,2),interpdata(:,3),'linear');
nonnan = ~isnan(interpdata(:,4));
depth  = scatteredInterpolant(interpdata(nonnan,1),interpdata(nonnan,2),interpdata(nonnan,4),'linear');

% Interpolate at SAMPLE Points
EPInodes = LVnodes(nEPI(:),:);
[~,MEPIn,TEPIn]  = LV_C2P(EPInodes(:,1),EPInodes(:,2),EPInodes(:,3),LVgeom.d);

% Interpolate
scartr = transmurality(TEPIn,MEPIn);
scardp = depth(TEPIn,MEPIn);


% VISUALIZE THE SCAR
% Manipulate to fit to Element data (data was interopolated to nodes)
% Each row is an element (corresponding to eEPI)
% Each column is the nodes making up that elements Epicardial surface
scartrans = reshape(scartr,numel(scartr)/4,4);
scardepth = reshape(scardp,numel(scartr)/4,4);

% Scar depth from Epi has no meaning when there is zero scar
% So we set any scartrans<=0 nodes to have nan scardepth
scardepth(scartrans<=0) = nan;

% Taking the mean of each row consolidates the nodal data to an individual
% element data value.
eEPItr = mean(scartrans,2);
eEPIwd = nanmean(scardepth,2); % mean that ignores nans

% Data projected on epicardium
FE_RenderPatchMesh(LVnodes,HEX(LVEpiMatrix(:),:),'data',eEPItr);
colormap(flipud(viridis()))

% How many elements deep is scar?
SCARlayers = round(eEPItr * LVRe);

% What layer does scar start (from epicardium)?
EPIstart = round(eEPIwd*LVRe) + 1;
%       # elems not scar from epi + 1 <- for where scar layer starts

% if scar is too thin (0 elements deep)- remove EPIstart value
EPIstart(SCARlayers==0) = nan;

% if SCARlayers and EPIstart don't add up to LVGEOM.LVLCR(3)+1 due to round error:
% Adjust SCARlayers
SCARlayers((SCARlayers+EPIstart)>(LVRe+1)) = ...
    SCARlayers((SCARlayers+EPIstart)>(LVRe+1)) - 1;

% scar mask
epl = LVLe*LVCe; % elements per layer
eEPI = LVEpiMatrix(:);

eSCAR = [];
for jz = 1:max(EPIstart)  
    outerscar = eEPI(EPIstart==jz) - epl*(jz-1);
    thruwallscar = SCARlayers(EPIstart==jz);
    thruwallscar(thruwallscar>(LVRe-jz+1)) = (LVRe-jz+1);
    
    for wh = 1:numel(outerscar)
        for by = 1:thruwallscar(wh)
            eSCAR = [eSCAR; outerscar(wh)- epl*(by-1)];
        end 
    end
end

H = FE_RenderPatchMesh(NODES,HEX,'elements',EPnodeIDX,'alpha',0,...
                       'edgealpha',0.15);
FE_RenderPatchMesh(NODES,HEX,'elements',eSCAR,'facecolor',[1 0 0],...
                   'handle',H);
    
% Remove scar HN pairs
HN = MODEL.HN;
HNscar = bsxfun(@or,ismember(HN(:,1),eSCAR),ismember(HN(:,2),eSCAR));
HN(HNscar,:) = [];

    % Store information
    MODEL.scar = eSCAR;
    MODEL.HN = HN;
    MODEL.EPnodeIDX = unique(HN(:));
else
    MODEL.scar = [];
end


%% Step 7: Set the velocities for the electrical model
% Set muscle fiber conduction velocity
% This value is arbitrary and can be scaled to match the QRS duration of
% the data. Changing this value only changes the absolute electrical
% activation times, but does not change relative activation of early and
% late regions.
v = 1;

% Set velocities (based on paper- Lee+ 2019)
% This sets the conduction velocities in 3D for all of the elements
% LV- left ventricle
% RV- right ventricle
%   M- myocardium
%   P- purkinje
%    f- fiber direction
%    c- cross-fiber direction (and also the radial direction)
LVMf = v;
LVMc = 0.4 * v;
LVPf = 6.0 * v;
LVPc = 0.4 * v;
RVMf = v;
RVMc = 0.4 * v;
RVPf = 6.0 * v;
RVPc = 0.4 * v;

% Assign all elements to have the RV conduction properties
VELO = repmat([RVMf RVMc RVMc],size(MODEL.HEXc,1),1);

% Reassign RV Purkinje
RVP = MODEL.elem.RVP;
VELO(RVP,:) = repmat([RVPf RVPc RVPc],numel(RVP),1);

% Reassign LV Purkinje
LVP = MODEL.elem.LVP;
VELO(LVP,:) = repmat([LVPf LVPc LVPc],numel(LVP),1);

% Reassign LV Muscle
LVM = MODEL.elem.LVM;
VELO(LVM,:) = repmat([LVMf LVMc LVMc],numel(LVM),1);

% Add the conduction velocities (VELO) to the MODEL struct variable
% (This is used in the pseudoECG calculation)
MODEL.VELO = VELO;


%% Step 8: Set which element(s) initiate the electrical propagation
% This matrix will be an n by 2 matrix that includes
%   elements that initiate propagation in column 1 
%   time at which the elements initate in column 2

STAR = [1 0];

% Visualize the element selected for STAR
SS = FE_RenderPatchMesh(MODEL.NODES,MODEL.HEX,'elements',MODEL.EPnodeIDX,'alpha',0,'edgealpha',0.25);
FE_RenderPatchMesh(MODEL.NODES,MODEL.HEX,'elements',STAR(:,1),...
                                         'facecolor',[1 0 0],...
                                         'handle',SS);


%% Step 9: Electrical Simulation
% Pull other variables from MODEL struct
Qfcr = MODEL.Qfcr; % Element rotation matrices to specify direction of fiber direction
HEXc = MODEL.HEXc; % Center point of each element
HN   = MODEL.HN;   % This specifies which elements are "touching" (aka Home-Neighbor pairs)

% The function ElecpropSPT simulates the electrical propagation
% Returns litTIME- an electrical activation time for all elements
%         litHOO- element index for who activated the element
% both of these variables including for scar elements and those not in the model
[litTIME,litHOO] = ElecpropSPT(HEXc,HN,Qfcr,VELO,STAR);

% Keep only the electrical activation times for elements in the electrical
% model
MODEL.litTIME = litTIME(MODEL.EPnodeIDX);

% Visualize the electrical activation
FE_RenderPatchMesh(MODEL.NODES,MODEL.HEX,...
                   'elements',MODEL.EPnodeIDX,...
                   'data',MODEL.litTIME);
colormap(flipud(plasma()))
