function MRIGEOM = MRI_register(varargin)
%MRI_register: Registers LGE and DENSE Short Axis locations on specified 
%CINE geometry surface fit. This is used so that we register the LGE and
%DENSE MRI-derived data onto the same CINE geometry. You can input either
%LGE and/or DENSE short axis stacks to register.
% MRIGEOM = MRI_register(varargin)
%   INPUTS:
%       If isempty(varargin)- prompt file selection for each input
%       varargin: 'CINESA'  - filepath/name for SA MRI
%       (optional)'LGESA'   - filepath/name for SA MRI
%       (optional)'DENSESA' - filepath/name for SA MRI 
%                             (include only slices used in DENSEsnalysis)
%                 'EPI'     - filepath/name for EPI.ipnode file generated
%                 	          CINE Bicubic Fit at specified time (frame 1)
%       (optional)'ENDO'    - filepath/name for ENDO.ipnode file
%                 'CINEPS'  - filepath/name for CINE ProcessedStack MRI
%       (optional)'plot'    - default (false), if true- visualize slice
%                             locations on surface geometry
%   CALCULATIONS:
%       STEP 1: Process Inputs
%       STEP 2: Process Short Axis MRI files & register Base-Apex distances
%       STEP 3: Load in CINE registration & surface fit
%       STEP 4: Solve for Surface Lambdas for the range of M given T
%           STEP 4a: Load in surfaces & pre-process
%           STEP 4b: Set up surface mesh for interpolation   
%           STEP 4c: Interpolate Lambda for all M,T pairs
%       STEP 5: Fitting LGE and DENSE locations data
%   OUTPUT:
%       MRIGEOM.CINE  \
%              .LGE    |- includes coordinates: L,M,T,x,y,z for each slice
%              .DENSE /
%              .focus
% 
% Based on test script Register_MRI.m
% Created by Thien-Khoi N. Phung (July 13, 2017)
% Edited by TNP to add DENSE (endo registration) (Oct. 23, 2017)

%% STEP 1: Process Inputs
plot_flag = false; % flag for visualizing
LGE_flag  = false;
DENSE_flag = false;
ENDO_flag = false;
if isempty(varargin) % no input files- prompt file selector
    % Identify directory with SHORT AXIS MRI segmentations
    [SAname,SApath] = uigetfile('*.mat','Pick CINE SA segmentation file');
    CINESA = [SApath SAname];
    
    [SAname,SApath] = uigetfile('*.mat','Pick LGE SA segmentation file',SApath);
    if SAname ~= 0
        LGESA = [SApath SAname];
        LGE_flag = true;
    end
    
    [SAname,SApath] = uigetfile('*.mat','Pick DENSE SA segmentation file',SApath);
    if SAname ~= 0
        DENSESA = [SApath SAname];
        DENSE_flag = true;
    end
    
    % Identify location for CINE ProcessedStack.m
    [SAname,SApath] = uigetfile('*.mat','Pick CINE ProcessedStack file',SApath);
    CINEPS = [SApath SAname];
    
    % Identify location for CINE EPI.ipnode
    [SAname,SApath] = uigetfile('*.ipnode','Pick CINE EPI.ipnode file',SApath);
    EPI = [SApath SAname];
    
    plot_flag = true;
    
elseif numel(varargin)/2 >= 4 % correct number of inputs
    % Parse through varargin
    for jz = 1:2:numel(varargin)
        switch varargin{jz}
            case 'CINESA'
                CINESA = varargin{jz+1};
			case 'LGESA'
				LGESA = varargin{jz+1};
                LGE_flag = true;
            case 'DENSESA'
				DENSESA = varargin{jz+1};
                DENSE_flag = true;
            case 'EPI'
				EPI = varargin{jz+1};
            case 'ENDO'
				ENDO = varargin{jz+1};
                ENDO_flag = true;
            case 'CINEPS'
				CINEPS = varargin{jz+1};
            case 'plot'
                plot_flag = varargin{jz+1};
            otherwise
                error('ERROR: Specified input is not valid.')
        end
    end
else
    error('ERROR: Number of inputs is wrong.')
end

%% STEP 2: Process Short Axis MRI files & register Base-Apex distances
% CINE
load(CINESA);

% Add KeptSlices variable to structure (this variable indicates the number
% of slices (1:setstruct(jz).ZSize) per imaging view (jz))
for jz = 1:length(setstruct)
    setstruct(jz).KeptSlices = [1:setstruct(jz).ZSize];
end

% Adjust KeptSlices variable to eliminate untraced images so that 
% epi/endo + scar traces will correctly register with each other
% (Note: scar is only traced on images with endo & sepi contours)
no_trace = sum(isnan(setstruct.EndoX(:, :,:))); % counts slices with no trace
delete_slices = no_trace ~= 0; % which have tracings
delete_slices = sum(delete_slices,2) == size(delete_slices,2);
setstruct.KeptSlices(:,squeeze(delete_slices)) = []; % removes those slice numbers
CKS = setstruct.KeptSlices';

% mm coordinates for most basal image
CINEpos = setstruct.ImagePosition(:);

% thickness between slices
CINEthick = setstruct.SliceThickness + setstruct.SliceGap; % mm

% BA axis
CINEBAaxis = -cross(setstruct.ImageOrientation(1:3),setstruct.ImageOrientation(4:6));

% Location of slices
CINElocs = ((repmat(CKS,1,3) - 1).*CINEthick) .* repmat(CINEBAaxis,numel(CKS),1);
CINElocs = repmat(CINEpos',numel(CKS),1) + CINElocs;

CINEstart = CINElocs(1,:)';

% Shifted location from 0 and normalized
CINEnorm = CKS'.*CINEthick;
CINEnorm = CINEnorm - min(CINEnorm); normfactor = max(CINEnorm);
CINEnorm = CINEnorm./normfactor;

% LGE
if LGE_flag
load(LGESA);

% Add KeptSlices variable to structure (this variable indicates the number
% of slices (1:setstruct(jz).ZSize) per imaging view (jz))
for jz = 1:length(setstruct)
    setstruct(jz).KeptSlices = [1:setstruct(jz).ZSize];
end

% Adjust KeptSlices variable to eliminate untraced images so that 
% epi/endo + scar traces will correctly register with each other
% (Note: scar is only traced on images with endo & sepi contours)
no_trace = sum(isnan(setstruct.EndoX(:, :,:))); % counts slices with no trace
delete_slices = no_trace ~= 0; % which have tracings
delete_slices = sum(delete_slices,2) == size(delete_slices,2);
setstruct.KeptSlices(:,squeeze(delete_slices)) = []; % removes those slice numbers
LKS = setstruct.KeptSlices';

% mm coordinates for most basal image
LGEpos = setstruct.ImagePosition(:);

% thickness between slices
LGEthick = setstruct.SliceThickness + setstruct.SliceGap; % mm

% BA axis
LGEBAaxis = -cross(setstruct.ImageOrientation(1:3),setstruct.ImageOrientation(4:6));

% Location of slices
LGElocs = ((repmat(LKS,1,3) - 1).*LGEthick) .* repmat(LGEBAaxis,numel(LKS),1);
LGElocs = repmat(LGEpos',numel(LKS),1) + LGElocs;

% Project point onto CINE BA axis
CINE2LGE = LGElocs(1,:)'-CINEstart;
LGEpt = dot(CINE2LGE,CINEBAaxis); % CINEBAaxis mag is 1

% normalized based on CINE
LGEnorm = (LKS-min(LKS))'.*LGEthick + LGEpt;
LGEnorm = LGEnorm./normfactor;
end

% DENSE
if DENSE_flag
load(DENSESA);

% Add KeptSlices variable to structure (this variable indicates the number
% of slices (1:setstruct(jz).ZSize) per imaging view (jz))
for jz = 1:length(setstruct)
    setstruct(jz).KeptSlices = [1:setstruct(jz).ZSize];
end

DKS = setstruct.KeptSlices';

% mm coordinates for most basal image
DENSEpos = setstruct.ImagePosition(:);

% thickness between slices
DENSEthick = setstruct.SliceThickness + setstruct.SliceGap; % mm

% BA axis
DENSEBAaxis = -cross(setstruct.ImageOrientation(1:3),setstruct.ImageOrientation(4:6));

% Location of slices
DENSElocs = ((repmat(DKS,1,3) - 1).*DENSEthick) .* repmat(DENSEBAaxis,numel(DKS),1);
DENSElocs = repmat(DENSEpos',numel(DKS),1) + DENSElocs;

% Project point onto CINE BA axis
CINE2DENSE = DENSElocs(1,:)'-CINEstart;
DENSEpt = dot(CINE2DENSE,CINEBAaxis); % CINEBAaxis mag is 1

% normalized based on CINE
DENSEnorm = (DKS-min(DKS))'.*DENSEthick + DENSEpt;
DENSEnorm = DENSEnorm./normfactor;
end

%% STEP 3: Load in CINE registration & surface fit
load(CINEPS);

landmarks = true;
[H,xyz] = MRI_SegmentRender(ProcessedStack,landmarks);
MRI_SurfaceRender(EPI,H);
if ENDO_flag
    MRI_SurfaceRender(ENDO,H);
end

if ~plot_flag
    close(H);
end

% Focus from the EPI.ipnode file
[EPIfit,focus]=MRI_IPNodeRead(EPI);
if ENDO_flag
    [ENDOfit,focus]=MRI_IPNodeRead(ENDO);
end

% Process xyz from MRI_SegmentRender into a grid


% Convert segmentation Cartesian->Prolate
% xyz matrix: going down rows is longitudinal
%             going left to right is circumferential
%             each page is X,Y,Z
% Cycle through each slice and create a matched-theta set of points
L = zeros(numel(xyz.x),49,1);
M = L;
T = L;
thetas = linspace(0,2*pi,50); thetas = thetas(1:end-1);
for jz = 1:numel(xyz.x)
    xseg = xyz.x{jz};
    yseg = xyz.y{jz};
    zseg = xyz.z{jz};
    
    [lraw,mraw,traw] = LV_C2P(xseg,yseg,zseg,focus);
    
    % Sort by theta
    LMT = sortrows([lraw mraw traw],3);
    
    % Pad in theta (circumferentially)
    nts = size(LMT,1);
    pLMT = [LMT((nts-2):nts,1:2) LMT((nts-2):end,3)-2*pi;
            LMT;
            LMT(1:3,1:2) LMT(1:3,3)+2*pi];
    
    % Interpolate thetas
    L(jz,:) = interp1(pLMT(:,3),pLMT(:,1),thetas,'linear');
    M(jz,:) = interp1(pLMT(:,3),pLMT(:,2),thetas,'linear');
    T(jz,:) = thetas;
end

[x,y,z] = LV_P2C(L,M,T,focus);

%% STEP 4: Solve for Surface Lambdas for the range of M given T
% Lambda(Mu,Theta)
% Based on MRI_FEMeshRender.m

% STEP 4a: Load in surfaces & pre-process
    % sort FE based on (lambda, mu, theta) & convert to radians
    EPIfit = sortrows(EPIfit,[3 2 1]);
    EPIfit(:,2:3) = EPIfit(:,2:3).*(pi/180);

    % determine size of bicubic mesh
    [m,n] = size(EPIfit);
    munum = length(find(EPIfit(:,3) == 0)); % # mu values
    thnum = m/munum;                         % # theta values
    
    % pad nodal mesh
    EPIfit(m+1:m+munum,:) = EPIfit(1:munum,:);
    EPIfit(m+1:m+munum,3) = 2*pi;
    thnum = thnum+1;   
    
    % reshape nodalmesh for INTERP function in STEP 3
    EPIdata  = reshape(EPIfit,munum,thnum,n);
    
% STEP 4b: Set up surface mesh for interpolation   
    Mvec = linspace(0,120,100).*(pi/180);
    Tvec = mean(T);
    
    [MUgrid,THgrid] = meshgrid(Mvec,Tvec);
        
    % set scaling derivative
    scale_der = 0; % not sure why this is 0 or what this does.. (TNP)
    
% STEP 4c: Interpolate Lambda for all M,T pairs

    EPIInterpSurf  = MRI_BiCubicInterp(MUgrid,THgrid,EPIdata,scale_der);
    % Output has dimension-3: lambda, mu, theta, & derivatives (look at
    % NodeData output written in MRI_BiCubicInterp)

    Ls = EPIInterpSurf(:,:,1);
    Ms = EPIInterpSurf(:,:,2); % Each row is one MU value
    Ts = EPIInterpSurf(:,:,3); % Each Column is one TH value


    [Xs,~,~] = LV_P2C(Ls,Ms,Ts,focus);
    
if ENDO_flag
    % STEP 4a: Load in surfaces & pre-process
        % sort FE based on (lambda, mu, theta) & convert to radians
        ENDOfit = sortrows(ENDOfit,[3 2 1]);
        ENDOfit(:,2:3) = ENDOfit(:,2:3).*(pi/180);

        % determine size of bicubic mesh
        [m,n] = size(ENDOfit);
        munum = length(find(ENDOfit(:,3) == 0)); % # mu values
        thnum = m/munum;                         % # theta values

        % pad nodal mesh
        ENDOfit(m+1:m+munum,:) = ENDOfit(1:munum,:);
        ENDOfit(m+1:m+munum,3) = 2*pi;
        thnum = thnum+1;   

        % reshape nodalmesh for INTERP function in STEP 3
        ENDOdata  = reshape(ENDOfit,munum,thnum,n);

    % STEP 4b: Set up surface mesh for interpolation   
        Mvec = linspace(0,120,100).*(pi/180);
        Tvec = mean(T);

        [MUgrid,THgrid] = meshgrid(Mvec,Tvec);

        % set scaling derivative
        scale_der = 0; % not sure why this is 0 or what this does.. (TNP)

    % STEP 4c: Interpolate Lambda for all M,T pairs

        ENDOInterpSurf  = MRI_BiCubicInterp(MUgrid,THgrid,ENDOdata,scale_der);
        % Output has dimension-3: lambda, mu, theta, & derivatives (look at
        % NodeData output written in MRI_BiCubicInterp)

        enLs = ENDOInterpSurf(:,:,1);
        enMs = ENDOInterpSurf(:,:,2); % Each row is one MU value
        enTs = ENDOInterpSurf(:,:,3); % Each Column is one TH value


        [enXs,~,~] = LV_P2C(enLs,enMs,enTs,focus);
end

%% STEP 5: Fitting LGE and DENSE locations data
% x y z are the cart coord for the CINE segs
% Lx Ly Lz are cart for LGE
% Dx Dy Dz are cart for DENSE
if LGE_flag
    Lx = zeros(numel(LGEnorm),size(x,2)); 
    Ly = Lx; 
    Lz = Lx; Ll = Lx; Lm = Lx; Lt = Lx;
end
if DENSE_flag
    Dx = zeros(numel(DENSEnorm),size(x,2));
    Dy = Dx;
    Dz = Dx; Dl = Dx; Dm = Dx; Dt = Dx;
end
if ENDO_flag
    enDx = zeros(numel(DENSEnorm),size(x,2));
    enDy = enDx;
    enDz = enDx; enDl = enDx; enDm = enDx; enDt = enDx;
end

for jz = 1:size(x,2)
    % x scaling using norms
    % x apex is max(x(:,1))
    % x basal is min(x(:,1))
    xscale = max(x(:,jz))-min(x(:,jz));
    CINEx  = CINEnorm.*xscale + min(x(:,jz));
    
    if LGE_flag
        LGEx   = LGEnorm.*xscale + min(x(:,jz));

        % fitting LGE and DENSE mu values
        LGEm = interp1(Xs(:,jz),Ms(:,jz),LGEx,'linear');

        % fitting LGE and DENSE lambda values
        LGEl = interp1(Xs(:,jz),Ls(:,jz),LGEx,'linear');

        LGEt = ones(numel(LGEm),1).*Ts(1,jz);

        [lx,ly,lz] = LV_P2C(LGEl(:),LGEm(:),LGEt(:),focus);

        Ll(:,jz) = LGEl(:); Lm(:,jz) = LGEm(:); Lt(:,jz) = LGEt(:);

        Lx(:,jz) = lx; Ly(:,jz) = ly; Lz(:,jz) = lz;
    end
    
    if DENSE_flag
        DENSEx = DENSEnorm.*xscale + min(x(:,jz));

        % fitting LGE and DENSE mu values
        DENSEm = interp1(Xs(:,jz),Ms(:,jz),DENSEx,'linear');

        % fitting LGE and DENSE lambda values
        DENSEl = interp1(Xs(:,jz),Ls(:,jz),DENSEx,'linear');

        DENSEt = ones(numel(DENSEm),1).*Ts(1,jz);

        [dx,dy,dz] = LV_P2C(DENSEl(:),DENSEm(:),DENSEt(:),focus);

        Dl(:,jz) = DENSEl(:); Dm(:,jz) = DENSEm(:); Dt(:,jz) = DENSEt(:);

        Dx(:,jz) = dx; Dy(:,jz) = dy; Dz(:,jz) = dz;
    end
    if DENSE_flag && ENDO_flag
        DENSEx = DENSEnorm.*xscale + min(x(:,jz));

        % fitting LGE and DENSE mu values
        enDENSEm = interp1(enXs(:,jz),enMs(:,jz),DENSEx,'linear');

        % fitting LGE and DENSE lambda values
        enDENSEl = interp1(enXs(:,jz),enLs(:,jz),DENSEx,'linear');

        enDENSEt = ones(numel(enDENSEm),1).*enTs(1,jz);

        [endx,endy,endz] = LV_P2C(enDENSEl(:),enDENSEm(:),enDENSEt(:),focus);

        enDl(:,jz) = enDENSEl(:); enDm(:,jz) = enDENSEm(:); enDt(:,jz) = enDENSEt(:);

        enDx(:,jz) = endx; enDy(:,jz) = endy; enDz(:,jz) = endz;
    end
end

% OUTPUT structure
MRIGEOM.focus = focus;
MRIGEOM.CINE.L = L;
MRIGEOM.CINE.M = M;
MRIGEOM.CINE.T = T;
MRIGEOM.CINE.x = x;
MRIGEOM.CINE.y = y;
MRIGEOM.CINE.z = z;
if LGE_flag
MRIGEOM.LGE.L = Ll;
MRIGEOM.LGE.M = Lm;
MRIGEOM.LGE.T = Lt;
MRIGEOM.LGE.x = Lx;
MRIGEOM.LGE.y = Ly;
MRIGEOM.LGE.z = Lz;
end
if DENSE_flag && ~ENDO_flag
MRIGEOM.DENSE.L = Dl;
MRIGEOM.DENSE.M = Dm;
MRIGEOM.DENSE.T = Dt;
MRIGEOM.DENSE.x = Dx;
MRIGEOM.DENSE.y = Dy;
MRIGEOM.DENSE.z = Dz;
end
if DENSE_flag && ENDO_flag
MRIGEOM.DENSE.EPI.L = Dl;
MRIGEOM.DENSE.EPI.M = Dm;
MRIGEOM.DENSE.EPI.T = Dt;
MRIGEOM.DENSE.EPI.x = Dx;
MRIGEOM.DENSE.EPI.y = Dy;
MRIGEOM.DENSE.EPI.z = Dz;
MRIGEOM.DENSE.ENDO.L = enDl;
MRIGEOM.DENSE.ENDO.M = enDm;
MRIGEOM.DENSE.ENDO.T = enDt;
MRIGEOM.DENSE.ENDO.x = enDx;
MRIGEOM.DENSE.ENDO.y = enDy;
MRIGEOM.DENSE.ENDO.z = enDz;
end

if plot_flag
    figure(H), hold on
    if LGE_flag
    plot3(Lx,Ly,Lz,'r<')
    end
    if DENSE_flag
    plot3(Dx,Dy,Dz,'bs')
    end
    
    figure, hold on
    plot(zeros(numel(CINEx),1).*xscale,CINEx,'ko-','LineWidth',2,'MarkerSize',10)
    if LGE_flag
    plot(0.1.*ones(numel(LGEx),1).*xscale,LGEx,'r<-','LineWidth',2,'MarkerSize',10)
    text(0.1.*xscale,LGEx(1),'LGE','FontWeight','b')
    end
    if DENSE_flag
    plot(0.2.*ones(numel(DENSEx),1).*xscale,DENSEx,'bs-','LineWidth',2,'MarkerSize',10)
    text(0.2.*xscale,DENSEx(1),'DENSE','FontWeight','b')
    end
    ylabel('X (mm)','FontWeight','b')
    text(0.*xscale,CINEx(1),'CINE','FontWeight','b')
    axis tight equal
    set(gca, 'YDir', 'reverse');
end