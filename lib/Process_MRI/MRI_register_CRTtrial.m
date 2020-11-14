function HEART = MRI_register_CRTtrial(varargin)
%MRI_register_CRTtrial: Registers MRI planes onto an epi surface of the
%heart. This is used to generate the virtual short-axis slices where the
%MRI planes are for the data sets with respect to the epicardial surface.
% HEART = MRI_register_CRTtrial(varargin)
%   INPUTS:
%       varargin: 'LGESA'   - filepath/name for SA MRI
%                 'DENSESA' - filepath/name for SA MRI 
%                 'EPI'     - filepath/name for EPI.ipnode file generated
%                 'LGEPS'   - filepath/name for CINE ProcessedStack MRI
%       (optional)'plot'    - default (false), if true- visualize slice
%                             locations on surface geometry
%   CALCULATIONS:
%       STEP 1: Process Inputs
%       STEP 2: Process Short Axis MRI files & register Base-Apex distances
%       STEP 3: Load in LGE registration & surface fit
%       STEP 4: Solve for Surface Lambdas for the range of M given T
%           STEP 4a: Load in surfaces & pre-process
%           STEP 4b: Set up surface mesh for interpolation   
%           STEP 4c: Interpolate Lambda for all M,T pairs
%       STEP 5: Fitting DENSE locations data
%   OUTPUT:
%       MRIGEOM.LGE    |- includes coordinates: L,M,T,x,y,z for each slice
%              .DENSE /
%              .focus
% 
% Based on MRI_register.m
% Created by Thien-Khoi N. Phung (November 13, 2017)
% Adding DENSE LA registration (Thien-Khoi N. Phung November 6, 2018)

%% STEP 1: Process Inputs
plot_flag = false; % flag for visualizing
DENSE_flag = false; DENSELA_flag = false;
for jz = 1:2:numel(varargin)
    switch varargin{jz}
        case 'LGESA'
            LGESA = varargin{jz+1};
        case 'LGEPS'
            LGEPS = varargin{jz+1};
        case 'DENSESA'
            DENSESA = varargin{jz+1};
            DENSE_flag = true;
        case 'DENSELA'
            DENSELA = varargin{jz+1};
            DENSELA_flag = true;
        case 'EPI'
            EPI = varargin{jz+1};
        case 'plot'
            plot_flag = varargin{jz+1};
        otherwise
            error('ERROR: Specified input is not valid.')
    end
end

%% STEP 2: Process Short Axis MRI files & register Base-Apex distances
% LGE
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

LGEstart = LGElocs(1,:)';

% Shifted location from 0 and normalized
LGEnorm = LKS'.*LGEthick;
LGEnorm = LGEnorm - min(LGEnorm); normfactor = max(LGEnorm);
LGEnorm = LGEnorm./normfactor;

% DENSE
if DENSE_flag
load(DENSESA);

% Add KeptSlices variable to structure (this variable indicates the number
% of slices (1:setstruct(jz).ZSize) per imaging view (jz))
for jz = 1:length(setstruct)
    setstruct(jz).KeptSlices = [1:setstruct(jz).ZSize];
end

% Changed to allow for unequal slice gaps (set using multiple stacks of SA
% loaded into Segment)
DENSElocslog = [];
for wh = 1:numel(setstruct)
    DKS = setstruct(wh).KeptSlices';

    % mm coordinates for most basal image
    DENSEpos = setstruct(wh).ImagePosition(:);

    % thickness between slices
    DENSEthick = setstruct(wh).SliceThickness + setstruct(wh).SliceGap; % mm

    % BA axis
    DENSEBAaxis = -cross(setstruct(wh).ImageOrientation(1:3),setstruct(wh).ImageOrientation(4:6));

    % Location of slices
    DENSElocs = ((repmat(DKS,1,3) - 1).*DENSEthick) .* repmat(DENSEBAaxis,numel(DKS),1);
    DENSElocs = repmat(DENSEpos',numel(DKS),1) + DENSElocs;
    DENSElocslog = [DENSElocslog; DENSElocs];
end
DENSElocs = DENSElocslog;

if ~DENSELA_flag
    % % Project point onto LGE BA axis
    % LGE2DENSE = DENSElocs(1,:)'-LGEstart;
    % DENSEpt = dot(LGE2DENSE,LGEBAaxis); % LGEBAaxis mag is 1
    % 
    % % normalized based on CINE
    % DENSEnorm = (DKS-min(DKS))'.*DENSEthick + DENSEpt;
    % DENSEnorm = DENSEnorm./normfactor;
    
    LGE2DENSE = DENSElocs - repmat(LGEstart',size(DENSElocs,1),1);
    DENSEpts = dot(LGE2DENSE',repmat(LGEBAaxis,size(DENSElocs,1),1)'); % LGEBAaxis mag is 1
    DENSEnorm = DENSEpts./normfactor;

else
    % Use DENSE LA points to determine DENSEnorm
    % First determine the relative scaling of the LGE points
    % Load ProcessedStack file
    load(LGEPS);
    
    % Access Base and Apex points
    BP = ProcessedStack.BasePoint;
    AP = ProcessedStack.ApexPoint;
    % Distance from LGE base to apex
    LGEbadist = norm(BP-AP);
    
    % Slice locations
    LGEsl = ProcessedStack.SliceCenters;
    
    % Determine base to apex normalized location
    normloc = zeros(size(LGEsl,1),1);
    for wh = 1:numel(normloc)
        normloc(wh) = dot((LGEsl(wh,:)-BP),(AP-BP)./LGEbadist);
    end
    normloc = sort(normloc);
    
    % Normalize by base-apex distance
    LGEnormba = normloc./LGEbadist;
    
    
    % Load DENSE LA pinpoints (apex-base)
    load(DENSELA)
    [ApexBasePts,~] = MRI_TransformLAWithPinPoints(setstruct);
    DENSE_AP = ApexBasePts(1,:);
    DENSE_BP = ApexBasePts(2,:);
    
    % B-A axis
    basetoapex = DENSE_AP-DENSE_BP;
    basetoapexdist = norm(basetoapex);
    basetoapexnorm = basetoapex./basetoapexdist;
    
    % Base to every point vectors
    basetopoints = DENSElocs - repmat(DENSE_BP,size(DENSElocs,1),1);
    
    % Determine dist of slice center to base poitn
    DENSEsl = zeros(size(DENSElocs,1),1);
    for wh = 1:numel(DENSEsl)
        DENSEsl(wh) = dot(basetopoints(wh,:),basetoapexnorm);
    end
    
    % Normalize by base-apex distance
    DENSEnormba = DENSEsl./basetoapexdist;
        
    % Renormalize with respect to LGE
    DENSEnorm = DENSEnormba - min(LGEnormba);
    DENSEnorm = DENSEnorm./(max(LGEnormba)-min(LGEnormba));
end
end

%% STEP 3: Load in LGE registration & surface fit
load(LGEPS);

landmarks = true;
[H,xyz] = MRI_SegmentRender(ProcessedStack,landmarks);
MRI_SurfaceRender(EPI,H);

if ~plot_flag
    close(H);
end

% Focus from the EPI.ipnode file
[EPIfit,focus]=MRI_IPNodeRead(EPI);

% Process xyz from MRI_SegmentRender into a grid


% Convert segmentation Cartesian->Prolate
% xyz matrix: going down rows is longitudinal
%             going left to right is circumferential
%             each page is X,Y,Z
% Cycle through each slice and create a matched-theta set of points
% L = zeros(numel(xyz.x),49,1);
% M = L;
% T = L;
L = [];
M = [];
T = [];
thetas = linspace(0,2*pi,50); thetas = thetas(1:end-1);
for jz = 1:numel(xyz.x)
    xseg = xyz.x{jz};
    yseg = xyz.y{jz};
    zseg = xyz.z{jz};
    
    if ~isnan(xseg)
    [lraw,mraw,traw] = LV_C2P(xseg,yseg,zseg,focus);
    
    % Sort by theta
    LMT = sortrows([lraw mraw traw],3);
    
    % Pad in theta (circumferentially)
    nts = size(LMT,1);
    pLMT = [LMT((nts-2):nts,1:2) LMT((nts-2):end,3)-2*pi;
            LMT;
            LMT(1:3,1:2) LMT(1:3,3)+2*pi];
    
    % Interpolate thetas
    Li = interp1(pLMT(:,3),pLMT(:,1),thetas,'linear');
    Mi = interp1(pLMT(:,3),pLMT(:,2),thetas,'linear');
    Ti = thetas;
    L = [L; Li];
    M = [M; Mi];
    T = [T; Ti];
    end
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

%% STEP 5: Fitting LGE and DENSE locations data
% x y z are the cart coord for the LGE segs
% Dx Dy Dz are cart for DENSE

if DENSE_flag
    Dx = zeros(numel(DENSEnorm),size(x,2));
    Dy = Dx;
    Dz = Dx; Dl = Dx; Dm = Dx; Dt = Dx;
end

for jz = 1:size(x,2)
    % x scaling using norms
    % x apex is max(x(:,1))
    % x basal is min(x(:,1))
    xscale = max(x(:,jz))-min(x(:,jz));
    LGEx  = LGEnorm.*xscale + min(x(:,jz));
    
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
end

% OUTPUT structure
HEART.focus = focus;
HEART.LGE.L = L;
HEART.LGE.M = M;
HEART.LGE.T = T;
HEART.LGE.x = x;
HEART.LGE.y = y;
HEART.LGE.z = z;
if DENSE_flag
HEART.DENSE.L = Dl;
HEART.DENSE.M = Dm;
HEART.DENSE.T = Dt;
HEART.DENSE.x = Dx;
HEART.DENSE.y = Dy;
HEART.DENSE.z = Dz;
end


if plot_flag
    figure(H), hold on
    plot3(x,y,z,'r<')
    if DENSE_flag
    plot3(Dx,Dy,Dz,'bs')
    end
end






