function GEOM = MRI_FEMeshRender(ENDO,EPI,LCR)
%MRI_FEMeshRender(ENDO,EPI,LCR): Interpolates a specified element density
%mesh based on the Endo and Epi surface fits from MRI_FitContours.m
% GEOM = MRI_FEMeshRender(ENDO,EPI,LCR)
%   INPUTS:
%       ENDO- endo surface .IPNODE file path
%       EPI-  epi surface .IPNODE file path
%       LCR-  FE density in long, circ, and radial directions
%   CALCULATIONS:
%       STEP 1: Load in surfaces & pre-process
%       STEP 2: Set up surface mesh for interpolation
%       STEP 3: Interpolate along new surface grid
%       STEP 4: Convert ENDO and EPI interp surfaces to cartesian
%       STEP 5: Interpolate nodes between ENDO and EPI  
%       STEP 6: Convert cartesian solution to prolate spheroidal coordinates
%       STEP 7: Rotate all matrices to follow vary in (theta,mu,lambda) for the
%               for the Geometry_Create library of code
%   OUTPUT:
%       GEOM- structure with cartesian coordinates (.x .y .z) and prolate 
%       sph. coord. (.lambda .mu .theta .d) for LV endo and epi (mm)
% 
% Created by Thien-Khoi N. Phung (February 22, 2017)

% STEP 1: Load in surfaces & pre-process
    % degrees to radians conversion
    rads=pi/180;

    % load in files
    [ENDOfit,focus] = MRI_IPNodeRead(ENDO);
    [EPIfit,~]      = MRI_IPNodeRead(EPI);
    
    % store focus in OUPUT
    GEOM.d = focus;

    % sort based on (lambda, mu, theta) & convert to radians
    ENDOfit = sortrows(ENDOfit,[3 2 1]);
    ENDOfit(:,2:3) = ENDOfit(:,2:3).*rads;
    EPIfit = sortrows(EPIfit,[3 2 1]);
    EPIfit(:,2:3) = EPIfit(:,2:3).*rads;
    
    % determine size of bicubic mesh
    [m,n] = size(ENDOfit);
    munum = length(find(ENDOfit(:,3) == 0)); % # mu values
    thnum = m/munum;                         % # theta values
    elenum = thnum*(munum - 1);              % # elements (2D patches)
    
    % pad nodal mesh
    ENDOfit(m+1:m+munum,:) = ENDOfit(1:munum,:);
    ENDOfit(m+1:m+munum,3) = 2*pi;
    EPIfit(m+1:m+munum,:) = EPIfit(1:munum,:);
    EPIfit(m+1:m+munum,3) = 2*pi;
    thnum = thnum+1;
    
    % reshape nodalmesh for INTERP function in STEP 3
    ENDOdata = reshape(ENDOfit,munum,thnum,n);
    EPIdata  = reshape(EPIfit,munum,thnum,n);

% STEP 2: Set up surface mesh for interpolation
    % pull unique MU values (do once, same for ENDO and EPI)
    MUvec = ENDOfit(1:munum,2); % vector of unique MUs
    
    % pull unique THETA values (do once, same for ENDO and EPI)
    sortbyMU = sortrows(ENDOfit,[2 3]); % sort data by MU values
    THvec = sortbyMU(1:thnum,3); % vector of unique THETAs
    
    % meshgrid of unique MU and THETA for OG nodes
    [MUgridOG,THgridOG] = meshgrid(MUvec,THvec);
    
    % new mesh to interpolate- based on FE mesh density
    MUvecFE = linspace(0,MUgridOG(end),LCR(1)+2); % # long elements + 2
                                                  % LCR(2) does not include
                                                  % apex cap node
    THvecFE = linspace(0,THgridOG(end),LCR(2)+1); % # circ elements
    THvecFE = THvecFE(2:end);
    [MUgridFE,THgridFE] = meshgrid(MUvecFE,THvecFE);
    
    % set scaling derivative
    scale_der = 0; % not sure why this is 0 or what this does.. (TNP)
    
% STEP 3: Interpolate along new surface grid
    ENDOInterpSurf = MRI_BiCubicInterp(MUgridFE,THgridFE,ENDOdata,scale_der);
    EPIInterpSurf  = MRI_BiCubicInterp(MUgridFE,THgridFE,EPIdata,scale_der);
    % Output has dimension-3: lambda, mu, theta, & derivatives (look at
    % NodeData output written in MRI_BiCubicInterp)

% STEP 4: Convert ENDO and EPI interp surfaces to cartesian
    lambda = ENDOInterpSurf(:,:,1);
    mu     = ENDOInterpSurf(:,:,2);
    theta  = ENDOInterpSurf(:,:,3);
    [ENDOx,ENDOy,ENDOz] = LV_P2C(lambda,mu,theta,focus);
    
    lambda = EPIInterpSurf(:,:,1);
    mu     = EPIInterpSurf(:,:,2);
    theta  = EPIInterpSurf(:,:,3);
    [EPIx,EPIy,EPIz] = LV_P2C(lambda,mu,theta,focus);
    
% STEP 5: Interpolate nodes between ENDO and EPI   
    % step size between ENDO and EPI for equally spaced radial nodes
    stepx = (EPIx - ENDOx)./LCR(3);
    stepy = (EPIy - ENDOy)./LCR(3);
    stepz = (EPIz - ENDOz)./LCR(3);
    
    % create stepnumber matrix
    steps = zeros(size(ENDOx,1),size(ENDOx,2),LCR(3)+1);
    for jz = 1:LCR(3)+1
        steps(:,:,jz) = jz - 1;
    end
    
    % multiply steps by stepnumber & store OUTPUT cartesian nodes
    GEOM.x = repmat(ENDOx,1,1,LCR(3)+1) + repmat(stepx,1,1,LCR(3)+1).*steps;
    GEOM.y = repmat(ENDOy,1,1,LCR(3)+1) + repmat(stepy,1,1,LCR(3)+1).*steps;
    GEOM.z = repmat(ENDOz,1,1,LCR(3)+1) + repmat(stepz,1,1,LCR(3)+1).*steps;
    
% STEP 6: Convert cartesian solution to prolate spheroidal coordinates
    [L,M,T] = LV_C2P(GEOM.x,GEOM.y,GEOM.z,focus);
    GEOM.lambda = L;
    GEOM.mu = M;
    GEOM.theta = T;

% STEP 7: Rotate all matrices to follow vary in (theta,mu,lambda) for the
% for the Geometry_Create library of code
    GEOM.x = rot90(GEOM.x,1);
    GEOM.y = rot90(GEOM.y,1);
    GEOM.z = rot90(GEOM.z,1);
    GEOM.lambda = rot90(GEOM.lambda,1);
    GEOM.mu     = rot90(GEOM.mu,1);
    GEOM.theta  = rot90(GEOM.theta,1);
    