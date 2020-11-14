 function [GEOM] = RV_geom(LVGEOM, LVNODES, RVESV, RVWTH, RVCE, RVRE)
%RV_geom: Creates cartesian coordinates for an RV based on a given LV (created
%from LV_geom) and a known RV cavity volume and wall thickness
% GEOM = RV_geom(LVGEOM, RVESV, RVWTH, RVCE, RVRE)
%   INPUT VARIABLE:
%       LVGEOM: structure that includes:
%               .d: focal length of LV
%               .lambda(end): LV epicardial lambda value
%       RVESV: RV end systolic volume (mL)
%       RVWTH: RV wall thickness (cm)
%       RVCE: number of circumferential elements
%       RVRE: number of radial elements
%       LVLCR: element dimensions of LV (Long, Circ, Rad)
%   CALCULATIONS:
%       1. Calculate alpha, scaling for distance between LV and RV SA centers
%       2. Calculate coordinates for RV
%       3. Calculate RV insertion nodes
%   OUTPUT VARIABLE:
%       GEOM: structure with cartesian coordinates for RV and other RV params
%
% Thien-Khoi N. Phung (September 15, 2016)

% Unpack LV element dimensions (for ease of reading code)
[Ce, Le, Re] = size(LVGEOM.x); % number of points [circ. long. rad.]
Le = Le - 2; Re = Re - 1; % correct to represent number of elements

% MAKE SURE LV CIRCUMFERENTIAL ELEMENTS IS MULTIPLE OF 12
if rem(Ce,12)~=0
    warning('LV Circumferential Element dimension not divisible by 12')
end

% Convert INPUT to units of mm
RVESV = RVESV*(1000); % mL to mm^3
RVWTH = RVWTH*10; % cm to mm

% STEP 1: Calculate alpha (RV centroid scaling factor)
alpharange = linspace(0,1.5,100); % test range for alpha
RVbloodorange = zeros(1,100);

% for each alpha value, calculate RV cavity volume
for by = 1:length(alpharange)
    alpha = alpharange(by);

    % Calculating RV cavity volume constants
    RVa = atan((sqrt(3)/2) / (1/2 - alpha));
    RVa = RVa + pi*(RVa<0);
    RV_lumpedconst = RVa*(3/4 + (1/2 - alpha)^2) - (pi/3) + (sqrt(3)/2)*alpha;

    % RV Blood Area vs x height of the model
    xmin = LVGEOM.d*cosh(LVGEOM.lambda(end))*cos(2*pi/3);
    xmax = LVGEOM.d*cosh(LVGEOM.lambda(end))*cos(0);
    x = linspace(xmin, xmax, 100);
    RVarea = (LVGEOM.d*sinh(LVGEOM.lambda(end)))^2 .* RV_lumpedconst .*...
             (sin(acos(x./(LVGEOM.d.*cosh(LVGEOM.lambda(end)))))).^2;

    % Integrate to get RV cavity volume
    RVbloodorange(by) = trapz(x,RVarea);
end

% Interpolate for alpha value
GEOM.alpha = interp1(RVbloodorange, alpharange, RVESV);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  VECTORIZE THIS STEP 2 AT SOME POINT  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 2: Calculate RV coordinates

% Create coordinates grid
% .x (circumferential, longitudinal, radial directions)
GEOM.x = zeros(RVCE+1, Le+1, RVRE+1);
GEOM.y = zeros(RVCE+1, Le+1, RVRE+1);
GEOM.z = zeros(RVCE+1, Le+1, RVRE+1);

% Store LV node indicies for insertions
GEOM.plusixty = zeros(Le+1,2); % RVENDO RVEPI insertions
GEOM.minsixty = zeros(Le+1,2);

% Cycle through each base-apex layer of nodes
for jz = 1:Le+1
    % Find LV epi node indices
    LVepi = (1:Ce) + (Ce*(Le+1) + 1)*(Re) + (jz-1)*Ce;

    % Set RV x coordinate to be same as LV epi
    RVx = LVNODES(LVepi(1),1);

    % Calculate RV centroid based on alpha
    zerodeg = fix(Ce/4); % 0 degree index
    RVcent = [LVNODES(LVepi(zerodeg),2) LVNODES(LVepi(zerodeg),3)].*GEOM.alpha;

    % Calculate RV radius
    sixtydeg = fix(Ce/12); % 60 degree index
    RVendorad = norm([LVNODES(LVepi(sixtydeg),2) LVNODES(LVepi(sixtydeg),3)]-RVcent);

    % RV insertion angle from RV centroid
    pRVa = atan((LVNODES(LVepi(sixtydeg),3)-RVcent(2)) /...
                (LVNODES(LVepi(sixtydeg),2)-RVcent(1)));
    pRVa = pRVa + pi*(pRVa<0);
    RVangs = linspace(-pRVa, pRVa, RVCE+1)';

    % Calculate nodes at endo and epi
    RVradial = repmat(linspace(RVendorad,RVendorad+RVWTH,RVRE+1),RVCE+1,1);
    RVy = RVradial.*repmat(cos(RVangs(1:end)),1,RVRE+1) + RVcent(1);
    RVz = RVradial.*repmat(sin(RVangs(1:end)),1,RVRE+1) + RVcent(2);

    % Store wall coordinates
    GEOM.x(:,jz,:) = RVx;
    GEOM.y(:,jz,:) = RVy;
    GEOM.z(:,jz,:) = RVz;

    % Store insertion indices
    GEOM.plusixty(jz,:) = [LVepi(sixtydeg) LVepi(sixtydeg-1)];
    GEOM.minsixty(jz,:) = [LVepi(fix(Ce*5/12)) LVepi(fix(Ce*5/12)+1)];
end

% STEP 3: Calculate RV insertion nodes
% Replace +60 and -60 degree nodes to attach RV to LV
% Anterior side
dxyz  = (LVNODES(GEOM.minsixty(:,2),:) - LVNODES(GEOM.minsixty(:,1),:))./RVRE;
AB    = repmat(dxyz,1,1,RVRE+1);
AB(:,:,1) = LVNODES(GEOM.minsixty(:,1),:);
AB = cumsum(AB,3);

GEOM.x(1,:,:) = AB(:,1,:);
GEOM.y(1,:,:) = AB(:,2,:);
GEOM.z(1,:,:) = AB(:,3,:);

% Posterior side
dxyz  = (LVNODES(GEOM.plusixty(:,2),:) - LVNODES(GEOM.plusixty(:,1),:))./RVRE;
AB    = repmat(dxyz,1,1,RVRE+1);
AB(:,:,1) = LVNODES(GEOM.plusixty(:,1),:);
AB = cumsum(AB,3);

GEOM.x(end,:,:) = AB(:,1,:);
GEOM.y(end,:,:) = AB(:,2,:);
GEOM.z(end,:,:) = AB(:,3,:);


