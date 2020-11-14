function [GEOM] = LV_geom(VOL,WTH,SI,m_density,t_density,n)
%LV_geom: Creates cartesian coordinates for an LV with known cavity volume,
%wall thickness, and sphericity index.
% XYZ = LV_geom(VOL,WTH,SI,m_density,t_density,n)
%   INPUT VARIABLE: (note: all from ES or ED)
%       VOL: cavity volume (mL)
%       WTH: wall thickness (cm)
%       SI: sphericity index, ratio between Base-Apex length and width at
%       midpoint of B-A length
%       m_density: number of rings of elements (Base to Apex) excludes apex
%       cap
%       t_density: number of elements per ring
%       n: number of elements thru wall
%   CALCULATIONS:
%       1. Calculate lambda_endo, len, from SI.
%       2. Calculate focal length, d, from VOL.
%       3. Calculate lambda_epi, lep, from WTH.
%       4. Calculate cartesian coordinates for LV.
%       5. Create nodes thru wall.
%   OUTPUT VARIABLE:
%       GEOM: structure with cartesian coordinates and prolate sph. coord. 
%       for LV endo and epi (mm)
%   NOTE: X is the Base-Apex axis
%
% Thien-Khoi N. Phung (April 26, 2016)
% Edited: Output mm nodes (TNP, May 12, 2016)
% Edited: Step 2 formulation incorrect (TNP, Aug 31, 2016)

% Convert INPUT to units of mm
VOL = VOL*(1000); % mL to mm^3
WTH = WTH*(10); % cm to mm

% STEP 1: Calculate lambda_endo, len, from SI.
len = atanh(1.5/(2*sin(acos(0.25))*SI));

% STEP 2: Calculate focal length, d, from VOL
    % OLD INCORRECT FORMULATION forgot power of 3 on cosh (TNP 08-31-16)
    % d = ((3*VOL)/(2*pi)*...
    %     (((1-cos(120*pi/180))*cosh(len) - (1-cos(120*pi/180)^3)*cosh(len))...
    %     -((1-cos(120*pi/180))*cosh(0)   - (1-cos(120*pi/180)^3)*cosh(0)))^(-1))^(1/3);
d = ((3*VOL)/(2*pi)/...
    (((1-cos(120*pi/180))*cosh(len)^3 - (1-cos(120*pi/180)^3)*cosh(len))...
    -((1-cos(120*pi/180))*cosh(0)^3   - (1-cos(120*pi/180)^3)*cosh(0))))^(1/3);

% STEP 3: Calculate lambda_epi, lep, from WTH.
lep = asinh(WTH/(d*sin(pi/2)) + sinh(len));

% STEP 4: Calculate cartesian coordinates for LV endo and epi.
m = linspace(0,120,m_density+2).*pi/180; % mu angles 0 to 120 degrees
t = linspace(0,360,t_density+1).*pi/180; t = t(2:end); % theta angles

[mg, tg] = meshgrid(m,t); % meshgrid for mu and theta values

% prolate spheroidal to cartesian
z_endo = d.*sinh(len).*sin(mg).*cos(tg);
y_endo = d.*sinh(len).*sin(mg).*sin(tg);
x_endo = d.*cosh(len).*cos(mg);

z_epi  = d.*sinh(lep).*sin(mg).*cos(tg);
y_epi  = d.*sinh(lep).*sin(mg).*sin(tg);
x_epi  = d.*cosh(lep).*cos(mg);

l_endo = repmat(len, length(t), length(m));
l_epi  = repmat(lep, length(t), length(m));

% Storing Variables for output
GEOM.x(:,:,1) = x_endo;
GEOM.y(:,:,1) = y_endo;
GEOM.z(:,:,1) = z_endo;
GEOM.lambda(:,:,1) = l_endo;
GEOM.mu(:,:,1)     = mg;
GEOM.theta(:,:,1)  = tg;
GEOM.d = d;

% STEP 5: Creating nodes to make n elements thru wall
% Calculate the spacing through lambda required
dl = (l_epi - l_endo)./n;

% Add new layers numbered 2 thru n
for j = 2:n
    l_layer = l_endo + dl.*(j-1);

    GEOM.z(:,:,j) = d.*sinh(l_layer).*sin(mg).*cos(tg);

    GEOM.y(:,:,j) = d.*sinh(l_layer).*sin(mg).*sin(tg);

    GEOM.x(:,:,j) = d.*cosh(l_layer).*cos(mg);

    GEOM.lambda(:,:,j) = l_layer;
    GEOM.mu(:,:,j)     = mg;
    GEOM.theta(:,:,j)  = tg;
end

% Storing epi variables
GEOM.x(:,:,n+1) = x_epi;
GEOM.y(:,:,n+1) = y_epi;
GEOM.z(:,:,n+1) = z_epi;
GEOM.lambda(:,:,n+1) = l_epi;
GEOM.mu(:,:,n+1)     = mg;
GEOM.theta(:,:,n+1)  = tg;
