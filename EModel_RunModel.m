% EModel_RunModel.m
% This script runs a single electrical model simulation and produces the
% resulting ECG.

% ***add "lib" folder to your directory

% Thien-Khoi N. Phung (October 7, 2020)

%% Load in the model geometry
% Select the subject to model (this is the folder name in the data folder)
sub = 'CRT006';

% Load the model geometry (this is the finite element mesh for the
% electrical simulation- this was generated from a different script)
% The variable MODEL is a struct which stores information about the
% geometry.
load(['data\' sub '\MODEL.mat'])


%% Set the velocities for the electrical models
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


%% Set which element(s) initiate the electrical propagation
% This matrix will be an n by 2 matrix that includes
%   elements that initiate propagation in column 1 
%   time at which the elements initate in column 2

STAR = [1 0];

% Visualize the element selected for STAR
SS = FE_RenderPatchMesh(MODEL.NODES,MODEL.HEX,'elements',MODEL.EPnodeIDX,'alpha',0,'edgealpha',0.25);
FE_RenderPatchMesh(MODEL.NODES,MODEL.HEX,'elements',STAR(:,1),...
                                         'facecolor',[1 0 0],...
                                         'handle',SS);

%% Pull other variables from MODEL struct
Qfcr = MODEL.Qfcr; % Element rotation matrices to specify direction of fiber direction
HEXc = MODEL.HEXc; % Center point of each element
HN   = MODEL.HN;   % This specifies which elements are "touching" (aka Home-Neighbor pairs)


%% Electrical Simulation
% The function ElecpropSPT simulates the electrical propagation
% Returns litTIME- an electrical activation time for all elements
%         litHOO- element index for who activated the element
% both of these variables including for scar elements and those not in the model
[litTIME,litHOO] = ElecpropSPT(HEXc,HN,Qfcr,VELO,STAR);

% Keep only the electrical activation times for elements in the electrical
% model
MODEL.litTIME = litTIME(MODEL.EPnodeIDX);


%% Visualize the electrical activation
FE_RenderPatchMesh(MODEL.NODES,MODEL.HEX,...
                   'elements',MODEL.EPnodeIDX,...
                   'data',MODEL.litTIME);
colormap(flipud(plasma()))


%% Pseudo-ECG simulation
% Now we take in the electrical activation times and predict what the ECG
% signals will look like at specified lead locations

% Calculate pseudo-Voltage at leads
pV = pseudoECG(MODEL); % Note that the times are already stored in MODEL struct

% pV contains signals at the "standard" 12 lead locations on the torso
% We need to do some post-processing to generate 12-lead ECG 
    % Calculate 12 Lead ECG
    % Lead Order (V1-V6 LA RA LL RL)
    % 12 Lead (I II III aVR aVL aVF V1-V6)
        % Lead I: RA (-) to LA (+) (Right Left, or lateral)
        pV12(1,:) = pV(7,:) - pV(8,:);
        % Lead II: RA (-) to LL (+) (Superior Inferior)
        pV12(2,:) = pV(9,:) - pV(8,:);
        % Lead III: LA (-) to LL (+) (Superior Inferior)
        pV12(3,:) = pV(9,:) - pV(7,:);

        % Lead aVR: RA (+) to [LA & LL] (-) (Rightward)
        pV12(4,:) = pV(8,:);% - pV12(3,:);
        % Lead aVL: LA (+) to [RA & LL] (-) (Leftward)
        pV12(5,:) = pV(7,:);% - pV12(2,:);
        % Lead aVF: LL (+) to [RA & LA] (-) (Inferior)
        pV12(6,:) = pV(9,:);% - pV12(1,:);

        pV12(7:12,:) = pV(1:6,:); 

    % Downsample
    pV12 = interp1(linspace(0,1,200),pV12',linspace(0,1,100));
    % Normalize the magnitudes
    pV12 = pV12./(max(max(pV12)) -  min(min(pV12))).*1.15;

% Visualize the 12 lead ECG
time = linspace(min(MODEL.litTIME),max(MODEL.litTIME),size(pV12,1));
ECG_plot(time,pV12)
