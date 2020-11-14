% EModel_RunModel.m
% This script sets up an electrical model and runs multiple simulations to
% compare ECG to data.

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


%% Load Data for optimization
% Here we load both the Baseline and LBBB data, but you can just choose one
% to optimize against.

% Baseline data
load(['data\' sub '\DATA_EP_bl.mat'])
% Interpolate and save BL data into model
MODEL.blt = linspace(DATA_EP.QRSthyme(1),DATA_EP.QRSthyme(end),100); % ECG time
MODEL.bl  = interp1(DATA_EP.QRSthyme,DATA_EP.QRS,MODEL.blt); % ECG signals

% LBBB data
load(['data\' sub '\DATA_EP_lbbb.mat'])
% Interpolate and save BL data into model
MODEL.lbbbt = linspace(DATA_EP.QRSthyme(1),DATA_EP.QRSthyme(end),100);
MODEL.lbbb  = interp1(DATA_EP.QRSthyme,DATA_EP.QRS,MODEL.blt);

% Clean up the workspace
clear DATA_EP


%% Find all possible electrical initation sites to simulate
% Star solution space
stars = [MODEL.elem.LVP(:); MODEL.elem.RVP(:)];
% Remove elements not in the electrical model (ie scar)
stars = stars(ismember(stars,MODEL.EPnodeIDX));
% stars includes all purkinje elements in the model


%% STEP 1. Run all of the possible models
% Simulate activation from every endocardial (Purkinje) element
% Calculate the 12 Lead ECG
% Quantify error against BL & LBBB data

% Pre-allocate (cause proper coding) error matrix (24 columns: 12 for BL,
% 12 for LBBB lead correlation coefficients)
n = size(stars,1);
solnspace = zeros(n,24);

% This is a brute force optimization- we run all of the possible models
% Depending on "n" this may take some time to run
% You can use a parfor loop to make it run "faster"
for jz = 1:n
    % Simulate model
    %   Assume velocity = 1 (will tune in next step)
    %   flag_ECG = true % calculate the leads
    [err,~,~] = runmodel(1,MODEL,[stars(jz) 0],true);
    solnspace(jz,:) = err; % Store the error
end


%% Step 2. Find best solution and tune velocity
% Find models with the highest correlation coefficient of all 12 leads
% For the baseline data
[~,bestBL] = max(mean(solnspace(:,1:12),2)); 
% For the LBBB data
[~,bestLBBB] = max(mean(solnspace(:,13:24),2));

% best____ is the index for the "star" that performed the best in matching
% the ECG shapes in the data

% Optimize conduction velocity to match the QRS
% Here we re-run the "best" models using many values for conduction
% velocity "v", and store the latest activation time "lat___"
vscale = linspace(0.4,1.6,40); % These is the velo range we will test

% Preallocate matrix for latest activation times
latbl = zeros(numel(vscale),1); % baseline simulations
latlbbb = latbl;                % lbbb simulations

% Cycle through each velocity value
for jz = 1:numel(vscale)
    % Simulate baseline model & store latest activation time 
    [~,LT] = runmodel(vscale(jz),MODEL,[stars(bestBL) 0],false);
    latbl(jz) = max(LT);
    
    % Simulate lbbb model & store latest activation time
    [~,LT] = runmodel(vscale(jz),MODEL,[stars(bestLBBB) 0],false);
    latlbbb(jz) = max(LT);
end

% Based on the data ECG, we will now interpolate the best velocity value
% that matches the latest activation time to the QRS.
% "ov___" is the velocity
ov_bl   = interp1(latbl,vscale,max(MODEL.blt),'linear','extrap');
ov_lbbb = interp1(latlbbb,vscale,max(MODEL.lbbbt),'linear','extrap');


%% STEP 3. Finally we simulate the best fit models
% Best Baseline Model
[~,litTIME_bl,PV_bl] = runmodel(ov_bl,MODEL,[stars(bestBL) 0],true);

FE_RenderPatchMesh(MODEL.NODES,MODEL.HEX,...
                   'elements',MODEL.EPnodeIDX,...
                   'data',litTIME_bl);
colormap(flipud(plasma()))

figure('WindowStyle','docked','numbertitle','off','Name','Baseline ECG')
for jz = 1:12
    subplot(3,4,jz)
    plot(MODEL.bl(:,jz),':','LineWidth',3)
    hold on
    plot(PV_bl(:,jz),'LineWidth',3)
end
legend('Data','Model')

% Best LBBB Model
[~,litTIME_lbbb,PV_lbbb] = runmodel(ov_lbbb,MODEL,[stars(bestLBBB) 0],true);

FE_RenderPatchMesh(MODEL.NODES,MODEL.HEX,...
                   'elements',MODEL.EPnodeIDX,...
                   'data',litTIME_lbbb);
colormap(flipud(plasma()))

figure('WindowStyle','docked','numbertitle','off','Name','LBBB ECG')
for jz = 1:12
    subplot(3,4,jz)
    plot(MODEL.lbbb(:,jz),':','LineWidth',3)
    hold on
    plot(PV_lbbb(:,jz),'LineWidth',3)
end
legend('Data','Model')


%% Run Model- output error metric
% We write a function here called "runmodel" to make our code more
% efficient. The code in this function is the same as in EModel_RunModel.m,
% but it also includes an error term to that determines if the simulation
% matches the data.
function [err,varargout] = runmodel(v,MODEL,STAR,flag_ECG)
% Set velocities (Lee+ 2019)
LVMf = v;
LVMc = 0.4 * v;
LVPf = 6.0 * v;
LVPc = 0.4 * v;
RVMf = v;
RVMc = 0.4 * v;
RVPf = 6.0 * v;
RVPc = 0.4 * v;

Qfcr = MODEL.Qfcr;

% Assign everything the RV properties
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

MODEL.VELO = VELO;

HEXc = MODEL.HEXc;
HN   = MODEL.HN;

% Run model
[litTIME,~] = ElecpropSPT(HEXc,HN,Qfcr,VELO,STAR);
MODEL.litTIME = litTIME(MODEL.EPnodeIDX);

% Export (optional) litTIME
varargout{1} = MODEL.litTIME;

if flag_ECG
    % Calculate pseudo-Voltage at leads
    pV = pseudoECG(MODEL);

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

    varargout{2} = pV12;

    % To calculate how the model ECG compares to data, we calculate the
    % correlation coefficient for each of the 12 leads against both the
    % baseline and lbbb data.
    % Data 
    databl = MODEL.bl;
    datalbbb = MODEL.lbbb;
    err = zeros(1,24);
    for jz = 1:12
        err(1,jz) = corr(pV12(:,jz),databl(:,jz));
        err(1,jz+12) = corr(pV12(:,jz),datalbbb(:,jz));
    end
else
    err = [];
end % End of flag_ECG
end % End of runmodel()