% EP_slab.m
% Toy model of electrical propagation through a 3D Slab
% layers of a transmural section of "myocardium"
%
% Thien-Khoi N. Phung (December 17, 2020)

%% Creating slab mesh
% Cube element length
res = 1; % mm

% Element density
Le = 10;
Ce = 10;
Re = 3;

% Make element nodes
Ln = (0:Le)*res;
Cn = (0:Ce)*res;
Rn = (0:Re)*res;

[x,y,z] = meshgrid(Ln, Cn, Rn);

NODES = [x(:) y(:) z(:)];

% Plot nodes to visualize numbering
figure('WindowStyle','docked','NumberTitle','off','name','Node Numbering')
hold on
for jz = 1:length(NODES)
    plot3(NODES(jz,1), NODES(jz,2), NODES(jz,3), 'k.')
    text(NODES(jz,1), NODES(jz,2), NODES(jz,3), num2str(jz))
end
axis equal
xlabel('Longitudinal'),ylabel('Circumferential'),zlabel('Radial')
view(3)

% Element connectivity
HEX = zeros(Le*Ce*Re, 8);
dad = 0; 
dadcheckLC = Ce:Ce:(Ce*Le*Re); % Check for next row
dadcheckR = (Ce*Le):(Ce*Le):(Ce*Le*Re); % Check for next layer
for jz = 1:(Le*Ce*Re)
    
    sd = jz + dad;
    HEX(jz,1:4) = [sd, sd+(Ce+1), (sd+1)+(Ce+1), sd+1];
    HEX(jz,5:8) = [sd, sd+(Ce+1), (sd+1)+(Ce+1), sd+1]+((Ce+1)*(Le+1));
    
    % Check if next element is next row or layer
    dad = dad + sum(dadcheckLC==jz) + (Ce+1)*sum(dadcheckR==jz);
    
end

% HEX centers & volumes
HEXc = LV_hexcent(NODES,HEX);
HEXv = LV_hexvol(NODES,HEX);

% Visualize Element Mesh
FE_RenderPatchMesh(NODES,HEX,'elements',1:(Le*Ce*Re),...
                        'alpha',0,'edgealpha',0.25,'title','Element Mesh');
axis on
xlabel('Longitudinal'),ylabel('Circumferential'),zlabel('Radial')
view([-45 30])

% Visualize Element Numbering
FE_RenderPatchMesh(NODES,HEX,'elements',1:(Le*Ce*Re),...
                        'alpha',0,'edgealpha',0.25,'title','Element Numbering');
hold on
for jz = 1:(Le*Ce*Re)
    text(HEXc(jz,1),HEXc(jz,2),HEXc(jz,3), num2str(jz))
end
axis on
xlabel('Longitudinal'),ylabel('Circumferential'),zlabel('Radial')
view([-45 30])


% Store model information
MODEL.LVLCR = [Le Ce Re]; % Mesh dimensions
MODEL.HEX   = HEX;   % Element connectivity matrix
MODEL.NODES = NODES; % Node coordinates
MODEL.HEXc  = HEXc;  % Element Centers
MODEL.HEXv  = HEXv;  % Element volumes
MODEL.EPnodeIDX = 1:(Le*Ce*Re); % Index of elements used in electrical model (all of them)


%% Assign FIBER
% lyrfib = linspace(90,90,Re).*pi/180; % all Longitudinal fibers
% lyrfib = linspace(0,0,Re).*pi/180;   % all Circumferential Fibers
lyrfib = linspace(60,-60,Re).*pi/180; % Endo to Epi gradient

% Assign by radial layer
Qfcr = zeros(3,3,Le*Ce*Re);
for jz = 1:Re
    FF = [sin(lyrfib(jz)); cos(lyrfib(jz)); 0];
    RR = [0; 0; 1];
    CC = cross(RR,FF);
    
    Qfcr(:,:,(1:Le*Ce)+(Le*Ce*(jz-1))) = repmat([FF CC RR],1,1,Le*Ce);
end

% Visualize Fibers in one transmural stack of elements
stack = 1:(Le*Ce):(Le*Ce*Re);
FE_RenderPatchMesh(NODES,HEX,'elements',stack,...
                        'alpha',0,'edgealpha',0.25,'title','Transmural Fibers');
hold on
quiver3(HEXc(stack,1),HEXc(stack,2),HEXc(stack,3),...
    squeeze(Qfcr(1,1,stack)),squeeze(Qfcr(2,1,stack)),squeeze(Qfcr(3,1,stack)),...
    0.25,'LineWidth',3,'ShowArrowHead','off');
axis on
xlabel('Longitudinal'),ylabel('Circumferential'),zlabel('Radial')
view([-45 30])

% Store model information
MODEL.Qfcr  = Qfcr;  % Store Rotation Matrix (for electrical model)


%% Indexing neighbors
% This section defines which elements are touching/connected
% Creating matrix with spatial elements laid out
el = reshape(1:(Le*Ce*Re),Le,Ce,Re);

% Pad with NaN on all borders
elpad = nan(Le+2, Ce+2, Re+ 2);
elpad(2:end-1, 2:end-1, 2:end-1) = el;

% Logical matrix (for searching neighbors)
L = ~isnan(elpad);

	nbr = zeros((Ce*Le*Re), 26);
    % Shift logical matrix across the neighbors                        C   L   R
    nbr(:,1) = elpad(L(    [2:end 1],     [2:end 1],     [2:end 1])); % [+1, +1, +1]
    nbr(:,2) = elpad(L(            :,     [2:end 1],     [2:end 1])); % [ 0, +1, +1]
    nbr(:,3) = elpad(L([end 1:end-1],     [2:end 1],     [2:end 1])); % [-1, +1, +1]

    nbr(:,4) = elpad(L(    [2:end 1],             :,     [2:end 1])); % [+1,  0, +1]
    nbr(:,5) = elpad(L(            :,             :,     [2:end 1])); % [ 0,  0, +1]
    nbr(:,6) = elpad(L([end 1:end-1],             :,     [2:end 1])); % [-1,  0, +1]

    nbr(:,7) = elpad(L(    [2:end 1], [end 1:end-1],     [2:end 1])); % [+1, -1, +1]
    nbr(:,8) = elpad(L(            :, [end 1:end-1],     [2:end 1])); % [ 0, -1, +1]
    nbr(:,9) = elpad(L([end 1:end-1], [end 1:end-1],     [2:end 1])); % [-1, -1, +1]

    nbr(:,10)= elpad(L(    [2:end 1],     [2:end 1],             :)); % [+1, +1,  0]
    nbr(:,11)= elpad(L(            :,     [2:end 1],             :)); % [ 0, +1,  0]
    nbr(:,12)= elpad(L([end 1:end-1],     [2:end 1],             :)); % [-1, +1,  0]

    nbr(:,13)= elpad(L(    [2:end 1],             :,             :)); % [+1,  0,  0]
    nbr(:,14)= elpad(L([end 1:end-1],             :,             :)); % [-1,  0,  0]

    nbr(:,15)= elpad(L(    [2:end 1], [end 1:end-1],             :)); % [+1, -1,  0]
    nbr(:,16)= elpad(L(            :, [end 1:end-1],             :)); % [ 0, -1,  0]
    nbr(:,17)= elpad(L([end 1:end-1], [end 1:end-1],             :)); % [-1, -1,  0]

    nbr(:,18)= elpad(L(    [2:end 1],     [2:end 1], [end 1:end-1])); % [+1, +1, -1]
    nbr(:,19)= elpad(L(            :,     [2:end 1], [end 1:end-1])); % [ 0, +1, -1]
    nbr(:,20)= elpad(L([end 1:end-1],     [2:end 1], [end 1:end-1])); % [-1, +1, -1]

    nbr(:,21)= elpad(L(    [2:end 1],             :, [end 1:end-1])); % [+1,  0, -1]
    nbr(:,22)= elpad(L(            :,             :, [end 1:end-1])); % [ 0,  0, -1]
    nbr(:,23)= elpad(L([end 1:end-1],             :, [end 1:end-1])); % [-1,  0, -1]

    nbr(:,24)= elpad(L(    [2:end 1], [end 1:end-1], [end 1:end-1])); % [+1, -1, -1]
    nbr(:,25)= elpad(L(            :, [end 1:end-1], [end 1:end-1])); % [ 0, -1, -1]
    nbr(:,26)= elpad(L([end 1:end-1], [end 1:end-1], [end 1:end-1])); % [-1, -1, -1]

    % Pair neighbors to house
    H = repmat(elpad(L(:)),26,1); % Houses
    N = nbr(:); % Neighbors 

    % Remove NaN neighbors
    nope = isnan(N);
    H(nope) = [];
    N(nope) = [];

% Unique Home Neighbor Pairs
HN = [H N];

% Store model information
MODEL.HN    = HN;    % Home-Neighbor pairs


%% Set Conduction Velocities
% Set muscle fiber conduction velocity
v = 1;

% Assign anisotropy of conduction velocity
LVMf = v;
LVMc = 0.4*v;

% Assign all elements to have the LV conduction velocities
VELO = repmat([LVMf LVMc LVMc],size(MODEL.HEXc,1),1);

% Add the conduction velocities (VELO) to the MODEL struct variable
MODEL.VELO = VELO;


%% Set which element(s) initiate the electrical propagation
% This matrix will be an n by 2 matrix that includes
%   elements that initiate propagation in column 1 
%   time at which the elements initate in column 2
STAR = [round(Ce/2) 0];

% Visualize the element selected for STAR
SS = FE_RenderPatchMesh(MODEL.NODES,MODEL.HEX,'elements',MODEL.EPnodeIDX,...
                        'alpha',0,'edgealpha',0.25,'title','Initiation Element');
FE_RenderPatchMesh(MODEL.NODES,MODEL.HEX,'elements',STAR(:,1),...
                                         'facecolor',[1 0 0],...
                                         'handle',SS);
xlabel('Longitudinal'),ylabel('Circumferential'),zlabel('Radial')
view([-45 30])
axis on


%% Electrical Simulation
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
                   'data',MODEL.litTIME,'title','Electrical Activation');
colormap(flipud(plasma()))
xlabel('Longitudinal'),ylabel('Circumferential'),zlabel('Radial')
view([-45 30])
axis on
