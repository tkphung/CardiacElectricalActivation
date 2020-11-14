function [H,varargout] = MRI_SegmentRender(ProcessedStack,landmarks,varargin)
%MRI_SegmentRender: takes processed MRI segmentation from
%(MRI_ProcessStack) to generate a plot of the contours
%   INPUT: ProcessedStack- from MRI_ProcessStack
%          varargin{1} can include figure handle to plot on top
%          landmarks- plot base-apex, septal inserts, mid septum
%   OUTPUT: H- figure handle
%           varargout- rotated X,Y,Z segmentation
% 
% Thien-Khoi N. Phung (February 16, 2017)
% TNP (June 28, 2017) added XYZ varargout

%% Figure window
if isempty(varargin)
    figure('WindowStyle','docked');
    hold on
else
    figure(varargin{1});
    hold on;
end

%% Process Data using same coordinates as MRI_FitContours
% Convert slices to cardiac coordinates
alldata_endo = ProcessedStack.alldata_endo;
alldata_epi = ProcessedStack.alldata_epi;
Endo = alldata_endo(:,1:3);
Epi = alldata_epi(:,1:3);
A = ProcessedStack.ApexPoint;
B = ProcessedStack.BasePoint;
Si = ProcessedStack.SeptumPoint;

    % First basis vector, subtract Base from Apex and divide by magnitude
    C = A-B;
    e1 = C ./ norm(C);

    % Calculate Origin location
    origin = B + (1/3)*C;

    Nber_points = length(Endo); % Assume confocal endo & epi

    % Second basis vector using NEW METHOD- plane intersects septal point &
    % e1
    D = Si(1,:) - origin;
    D2 = D - dot(D,e1)*e1;
    e2 = D2 ./ norm(D2);

    % Third basis vector
    E = cross(e1,e2);
    e3 = E ./ norm(E);

    % Transformation matrix
    Transform = [e1; e2; -e3];

    % Substract origin and transform (rotate) data
    DataEndo = (Endo - repmat(origin,Nber_points,1))*Transform';
    DataEpi = (Epi - repmat(origin,Nber_points,1))*Transform';
    A_transf = (A - origin)*Transform';
    B_transf = (B - origin)*Transform';   
    S_transf  = (Si(1,:) - origin)*Transform';% Septal midpoint
    S_transf1 = (Si(2,:) - origin)*Transform';% Septal intersection 1
    S_transf2 = (Si(3,:) - origin)*Transform';% Septal intersection 2
    % DataScar = [Data alldata(:,4)];
    
    % X Y Z coordinates (cartesian) and Scar Transmurality
    DataEndo = [DataEndo alldata_endo(:,4)];
    DataEpi = [DataEpi alldata_epi(:,4)];

%% Plot around the segmentations
bins = unique(alldata_endo(:,3));
for jz = 1:numel(bins)
    tracing = alldata_endo(:,3) == bins(jz);
    tracing = find(tracing);
    tracing = [tracing; tracing(1)];
    % Plot endocardium
    x = DataEndo(tracing,1);
    y = DataEndo(tracing,2);
    z = DataEndo(tracing,3);

    plot3(x,y,z,'.-','color',[253,141,60]./255,'LineWidth',2)

    % Plot epicardium
    x = DataEpi(tracing,1);
    y = DataEpi(tracing,2);
    z = DataEpi(tracing,3);

    rotdata.x{jz} = x(1:end-1);
    rotdata.y{jz} = y(1:end-1);
    rotdata.z{jz} = z(1:end-1);
    
    plot3(x,y,z,'.-','color',[37,52,148]./255,'LineWidth',2)
end

varargout{1} = rotdata;

%% Plot landmarks
if landmarks
    AB = [A_transf; B_transf];
    Si = [S_transf1;S_transf2];
    S  = S_transf;
    
    % plot apex-base
    plot3(AB(:,1),AB(:,2),AB(:,3),'k.-','LineWidth',3)
    
    % plot septal stuff
    plot3(Si(:,1),Si(:,2),Si(:,3),'bo','MarkerSize',8,'MarkerFaceColor','b')
    plot3(S(1),S(2),S(3),'ro','MarkerSize',8,'MarkerFaceColor','r')
end

%% Clean up figure and send back handle.
axis equal;
xlabel('X'); 
ylabel('Y'); 
zlabel('Z');

H=gcf;

