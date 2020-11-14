function [H,varargout] = MRI_SegmentRenderScar(ProcessedStack,varargin)
%MRI_SegmentRenderScar: takes processed MRI segmentation from
%(MRI_ProcessStack) to generate a plot of Scar Segmentation
%   INPUT: ProcessedStack- from MRI_ProcessStack
%   OUTPUT: H- figure handle
%           varargout- rotated cartesian coordinates of scar
% 
% Thien-Khoi N. Phung (April 13, 2017)

%% Figure window
if isempty(varargin)
    H = figure('WindowStyle','docked');
    hold on
else
    figure(varargin{1});
    hold on;
end

%% Process SCAR countours using same coordinates as MRI_FitContours
% Convert scar to cardiac coordinates
% Assume ProcessedStack.Scar is 1 by n
scar = [];
for jz = 1:length(ProcessedStack.Scar)
    if ~isempty(ProcessedStack.Scar{jz})
        scar = [scar; ProcessedStack.Scar{jz}];
    end
end

A = ProcessedStack.ApexPoint;
B = ProcessedStack.BasePoint;
Si = ProcessedStack.SeptumPoint;

    % First basis vector, subtract Base from Apex and divide by magnitude
    C = A-B;
    e1 = C ./ norm(C);

    % Calculate Origin location
    origin = B + (1/3)*C;

    Nber_points = length(scar); % Assume confocal endo & epi

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
    DataScar = (scar - repmat(origin,Nber_points,1))*Transform';
    varargout{1} = DataScar;
    
%% Plot scar
x = DataScar(:,1);
y = DataScar(:,2);
z = DataScar(:,3);
plot3(x,y,z,'.','color',[136,86,167]./255,'MarkerSize',12);


