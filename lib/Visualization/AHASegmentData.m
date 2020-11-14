function [AHAseg,varargout] = AHASegmentData(LOCS,DATA,focus,dataname,plotflag,varargin)
%AHASegmentData: Takes a data from an LV model and maps it onto a 16
%segment AHA visualization
% AHASeg = AHASegmentData(LOCS,DATA,focus,varargin)
%   INPUTS:
%       LOCS- cartesian points in cardiac coordinates
%       DATA- data associated with each point
%       focus- prolate spheroidal coord focus 
%       plotflag- produce plot? (true or false) 
%       varargin- axesflag
%   OUTPUT:
%       AHASeg
%       varargout- region for each element
%                  timing for each element
% 
% Created by Thien-Khoi N. Phung (June 21, 2019)
axes_flag      = false;
if ~isempty(varargin)
    for jz = 1:2:numel(varargin)
        switch varargin{jz}
            case 'axes' % plot data colormap on FE model
                axes_flag = true;
                axeshandle = varargin{jz+1};
            otherwise
            error('ERROR: Check your varargins.')
        end
    end
end
%% Convert Data from Cartesian to PSC
[L,M,T] = LV_C2P(LOCS(:,1),LOCS(:,2),LOCS(:,3),focus);

% NOTE: theta (T) starts at mid-septum 0 rad and follows SPLAS to 2pi rads
theta = T;
% Reverse to SALPS
theta = -theta + 2*pi;

rad = M.*(180/pi);

%% Sort data in segments - data is in SALPS with theta = 0 at septum
% Number of segments
nSegments = 16;

% Construct topology, detect which of the points are contained in each
% segment (see Pim's NB Growth>UVA dog study>Infarct Sizes
AHAseg(1).topo  = find((theta >= 1/3*pi) & (theta < 2/3*pi) & (rad >= 86));
AHAseg(2).topo  = find((theta >= 0)      & (theta < 1/3*pi) & (rad >= 86));
AHAseg(3).topo  = find((theta >= 5/3*pi) & (theta < 6/3*pi) & (rad >= 86));
AHAseg(4).topo  = find((theta >= 4/3*pi) & (theta < 5/3*pi) & (rad >= 86));
AHAseg(5).topo  = find((theta >= 3/3*pi) & (theta < 4/3*pi) & (rad >= 86));
AHAseg(6).topo  = find((theta >= 2/3*pi) & (theta < 3/3*pi) & (rad >= 86));

AHAseg(7).topo  = find((theta >= 1/3*pi) & (theta < 2/3*pi) & (rad >= 45) & (rad < 86));
AHAseg(8).topo  = find((theta >= 0)      & (theta < 1/3*pi) & (rad >= 45) & (rad < 86));
AHAseg(9).topo  = find((theta >= 5/3*pi) & (theta < 6/3*pi) & (rad >= 45) & (rad < 86));
AHAseg(10).topo  = find((theta >= 4/3*pi) & (theta < 5/3*pi) & (rad >= 45) & (rad < 86));
AHAseg(11).topo  = find((theta >= 3/3*pi) & (theta < 4/3*pi) & (rad >= 45) & (rad < 86));
AHAseg(12).topo  = find((theta >= 2/3*pi) & (theta < 3/3*pi) & (rad >= 45) & (rad < 86));

AHAseg(13).topo  = find((theta >= 1/4*pi) & (theta < 3/4*pi) & (rad < 45));
AHAseg(15).topo  = find((theta >= 5/4*pi) & (theta < 7/4*pi) & (rad < 45));
AHAseg(16).topo  = find((theta >= 3/4*pi) & (theta < 5/4*pi) & (rad < 45));
AHAseg(14).topo  = find((  ((theta >= 7/4*pi) & (theta < 8/4*pi)) | ...
                            ((theta >= 0)      & (theta < 1/4*pi)) )...
                         & (rad < 45));
                     
% Extract coordinates and scar transmurality for each segment                    
for jz = 1:nSegments

    % Coordinates in polar and cartesian systems
    AHAseg(jz).theta = theta(AHAseg(jz).topo);
    AHAseg(jz).rad = rad(AHAseg(jz).topo);
    [AHAseg(jz).x, AHAseg(jz).y] = pol2cart(AHAseg(jz).theta, AHAseg(jz).rad);
    
    % Data for each data point
    AHAseg(jz).data = DATA(AHAseg(jz).topo);
    AHAseg(jz).dataAvg = nanmean(AHAseg(jz).data);
    AHAseg(jz).dataStd = nanstd(AHAseg(jz).data);
end

%% Other things for plotting
AHAseg(1).thetaLim =  [1/3*pi 2/3*pi];      AHAseg(1).radLim =  [90 120];
AHAseg(2).thetaLim =  [0      1/3*pi];      AHAseg(6).radLim =  [90 120];
AHAseg(3).thetaLim =  [5/3*pi 6/3*pi];      AHAseg(5).radLim =  [90 120];
AHAseg(4).thetaLim =  [4/3*pi 5/3*pi];      AHAseg(4).radLim =  [90 120];
AHAseg(5).thetaLim =  [3/3*pi 4/3*pi];      AHAseg(3).radLim =  [90 120];
AHAseg(6).thetaLim =  [2/3*pi 3/3*pi];      AHAseg(2).radLim =  [90 120];

AHAseg(7).thetaLim  =  [1/3*pi 2/3*pi];      AHAseg(7).radLim =  [55 90];
AHAseg(8).thetaLim  = [0      1/3*pi];      AHAseg(12).radLim = [55 90];
AHAseg(9).thetaLim  = [5/3*pi 6/3*pi];      AHAseg(11).radLim = [55 90];
AHAseg(10).thetaLim = [4/3*pi 5/3*pi];      AHAseg(10).radLim = [55 90];
AHAseg(11).thetaLim =  [3/3*pi 4/3*pi];      AHAseg(9).radLim =  [55 90];
AHAseg(12).thetaLim =  [2/3*pi 3/3*pi];      AHAseg(8).radLim =  [55 90];

AHAseg(13).thetaLim = [1/4*pi 3/4*pi];      AHAseg(13).radLim = [15 55];
AHAseg(14).thetaLim = [7/4*pi 1/4*pi];      AHAseg(16).radLim = [15 55];
AHAseg(15).thetaLim = [5/4*pi 7/4*pi];      AHAseg(15).radLim = [15 55];
AHAseg(16).thetaLim = [3/4*pi 5/4*pi];      AHAseg(14).radLim = [15 55];

[contours.X1, contours.Y1] =pol2cart(0:.01:2*pi,120);
[contours.X2, contours.Y2] =pol2cart(0:.01:2*pi,90);
[contours.X3, contours.Y3] =pol2cart(0:.01:2*pi,55);
[contours.X0, contours.Y0] =pol2cart(0:.01:2*pi,15);
[contours.X4, contours.Y4] =pol2cart(0*pi/180,55:1:120);
[contours.X5, contours.Y5] =pol2cart(300*pi/180,55:1:120);
[contours.X6, contours.Y6] =pol2cart(240*pi/180,55:1:120);
[contours.X7, contours.Y7] =pol2cart(180*pi/180,55:1:120);
[contours.X8, contours.Y8] =pol2cart(120*pi/180,55:1:120);
[contours.X9, contours.Y9] =pol2cart(60*pi/180,55:1:120);
[contours.X10,contours.Y10]=pol2cart(315*pi/180,15:1:55);
[contours.X11,contours.Y11]=pol2cart(225*pi/180,15:1:55);
[contours.X12,contours.Y12]=pol2cart(135*pi/180,15:1:55);
[contours.X13,contours.Y13]=pol2cart(45*pi/180,15:1:55);
    
    
% Segment patches for plotting
for i = 1:nSegments

    if (i~=14)
        % Theta
        AHAseg(i).thetaPatch = [linspace(AHAseg(i).thetaLim(1), AHAseg(i).thetaLim(2), 100)...
                flip(linspace(AHAseg(i).thetaLim(1), AHAseg(i).thetaLim(2), 100))];
        
        % Center
        [AHAseg(i).xCenter, AHAseg(i).yCenter] = pol2cart(mean(AHAseg(i).thetaLim), mean(AHAseg(i).radLim));
    % Segment 14 crosses over theta = 0
    else
        AHAseg(i).thetaPatch = [linspace(AHAseg(i).thetaLim(1), 2*pi, 50) linspace(0, AHAseg(i).thetaLim(2), 50)...
                flip([linspace(AHAseg(i).thetaLim(1), 2*pi, 50) linspace(0, AHAseg(i).thetaLim(2), 50)])];
        [AHAseg(i).xCenter,AHAseg(i).yCenter] = pol2cart(0, mean(AHAseg(i).radLim));    
    end
    AHAseg(i).radPatch = [ones(1,100)*AHAseg(i).radLim(1) ones(1,100)*AHAseg(i).radLim(2)];    
    [AHAseg(i).xPatch, AHAseg(i).yPatch] = pol2cart(AHAseg(i).thetaPatch, AHAseg(i).radPatch);
    
end

%% Bullseye plot 
c = zeros(nSegments,1);
for i = 1:nSegments; c(i) = AHAseg(i).dataAvg;  end

if plotflag
    if axes_flag
        plotBE(AHAseg,c,contours,[0 round(max(DATA))],plasma(),dataname,axeshandle)
    else
        plotBE(AHAseg,c,contours,[0 round(max(DATA))],plasma(),dataname)
    end
end

%% Outputs
reg = DATA;
for jz = 1:16
    reg(AHAseg(jz).topo) = jz;
end

varargout{1} = reg;
varargout{2} = c;

%%
function plotBE(segments,c,contours,cLim,cMap,cBarTitle,varargin)
                             
txtFormat = '%1.0f';

af = false;
if ~isempty(varargin)
   af = true;
   ahandle = varargin{1};
end

%% Bullseye plot
if af
    axes(ahandle);
    hold on
else
h = figure('Name',cBarTitle,'NumberTitle','off','WindowStyle','docked'); hold on
end

% Patches - flip x axis as theta = 0 is in the septum (on the right)
for biv = 1:length(c)
    patch(-segments(biv).xPatch, segments(biv).yPatch, c(biv,:), 'EdgeColor', 'None')    
end

% Contour lines
plot(contours.X1,contours.Y1, 'k', 'LineWidth', 3);
plot(contours.X2,contours.Y2, 'k', 'LineWidth', 3);
plot(contours.X3,contours.Y3, 'k', 'LineWidth', 3);
plot(contours.X0,contours.Y0, 'k', 'LineWidth', 3);
plot(contours.X4,contours.Y4, 'k', 'LineWidth', 3);
plot(contours.X5,contours.Y5, 'k', 'LineWidth', 3);
plot(contours.X6,contours.Y6, 'k', 'LineWidth', 3);
plot(contours.X7,contours.Y7, 'k', 'LineWidth', 3);
plot(contours.X8,contours.Y8, 'k', 'LineWidth', 3);
plot(contours.X9,contours.Y9, 'k', 'LineWidth', 3);
plot(contours.X10,contours.Y10, 'k', 'LineWidth', 3);
plot(contours.X11,contours.Y11, 'k', 'LineWidth', 3);
plot(contours.X12,contours.Y12, 'k', 'LineWidth', 3);
plot(contours.X13,contours.Y13, 'k', 'LineWidth', 3);

text(0,128,'Anterior','fontsize',12,'fontweight','b','HorizontalAlignment','center')
text(-128,0,'Septal','fontsize',12,'Rotation',90,'fontweight','b','HorizontalAlignment','center')
text(128,0,'Lateral','fontsize',12,'Rotation',270,'fontweight','b','HorizontalAlignment','center')
text(0,-128,'Posterior','fontsize',12,'fontweight','b','HorizontalAlignment','center')
     
set(gca, 'LineWidth', 3)

shading flat;
colormap(cMap)

% axis([-120 120 -120 120]);
axis off equal;
set(gcf,'Color','white')

% % Set figure square
% pos = get(gcf,'Position');
% set(gcf, 'Position', [pos(1) pos(2) pos(3) pos(3)])


for wh = 1:length(c)
    
    % Print segments numbers if no patch data color coding
    % Else, place data value in center of each patch
        text(-segments(wh).xCenter, segments(wh).yCenter,num2str(c(wh),txtFormat),...
                'FontSize', 14, 'Color', [0 0 0],'HorizontalAlignment','center')
end
end
end