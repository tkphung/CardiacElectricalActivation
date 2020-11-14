function H = MC_HammerMap_WallThickness(STACK,varargin)
%MC_HammerMap_Scar: takes processed MRI segmentation from
%(MRI_ProcessStack) to generate a SALPS oriented Hammer Map
%   INPUT: STACK- from MRI_ProcessStack
%          VARARGIN:
%          		'title' - plot title
%               'surface' - 'endo' (auto) or 'epi'
%               'windowstyle' - 'docked' (auto)
%               'segment' - plot SA segs on Hammer
%               'savefig' - input Save Path
%               'bullseye' - false (auto) plot bullseye plot
%               'savename' - filename for saving
%               'clim' - color axis limits [min max]
%     STEP 1: Parse data points for Endo or Epi
%     STEP 2: Rotate into X (base-axis) Y (lateral-septal) Z (anterior-posterior)
%     STEP 3: Cartesian to Prolate
%     STEP 4: Create Hammer Grid (in prolate) & Fit SCAR data
%     STEP 5: Render Hammer Map
%     STEP 6: (OPTIONAL)- Bullseye Plot
%   OUTPUT: H- figure handle
%
% This particular function is specifically formatted for
% making a presentation-ready scar transmurality map.
% Thien-Khoi N. Phung (April 13, 2017)
% Edited from MC_HammerMap_Scar.m by TNP on May the Fourth be with you 2017
% This version maps WallThickness instead of Scar

% Deal with VARARGIN
plottitle = 'Wall Thickness Map';
windowstyle = 'docked';
surface = 'endo';
segment_flag = false;
savefig_flag = false;
bullseye_flag = false;
savename_flag = false;
clim_flag = false;
if ~isempty(varargin)
    for jz = 1:2:numel(varargin)
        switch varargin{jz}
            case 'title' % title figure window
                plottitle = varargin{jz+1};
			case 'surface'
				surface = varargin{jz+1};
            case 'windowstyle'
				windowstyle = varargin{jz+1};
            case 'segment'
				segment_flag = varargin{jz+1};
            case 'savefig'
				savefig_flag = true;
                savepath = varargin{jz+1};
            case 'bullseye'
				bullseye_flag = varargin{jz+1};
            case 'savename'
                savename_flag = true;
                savename = varargin{jz+1};
            case 'clim'
                clim_flag = true;
                clim = varargin{jz+1};
            otherwise
                error('ERROR: Check your varargins.')
        end
    end
end

%% STEP 1: Parse data points for Endo or Epi
switch surface
	case 'endo'
		data = STACK.alldata_endo;
	case 'epi'
        data = STACK.alldata_epi;
	otherwise
		error('ERROR: Specify endo or epi surface.')
end

% Cartesian XYZ 
XYZ = data(:,1:3);
WTH = data(:,5); % wall thickness (mm)

%% STEP 2: Rotate into X (base-axis) Y (lateral-septal) Z (anterior-posterior)
A = STACK.ApexPoint;
B = STACK.BasePoint;
S = STACK.SeptumPoint;

% First basis vector, subtract Base from Apex and divide by magnitude
C = A-B;
e1 = C ./ norm(C);

% Calculate Origin location
origin = B + (1/3)*C;

% Calculate focus length from norm(C)
Focus  = (2*norm(C)/3)/ cosh(1);
Nber_points = length(XYZ); % Assume confocal endo & epi

% Second basis vector using NEW METHOD- plane intersects septal point &
% e1
D = S(1,:) - origin;
D2 = D - dot(D,e1)*e1;
e2 = D2 ./ norm(D2);

% Third basis vector
E = cross(e1,e2);
e3 = E ./ norm(E);

% Transformation matrix
Transform = [e1; e2; -e3];

% Substract origin and transform (rotate) data
rotXYZ = (XYZ - repmat(origin,Nber_points,1))*Transform';

%% STEP 3: Cartesian to Prolate
[~,mu,the] = LV_C2P(rotXYZ(:,1),rotXYZ(:,2),rotXYZ(:,3),Focus);

%% STEP 4: Create Hammer Grid (in prolate) & Fit SCAR data
% Grid for Hammer Map (mu and theta values)
Ht = (0:1:359)*pi/180;
Hm = (0:1:120)*pi/180;
[Hthe, Hmu] =  meshgrid(Ht,Hm);

% Fit data for SCAR transmurality to the other mesh
%                known data -->  Hammer Grid
Hwth = griddata(the,mu,WTH(:,1),Hthe,Hmu);
Hwth(isnan(Hwth)) = 0;

%% STEP 4: Transform into Hammer Coordinates
% Hammer transformation
k = (1+cos(Hmu-(pi/2)).*cos((Hthe-pi)./2)).^(-0.5);
HamX = -k.*(cos(Hmu-(pi/2)).*sin((Hthe-pi)./2));
HamY = k.*(sin(Hmu-(pi/2)));

% (OPTIONAL)- Check where Hammer grid points are for segmentation
if segment_flag
    k = (1+cos(mu-(pi/2)).*cos((the-pi)./2)).^(-0.5);
    HamsegX = -k.*(cos(mu-(pi/2)).*sin((the-pi)./2));
    HamsegY = k.*(sin(mu-(pi/2)));
end

%% STEP 5: Render Hammer Map
%--------------------------------------------------------------------------
% NOTE: THIS IS FOR HAMMER MAP WITH SALPS ORIENTATION (left to right:
% septum, anterior, lateral, posterior, septum). 
%--------------------------------------------------------------------------

% for correct coloring, make areas with no scar =-0.25 instead of 0
% Hwth(Hwth==0)=-0.25; 
Hwth(Hwth==0) = NaN;

% Divisions for contour lines
contlines=linspace(min(min(Hwth)),max(max(Hwth)),7);

% Plot SCAR
H = figure('Name',plottitle,'NumberTitle','off','WindowStyle',windowstyle); hold on
[C, h] = contourf(HamX, HamY, Hwth, contlines);
    set(h,'LineWidth',2)
%   colormap(repmat((1:-0.1:0)',1,3));
	if clim_flag
        caxis(clim)
    else
        caxis([5 20])
    end
	axis equal off
    set(h,'LineWidth',2)
	set(gcf,'Color','white')
    text(-.9,.55,'Septum','fontsize',12,'fontweight','b')
    text(-.48,.45,'Anterior','fontsize',12,'fontweight','b')
    text(-.05,.42,'Lateral','fontsize',12,'fontweight','b')
    text(.28,.48,'Posterior','fontsize',12,'fontweight','b')
    text(.78,.55,'Septum','fontsize',12,'fontweight','b')
    text(-.7,-.95,'Apex','fontsize',12,'fontweight','b')


    hold on;
    plot(HamX(120,1:360), HamY(120,1:360), 'k','LineWidth',2);
    plot(HamX(45,1:360), HamY(45,1:360), 'k','LineWidth',2);
    plot(HamX(86,1:360), HamY(86,1:360), 'k','LineWidth',2);
    plot(HamX(1:120,360), HamY(1:120,360), 'k','LineWidth',2);
    plot(HamX(45:120,300), HamY(45:120,300), 'k','LineWidth',2);
    plot(HamX(45:120,240), HamY(45:120,240), 'k','LineWidth',2);
    plot(HamX(45:120,180), HamY(45:120,180), 'k','LineWidth',2);
    plot(HamX(45:120,120), HamY(45:120,120), 'k','LineWidth',2);
    plot(HamX(45:120,60), HamY(45:120,60), 'k','LineWidth',2);
    plot(HamX(1:120,1), HamY(1:120,1), 'k','LineWidth',2);
    plot(HamX(1:45,315), HamY(1:45,315), 'k','LineWidth',2);
    plot(HamX(1:45,225), HamY(1:45,225), 'k','LineWidth',2);
    plot(HamX(1:45,135), HamY(1:45,135), 'k','LineWidth',2);
    plot(HamX(1:45,45), HamY(1:45,45), 'k','LineWidth',2);

colorbar

% (OPTIONAL)- Plot segment contours
if segment_flag
    plot(HamsegX,HamsegY,'k.','MarkerSize',15)
end

if savefig_flag
    if savename_flag
        saveas(gcf,[savepath savename '.fig']);
        cdhome = pwd;
        cd(savepath);
        print(savename,'-dtiffn')
        cd(cdhome);
    else
        saveas(gcf,[savepath 'WTH hammer map ' date '.fig']);
    end
end

%% STEP 6: (OPTIONAL)- Bullseye Plot
if bullseye_flag
%Bullseye Map
rho=linspace(0,120,100);
angl=linspace(0,2*pi,100);
[thetapolar,radpolar]=meshgrid(angl,rho);
BullseyeTheta=[0:1:360];
BullseyeMu=[0:1:120]';
Hwth(:,361)=Hwth(:,1);
Z_subtract_bullseye=[Hwth(:,181:361) Hwth(:,1:180)];
A=interp2((BullseyeTheta*pi/180),BullseyeMu,Z_subtract_bullseye,angl,rho');
[xpolar,ypolar,zpolar]=pol2cart(thetapolar,radpolar,A);
figure('Name','Bullseye Thickness','NumberTitle','off','WindowStyle',windowstyle); hold on
polar([0 2*pi],[0 100]);
[C1,h1]=contourf(xpolar,ypolar,zpolar,contlines);
shading flat;
axis off; axis([-120 120 -120 120]);
contour(xpolar,ypolar,zpolar,contlines,'k','LineWidth',2);
% colormap(repmat((1:-0.1:0)',1,3));
if clim_flag
    caxis(clim)
else
    caxis([5 20])
end
axis equal;
[X1,Y1]=pol2cart([0:.01:2*pi],120);%draw circle
plot(X1,Y1,'k','LineWidth',2);%draw circle
set(h1,'LineWidth',2)
set(gcf,'Color','white')
text(0,128,'Anterior','fontsize',12,'fontweight','b','HorizontalAlignment','center')
text(-128,0,'Septum','fontsize',12,'Rotation',90,'fontweight','b','HorizontalAlignment','center')
text(128,0,'Lateral','fontsize',12,'Rotation',270,'fontweight','b','HorizontalAlignment','center')
text(0,-128,'Posterior','fontsize',12,'fontweight','b','HorizontalAlignment','center')
    [X2,Y2]=pol2cart([0:.01:2*pi],86);
    [X3,Y3]=pol2cart([0:.01:2*pi],45);
    [X4,Y4]=pol2cart(0*pi/180,45:1:120);
    [X5,Y5]=pol2cart(300*pi/180,45:1:120);
    [X6,Y6]=pol2cart(241*pi/180,45:1:120);
    [X7,Y7]=pol2cart(180*pi/180,45:1:120);
    [X8,Y8]=pol2cart(120*pi/180,45:1:120);
    [X9,Y9]=pol2cart(60*pi/180,45:1:120);
    [X10,Y10]=pol2cart(315*pi/180,1:1:45);
    [X11,Y11]=pol2cart(225*pi/180,1:1:45);
    [X12,Y12]=pol2cart(135*pi/180,1:1:45);
    [X13,Y13]=pol2cart(45*pi/180,1:1:45);
    plot(X2,Y2,'k','LineWidth',2);plot(X3,Y3,'k','LineWidth',2);plot(X4,Y4,'k','LineWidth',2);
    plot(X5,Y5,'k','LineWidth',2);plot(X6,Y6,'k','LineWidth',2);plot(X7,Y7,'k','LineWidth',2);
    plot(X8,Y8,'k','LineWidth',2);plot(X9,Y9,'k','LineWidth',2);plot(X10,Y10,'k','LineWidth',2);
    plot(X11,Y11,'k','LineWidth',2);plot(X12,Y12,'k','LineWidth',2);plot(X13,Y13,'k','LineWidth',2);

if savefig_flag
    if savename_flag
        saveas(gcf,[savepath savename 'bullseye' '.fig']);
    else
        saveas(gcf,[savepath 'WTH bullseye map ' date '.fig']);
    end
end

end