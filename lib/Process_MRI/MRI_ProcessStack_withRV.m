function ProcessedStack = MRI_ProcessStack_withRV(SA,LA,varargin)
%MRI_ProcessStack_withRV: Takes in segmented MRI (using Segment) in a short axis
%and long axis view and shifts & rotates it into a ProcessedStack that is
%formatted for fitting to a FE mesh
% ProcessedStack = MRI_ProcessStack(SA,LA,varagin)
%   INPUT VARIABLE:
%       SA- short axis path+name; includes endo/epi contour & RV insert
%       pinpoints (and optional- scar contours)
%       LA- long axis path+name; includes base-apex pinpoints
%       varargin- 'scar' if scar contour is present
%                 'time' to specify time frame
%   CALCULATIONS:
%		Preprocess Flags
%		STEP 1: Load in Short Axis Stack
%		STEP 2: Process KeptSlices Variable
%				Add KeptSlices variable to structure (this variable indicates the number
%				of slices (1:setstruct(jz).ZSize) per imaging view (jz))
%				Adjust KeptSlices variable to eliminate untraced images so that 
%				epi/endo + scar traces will correctly register with each other
%				(Note: scar is only traced on images with endo & sepi contours)
%		STEP 3: Identify RV insertion slice
%		STEP 4: Find endocardial mid-septal point
%			4a: Find midpoint of the line connecting the RV insertion points
%			4b: Find slope of line perpendicular to line connecting RV insertion points
%			4c: Find intersection of endocardium with perpendicular line projected
%				from midpoint that is on the septum (automated by TNP)
%				Shift slice to MidPt as the center of endocardium
%				Convert to polar coordinates
%				Find theta for perpendicular slope (perpslope)
%				Interpolate the radii for the two angles
%				Choose the smaller radius to designate the septal point
%				Calculate mid-septal endocardial coordinate and shift back to original
%				origin
%				Store those points in the original setstruct
%       STEP CORRECTION pre-5: Elimate the tracings in other time frames
%		STEP 5: Rotate countour stack, RV insertion, and mid-septal points
%				rotate stack + RV insertion + mid-septal points
%				Rotate SA Epicardial contours
%		STEP 6: (OPTIONAL-Flag) Restructure Scar Variable Information
%				Restructure Scar Variable to be similar to Endo + Epi
%				Rotate Scar Contours  
%		STEP 7: Load in Long Axis Stack
%		STEP 8: Rotate Endocardial Apex + Base Points from LA Image
%			8a: Transform Apex and Base Points
%			8b: Rotate stack back so that slices are in XY plane and Z is vertical
%			8c: Define Base-Apex Axis, Base point - Apex point
%			8d: Find m value of each slice based on Z location of slice
%			8e: Use m values to find the x+y coordinates where the slice
%				intersects the BA axis
%			8f: Shift center of slice to correspond to base-apex axis THIS
%			STEP CHANGED FOR RV Fit
%			8g: Shift septal point by the same amount the slice was shifted by to align
%				with the AB axis
%				Find new number of slice with septal point, now that untraced basal
%				slices might have been removed
%		STEP 9: Sort Shifted + Rotated Slices into a single structure
%       STEP 10: Calculate Wall Thickness (And optional scar transmurality)
%       STEP 11: Construct alldata variable for fitting code
%   OUTPUT VARIABLE:
%       ProcessedStack- one structure containing contours (epi, endo, and
%       (optional) scar) rotated into MRI coordinates (3D)
% 
% Thien-Khoi N. Phung (started January 18, 2017)
% CHANGED TNP March 17, 2017
% CHANGED TNP March 29, 2017- adapted from MRI_ProcessStack.m to include RV
% fitting. The centering of each stack is based on the endocardial centroid
% instead of the epicardial centroid. The endocardial segmentation and
% pinpoints should be the same for all of the 3 segmentation files with the
% epicardial tracing including from LVepi, RVendo, to RVepi.
% NOTE: This code should work the same as the MRI_ProcessStack EXCEPT that
% is used the ENDO centroid INSTEAD of the EPI centroid.
% Edited by TNP (October 10, 2017)- updated the TimeID_flag from
% MRI_ProcessStack.m 

%% Preprocess Flags
scar_flag = false;
timeID_flag = false;
if ~isempty(varargin)
    for jz = 1:2:numel(varargin)
        switch varargin{jz}
            case 'scar' % title figure window
                scar_flag = varargin{jz+1};
			case 'time'
                timeID_flag = true;
				timeID = varargin{jz+1};
            otherwise
                error('ERROR: Check your varargins.')
        end
    end
end

%% STEP 1: Load in Short Axis Stack
load(SA);

%% STEP 2: Process KeptSlices Variable
% Add KeptSlices variable to structure (this variable indicates the number
% of slices (1:setstruct(jz).ZSize) per imaging view (jz))
for jz = 1:length(setstruct)
    setstruct(jz).KeptSlices = [1:setstruct(jz).ZSize];
end

% Adjust KeptSlices variable to eliminate untraced images so that 
% epi/endo + scar traces will correctly register with each other
% (Note: scar is only traced on images with endo & sepi contours)
no_trace = sum(isnan(setstruct.EndoX(:, :,:))); % counts slices with no trace
delete_slices = no_trace ~= 0; % which have tracings
delete_slices = sum(delete_slices,2) == size(delete_slices,2);
setstruct.KeptSlices(:,squeeze(delete_slices)) = []; % removes those slice numbers

% Determine which time index (column) was traced
if ~timeID_flag
    timeID = find(no_trace(:,:,setstruct.KeptSlices(1)) == 0);
end

%% STEP 3: Identify RV insertion slice
for jz = setstruct.KeptSlices
    convert = cell2mat(setstruct.EndoPinX(:,jz));
    contents(jz) = sum(convert);
end
septal_slice = find(contents);

% Check for if more than 1 time frame is segmented in SA MRI
% Use the time frame that the RV insertions were pinpointed
if numel(timeID)>1 && ~timeID_flag
    [timeID,~] = find(~cellfun('isempty',setstruct.EndoPinX));
end

%% STEP 4: Find endocardial mid-septal point
% 4a: Find midpoint of the line connecting the RV insertion points
RVInsertionPts = [setstruct.EndoPinX{timeID,septal_slice},setstruct.EndoPinY{timeID,septal_slice}];
MidPt = mean(RVInsertionPts);

% 4b: Find slope of line perpendicular to line connecting RV insertion points
slope = (RVInsertionPts(2,2)-RVInsertionPts(1,2))/(RVInsertionPts(2,1)-RVInsertionPts(1,1));
perpslope = -1/slope;

% 4c: Find intersection of endocardium with perpendicular line projected
% from midpoint that is on the septum (automated by TNP)
% Shift slice to MidPt as the center of endocardium
    slice = [setstruct.EndoX(:,timeID,septal_slice),setstruct.EndoY(:,timeID,septal_slice)];
    slice = slice-repmat(MidPt,size(slice,1),1);
% Convert to polar coordinates
    [TH,R] = cart2pol(slice(:,1),slice(:,2));
    TH(TH<0) = TH(TH<0) + 2*pi;
    TH = TH(2:end); % beginning and end points are the same
    R = R(2:end);
% Find theta for perpendicular slope (perpslope)
    dotprod =  dot([1 1*perpslope],[1 0]);
    th1 = acos(dotprod/norm([1 1*perpslope]));
    th2 = th1 + pi; % 180 degree rotation
% Interpolate the radii for the two angles
    r1 = interp1(TH,R,th1,'linear');
    r2 = interp1(TH,R,th2,'linear');
% Choose the smaller radius to designate the septal point
    r = r1*(r1<r2) + r2*(r2<r1);
    th = th1*(r1<r2) + th2*(r2<r1);
% Calculate mid-septal endocardial coordinate and shift back to original
% origin
    [septx,septy] = pol2cart(th,r);
    septx = septx + MidPt(1);
    septy = septy + MidPt(2);
% Store those points in the original setstruct
setstruct.EndoPinX{timeID,septal_slice} = [setstruct.EndoPinX{timeID,septal_slice};septx];
setstruct.EndoPinY{timeID,septal_slice} = [setstruct.EndoPinY{timeID,septal_slice};septy];

%% STEP CORRECTION pre-5: Elimate the tracings in other time frames
% CORRECTION MADE October 10, 2017 by TNP
% Replace the non-timeID columns in KeptSlices with NaN
nonTIME = 1:size(setstruct.EndoX,2); nonTIME = nonTIME~=timeID;

setstruct.EndoX(:,nonTIME,setstruct.KeptSlices) = NaN;
setstruct.EndoY(:,nonTIME,setstruct.KeptSlices) = NaN;
setstruct.EpiX(:,nonTIME,setstruct.KeptSlices) = NaN;
setstruct.EpiY(:,nonTIME,setstruct.KeptSlices) = NaN;


%% STEP 5: Rotate countour stack, RV insertion, and mid-septal points
% rotate stack + RV insertion + mid-septal points
[Cxyz_SAEndo,RVInsertionPts,~] = MRI_RotateEndoStack(setstruct);

% Rotate SA Epicardial contours
[Cxyz_SAEpi,~] = MRI_RotateEpiStack(setstruct);

SAstack = setstruct;
clear im info preview setstruct

%% STEP 6: (OPTIONAL-Flag) Restructure Scar Variable Information
if scar_flag
    % Restructure Scar Variable to be similar to Endo + Epi
    max_scar = max(sum(sum(SAstack.Scar.Manual,1),2));
    SAstack.Scar.ScarX = zeros(max_scar,1,SAstack.ZSize);
    SAstack.Scar.ScarY = zeros(max_scar,1,SAstack.ZSize);
    for i = 1:SAstack.ZSize
        [x,y] = find(SAstack.Scar.Manual(:,:,i) == 1);
        SAstack.Scar.ScarX(1:length(x),1,i) = x;
        SAstack.Scar.ScarY(1:length(y),1,i) = y;
    end

    %Rotate Scar Contours
    [Cxyz_SAScar,~] = MRI_RotateScarStack(SAstack);
    for j = 1:SAstack.ZSize
        slice_ind = Cxyz_SAScar(:,5)==j;
        slice = Cxyz_SAScar(slice_ind,:);
        xmode = mode(slice(:,1));
        mode_idx = Cxyz_SAScar(:,1) == xmode;
        Cxyz_SAScar(mode_idx,:) = [];
    end
end

%%% PLOT OPTION
% figure
% xlabel('x');ylabel('y');zlabel('z');
% plot3(Cxyz_SAEndo(:,1),Cxyz_SAEndo(:,2),Cxyz_SAEndo(:,3),'k.','MarkerSize',5)
% hold on
% plot3(Cxyz_SAEpi(:,1),Cxyz_SAEpi(:,2),Cxyz_SAEpi(:,3),'k.','MarkerSize',5)
% hold on
% plot3(RVInsertionPts(1:2,1),RVInsertionPts(1:2,2),RVInsertionPts(1:2,3),'g*','MarkerSize',10)
% plot3(RVInsertionPts(3,1),RVInsertionPts(3,2),RVInsertionPts(3,3),'b*','MarkerSize',10)
% axis equal
% hold on
% plot3(Cxyz_SAScar(:,1),Cxyz_SAScar(:,2),Cxyz_SAScar(:,3),'m.','MarkerSize',5)

%% STEP 7: Load in Long Axis Stack
load(LA);

%% STEP 8: Rotate Endocardial Apex + Base Points from LA Image
% 8a: Transform Apex and Base Points
if timeID_flag % IF TIME IS SPECIFIED - edit TNP 05-09-17 for CINE MRI processing
    [ApexBasePts,~] = MRI_TransformLAWithPinPoints(setstruct,'time',timeID);
else
    [ApexBasePts,~] = MRI_TransformLAWithPinPoints(setstruct);
end

A = ApexBasePts(end-1,:);
B = ApexBasePts(end,:);
% S=RVInsertionPts(3,:);
S = RVInsertionPts([3,1,2],:); %CHANGED to INCLUDE MID SEPTUM & INSERTIONS

%%% PLOT OPTION
% plot3(ApexBasePts(:,1),ApexBasePts(:,2),ApexBasePts(:,3),'b*','MarkerSize',10)
% xlabel('x');ylabel('y');zlabel('z');

% 8b: Rotate stack back so that slices are in XY plane and Z is vertical
x_image_orientation = SAstack.ImageOrientation(4:6);
y_image_orientation = SAstack.ImageOrientation(1:3);
z_image_orientation = cross(y_image_orientation,x_image_orientation);
M = [x_image_orientation(:), y_image_orientation(:), z_image_orientation(:)]';
transform_endo = [M*Cxyz_SAEndo(:,1:3)']';
transform_epi = [M*Cxyz_SAEpi(:,1:3)']';
if scar_flag
    transform_scar = [M*Cxyz_SAScar(:,1:3)']';
end
transform_ABSpts = [M*[A;B;S]']';

% 8c: Define Base-Apex Axis, Base point - Apex point
r = transform_ABSpts(2,:)-transform_ABSpts(1,:);
m = 0:0.01:1;
x = transform_ABSpts(1,1)+r(1,1).*m;
y = transform_ABSpts(1,2)+r(1,2).*m;
z = transform_ABSpts(1,3)+r(1,3).*m;

% 8d: Find m value of each slice based on Z location of slice
count=0;
for i = SAstack.KeptSlices
    count=count+1;
    idx = find(Cxyz_SAEndo(:,5)==i,1);
    z_loc(count) = transform_endo(idx,3);
end
m_val_slices = ((z_loc-transform_ABSpts(1,3))./r(1,3))';

% 8e: Use m values to find the x+y coordinates where the slice
% intersects the BA axis
BAaxis_intersect_pts = [transform_ABSpts(1,1)+r(1,1).*m_val_slices,transform_ABSpts(1,2)+r(1,2).*m_val_slices,z_loc'];

% 8f: Shift center of slice to correspond to base-apex axis
% CHANGED 3-17-2017 TNP: Shift each slice using ENDO centroid
count=0;
for i = SAstack.KeptSlices
    count = count+1;
    slice_idx_endo = find(Cxyz_SAEndo(:,5) == i);
    slice_idx_epi = find(Cxyz_SAEpi(:,5) == i);
    if scar_flag
        slice_idx_scar = find(Cxyz_SAScar(:,5) == i);
    end
    slice_center(i,:) = mean(transform_endo(slice_idx_endo,:)); % CHANGED TO _endo FROM _epi
    axis_center(i,:) = BAaxis_intersect_pts(count,:);
    center_axis_diff(i,:) = slice_center(i,:)-axis_center(i,:);
    Endo_shifted(slice_idx_endo,:) = [transform_endo(slice_idx_endo,1)-center_axis_diff(i,1),transform_endo(slice_idx_endo,2)-center_axis_diff(i,2),transform_endo(slice_idx_endo,3)];
    Epi_shifted(slice_idx_epi,:) = [transform_epi(slice_idx_epi,1)-center_axis_diff(i,1),transform_epi(slice_idx_epi,2)-center_axis_diff(i,2),transform_epi(slice_idx_epi,3)];
    if scar_flag
        Scar_shifted(slice_idx_scar,:) = [transform_scar(slice_idx_scar,1)-center_axis_diff(i,1),transform_scar(slice_idx_scar,2)-center_axis_diff(i,2),transform_scar(slice_idx_scar,3)];
    end
end

% 8g: Shift septal point by the same amount the slice was shifted by to align
% with the AB axis
% Find new number of slice with septal point, now that untraced basal
% slices might have been removed
stack_shift = Cxyz_SAEndo(1,5)-1;
septal_slice_new = septal_slice-stack_shift;
% ABSpts_shifted = transform_ABSpts- [0,0,0;0,0,0;center_axis_diff(septal_slice_new,1),center_axis_diff(septal_slice_new,2),0];
ABSpts_shifted = transform_ABSpts- [0,0,0;0,0,0;center_axis_diff(septal_slice_new,1),center_axis_diff(septal_slice_new,2),0;...
    center_axis_diff(septal_slice_new,1),center_axis_diff(septal_slice_new,2),0;center_axis_diff(septal_slice_new,1),center_axis_diff(septal_slice_new,2),0]; %CHANGED
Sslice = [ABSpts_shifted(3,:)-axis_center(septal_slice_new,:)];

%%% PLOT OPTION
% figure
% plot3(transform_ABSpts(:,1),transform_ABSpts(:,2),transform_ABSpts(:,3),'b*','MarkerSize',10)
% hold on
% plot3(x,y,z,'c')
% hold on
% plot3(Endo_shifted(:,1),Endo_shifted(:,2),Endo_shifted(:,3),'k.','MarkerSize',5)
% hold on
% plot3(Epi_shifted(:,1),Epi_shifted(:,2),Epi_shifted(:,3),'k.','MarkerSize',5)
% % hold on
% % plot3(Scar_shifted(:,1),Scar_shifted(:,2),Scar_shifted(:,3),'mo','MarkerSize',5)
% xlabel('x'); ylabel('y'); zlabel('z');
% axis equal
% title('Slices Shifted so Slice X-Y Center (Calculated from ENDO Contours) Lies on Apex-Base Axis');

%% STEP 9: Sort Shifted + Rotated Slices into a single structure
ProcessedStack={};
ProcessedStack.Slices = SAstack.KeptSlices;
ProcessedStack.ApexPoint = ABSpts_shifted(1,:);
ProcessedStack.BasePoint = ABSpts_shifted(2,:);
% ProcessedStack.SeptumPoint = ABSpts_shifted(3,:);
ProcessedStack.SeptumPoint = ABSpts_shifted([3,4,5],:); % Stores midseptum and septal intercepts
ProcessedStack.SliceCenters = axis_center;
for jz = 1:max(Cxyz_SAEndo(:,5))
    slice_idx_endo = find(Cxyz_SAEndo(:,5) == jz);
    slice_idx_epi = find(Cxyz_SAEpi(:,5) == jz);
    ProcessedStack.Endo{jz} = Endo_shifted(slice_idx_endo,:);
    ProcessedStack.Epi{jz} = Epi_shifted(slice_idx_epi,:);
    if scar_flag
        slice_idx_scar = find(Cxyz_SAScar(:,5) == jz);
        ProcessedStack.Scar{jz} = Scar_shifted(slice_idx_scar,:);
    end
end
ProcessedStack.Xres = SAstack.ResolutionX;
ProcessedStack.Yres = SAstack.ResolutionY;
ProcessedStack.SliceThickness = SAstack.SliceThickness;
ProcessedStack.SliceGap = SAstack.SliceGap;
if scar_flag
    ProcessedStack.WallMask = SAstack.Scar.MyocardMask;
    ProcessedStack.ScarMask = SAstack.Scar.Result;
end
ProcessedStack.center_axis_diff = center_axis_diff; % Saved for MRI_ProcessCoronaries.m

%% STEP 10: Calculate Wall Thickness (And optional scar transmurality)
% ProcessedStack = MRI_WallThicknessScarTransmurality(ProcessedStack,scar_flag);
ProcessedStack = MRI_WallThicknessScarTransmurality_withRV(ProcessedStack,scar_flag);

%% STEP 11: Construct alldata variable for fitting code
slices_all= ProcessedStack.Slices;
alldata_endo=[];
alldata_epi=[];

for i = slices_all
    endopol_all = ProcessedStack.EndoPolar{i};
    endopol = [endopol_all(:,2),endopol_all(:,4)];
    [xendocart,yendocart] = pol2cart(endopol(:,1),endopol(:,2));
    endoz_all = ProcessedStack.Endo{i};
    endoz = repmat(endoz_all(1,3),length(xendocart),1);
    epipol_all = ProcessedStack.EpiPolar{i};
    epipol = [epipol_all(:,2),epipol_all(:,4)];
    [xepicart,yepicart] = pol2cart(epipol(:,1),epipol(:,2));
    epiz_all = ProcessedStack.Epi{i};
    epiz = repmat(epiz_all(1,3),length(xepicart),1);
    
    % Shift data in x-y to reverse polar origin shift + line back up with original contours
    % CHANGED TNP 03-17-2017 to use ENDO center
    polar_center = [mean(xendocart),mean(yendocart)];
    center_diff = [polar_center(1,1)-ProcessedStack.SliceCenters(i,1),polar_center(1,2)-ProcessedStack.SliceCenters(i,2)];
    xendocart = xendocart-center_diff(1,1);
    yendocart = yendocart-center_diff(1,2);
    xepicart = xepicart-center_diff(1,1);
    yepicart = yepicart-center_diff(1,2);
    
    transmurality_all = ProcessedStack.TransmuralityPolar{i};
    transmurality = transmurality_all(:,2);
    wall_all = ProcessedStack.WallThicknessPolar{i};
    wall = wall_all(:,2);
    
    % alldata matrix contains columns of x, y, z, percent scar, and wall thickness
    alldata_endo_slice = [xendocart yendocart endoz transmurality wall];
    alldata_epi_slice = [xepicart yepicart endoz transmurality wall];
    alldata_endo=[alldata_endo; alldata_endo_slice];
    alldata_epi=[alldata_epi; alldata_epi_slice];
    ProcessedStack.alldata_endo= alldata_endo;
    ProcessedStack.alldata_epi= alldata_epi;
end
ProcessedStack.AvgWallThickness = mean(alldata_endo(:,5)); %mm