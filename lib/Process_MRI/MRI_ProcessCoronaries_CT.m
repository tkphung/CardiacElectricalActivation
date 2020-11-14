function CoronariesXYZ = MRI_ProcessCoronaries_CT(filename,LA,SA)
%MRI_ProcessCoronaries: Loads in ROI segmentation for coronaries and
%rotates them into MRI coordinates (torso coordinates) that registers with
%the ProcessedStack coordinates from MRI_ProcessStack.m
% CoronariesXYZ = MRI_ProcessCoronaries(filename)
%   INPUT VARIABLE:
%       filename- path+filename of the coronaries segmentation (from
%       Segment)
%       ProcessedStack- from MRI_ProcessStack to get center_axis_diff
%       LGESA- filename of the matched SA.mat used in processedstack
%   CALCULATION:
%       Load file
%       Parse ROI Coordinates
%       Parse coordinate tranformations
%       Transformations to Human Coordinates
%       Rotate back to MRI coordinates
%       Shift center based on ProcessedStack centroids
%   OUTPUT VARIABLE:
%       CoronariesXYZ- MRI coordinates of adjusted coronaries
%
% Thien-Khoi N. Phung (February 8, 2017)

load(filename);

x= []; y = []; z= [];
% Parse ROI coordinates
for jz = 1:length(setstruct.Roi)
    x = [x; setstruct.Roi(jz).X];
    y = [y; setstruct.Roi(jz).Y];
    z = [z; repmat(setstruct.Roi(jz).Z,80,1)];    
end

% Parse coordinate transformations
    xres = setstruct.ResolutionX;
    yres = setstruct.ResolutionY;
    imagepos = setstruct.ImagePosition;
    xorient = setstruct.ImageOrientation(4:6);
    yorient = setstruct.ImageOrientation(1:3);
    zorient = cross(yorient,xorient);
    zoffset = (z-1);
    z = -zoffset;

% Transformations to 3D Imaging Coordinates for the Coronaries
    To = eye(4);
        To(1:3,4) = [-1 -1 0];
    S = eye(4);
        S(1,1) = xres;
        S(2,2) = yres;
        S(3,3) = setstruct.SliceThickness + setstruct.SliceGap;
    R = eye(4);
        R(1:3,1:3) = [xorient(:), yorient(:), zorient(:)];
    Tipp = eye(4);
        Tipp(1:3,4) = imagepos';
    M = Tipp*R*S*To;

    xyz = M * [x y z ones(length(x),1)]';

    CoronariesXYZ = xyz(1:3,:)';

        % Rotate to 3D CT coordinates for the specified LGESA file
        % load(SA)
        %     x_image_orientation = setstruct.ImageOrientation(4:6);
        %     y_image_orientation = setstruct.ImageOrientation(1:3);
        %     z_image_orientation = cross(y_image_orientation,x_image_orientation);
        %     M = [x_image_orientation(:), y_image_orientation(:), z_image_orientation(:)]';
        %     transform_cor = [M*xyz(:,1:3)']';
        % 
        %     % figure(1), hold on
        %     % plot3(xyz(:,1),xyz(:,2),xyz(:,3),'r.')
        %     % axis equal
        % 
        % % Shift center
        % cor_shifter = mean(ProcessedStack.center_axis_diff);
        % CoronariesXYZ = transform_cor - cor_shifter;
        % 
        %     % figure(2), hold on
        %     % plot3(cor_shifted(:,1),cor_shifted(:,2),cor_shifted(:,3),'r.')
        %     % axis equal
    
% Make Transformation matrix and origin (3D MRI to Cardiac Coordinates)
%% Septal Point
load(SA)
    % STEP 2: Process KeptSlices Variable
    % Add KeptSlices variable to structure (this variable indicates the number
    % of slices (1:setstruct(jz).ZSize) per imaging view (jz))
    for jz = 1:length(setstruct)
        setstruct(jz).KeptSlices = [1:setstruct(jz).ZSize];
    end

    % Adjust KeptSlices variable to eliminate untraced images so that 
    % epi/endo + scar traces will correctly register with each other
    % (Note: scar is only traced on images with endo & epi contours)
    no_trace = sum(isnan(setstruct.EndoX(:, :,:))); % counts slices with no trace
    delete_slices = no_trace ~= 0; % which have tracings
    delete_slices = sum(delete_slices,2) == size(delete_slices,2);
    setstruct.KeptSlices(:,squeeze(delete_slices)) = []; % removes those slice numbers

    % Determine which time index (column) was traced
        timeID = find(no_trace(:,:,setstruct.KeptSlices(1)) == 0);

    % STEP 3: Identify RV insertion slice
    for jz = setstruct.KeptSlices
        convert = cell2mat(setstruct.EndoPinX(:,jz));
        contents(jz) = sum(convert);
    end
    septal_slice = find(contents);

    % Check for if more than 1 time frame is segmented in SA MRI
    % Use the time frame that the RV insertions were pinpointed
    if numel(timeID)>1
        [timeID,~] = find(~cellfun('isempty',setstruct.EndoPinX));
    end

    % STEP 4: Find endocardial mid-septal point
    % 4a: Find midpoint of the line connecting the RV insertion points
    RVInsertionPts = [setstruct.EndoPinX{timeID,septal_slice},setstruct.EndoPinY{timeID,septal_slice}];
if size(RVInsertionPts,1) > 1    

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
elseif size(RVInsertionPts,1) == 1
    septx = setstruct.EndoPinX{timeID,septal_slice};
    septy = setstruct.EndoPinY{timeID,septal_slice};
end
    % Store those points in the originals etstruct
    setstruct.EndoPinX{timeID,septal_slice} = [setstruct.EndoPinX{timeID,septal_slice};septx];
    setstruct.EndoPinY{timeID,septal_slice} = [setstruct.EndoPinY{timeID,septal_slice};septy];

    % STEP CORRECTION pre-5: Elimate the tracings in other time frames
    % CORRECTION MADE JULY 18, 2017 by TNP
    % Replace the non-timeID columns in KeptSlices with NaN
    nonTIME = 1:size(setstruct.EndoX,2); nonTIME = nonTIME~=timeID;

    setstruct.EndoX(:,nonTIME,setstruct.KeptSlices) = NaN;
    setstruct.EndoY(:,nonTIME,setstruct.KeptSlices) = NaN;
    setstruct.EpiX(:,nonTIME,setstruct.KeptSlices) = NaN;
    setstruct.EpiY(:,nonTIME,setstruct.KeptSlices) = NaN;

    % STEP 5: Rotate countour stack, RV insertion, and mid-septal points
    % rotate stack + RV insertion + mid-septal points
    [~,RVInsertionPts,~] = MRI_RotateEndoStack(setstruct);
    
    Si = RVInsertionPts(end,:);
%% Apex and Base Points
load(LA)   
[ApexBasePts,~] = MRI_TransformLAWithPinPoints(setstruct); 
    
    A  = ApexBasePts(1,:);
    B  = ApexBasePts(2,:);
    
%%
    % First basis vector, subtract Base from Apex and divide by magnitude
    C = A-B;
    e1 = C ./ norm(C);

    % Calculate Origin location
    origin = B + (1/3)*C;

    % Second basis vector using NEW METHOD- plane intersects septal point &
    % e1
    D = Si(1,:) - origin;
    D2 = D - dot(D,e1)*e1;
    e2 = D2 ./ norm(D2);

    % Third basis vector
    E = cross(e1,e2);
    e3 = E ./ norm(E);

    % Transformation matrix (3D MRI to Cardiac Coordinates)
    Transform = [e1; e2; -e3];
    
CoronariesXYZ = (CoronariesXYZ - repmat(origin,length(CoronariesXYZ),1))*Transform';
