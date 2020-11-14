function CoronariesXYZ = MRI_ProcessCoronaries(filename,ProcessedStack,LGESA)
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
load(ProcessedStack);

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

% Transformations to 3D MRI for the Coronaries
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

    xyz = xyz(1:3,:)';

% Rotate to 3D MRI coordinates for the specified LGESA file
load(LGESA)
    x_image_orientation = setstruct.ImageOrientation(4:6);
    y_image_orientation = setstruct.ImageOrientation(1:3);
    z_image_orientation = cross(y_image_orientation,x_image_orientation);
    M = [x_image_orientation(:), y_image_orientation(:), z_image_orientation(:)]';
    transform_cor = [M*xyz(:,1:3)']';

    % figure(1), hold on
    % plot3(xyz(:,1),xyz(:,2),xyz(:,3),'r.')
    % axis equal

% Shift center
cor_shifter = mean(ProcessedStack.center_axis_diff);
CoronariesXYZ = transform_cor - cor_shifter;

    % figure(2), hold on
    % plot3(cor_shifted(:,1),cor_shifted(:,2),cor_shifted(:,3),'r.')
    % axis equal
    
% Make Transformation matrix and origin (3D MRI to Cardiac Coordinates)
    % Pull Pinpoints
    A  = ProcessedStack.ApexPoint;
    B  = ProcessedStack.BasePoint;
    Si = ProcessedStack.SeptumPoint;

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
