function [xyz_pts,pinpts3d,M] = MRI_TransformEndoStackWithPinPoints(imgstruct,varargin)
% Name changed by TK Phung (January 18, 2017)
% This function rotates all of the segmentations into the 3D MRI
% coordinates based on the imgstruct information (from Segment) 
% varargin indicates the slice number to rotate (otherwise it rotates them
% all)
% Previous function name: TransformSAEndoStackWithPinPoints.m

if ~isempty(varargin)
    slicenumber = varargin{1};
    x_pix = imgstruct.EndoX(:,:,slicenumber);
    y_pix = imgstruct.EndoY(:,:,slicenumber);
else
    slicenumber = 1;
    x_pix = imgstruct.EndoX;
    y_pix = imgstruct.EndoY;
end

x_pix_round = round(x_pix);
y_pix_round = round(y_pix);

perim_length = zeros(size(x_pix_round,2),1);
xy_pts = cell(size(x_pix_round,2),1);

for idx=1:size(x_pix_round,2)
    
    if any(isnan(x_pix_round(:,idx)))
        xy_pts{idx} = [NaN NaN NaN];
        continue
    end
    
    xy_pix_round = [x_pix_round(:,idx), y_pix_round(:,idx)];
    perim_length(idx) = size(unique(xy_pix_round,'rows'),1);
    
    perim_pts = linspace(0,1,size(x_pix,1)+1);
    interp_perim_pts = linspace(0,1,perim_length(idx)+1);
    perim_xy_pts = [[x_pix(:,idx), y_pix(:,idx)]; [x_pix(1,idx), y_pix(1,idx)]];
    interp_xy_pts = interp1(perim_pts,perim_xy_pts,interp_perim_pts,'spline');
    
    xy_pts{idx} = interp_xy_pts(1:end-1,:);
    
end


x_resolution = imgstruct.ResolutionX;
y_resolution = imgstruct.ResolutionY;

image_position = imgstruct.ImagePosition;

x_image_orientation = imgstruct.ImageOrientation(4:6);
y_image_orientation = imgstruct.ImageOrientation(1:3);
z_image_orientation = cross(y_image_orientation,x_image_orientation);

slice_thickness = imgstruct.SliceThickness;
slice_gap = imgstruct.SliceGap;
z_offset = (slicenumber - 1);

xyz_pts = xy_pts;
for idx=1:size(x_pix_round,2)
    z_pix = -z_offset * ones(perim_length(idx),1);
    if ~any(isnan(xy_pts{idx}(:)))
        xyz_pts{idx} = [xy_pts{idx}, z_pix];
    end
end

To = eye(4);
To(1:3,4) = [-1 -1 0]';

S = eye(4);
S(1,1) = x_resolution;
S(2,2) = y_resolution;
S(3,3) = slice_thickness+slice_gap; 

R = eye(4);
R(1:3,1:3) = [x_image_orientation(:), y_image_orientation(:), z_image_orientation(:)];

Tipp = eye(4);
Tipp(1:3,4) = image_position';

M = Tipp * R * S * To;

for idx=1:size(x_pix_round,2)
    if ~any(isnan(xyz_pts{idx}(:)))
        try
        X = ( M * [xyz_pts{idx}, ones(size(xyz_pts{idx},1),1)]')';
        catch ex
            keyboard
        end
        X = X(:,1:3);
        xyz_pts{idx} = X;
    end
end

[timeID,PPSlice] = find(~cellfun(@isempty,imgstruct.EndoPinX));

% Added by TNP May 9, 2017- Specifically for processing CINE that has
% pinpts on all time frames- choose frame that has 3 pinpoints (from
% ProcessedStack), OG has 2 pinpoints- 3rd comes from the midpoint
% selection in MRI_ProcessStack
for jz = 1:numel(timeID)
    if size(imgstruct.EndoPinX{timeID(jz),PPSlice(jz)},1)==3
        timeID = timeID(jz);
        PPSlice = PPSlice(jz);
        break
    end
end

z_offset_PP = (PPSlice - 1);
x_pinpts = imgstruct.EndoPinX{timeID,PPSlice};
y_pinpts = imgstruct.EndoPinY{timeID,PPSlice};
pinpts_round = [round(x_pinpts),round(y_pinpts)];
pinpts = {pinpts_round};
z_pix = -z_offset_PP * ones(length(x_pinpts),1);
pinpts3d{1} = [pinpts{1}, z_pix];

PP = ( M * [pinpts3d{1}, ones(size(pinpts3d{1},1),1)]')';
PP = PP(:,1:3);
pinpts3d{1} = PP;

end
