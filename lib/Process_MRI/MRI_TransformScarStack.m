function [xyz_pts,M] = MRI_TransformScarStack(imgstruct,varargin)
% Name changed by TK Phung (January 18, 2017)
% This function rotates all of the segmentations into the 3D MRI
% coordinates based on the imgstruct information (from Segment) 
% varargin indicates the slice number to rotate (otherwise it rotates them
% all)
% Previous function name: TransformSAScarStack.m

if ~isempty(varargin)
    slicenumber = varargin{1};
    x_pix = imgstruct.Scar.ScarX(:,1,slicenumber);
    y_pix = imgstruct.Scar.ScarY(:,1,slicenumber);
else
    slicenumber = 1;
    x_pix = imgstruct.Scar.ScarX(:,1,slicenumber);
    y_pix = imgstruct.Scar.ScarY(:,1,slicenumber);
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


if isnan(sum(sum(imgstruct.Scar.ScarX(:,:,slicenumber)))) < 1
    
    x_pix_round = round(x_pix);
    y_pix_round = round(y_pix);

    perim_length = zeros(size(x_pix_round,2),1);
    xy_pts = cell(size(x_pix_round,2),1);

    for idx=1:size(x_pix_round,2)

    %     if any(isnan(x_pix_round(:,idx)))
    %         xy_pts{idx} = [NaN NaN NaN];
    %         continue
    %     end

        xy_pix_round = [x_pix_round(:,idx), y_pix_round(:,idx)];
        perim_length(idx) = size(xy_pix_round,1);

        perim_pts = linspace(0,1,size(x_pix,1)+1);
        interp_perim_pts = linspace(0,1,perim_length(idx)+1);
        perim_xy_pts = [[x_pix(:,idx), y_pix(:,idx)]; [x_pix(1,idx), y_pix(1,idx)]];
        interp_xy_pts = interp1(perim_pts,perim_xy_pts,interp_perim_pts,'spline');

        xy_pts{idx} = interp_xy_pts(1:end-1,:);

    end

    xyz_pts = xy_pts;
    for idx=1:size(x_pix_round,2)
        z_pix = -z_offset * ones(perim_length(idx),1);
    %     if ~any(isnan(xy_pts{idx}(:)))
            xyz_pts{idx} = [xy_pts{idx}, z_pix];
    %     end
    end

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
else
    
    xyz_pts{1} = [NaN,NaN,NaN];
end

end
