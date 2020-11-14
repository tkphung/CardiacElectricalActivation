function [pinpts3d,M] = MRI_TransformLAWithPinPoints(imgstruct,varargin)
% Name changed by TK Phung (January 18, 2017)
% This function rotates all of the segmentations into the 3D MRI
% coordinates based on the imgstruct information (from Segment) 
% varargin indicates the slice number to rotate (otherwise it rotates them
% all)
% Previous function name: TransformLAWithPinPoints.m
% Edited by TNP May 9, 2017 to add varargin:
% 'time' Time Frame
timeID_flag = false;
if ~isempty(varargin)
    for jz = 1:2:numel(varargin)
        switch varargin{jz}
			case 'time'
                timeID_flag = true;
				timeID = varargin{jz+1};
            otherwise
                error('ERROR: Check your varargins.')
        end
    end
end

x_resolution = imgstruct.ResolutionX;
y_resolution = imgstruct.ResolutionY;

image_position = imgstruct.ImagePosition;

x_image_orientation = imgstruct.ImageOrientation(4:6);
y_image_orientation = imgstruct.ImageOrientation(1:3);
z_image_orientation = cross(y_image_orientation,x_image_orientation);

slice_thickness = imgstruct.SliceThickness;
slice_gap = imgstruct.SliceGap;
z_offset = 0;

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

if timeID_flag % IF TIME IS SPECIFIED - edit TNP 05-09-17 for CINE MRI processing
    PPslice = 1; % ASSUMING LA only has ONE slice
    % timeID was set above with the timeID_flag
else
    [timeID,PPslice] = find(~cellfun(@isempty,imgstruct.EndoPinX));
end

z_offset_PP = (PPslice - 1)';
% ADDEDED 01-15-18 by TNP 
pinpts_round = [];
for db = 1:numel(PPslice)
    x_pinpts = imgstruct.EndoPinX{timeID,PPslice(db)};
    y_pinpts = imgstruct.EndoPinY{timeID,PPslice(db)};
    pinpts_round = [pinpts_round; round(x_pinpts),round(y_pinpts)];
end
% end edit TNP
pinpts = {pinpts_round};
z_pix = -z_offset_PP .* ones(size(pinpts{1},1),1);
pinpts3d{1} = [pinpts{1}, z_pix];

PP = ( M * [pinpts3d{1}, ones(size(pinpts3d{1},1),1)]')';
PP = PP(:,1:3);
pinpts3d{1} = PP;
pinpts3d = cell2mat(pinpts3d);
end
