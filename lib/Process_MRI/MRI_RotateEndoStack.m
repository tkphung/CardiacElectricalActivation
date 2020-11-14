function [Cxyz,Pd,M,varargout] = MRI_RotateEndoStack(setstruct,varargin)
%MRI_RotateEndoStack: Rotates endocardial contours, rv insertions, and
%mid-septal points based on septum location
% [Cxyz,Pd,M,varargout] = MRI_RotateEndoStack(setstruct,varargin)
%   INPUT VARIABLE:
%       setstruct- structure of LV segmentation (short axis stack)
%       varargin- 'axial' if axial stack used instead of short axis
%   CALCULATIONS:
%       LEGACY function from Katie Parker & Sam Clarke that rotates endo
%       stack
%   OUTPUT VARIABLE:
%       Cxyz- cartesian coordinates of endo [x y z timeindx sliceindx]
%       Pd- RV insertion points [x y z]
%       M- Rotation Matrix from MRI_TransformEndoStackWithPinPoints()
%       varargout- heart rate
% 
% Thien-Khoi N. Phung (Based on PlotSAEndoStack.m- January 18, 2017)

%% Assume image set is Short Axis view (unless varargin indicated- axial)
axial_flag = false;
if ~isempty(varargin)
    if strcmp(varargin{1},'axial')
        axial_flag = true;
    end   
end

%% Cycle through contours and rotate
imaging_planes = length(setstruct);
slice_counter = 1;

Cxyz = [];
heartrate = [];

for idx=1:(imaging_planes)

    imgstruct = setstruct(idx);
    slice_labels = imgstruct.KeptSlices;
    hr = imgstruct.HeartRate;
    
    if axial_flag
        if (abs(abs(imgstruct.ImageOrientation(1))-1) < 1e-7) & ...
                (abs(abs(imgstruct.ImageOrientation(5))-1) < 1e-7)
            % Do nothing here
        else
            continue
        end
    end

    time_indicies = (1:imgstruct.TSize);
    
    if length(slice_labels) > 1
        for jdx=slice_labels
       
            [Xs,Ps,M] = MRI_TransformEndoStackWithPinPoints(imgstruct,jdx);
            Xd = cell2mat(Xs);
            Pd = cell2mat(Ps);

            if all(isnan(Xd(:)))
                continue
            end
            
            for kdx=1:imgstruct.TSize
                Xd = Xs{kdx};
                
                if all(isnan(Xd(:)))
                    continue
                end
                
                slice_indicies = ones(size(Xd,1),1) * jdx;
                Cxyz = [Cxyz; ...
                    [Xd, ones(size(Xd,1),1)*time_indicies(kdx), slice_indicies] ];
            end
            
            heartrate = [heartrate; [hr, slice_counter]];
            slice_counter = slice_counter + 1;
        end
    else
        [Xs,Ps,M] = MRI_TransformEndoStackWithPinPoints(imgstruct);
        Xd = cell2mat(Xs);
        Pd = cell2mat(Ps);        

        if all(isnan(Xd(:)))
            continue
        end

        for kdx=1:imgstruct.TSize
            Xd = Xs{kdx};
            
            if all(isnan(Xd(:)))
                continue
            end

            slice_indicies = ones(size(Xd,1),1) * slice_counter;
            Cxyz = [Cxyz; ...
                [Xd, ones(size(Xd,1),1)*time_indicies(kdx), slice_indicies] ];
        end

        heartrate = [heartrate; [hr, slice_counter]];
        slice_counter = slice_counter + 1; 
    end
   
end
varargout{1} = heartrate;
end