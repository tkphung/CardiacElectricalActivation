function MRI_FitContours(ProcessedStack,SApath,meshd,varargin)
%MRI_FitContours: fits contours to a bicubic surface mesh for the
%endocardium and epicardium saving two IPNODE files
% OUTPUT = MRI_FitContours(INTPUT)
%   INPUT VARIABLE:
%       ProcessedStack from MRI_ProcessStack.m
%       SApath being the save path for the IPNODE files
%       meshd string indicating mesh density
%       VARARGIN: 'epiname', 'endoname'
%   CALCULATION:
%       Adapted from SAC_3DScarMap_CRT102_Coronaries.m script
%   OUTPUT VARIABLE:
%       Saves two IPNODE files in SApath folder for Endo & Epi
%
% Thien-Khoi N. Phung (February 10, 2017)
% CHANGED March 23, 2017 TNP: Include Names for IPNODE files

% VARARGIN Processing
ENDONAME_flag = false;
EPINAME_flag  = false;

if ~isempty(varargin)
    for jz = 1:2:numel(varargin)
        switch varargin{jz}
            case 'endoname' % rename endo IPNODE file
                ENDONAME_flag = true;
                ENDONAME = varargin{jz+1};
            case 'epiname' % rename epi IPNODE file
                EPINAME_flag = true;
                EPINAME = varargin{jz+1};
            otherwise
                error('ERROR: Check your varargins.')
        end
    end
end

% Flags to write ipdata and ipnode files
% writeipdata = false; % Bicubic fit does not use ipdata
writeipnode = true;

% File paths and names for Endo and Epi fit data
% name_ipdata_endo = [SApath 'Endo.IPDATA']; %IPDATA NOT USED
% name_ipdata_epi = [SApath 'Epi.IPDATA'];
if ENDONAME_flag
    NameIPNODEendo = [SApath ENDONAME '.IPNODE'];
else
    NameIPNODEendo = [SApath 'Endo.IPNODE'];
end
if EPINAME_flag
    NameIPNODEepi = [SApath EPINAME '.IPNODE'];
else
    NameIPNODEepi = [SApath 'Epi.IPNODE'];
end

% Convert slices to cardiac coordinates
alldata_endo = ProcessedStack.alldata_endo;
alldata_epi = ProcessedStack.alldata_epi;
Endo = alldata_endo(:,1:3);
Epi = alldata_epi(:,1:3);
A = ProcessedStack.ApexPoint;
B = ProcessedStack.BasePoint;
S = ProcessedStack.SeptumPoint;

    % First basis vector, subtract Base from Apex and divide by magnitude
    C = A-B;
    e1 = C ./ norm(C);

    % Calculate Origin location
    origin = B + (1/3)*C;

    % Calculate focus length from norm(C)
    Focus  = (2*norm(C)/3)/ cosh(1);
    Nber_points = length(Endo); % Assume confocal endo & epi

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
    DataEndo = (Endo - repmat(origin,Nber_points,1))*Transform';
    DataEpi = (Epi - repmat(origin,Nber_points,1))*Transform';
    % A_transf = (A - origin)*Transform';
    % B_transf = (B - origin)*Transform';   
    % S_transf  = (S(1,:) - origin)*Transform';% Septal midpoint
    % S_transf1 = (S(2,:) - origin)*Transform';% Septal intersection 1
    % S_transf2 = (S(3,:) - origin)*Transform';% Septal intersection 2
    % DataScar = [Data alldata(:,4)];
    
    % X Y Z coordinates (cartesian) and Scar Transmurality
    DataEndo = [DataEndo alldata_endo(:,4)];
    DataEpi = [DataEpi alldata_epi(:,4)];
    
%% IPDATA NOT USED IN RENDERING- Commented out instead of deleted
% % This Section may be obsolete because fitting codes do not use IPDATA
% % information- it takes raw coordinate data
% % Add point numbering and non-uniform weighting for .IPDATA
% Points_indexes = (1:(Nber_points))';
% 
%     % For Long axis rotational tracing want to keep weighting factor of 1
%     Weight  = ones(Nber_points,1);
% 
%     Data_all_endo = [Points_indexes DataEndo(:,1) DataEndo(:,2) DataEndo(:,3) Weight Weight Weight];
%     Data_all_epi  = [Points_indexes DataEpi(:,1) DataEpi(:,2) DataEpi(:,3) Weight Weight Weight];
% 
% % Write transformed data points to .IPDATA file
% 
% if writeipdata
%      % ENDO
%     [fid,~] = fopen(name_ipdata_endo,'w');
%     fcomment  = strcat(name_ipdata_endo,' data from segmented MRI, focus=',num2str(Focus));
%     fprintf(fid,'%s\n',fcomment);
%     for jz=1:size(Data_all_endo,1)
%         fprintf(fid,'%4.0f %11.6f %11.6f %11.6f %8.4f %8.4f %8.4f\n',Data_all_endo(jz,:));
%     end
%     fclose(fid);
%     
%     % EPI
%     [fid,~] = fopen(name_ipdata_epi,'w');
%     fcomment  = strcat(name_ipdata_epi,' data from segmented MRI, focus=',num2str(Focus));
%     fprintf(fid,'%s\n',fcomment);
%     for jz=1:size(Data_all_epi,1)
%         fprintf(fid,'%4.0f %11.6f %11.6f %11.6f %8.4f %8.4f %8.4f\n',Data_all_epi(jz,:));
%     end
%     fclose(fid);
% end

%% Fit segmented data- OUTPUT IPNODE FILES
% options: '4x2' '4x4' '4x8'
[~,errorendo] = MRI_Bicubic_fit_data(DataEndo(:,1:3),Focus,writeipnode,NameIPNODEendo,meshd);
[~,errorepi] = MRI_Bicubic_fit_data(DataEpi(:,1:3),Focus,writeipnode,NameIPNODEepi,meshd);
