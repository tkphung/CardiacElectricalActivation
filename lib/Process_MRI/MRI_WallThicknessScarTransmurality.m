function ProcessedStack = MRI_WallThicknessScarTransmurality(ProcessedStack,scar_flag)
%MRI_ScarTransmurality: Takes formatted ProcessedStack from
%MRI_ProcessStack that should include scar tracings and calculates teh
%percent wall thickness
% ProcessedStack = MRI_ScarTransmurality(ProcessedStack)
%   INPUT VARIABLE:
%       ProcessedStack- ProcessedStack from MRI_ProcessStack([],[],'scar')
%       scar_flag- is there scar segmented?
%   CALCULATION:
%       LEGACY function- directly from SAC_3DScarMap_.m series
%   OUTPUT VARIABLE:
%       ProcessedStack- with WallThicknessPolar, TransmuralityPolar, and
%       CircExtendDeg fields
%
% Adapted by Thien-Khoi N. Phung from SAC_3DScarMap_CRT102_Coronaries.m
% code (February 7, 2017)

angles = linspace(-pi,pi,50);
for i = ProcessedStack.Slices;
    endo = cell2mat(ProcessedStack.Endo(1,i));
    epi = cell2mat(ProcessedStack.Epi(1,i)); 
    if scar_flag
        scar = cell2mat(ProcessedStack.Scar(1,i));
    end
    
    % Convert Endo, Epi, Scar Contours to Polar Coordinates
    endo_shift = [endo(:,1)-mean(epi(:,1)),endo(:,2)-mean(epi(:,2)),endo(:,3)];
    epi_shift = [epi(:,1)-mean(epi(:,1)),epi(:,2)-mean(epi(:,2)),epi(:,3)];
    if scar_flag
        scar_shift = [scar(:,1)-mean(epi(:,1)),scar(:,2)-mean(epi(:,2)),scar(:,3)];
    end
    
    % Convert Endo, Epi, Scar Contours to Polar Coordinates
    [theta_endo,rho_endo] = cart2pol(endo_shift(:,1),endo_shift(:,2));
    [theta_epi,rho_epi] = cart2pol(epi_shift(:,1),epi_shift(:,2));
    if scar_flag
        [theta_scar,rho_scar] = cart2pol(scar_shift(:,1),scar_shift(:,2));
    end
    
    % figure
    % polar([0 2*pi],[0 60])
    % hold on
    % polar(theta_endo,rho_endo)
    % hold on
    % polar(theta_epi,rho_epi)
    % hold on
    % polar(theta_scar,rho_scar,'m.')

    % Bin contours by angle
    for j = 1:(length(angles)-1)
        range = [angles(j),angles(j+1)];
        endo_idx = find(range(1) <= theta_endo & theta_endo < range(2));
        epi_idx = find(range(1) <= theta_epi & theta_epi < range(2));
        if scar_flag
            scar_idx = find(range(1) <= theta_scar & theta_scar < range(2));
        end
        if isempty(endo_idx)
            endo_binned(j,1:4) = [range(1),mean(range),range(2),NaN];
        else
            endo_binned(j,1:4) = [range(1),mean(range),range(2),mean(rho_endo(endo_idx))];
        end
        if isempty(epi_idx)
            epi_binned(j,1:4) = [range(1),mean(range),range(2),NaN];
        else
            epi_binned(j,1:4) = [range(1),mean(range),range(2),mean(rho_epi(epi_idx))];
        end
        if scar_flag
            if isempty(scar_idx)
                scar_binned(j,1:5) = [range(1),mean(range),range(2),NaN,NaN];
            else
                scar_binned(j,1:5) = [range(1),mean(range),range(2),min(rho_scar(scar_idx)),max(rho_scar(scar_idx))];
            end
        end
    end
    
    epi_binned(isnan(endo_binned(:,4)),:) = [];
    if scar_flag
        scar_binned(isnan(endo_binned(:,4)),:) = [];
    end
    endo_binned(isnan(endo_binned(:,4)),:) = [];
    
    ProcessedStack.EndoPolar{i} = endo_binned;
    ProcessedStack.EpiPolar{i} = epi_binned;
    if scar_flag
        ProcessedStack.ScarPolar{i} = scar_binned;
    end
    
    %Calculate wall thickness + transmurality
    wall_thickness = zeros(length(endo_binned),2); % CHANGED BY TNP 04-03-2017 to preallocate
    transmurality = wall_thickness;
    for j = 1:length(endo_binned);
        wall_thickness(j,1:2) = [endo_binned(j,2),epi_binned(j,4)-endo_binned(j,4)];
    end
    
    if scar_flag
        for j = 1:length(endo_binned);
            transmurality(j,1:2) = [scar_binned(j,2),(scar_binned(j,5)-scar_binned(j,4))/wall_thickness(j,2)];
        end
    else
        for j = 1:length(endo_binned);
            transmurality(j,1:2) = [NaN NaN];
        end
    end
    transmurality(isnan(transmurality)) = 0;
    
    ProcessedStack.WallThicknessPolar{i} = wall_thickness;
    ProcessedStack.TransmuralityPolar{i} = transmurality;

    %Calculate circumferential extent of each slice
    if scar_flag
        scar_angles = find(transmurality(:,2) ~= 0);
        if isempty(scar_angles)
            circ_extent(i) = 0;
        else
            % Multiply the number of bins by the size of each bin, in rad
            circ_extent(i) = length(scar_angles)*((2*pi)/50);
        end
        ProcessedStack.CircExtentDeg = circ_extent*(180/pi);
    end
    % polar(endo_binned(:,2),endo_binned(:,4),'c.')
    % hold on
    % polar(epi_binned(:,2),epi_binned(:,4),'g.')   
end
