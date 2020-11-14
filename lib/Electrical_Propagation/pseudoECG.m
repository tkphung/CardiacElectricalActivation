function [pV] = pseudoECG(MODEL)
%[pV,varargout] = pseudoECG(MODEL) calculates the pseudoECG based on
%(Roberts et al 1979; Miller et al 1978)
% INPUTS:
%         MODEL must inlcude: EPnodeIDX (elements used in model)
%                             litTIME (element activation time)
%                             Qfcr (transformation matrix)
%                             VELO (element velocities)
%                             HN
%                             HEXc
%                             HEXv
%                             leads.xyz
% STEPS: Calculate pseudo-Voltage at leads
%         1. Create TMP
%         2. Calculate intracellular conductivity matrix
%         3. Calculate spatial gradient of TMP = Electric Field (mV/mm = V/m)
%         4. Calculate element volumes
%         5. Calculate pVoltage at leads
% OUTPUTS:
%         pV is the lead voltages
%
% Thien-Khoi N. Phung (May 31, 2019)

% 1. Create TMP
        % Elements Used
        ele = MODEL.EPnodeIDX;

        % Assign transmembrane potential through time (mV)
        litTIME = MODEL.litTIME;

        % Solve over n equally spaced time points
        nt = 200;
        time = linspace(min(litTIME),max(litTIME),nt);
        state = repmat(time,numel(ele),1);

        % Determine of element (rows) are resting (1) or depolarized (0)
        state = state<repmat(litTIME(:),1,nt);

        % Assign TMP at each time point
        % -90 mV if element is resting
        %   0 mV if element is depolarized
        tmp = state.*(-90);

% 2. Intracellular Conductivity Tensor (ohm meter)^-1
        % Body Coordinate Matrix (columns are vectors)
        Qfcr = MODEL.Qfcr(:,:,ele);

        % Create Diagonal Conductivity Tensors for each element
        % Set the resistance constants (ohm meter)
        rit = 38;  % intracellular transverse
        rot = 7.5; % extra         transverse
        rof = 4.5; % extra         fiber

        % Solve for resistance intracellular fiber
        % Roberts Circ Res 1979
        vt2vf = MODEL.VELO(ele,2)./MODEL.VELO(ele,1);
        rif = vt2vf.^2 * (rit+rot) - rof;

        Dlocal = zeros(3,3,numel(ele));
        Dlocal(1,1,:) = 1/rif; % fiber [1/(ohm meter)]
        Dlocal(2,2,:) = 1/rit; % transverse
        Dlocal(3,3,:) = 1/rit; % transverse

        % Convert to body coordinates
        Dbody = Dlocal;
        for jz = 1:numel(ele)
            Dbody(:,:,jz) = Qfcr(:,:,jz)*Dlocal(:,:,jz)*Qfcr(:,:,jz)';
        end

% 3. Calculate spatial gradient of TMP = Electric Field (mV/mm = V/m)
        % Define unique local neighbors for each element
        HN = MODEL.HN;
        % Renumber from 1:numel(ele) instead of skipping element indices
        [~,HN] = ismember(HN,ele);
        % Repeat [HN; NH] definition (HNNH is redundant)
        HNNH = [HN; HN(:,[2 1])];

        % Sort rows to make life easier
        HNNH = sortrows(HNNH);

        % Approximating spatial gradient using HNNH
        % Element centers
        HEXc = MODEL.HEXc(ele,:);
        % Unit vector from Home (col 1) to Neighbor (col 2)
        uvec = HEXc(HNNH(:,2),:) - HEXc(HNNH(:,1),:);
        duv = sqrt(sum(uvec.^2,2)); % distance between home and neighbor
        uv = uvec./repmat(duv,1,3);

        % At each time point calculate the spatial gradient between HN
        % Calculate "forward" difference in TMPs for all HNs
        %           (neighbor - home) for all matched time points
        %           divide by the distance between the h-n
        gradHNNH = (tmp(HNNH(:,2),:) - tmp(HNNH(:,1),:))./repmat(duv,1,nt);

        % Multiply the gradHNNH by the unit vector (V/m)
        gradX = gradHNNH.*repmat(uv(:,1),1,size(tmp,2));
        gradY = gradHNNH.*repmat(uv(:,2),1,size(tmp,2));
        gradZ = gradHNNH.*repmat(uv(:,3),1,size(tmp,2));

        % Number of neighbors for each element
        [n,~] = histcounts(HNNH(:,1),numel(ele));

        % Split matrices into cells (each cell is an element)
        gradXc = mat2cell(gradX,n,size(gradX,2));
        gradYc = mat2cell(gradY,n,size(gradX,2));
        gradZc = mat2cell(gradZ,n,size(gradX,2));

        % Add each elements' spatial gradient components
        gradXc = cellfun(@(x) sum(x,1),gradXc,'UniformOutput',0);
        gradYc = cellfun(@(x) sum(x,1),gradYc,'UniformOutput',0);
        gradZc = cellfun(@(x) sum(x,1),gradZc,'UniformOutput',0);

        % Cell 2 Matrix: each row is an element, each col is a time point
        gradXe = cell2mat(gradXc);
        gradYe = cell2mat(gradYc);
        gradZe = cell2mat(gradZc);

        % Multiply intracellular conductivity into the gradient
        gX = repmat(squeeze(Dbody(1,1,:)),1,nt).*gradXe + ...
             repmat(squeeze(Dbody(1,2,:)),1,nt).*gradYe + ... 
             repmat(squeeze(Dbody(1,3,:)),1,nt).*gradZe;
        gY = repmat(squeeze(Dbody(2,1,:)),1,nt).*gradXe + ...
             repmat(squeeze(Dbody(2,2,:)),1,nt).*gradYe + ...
             repmat(squeeze(Dbody(2,3,:)),1,nt).*gradZe;
        gZ = repmat(squeeze(Dbody(3,1,:)),1,nt).*gradXe + ...
             repmat(squeeze(Dbody(3,2,:)),1,nt).*gradYe + ...
             repmat(squeeze(Dbody(3,3,:)),1,nt).*gradZe;


% 4. Calculate element volumes
        HEXv = MODEL.HEXv(ele,:); % cubic mm

% 5. Calculate pVoltage at leads
% Lead positions
lead = MODEL.leads.xyz;

pV = zeros(size(lead,1),nt);
for ld = 1:size(lead,1)  

    % Distance between the leads and element centers (dipole origins)
    r = sqrt(sum((HEXc - lead(ld,:)).^2,2));

    % [(x-a)/r^3 (y-b)/r^3 (z-c)/r^3]
    gr = (HEXc - lead(ld,:))./(repmat(r,1,3).^3);

    % multiply gr by the (gradient of Vm * intracellular conductivity)
    pVx = repmat(gr(:,1),1,nt).*gX;
    pVy = repmat(gr(:,2),1,nt).*gY;
    pVz = repmat(gr(:,3),1,nt).*gZ;

    % integrate pVx, pVy, pVz
    pV(ld,:) = sum((pVx + pVy + pVz).*repmat(HEXv,1,nt),1); 
end
