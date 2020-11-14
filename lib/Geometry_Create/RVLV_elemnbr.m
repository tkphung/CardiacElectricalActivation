function [E, varargout] = RVLV_elemnbr(LVLCR, RVLCR)
%RVLV_elemnbr: returns a cell array that includes the element
%index for each element's neighboring elements in the RV & LV EP model
% E = RVLV_elemnbr(LVLCR, RVLCR)
%   INPUT VARIABLES:
%       LVLCR: element dim of LV Long Circ Rad
%       RVLCR: element dim of RV Long Circ Rad
%       ASSUMPTION: RV insertion is 1 element wide on LV outside of 120 deg
%                   septal wall
%   OUTPUT VARIABLES:
%       E: cell indexed by elem number, contains indices of nbr elements
%       vargout: Matrix of unique home-neighbor pairs (no forward-backward
%                repeats)
%
% adapted from Wolfgang Schwanghart function "ixneighbors"
% Thien-Khoi N. Phung (September 26, 2016)

%% FIND ALL OF THE NEIGHBORS
% Use existing LV_elemnbrs to find LV neighbors
[~, LVHN] = LV_elemnbr(LVLCR(1), LVLCR(2), LVLCR(3));

% Find RV neighbors
Le = RVLCR(1); Ce = RVLCR(2); Re = RVLCR(3);
    % Create matrix with elem indices in the correct spatial position
    egrid = reshape(1:(Ce*Le*Re), Ce, Le, Re);

    % Pad with NaN on all borders
    eg = nan(Ce+2, Le+2, Re+ 2);
    eg(2:end-1, 2:end-1, 2:end-1) = egrid;

    % Logical matrix (for searching neighbors)
    L = ~isnan(eg);

    nbr = zeros((Ce*Le*Re), 26);
    % Shift logical matrix across the neighbors                        C   L   R
    nbr(:,1) = eg(L(    [2:end 1],     [2:end 1],     [2:end 1])); % [+1, +1, +1]
    nbr(:,2) = eg(L(            :,     [2:end 1],     [2:end 1])); % [ 0, +1, +1]
    nbr(:,3) = eg(L([end 1:end-1],     [2:end 1],     [2:end 1])); % [-1, +1, +1]

    nbr(:,4) = eg(L(    [2:end 1],             :,     [2:end 1])); % [+1,  0, +1]
    nbr(:,5) = eg(L(            :,             :,     [2:end 1])); % [ 0,  0, +1]
    nbr(:,6) = eg(L([end 1:end-1],             :,     [2:end 1])); % [-1,  0, +1]

    nbr(:,7) = eg(L(    [2:end 1], [end 1:end-1],     [2:end 1])); % [+1, -1, +1]
    nbr(:,8) = eg(L(            :, [end 1:end-1],     [2:end 1])); % [ 0, -1, +1]
    nbr(:,9) = eg(L([end 1:end-1], [end 1:end-1],     [2:end 1])); % [-1, -1, +1]

    nbr(:,10)= eg(L(    [2:end 1],     [2:end 1],             :)); % [+1, +1,  0]
    nbr(:,11)= eg(L(            :,     [2:end 1],             :)); % [ 0, +1,  0]
    nbr(:,12)= eg(L([end 1:end-1],     [2:end 1],             :)); % [-1, +1,  0]

    nbr(:,13)= eg(L(    [2:end 1],             :,             :)); % [+1,  0,  0]
    nbr(:,14)= eg(L([end 1:end-1],             :,             :)); % [-1,  0,  0]

    nbr(:,15)= eg(L(    [2:end 1], [end 1:end-1],             :)); % [+1, -1,  0]
    nbr(:,16)= eg(L(            :, [end 1:end-1],             :)); % [ 0, -1,  0]
    nbr(:,17)= eg(L([end 1:end-1], [end 1:end-1],             :)); % [-1, -1,  0]

    nbr(:,18)= eg(L(    [2:end 1],     [2:end 1], [end 1:end-1])); % [+1, +1, -1]
    nbr(:,19)= eg(L(            :,     [2:end 1], [end 1:end-1])); % [ 0, +1, -1]
    nbr(:,20)= eg(L([end 1:end-1],     [2:end 1], [end 1:end-1])); % [-1, +1, -1]

    nbr(:,21)= eg(L(    [2:end 1],             :, [end 1:end-1])); % [+1,  0, -1]
    nbr(:,22)= eg(L(            :,             :, [end 1:end-1])); % [ 0,  0, -1]
    nbr(:,23)= eg(L([end 1:end-1],             :, [end 1:end-1])); % [-1,  0, -1]

    nbr(:,24)= eg(L(    [2:end 1], [end 1:end-1], [end 1:end-1])); % [+1, -1, -1]
    nbr(:,25)= eg(L(            :, [end 1:end-1], [end 1:end-1])); % [ 0, -1, -1]
    nbr(:,26)= eg(L([end 1:end-1], [end 1:end-1], [end 1:end-1])); % [-1, -1, -1]

    % Pair neighbors to house
    H = repmat(eg(L(:)),26,1); % Houses
    N = nbr(:); % Neighbors

    % Remove NaN neighbors
    nope = isnan(N);
    H(nope) = [];
    N(nope) = [];

    % Home Neighbor Pairs
    RVHN = [H N];
    RVHN = sortrows(RVHN(RVHN(:,1)<RVHN(:,2),:));

    % Index RV & LV elements together
    RVHN = RVHN + max(max(LVHN));

% Combine LV and RV neighbors
RVLVHN = [LVHN; RVHN];

% Add RV-LV Boundary neighbors
    % LV Apex, RV Apex, RV +60 insert, RV -60 insert
    % LV APEX (all circ elements for each radial layer are connected)
    for jz = 1:LVLCR(3) % combing through radial layers
        % LV elements around apex
        LVapex = [(LVLCR(2)*(LVLCR(1)-1) + 1):(LVLCR(2)*LVLCR(1))] + ...
                 (LVLCR(2)*LVLCR(1)*(jz-1));
        for by = 1:LVLCR(2)
            LVapexhome = LVapex(by);
            LVapexnbrs = LVapex((by+1):end)';
            % Appending radial neighbors
            RVLVHN = [RVLVHN;...
                      repmat(LVapexhome,length(LVapexnbrs),1) LVapexnbrs];
        end
    end

    % RV Apex (connect RV apex endo to LV apex epi)
    RVapexendo = (RVLCR(2)*(RVLCR(1)-1) + 1):(RVLCR(2)*RVLCR(1));
    LVseptepi  = LVLCR(1)*LVLCR(2)*(LVLCR(3)-1) +...
                 LVLCR(2)*(LVLCR(1)-1) +...
                 [(fix(LVLCR(2)/12)-1):(LVLCR(2)*5/12)];
    for jz = 1:length(LVseptepi)
        RVLVHN = [RVLVHN;...
                  repmat(LVseptepi(jz),length(RVapexendo),1) RVapexendo'];
    end

    % RV +60 insert
    RVplus_base = LVLCR(1)*LVLCR(2)*(LVLCR(3)-1) + (fix(LVLCR(2)/12) - 1);
    RVplus = RVplus_base:LVLCR(2):LVLCR(1)*LVLCR(2)*LVLCR(3);
    for jz = 1:length(RVplus)
        RVinsert = max(max(LVHN)) +...
                  [1:(RVLCR(1)*RVLCR(2)):(RVLCR(1)*RVLCR(2)*RVLCR(3))] +...
                  RVLCR(2)*(jz-1);
        RVLVHN = [RVLVHN;...
                  repmat(RVplus(jz),RVLCR(3),1) RVinsert'];
    end
    
    % RV -60 insert
    RVmin_base = LVLCR(1)*LVLCR(2)*(LVLCR(3)-1) + fix(LVLCR(2)*5/12);
    RVmin = RVmin_base:LVLCR(2):LVLCR(1)*LVLCR(2)*LVLCR(3);
    for jz = 1:length(RVmin)
        RVinsert = max(max(LVHN)) +...
           [RVLCR(2):(RVLCR(1)*RVLCR(2)):(RVLCR(1)*RVLCR(2)*RVLCR(3))] +...
           RVLCR(2)*(jz-1);
        RVLVHN = [RVLVHN;...
                  repmat(RVmin(jz),RVLCR(3),1) RVinsert'];
    end

% Sort Neighbors
HN = unique(sortrows(RVLVHN(RVLVHN(:,1)<RVLVHN(:,2),:)),'rows');
varargout{1} = HN;
    
for el = 1:max(max(RVLVHN))
    E{el} = HN(HN(:,1)==el, 2);
end

    