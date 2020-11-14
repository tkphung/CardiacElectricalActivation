function [E,varargout] = LV_elemnbr(Le, Ce, Re)
%LV_elemnbr: returns a cell array that includes the element
%index for each element's neighboring elements (up to 26)
% [E] = LV_elemnbr(Le, Ce, Re)
%   INPUT VARIABLES:
%       Le: number of long. elements excluding apex cap elem.
%       Ce: number of circ. elements around each ring
%       Re: number of rad. elements through wall
%   OUTPUT VARIABLES:
%       E: cell indexed by elem number, contains indices of nbr elements
%       vargout: Matrix of unique home-neighbor pairs
%
% adapted from Wolfgang Schwanghart function "ixneighbors"
% Thien-Khoi N. Phung (April 27, 2016)
%                     (June 6, 2016- added vargout)

% Create matrix with elem indices in the correct spatial position
egrid = reshape(1:(Ce*Le*Re), Ce, Le, Re);

% Pad with NaN on all borders
eg = nan(Ce+2, Le+2, Re+ 2);
eg(2:end-1, 2:end-1, 2:end-1) = egrid;

% Logical matrix (for searching neighbors)
L = ~isnan(eg);

% Pad circumferentially
eg(1,:,:)   = eg(end-1,:,:);
eg(end,:,:) = eg(2,:,:);

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

% Sort into structure
for el = 1:(Ce*Le*Re)
    E{el} = N(H==el);
end

% Variable output arguments
% Unique Home Neighbor Pairs
HN = [H N];
HN = sortrows(HN(HN(:,1)<HN(:,2),:));
varargout{1} = HN;


