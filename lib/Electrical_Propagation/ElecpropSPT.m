function [litTIME,litHOO,varargout] = ElecpropSPT(HEXC,HN,QFCR,CV,STARLIT)
%ElecpropSPT: Calculates onset times for electrical activation for each
%element given the fiber direction, conduction CVcities, and seed
%element(s)
%   [litTIME,litHOO] = ElecpropSPT(HEXC,HN,QFCR,CV,STARLIT)
%   INPUTS:
%       HEXC- nodes (FE-element centroids)
%       HN- edges (home-neighbor pairs)- *unique pairs
%       QFCR- local coordinate system
%       CV: conduction CVcity for elements (F, CF, R)
%       STARLIT: element indicies initially activated (can be matrix)
%                if size(STARLIT)>[1 1], col 1 is index, col 2 is time
%   CALCULATIONS: DIJKSTRA's ALGORITHM
%       1. Calculate electrical propagation dt "space"
%		2. Algorithm for EP assignment timing calculation
%   OUTPUTS:
%       litTIME: LIT times for all elements
%       litHOO: who lit which node
%       varargout: G- the tree
%
% Adapted from RVLV_elecprop
% Thien-Khoi Phung (October 7, 2016)
% Updated (October 12, 2016- TNP) vectorized algorithm
% Changed (October 17, 2016- TNP) forwards and backwards edge differences
% Changed to RVLV_elecpropSPT using MATLAB SPT function (October 24,
% 2016-TNP)
% Changes (December 13, 2016- TNP) allows for input of multiple seed
% elements during initiation of active contraction
% ADAPTED from RVLV_elecpropSPT for a more generalized function to use with
% MRI-derived electrical models (TNP- March 31, 2017)
% Changed (May 14, 2017- Happy Mother's Day) add varagout for the GRAPH, G

%%% 1. CALCULATE EP dt SPACE
    % FORWARD SPACE
    hm = HN(:,1);
    nbr = HN(:,2);

    % Vector for home to neighbor
    rv = HEXC(nbr,:) - HEXC(hm,:);
    R = (rv(:,1).^2 + rv(:,2).^2 + rv(:,3).^2).^(1/2);
    rv = rv./repmat(R,1,3);

    % Fiber
    fib = squeeze(QFCR(:,1,hm))';
    c1 = CV(hm,1); % c1

    % Cross-fiber
    cro = squeeze(QFCR(:,2,hm))';
    c2 = CV(hm,2); % c2

    % Radial
    rad = squeeze(QFCR(:,3,hm))';
    c3 = CV(hm,3); % c3

    % Calculate dt
    dtfwd = R.*((dot(rv,fib,2)./c1).^2 + (dot(rv,cro,2)./c2).^2 ...
                + (dot(rv,rad,2)./c3).^2).^(1/2);
            
    % BACKWARD SPACE
    hm = HN(:,2);
    nbr = HN(:,1);

    % Vector for home to neighbor
    rv = HEXC(nbr,:) - HEXC(hm,:);
    R = (rv(:,1).^2 + rv(:,2).^2 + rv(:,3).^2).^(1/2);
    rv = rv./repmat(R,1,3);

    % Fiber
    fib = squeeze(QFCR(:,1,hm))';
    c1 = CV(hm,1); % c1

    % Cross-fiber
    cro = squeeze(QFCR(:,2,hm))';
    c2 = CV(hm,2); % c2

    % Radial
    rad = squeeze(QFCR(:,3,hm))';
    c3 = CV(hm,3); % c3

    % Calculate dt
    dtbkwd = R.*((dot(rv,fib,2)./c1).^2 + (dot(rv,cro,2)./c2).^2 ...
                + (dot(rv,rad,2)./c3).^2).^(1/2);
    
    % CHOOSE AVERAGE TO STORE OF FORWARD-BACKWARD EDGES
    dt = (dtbkwd+dtfwd)./2;
    
%%% 2. ALGORITHM FOR EP
    % Create GRAPH of model
    % HN should be unique edges
    edges = length(dt);
    G = graph(HN(1:edges, 1), HN(1:edges,2), dt(1:edges));
    varargout{1} = G;
    
    % Solve SPT problem
    if numel(STARLIT) == 1
        [litHOO, litTIME] = shortestpathtree(G,STARLIT,'Method','positive',...
                            'OutputForm','vector');
    else
        [litHOO, litTIME] = shortestpathtree(G,STARLIT(1,1),'Method','positive',...
                            'OutputForm','vector');
        % Add time delay @STARLIT
        litTIME = litTIME + STARLIT(1,2); 
        for jz = 2:size(STARLIT,1)
            [HOO, TIME] = shortestpathtree(G,STARLIT(jz,1),'Method','positive',...
                            'OutputForm','vector');
            TIME = TIME + STARLIT(jz,2);
            % Replace slower times & neighbors from other parts of SPTs
            litHOO(TIME<litTIME) = HOO(TIME<litTIME);
            litTIME(TIME<litTIME) = TIME(TIME<litTIME);
        end
    end
    
    litHOO(STARLIT(:,1)) = STARLIT(:,1);
        
end