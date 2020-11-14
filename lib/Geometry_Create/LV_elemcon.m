function [HEX, PENT] = LV_elemcon(Le, Ce, Re)
%LV_HEXelemcon: creates element connnectivity matrix for given number of
%elements in the longitudinal (excluding cap), circumferential, and radial
%directions
% EC = LV_HEXelemcon(Le, Ce, Re)
%   INPUT VARIABLES:
%       Le: number of long. elements excluding apex cap elem.
%       Ce: number of circ. elements around each ring
%       Re: number of rad. elements through wall
%   OUTPUT VARIABLE:
%       HEX: matrix with each row being one elements connectivity with
%       respect to nodes defined by LV_nodenumb
%		PENT: apex cap elements
% 
% Thien-Khoi N. Phung (April 19, 2016)

NPL = (Le+1)*Ce + 1; % number of nodes per layer
EPL = Le*Ce;         % number of elements per layer (excludes apex cap)

% Node connectivity for HEX elements
HEX = zeros(EPL*Re, 8);
for L = 1:Re
	for n = 1:EPL
        Nn = (EPL*(L-1) + 1):(EPL*L); % elem number
        Sn = (NPL)*(L-1)+n; %seed node
		if mod(n,Ce)~= 0
			HEX(Nn(n),:) = [Sn   , Sn+1   , Sn+NPL+1   , Sn+NPL,...
                            Sn+Ce, Sn+Ce+1, Sn+NPL+Ce+1, Sn+NPL+Ce];
		else
			HEX(Nn(n),:) = [Sn   , Sn+1-Ce, Sn+NPL+1-Ce, Sn+NPL,...
                            Sn+Ce, Sn+1   , Sn+NPL+1   , Sn+NPL+Ce];
		end
	end
end

% Node connectivity for PENT elements in apex cap
PENT = zeros(Ce*Re, 6);
for L = 1:Re
	ENum = (Ce*(L-1) + 1):(Ce*L);
	for n = 1:Ce
		if mod(n,Ce)~= 0
			PENT(ENum(n),:) = [NPL*(L+1), NPL*(L+1)-Ce+(n-1), NPL*(L+1)-Ce+n,...
                               NPL*L    , NPL*L-Ce+(n-1)    , NPL*L-Ce+n];
		else
			PENT(ENum(n),:) = [NPL*(L+1), NPL*(L+1)-Ce+(n-1), NPL*(L+1)-Ce,...
						       NPL*L    , NPL*L-Ce+(n-1)    , NPL*L-Ce];
		end
	end
end


