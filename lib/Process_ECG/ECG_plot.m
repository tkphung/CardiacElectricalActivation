function H = ECG_plot(time,leads,varargin)
% H = ECG_plot(time,leads,varargin) plots the 12 lead ECG data in a 
% standard format leads matrix should be n by 12
% INPUT: time- numel(time) == size(leads,1)
%        leads- n by 12 matrix
%        varargin: 'handle' followed by plot handle
% OUTPUT: H- figure handle
% 
% TNPhung April 4, 2017

% Deal with VARARGIN
handle_flag = false;
if ~isempty(varargin)
    for jz = 1:2:numel(varargin)
        switch varargin{jz}
            case 'handle' % plot on existing figure
                handle_flag = true;
                handle = varargin{jz+1};
            otherwise
                error('ERROR: Check your varargins.')
        end
    end
end

% 12-Lead Plot
if handle_flag
    figure(handle), hold on
else
H = figure('numbertitle','off','name','12-Lead Plots');
end

subplot(3,4,1); hold on
plot(time,(leads(:,1)),'LineWidth',2);
title('I'); xlim([min(time) max(time)])
subplot(3,4,5); hold on
plot(time,(leads(:,2)),'LineWidth',2);
title('II'); xlim([min(time) max(time)])
subplot(3,4,9); hold on
plot(time,(leads(:,3)),'LineWidth',2);
title('III'); xlim([min(time) max(time)])
subplot(3,4,2); hold on
plot(time,(leads(:,4)),'LineWidth',2);
title('aVR'); xlim([min(time) max(time)])
subplot(3,4,6); hold on
plot(time,(leads(:,5)),'LineWidth',2);
title('aVL'); xlim([min(time) max(time)])
subplot(3,4,10); hold on
plot(time,(leads(:,6)),'LineWidth',2);
title('aVF'); xlim([min(time) max(time)])
subplot(3,4,3); hold on
plot(time,(leads(:,7)),'LineWidth',2);
title('V1'); xlim([min(time) max(time)])
subplot(3,4,7); hold on
plot(time,(leads(:,8)),'LineWidth',2);
title('V2'); xlim([min(time) max(time)])
subplot(3,4,11); hold on
plot(time,(leads(:,9)),'LineWidth',2);
title('V3'); xlim([min(time) max(time)])
subplot(3,4,4); hold on
plot(time,(leads(:,10)),'LineWidth',2);
title('V4'); xlim([min(time) max(time)])
subplot(3,4,8); hold on
plot(time,(leads(:,11)),'LineWidth',2);
title('V5'); xlim([min(time) max(time)])
subplot(3,4,12); hold on
plot(time,(leads(:,12)),'LineWidth',2);
title('V6'); xlim([min(time) max(time)])