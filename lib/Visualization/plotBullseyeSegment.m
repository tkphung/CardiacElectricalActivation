function plotBullseyeSegment(segments, c, contours, cLim, cMap, cBarTitle,...
                             shadingSwitch)
                             
txtFormat = '%1.0f';

%% Bullseye plot

h = figure('Visible', 'on'); hold on

% Patches - flip x axis as theta = 0 is in the septum (on the right)
for i = 1:length(c)
    patch(-segments(i).xPatch, segments(i).yPatch, c(i,:), 'EdgeColor', 'None')    
end

% Contour lines
plot(contours.X1,contours.Y1, 'k', 'LineWidth', 3);
plot(contours.X2,contours.Y2, 'k', 'LineWidth', 3);
plot(contours.X3,contours.Y3, 'k', 'LineWidth', 3);
plot(contours.X0,contours.Y0, 'k', 'LineWidth', 3);
plot(contours.X4,contours.Y4, 'k', 'LineWidth', 3);
plot(contours.X5,contours.Y5, 'k', 'LineWidth', 3);
plot(contours.X6,contours.Y6, 'k', 'LineWidth', 3);
plot(contours.X7,contours.Y7, 'k', 'LineWidth', 3);
plot(contours.X8,contours.Y8, 'k', 'LineWidth', 3);
plot(contours.X9,contours.Y9, 'k', 'LineWidth', 3);
plot(contours.X10,contours.Y10, 'k', 'LineWidth', 3);
plot(contours.X11,contours.Y11, 'k', 'LineWidth', 3);
plot(contours.X12,contours.Y12, 'k', 'LineWidth', 3);
plot(contours.X13,contours.Y13, 'k', 'LineWidth', 3);


text(0,128,'Anterior','fontsize',12,'fontweight','b','HorizontalAlignment','center')
text(-128,0,'Septal','fontsize',12,'Rotation',90,'fontweight','b','HorizontalAlignment','center')
text(128,0,'Lateral','fontsize',12,'Rotation',270,'fontweight','b','HorizontalAlignment','center')
text(0,-128,'Posterior','fontsize',12,'fontweight','b','HorizontalAlignment','center')
     
set(gca, 'LineWidth', 3)


if shadingSwitch

    shading flat;
    colormap(cMap)
    
    if ~isempty(cLim)
        caxis(cLim)
    end

end

axis([-120 120 -120 120]);
axis off equal;
set(gcf,'Color','white')

% Set figure square
pos = get(gcf,'Position');
set(gcf, 'Position', [pos(1) pos(2) pos(3) pos(3)])

% removePadding
% fixPaperSize

% print(h, '-dpdf', fName)

for i = 1:length(c)
    
    % Print segments numbers if no patch data color coding
    if ~shadingSwitch
        text(-segments(i).xCenter, segments(i).yCenter, num2str(i, txtFormat),...
                'FontSize', 14, 'Color', [0 0 0],'HorizontalAlignment','center')
    % Else, place data value in center of each patch
    else
        text(-segments(i).xCenter, segments(i).yCenter, num2str(c(i), txtFormat),...
                'FontSize', 14, 'Color', [0 0 0],'HorizontalAlignment','center')

    end
end

% if ~isempty(fName)
    % print(h, '-dpdf', strcat(fName, 'Txt'))
% end
% close(h)


%% Print colorbar in a seperate figure

% if (~isempty(cLim) && shadingSwitch)
% 
%     h = figure('Visible', 'Off'); hold on
% 
%     % Set figure square
%     pos = get(gcf,'Position');
%     set(gcf, 'Position', [pos(1:3) pos(4)/4])
% 
%     hcb = colorbar('Orientation', 'Horizontal', 'position',[0.2 0.3 0.6 .4],...
%              'LineWidth', 3, 'FontSize', 20, 'TickLength', 0, 'Ticks', [cLim(1) mean(cLim) cLim(end)]);
% 
%     caxis(cLim)     
%     axis off
%     set(gcf,'Color','white')
% 
%     set(get(hcb,'Title'),'String',cBarTitle)
% 
%     fixPaperSize
% 
%     if ~isempty(fName)
%         print(h, '-dpdf', strcat(fName, 'Cbar'))
%     end
% end