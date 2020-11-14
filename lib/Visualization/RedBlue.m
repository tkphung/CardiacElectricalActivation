function rb = RedBlue(d)
%rb = RedBlue(d) takes in a scalar value (number of color divisions) and
%creates the divergent colormap. Best to use an odd number to have central
%value.
%OR if d input is a vector with values between 0 and 1, we will interpolate
%at those given values.
%Based on colorbrewer2.org divergent color scheme

% Original color scheme from  
og = [103,0,31    % red
      178,24,43
      214,96,77
      244,165,130
      253,219,199
      247,247,247 % Neutral
      209,229,240
      146,197,222
      67,147,195
      33,102,172
      5,48,97]./255;   % Blue
  
% % Plotting original RGB
% figure, hold on
% plot(linspace(0,1,11),og(:,1),'r','LineWidth',2)
% plot(linspace(0,1,11),og(:,2),'g','LineWidth',2)
% plot(linspace(0,1,11),og(:,3),'b','LineWidth',2)

if numel(d) == 1
    % Red to White
    % Number of divisions
    divs = ceil(d/2);
    R2W = interp1(linspace(0,1,6)',og(1:6,:),linspace(0,1,divs)','pchip');

    % White to Blue
    divs2 = d - divs + 1;
    W2B = interp1(linspace(0,1,6)',og(6:11,:),linspace(0,1,divs2)','pchip');

    % Combine them
    rb = [R2W; W2B(2:end,:)];

elseif numel(d)>1
    rb = interp1(linspace(0,1,11)',og,d,'pchip');
end



% % Plotting New RGB
% hold on
% plot(linspace(0,1,d),rb(:,1),'ro')
% plot(linspace(0,1,d),rb(:,2),'go')
% plot(linspace(0,1,d),rb(:,3),'bo')



