function rb = BlueYellowRed(d)
%rb = BlueYellowRed(d)) takes in a scalar value (number of color divisions) and
%creates the divergent colormap. Best to use an odd number to have central
%value. Blue-to-Yellow-to-Red
%OR if d input is a vector with values between 0 and 1, we will interpolate
%at those given values.
%Based on colorbrewer2.org divergent color scheme

% Original color scheme from  
og = [158,1,66
      213,62,79
      244,109,67
      253,174,97
      254,224,139
      255,255,191
      230,245,152
      171,221,164
      102,194,165
      50,136,189
      94,79,162]./255;

og = flipud(og);
	  
% % Plotting original RGB
% figure, hold on
% plot(linspace(0,1,11),og(:,1),'r','LineWidth',2)
% plot(linspace(0,1,11),og(:,2),'g','LineWidth',2)
% plot(linspace(0,1,11),og(:,3),'b','LineWidth',2)

% Red to White
% Number of divisions
divs = d;
if numel(d) == 1
rb = interp1(linspace(0,1,11)',og,linspace(0,1,divs)','pchip');
elseif numel(d)>1
rb = interp1(linspace(0,1,11)',og,d,'pchip');
end

% % Plotting New RGB
% hold on
% plot(linspace(0,1,d),rb(:,1),'ro')
% plot(linspace(0,1,d),rb(:,2),'go')
% plot(linspace(0,1,d),rb(:,3),'bo')



