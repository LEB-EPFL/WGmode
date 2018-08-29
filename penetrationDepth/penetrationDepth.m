%------------------------------------------------------------------%
%                         penetrationDepth                         %
%------------------------------------------------------------------%

%% PMMA/cytop

lambda = 0.647;  % wavelength [um]
n1 = 1.488;    % refractive index core PMMA

% H20 Solution/cytop
n2 = 1.40;    % refractive index cladding cytop


thetaH20 = (73:0.1:89); % angle of incidence at the interface core/cladding

% Penetration depth in H20
dH20 = lambda ./ (4.*pi) .* (n1.^2.*sind(thetaH20).^2 - n2.^2).^(-0.5) ;

% Core change ref Idx
n1Vect = (1.49: 0.0001 : 1.5);    
%n1Vect = (1.57: 0.0001 : 1.6);    

figure,
theta = (73: 1 :73); % angle of incidence at the interface core/cladding
for angIdx = 1: length(theta)
    % Penetration depth in SiN
    dRIdx(:, angIdx) = lambda ./ (4.*pi) .* (n1Vect.^2.*sind(theta(angIdx)).^2 - n2.^2).^(-0.5) ;

    % Create a plot with 2 x axes using the plotxx function
    % [ax, h1, h2, h3] = plotyy( dH20,  thetaH20,  dGly, thetaGly, dRIdx, n2Vect, 'plot');

    % Create a plot with 2 x axes using the plotxx function
    [ax, h1, h2] = plotyy(  dH20, thetaH20,  dRIdx(:, angIdx),  n1Vect,'plot');
    hold on
end
%style the plot
set(h1,'Color','b','LineWidth',1);
set(h2,'Color','k','LineWidth',1);
% set(h3,'Color','r','LineWidth',5);
set(ax(1),'XColor','b');
set(ax(2),'XColor','k');
% set(ax(3),'XColor','r');

% Use the axis handles to set the labels of the y axes
set(get(ax(1), 'Ylabel'), 'String', 'Incident Angle (degree)');
set(get(ax(2), 'Ylabel'), 'String', 'Refractive Index');
% set(get(ax(3), 'Xlabel'), 'String', 'Refractive Index');
hold on

% xlabel ('Incident angle [degree]')
xlabel ('Penetration depth [um]')
legend ('Changing incident Angle (n1(PMMA) = 1.49; n2Media = 1.42)', 'Changing Refractive Index (theta = 70 degree)')
title('Changing refractive index PMMA - lambda = 561')

%% Si2N3

lambda = 0.647;  % wavelength [um]
n1 = 2.038;    % refractive index core Si2N3

% H20 Solution
% n2 = 1.42;    % refractive index cladding gel close to water
n2 = 1.5;    % refractive index cladding SiO2


% thetaH20 = (44.5:0.1:89); % angle of incidence at the interface core/cladding
thetaH20 = (48:0.1:89); % angle of incidence at the interface core/cladding

% Penetration depth in H20
dH20 = lambda ./ (4.*pi) .* (n1.^2.*sind(thetaH20).^2 - n2.^2).^(-0.5) ;

% SiN solution change ref Idx
n1Vect = (2.038: 0.0001 : 2.05);    % refractive index cladding glycerol

figure,
% theta = (44.5: 1 :44.5); % angle of incidence at the interface core/cladding
theta = (49: 1 :49); % angle of incidence at the interface core/cladding
for angIdx = 1: length(theta)

    % Penetration depth in SiN
    dRIdx(:, angIdx) = lambda ./ (4.*pi) .* (n1Vect.^2.*sind(theta(angIdx)).^2 - n2.^2).^(-0.5) ;

    % Create a plot with 2 x axes using the plotxx function
    % [ax, h1, h2, h3] = plotyy( dH20,  thetaH20,  dGly, thetaGly, dRIdx, n2Vect, 'plot');

    % Create a plot with 2 x axes using the plotxx function
    [ax, h1, h2] = plotyy(  dH20, thetaH20,  dRIdx(:, angIdx),  n1Vect,'plot');
    hold on
end
%style the plot
set(h1,'Color','b','LineWidth',1);
set(h2,'Color','k','LineWidth',1);
% set(h3,'Color','r','LineWidth',5);
set(ax(1),'XColor','b');
set(ax(2),'XColor','k');
% set(ax(3),'XColor','r');

% Use the axis handles to set the labels of the y axes
set(get(ax(1), 'Ylabel'), 'String', 'Incident Angle (degree)');
set(get(ax(2), 'Ylabel'), 'String', 'Refractive Index');
% set(get(ax(3), 'Xlabel'), 'String', 'Refractive Index');
hold on

% xlabel ('Incident angle [degree]')
xlabel ('Penetration depth [um]')
legend ('Changing incident Angle (n1 (Si3N4) = 2.038; n2Media = 1.4)', 'Changing Refractive Index (theta = 43.4 degree)')
title('Changing refractive index of core Si3N4 - lambda = 0.641')



%% Checking Axial res via multiangle TIRF paper - changing solution refractive index

% lambda = 0.561;  % wavelength [um]
% lambda = 0.488;  % wavelength [um]
lambda = 0.647;  % wavelength [um]
n1 = 1.5168;    % refractive index core glass

% H20 Solution
n2 = 1.3333;    % refractive index cladding gel close to water
thetaH20 = (62:0.1:70); % angle of incidence at the interface core/cladding 

% Penetration depth in H20
dH20 = lambda ./ (4.*pi) .* (n1.^2.*sind(thetaH20).^2 - n2.^2).^(-0.5) ;

% Gly solution
n2 = 1.42;    % refractive index cladding glycerol
thetaGly = (69.8:0.1:78); % angle of incidence at the interface core/cladding

% Penetration depth in Gly
dGly = lambda ./ (4.*pi) .* (n1.^2.*sind(thetaGly).^2 - n2.^2).^(-0.5) ;

% Gly solution change ref Idx
n2Vect = (1.42: 0.0001 : 1.425);    % refractive index cladding glycerol
theta = 70; % angle of incidence at the interface core/cladding

% Penetration depth in Gly
dRIdx = lambda ./ (4.*pi) .* (n1.^2.*sind(theta).^2 - n2Vect.^2).^(-0.5) ;

figure,
% Create a plot with 2 x axes using the plotxx function
% [ax, h1, h2, h3] = plotyy( dH20,  thetaH20,  dGly, thetaGly, dRIdx, n2Vect, 'plot');

% Create a plot with 2 x axes using the plotxx function
[ax, h1, h2] = plotyy(  dGly, thetaGly,  dRIdx,  n2Vect,'plot');

%style the plot
set(h1,'Color','b','LineWidth',1);
set(h2,'Color','k','LineWidth',1);
% set(h3,'Color','r','LineWidth',5);
set(ax(1),'XColor','b');
set(ax(2),'XColor','k');
% set(ax(3),'XColor','r');

% Use the axis handles to set the labels of the y axes
set(get(ax(1), 'Ylabel'), 'String', 'Incident Angle (degree)');
set(get(ax(2), 'Ylabel'), 'String', 'Refractive Index');
% set(get(ax(3), 'Xlabel'), 'String', 'Refractive Index');
hold on

% xlabel ('Incident angle [degree]')
xlabel ('Penetration depth [um]')
legend ('Changing incident Angle (n2Gly = 1.42; n1Glass = 1.517)', 'Changing Refractive Index (theta = 70 degree)')

%% Checking Axial res via multiangle TIRF paper - changing ref index by changing n1

lambda = 0.561;  % wavelength [um]
% lambda = 0.488;  % wavelength [um]
n1 = 1.5168;    % refractive index core glass

% H20 Solution
% n2 = 1.3333;    % refractive index cladding gel close to water
% thetaH20 = (62:0.1:70); % angle of incidence at the interface core/cladding 

% Penetration depth in H20
% dH20 = lambda ./ (4.*pi) .* (n1.^2.*sind(thetaH20).^2 - n2.^2).^(-0.5) ;

% Gly solution
n2 = 1.42;    % refractive index cladding glycerol
thetaGly = (69.8:0.1:78); % angle of incidence at the interface core/cladding

% Penetration depth in Gly
dGly = lambda ./ (4.*pi) .* (n1.^2.*sind(thetaGly).^2 - n2.^2).^(-0.5) ;

% Core change ref Idx
n1Vect = (1.514: 0.0001 : 1.53);    % refractive index core
theta = 70; % angle of incidence at the interface core/cladding

% Penetration depth in Gly
dRIdx = lambda ./ (4.*pi) .* (n1Vect.^2.*sind(theta).^2 - n2.^2).^(-0.5) ;

figure,
% Create a plot with 2 x axes using the plotxx function
% [ax, h1, h2, h3] = plotyy( dH20,  thetaH20,  dGly, thetaGly, dRIdx, n2Vect, 'plot');

% Create a plot with 2 x axes using the plotxx function
[ax, h1, h2] = plotyy(  dGly, thetaGly,  dRIdx,  n1Vect,'plot');

%style the plot
set(h1,'Color','b','LineWidth',1);
set(h2,'Color','k','LineWidth',1);
% set(h3,'Color','r','LineWidth',5);
set(ax(1),'XColor','b');
set(ax(2),'XColor','k');
% set(ax(3),'XColor','r');

% Use the axis handles to set the labels of the y axes
set(get(ax(1), 'Ylabel'), 'String', 'Incident Angle (degree)');
set(get(ax(2), 'Ylabel'), 'String', 'Refractive Index');
% set(get(ax(3), 'Xlabel'), 'String', 'Refractive Index');
hold on

% xlabel ('Incident angle [degree]')
xlabel ('Penetration depth [um]')
legend ('Changing incident Angle (n2Gly = 1.42; n1Glass = 1.517)', 'Changing Refractive Index (theta = 70 degree)')

title('Changing refractive index of core glass - lambda = 0.561')