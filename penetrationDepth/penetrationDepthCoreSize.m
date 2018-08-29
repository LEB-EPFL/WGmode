%-------------------------------------------------------------------%
%                    Penetration depth vs core size                 %
%-------------------------------------------------------------------%

%% WG given parameters mode m = 0
%corDim = 170;     % um 2a 2
n1 = 1.488;       % core pmma 1.488
n0 = 1.33;       % cladding 1.42 cytop 1.33
ns = 1.33;      % substrate 1.42
lambda = 0.647;   % um

gamma = 0;
k = 2*pi/lambda;
%gamma = (ns^2 - n0^2)/(n1^2 - ns^2);

corDim = (0.100:0.010:0.500); % [um]

v = k.*corDim./2.*sqrt(n1^2 - ns^2);

% DEBUG
% bx = (0:0.01:1);
% by = atan(sqrt(bx./(1 - bx))) +  ...
%                     atan(sqrt((bx + gamma)./(1 - bx))) - ...
%                     2.*v.*sqrt(1 - bx);                
% figure, plot(bx,by), hold on % END debug       

%b = zeros(1:length(v));
for vIdx = 1:length(v)
% Define the dispersion b-v function(v,b,m,gamma)
myFunc = @(br) atan(sqrt(br./(1 - br))) +  ...
                    atan(sqrt((br + gamma)./(1 - br))) - ...
                    2.*v(vIdx).*sqrt(1 - br);
% Get the b number from f(v,b,m,gamma) function roots               
[b(vIdx), byRoot] = findCFRootBisection(myFunc, [0, 1], 0.01, 10^8);   

end
% Get the effective index
ne = sqrt(b.*(n1^2 - ns^2) + ns^2);

beta = ne.*k;
u = v.*sqrt(1 - b);
w =  v.*sqrt(b);
wP = v.*sqrt(b + gamma);

% transverse propagation constants
ki = u./(corDim/2);        %sqrt(k^2*n1^2 - betak^2);
sigma = w./(corDim/2);     %sqrt(betak^2 - k^2*n0^2);
xi = wP./(corDim/2);       %sqrt(betak^2 - k^2*ns^2);

% phase
phi = 0.5.*atan(w./u) - 0.5.*atan(wP./u);

% power confinment 
% powerConf = (1 + (sin(u + phi))^2/(2*w) + (sin(u - phi))^2/(2*wP)) / ...
%                 (1 + 1/(2*w) + 1/(2*wP) );
% effective width
he = corDim + 1./sigma + 1./xi;

% Core size solutions
%dSigma = 1./sigma;
dSigma = lambda./(2.*pi)./sqrt(ne.^2 - n0^2); 
% Gly solution
n2 = 1.33;    % refractive index cladding glycerol 1.42
%thetaGly = (69.8:0.1:78); % angle of incidence at the interface core/cladding
thetaGly = (64:0.1:90); % angle of incidence at the interface core/cladding

% Get the effective index from angle
neAngle = sind(thetaGly).*n1;

% Penetration depth in Gly
dGly = lambda ./ (4.*pi) .* (n1.^2.*sind(thetaGly).^2 - n2.^2).^(-0.5) ;

figure,
% Create a plot with 2 x axes using the plotxx function
% [ax, h1, h2, h3] = plotyy( dH20,  thetaH20,  dGly, thetaGly, dRIdx, n2Vect, 'plot');

% Create a plot with 2 x axes using the plotxx function
[ax, h1, h2] = plotyy(  dGly, thetaGly,  dSigma,  corDim,'plot');

%style the plot
set(h1,'Color','b','LineWidth',1);
set(h2,'Color','k','LineWidth',1);
% set(h3,'Color','r','LineWidth',5);
set(ax(1),'XColor','b');
set(ax(2),'XColor','k');
% set(ax(3),'XColor','r');

% Use the axis handles to set the labels of the y axes
set(get(ax(1), 'Ylabel'), 'String', 'Incident Angle (degree)');
set(get(ax(2), 'Ylabel'), 'String', 'Core size');
% set(get(ax(3), 'Xlabel'), 'String', 'Refractive Index');
hold on

% xlabel ('Incident angle [degree]')
xlabel ('Penetration depth [um]')
legend (['Changing incident Angle (n2Sample = ' num2str(n2) '; n1Core = '  num2str(n1)], 'Changing Core Size')

figure(4),
plot( thetaGly, dGly, 'k--')
xlabel('Angle [degree]')
ylabel('Penetration depth [um]')
title(['Changing incident Angle (n2Sample = ' num2str(n2) '; n1Core =' num2str(n1)])
hold on

figure(5);
plot( corDim, dSigma, 'k')
xlabel('Core Dimention [um]')
ylabel('Penetration depth [um]')
title(['Changing core size (n2Sample =' num2str(n2) '; n1Core =' num2str(n1) ])
hold on

figure(6);
plot( corDim, ne, 'k')
xlabel('Core Dimention [um]')
ylabel('Effective index')
title(['Changing core size (n2Sample =' num2str(n2) '; n1Core =' num2str(n1) ])
hold on

figure(7);
plot( thetaGly, neAngle, 'k')
xlabel('Core Dimention [um]')
ylabel('Effective index')
title(['Changing incident Angle (n2Sample =' num2str(n2) '; n1Core =' num2str(n1) ])
hold on



%% Si3N4
%corDim = 170;     % um 2a 2
n1 = 2.038;       % core pmma 1.488
n0 = 1.33;       % cladding 1.42 cytop 1.33
ns = 1.52;      % substrate 1.42
lambda = 0.647;   % um

gamma = 0;
k = 2*pi/lambda;
%gamma = (ns^2 - n0^2)/(n1^2 - ns^2);

corDim = (0.100:0.010:2); % [um]
v = k.*corDim./2.*sqrt(n1^2 - ns^2);

% DEBUG
% bx = (0:0.01:1);
% by = atan(sqrt(bx./(1 - bx))) +  ...
%                     atan(sqrt((bx + gamma)./(1 - bx))) - ...
%                     2.*v.*sqrt(1 - bx);                
% figure, plot(bx,by), hold on % END debug       

figure
%b = zeros(1:length(v));
for vIdx = 1:length(v)
% Define the dispersion b-v function(v,b,m,gamma)
myFunc = @(br) atan(sqrt(br./(1 - br))) +  ...
                    atan(sqrt((br + gamma)./(1 - br))) - ...
                    2.*v(vIdx).*sqrt(1 - br);
% Get the b number from f(v,b,m,gamma) function roots               
[b(vIdx), byRoot] = findCFRootBisection(myFunc, [0, 1], 0.01, 10^8);   

end
% Get the effective index
ne = sqrt(b.*(n1^2 - ns^2) + ns^2);

beta = ne.*k;
u = v.*sqrt(1 - b);
w =  v.*sqrt(b);
wP = v.*sqrt(b + gamma);

% transverse propagation constants
ki = u./(corDim/2);        %sqrt(k^2*n1^2 - betak^2);
sigma = w./(corDim/2);     %sqrt(betak^2 - k^2*n0^2);
xi = wP./(corDim/2);       %sqrt(betak^2 - k^2*ns^2);

% phase
phi = 0.5.*atan(w./u) - 0.5.*atan(wP./u);

% power confinment 
% powerConf = (1 + (sin(u + phi))^2/(2*w) + (sin(u - phi))^2/(2*wP)) / ...
%                 (1 + 1/(2*w) + 1/(2*wP) );
% effective width
he = corDim + 1./sigma + 1./xi;

% Core size solutions
%dSigma = 1./sigma;
dSigma = lambda./(2.*pi)./sqrt(ne.^2 - n0^2); 

% Critical angle
thetaCoCritical = n1*sqrt((n1^2-n0^2)/n1^2)*180/pi; % in degree
phiCritical = asin(sind(thetaCritical)/n1)*180/pi; % in degree
thetaInCritical = 90 - phiCritical;

% Changing angle
thetaIn = (thetaInCritical:0.1:90); % angle of incidence at the interface core/cladding
phi = (0:0.1:phiCritical); % angle between ray direction inside the core and tranmission direction (beta direction)
thetaCo = (0:0.1:thetaCoCritical); % angle between the coupling ray and the transmission


% Get the effective index from angle
% with angle equal to thetaIn
%neAngle = cosd(thetaIn).*n1;
% with angle equal to phi
%neAngle = sind(phi).*n1;
% with angle equal to thetaCo
neAngle = sqrt(n1.^2-sind(thetaCo).^2);

% Penetration depth 
%dAngle = lambda ./ (2.*pi) .* (n1.^2.*sind(thatIn).^2 - n0.^2).^(-0.5); 
%dAngle = lambda ./ (2.*pi) .* (n1.^2.*cosd(phi).^2 - n0.^2).^(-0.5);
dAngle = lambda ./ (2.*pi) .* (n1.^2 - sind(thetaCo).^2 - n0.^2).^(-0.5);

theta = thetaCo;
figure(5)
plot( corDim, dSigma, 'b')
xlabel('Core Dimention [um]')
ylabel('Penetration depth [um]')
title('Changing core size (n2Sample = 1.33; n1Core = 1.488/2.038')

figure(4),
plot( theta, dAngle, 'b')
hold on
plot(repmat(thetaCoCritical, 1, length(dAngle)), dAngle, 'r');
xlabel('Angle [degree]')
ylabel('Penetration depth [um]')
title(['Changing incident Angle (n2Sample = ' num2str(n0) '; n1Core =' num2str(n1)])
legend('angle', 'critical angle')
hold off

figure(6);
plot( corDim, ne, 'b')
xlabel('Core Dimention [um]')
ylabel('Effective index')
title(['Changing core size (n2Sample =' num2str(n0) '; n1Core =' num2str(n1) ])

figure(7);
plot( theta, neAngle, 'b')
xlabel('Angle [degree]')
ylabel('Effective index')
title(['Changing incident Angle (n2Sample =' num2str(n0) '; n1Core =' num2str(n1) ])


%% Glass
%corDim = 170;     % um 2a 2
n1 = 1.515;       % core pmma 1.488
n0 = 1.33;       % cladding 1.42 cytop 1.33
ns = 1.33;      % substrate 1.42
lambda = 0.647;   % um

gamma = 0;
k = 2*pi/lambda;
%gamma = (ns^2 - n0^2)/(n1^2 - ns^2);

corDim = (0.100:0.010:0.500); % [um]
v = k.*corDim./2.*sqrt(n1^2 - ns^2);

% DEBUG
% bx = (0:0.01:1);
% by = atan(sqrt(bx./(1 - bx))) +  ...
%                     atan(sqrt((bx + gamma)./(1 - bx))) - ...
%                     2.*v.*sqrt(1 - bx);                
% figure, plot(bx,by), hold on % END debug       

figure
%b = zeros(1:length(v));
for vIdx = 1:length(v)
% Define the dispersion b-v function(v,b,m,gamma)
myFunc = @(br) atan(sqrt(br./(1 - br))) +  ...
                    atan(sqrt((br + gamma)./(1 - br))) - ...
                    2.*v(vIdx).*sqrt(1 - br);
% Get the b number from f(v,b,m,gamma) function roots               
[b(vIdx), byRoot] = findCFRootBisection(myFunc, [0, 1], 0.01, 10^8);   

end
% Get the effective index
ne = sqrt(b.*(n1^2 - ns^2) + ns^2);

beta = ne.*k;
u = v.*sqrt(1 - b);
w =  v.*sqrt(b);
wP = v.*sqrt(b + gamma);

% transverse propagation constants
ki = u./(corDim/2);        %sqrt(k^2*n1^2 - betak^2);
sigma = w./(corDim/2);     %sqrt(betak^2 - k^2*n0^2);
xi = wP./(corDim/2);       %sqrt(betak^2 - k^2*ns^2);

% phase
phi = 0.5.*atan(w./u) - 0.5.*atan(wP./u);

% power confinment 
% powerConf = (1 + (sin(u + phi))^2/(2*w) + (sin(u - phi))^2/(2*wP)) / ...
%                 (1 + 1/(2*w) + 1/(2*wP) );
% effective width
he = corDim + 1./sigma + 1./xi;

% Core size solutions
%dSigma = 1./sigma;
dSigma = lambda./(2.*pi)./sqrt(ne.^2 - n0^2); 

% Changing angle
theta = (64:0.1:90); % angle of incidence at the interface core/cladding

% Get the effective index from angle
neAngle = sind(theta).*n1;

% Penetration depth 
dAngle = lambda ./ (4.*pi) .* (n1.^2.*sind(theta).^2 - n2.^2).^(-0.5) ;

figure(5)
plot( corDim, dSigma, 'c')
xlabel('Core Dimention [um]')
ylabel('Penetration depth [um]')
title('Changing core size (n2Sample = 1.33; n1Core = 1.488/1.515/2.038')
legend('PMMA core', 'Si3N4 core', 'Glass')

figure(4),
plot( theta, dAngle, 'c')
xlabel('Angle [degree]')
ylabel('Penetration depth [um]')
title(['Changing incident Angle (n2Sample = ' num2str(n2) '; n1Core =' num2str(n1)])
legend('PMMA core', 'Si3N4 core', 'Glass')

figure(6);
plot( corDim, ne, 'c')
xlabel('Core Dimention [um]')
ylabel('Effective index')
title(['Changing core size (n2Sample =' num2str(n2) '; n1Core =' num2str(n1) ])
legend('PMMA core', 'Si3N4 core', 'Glass')

figure(7);
plot( theta, neAngle, 'c')
xlabel('Core Dimention [um]')
ylabel('Effective index')
title(['Changing incident Angle (n2Sample =' num2str(n2) '; n1Core =' num2str(n1) ])
legend('PMMA core', 'Si3N4 core', 'Glass')

