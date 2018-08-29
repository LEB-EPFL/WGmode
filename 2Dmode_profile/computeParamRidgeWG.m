
%---------------------------------------------------------------------%
%           Ridge waveguide: "effective index" approximation
%---------------------------------------------------------------------%
%
%                       ______________________w__
%    __________________|                         |__________________no
%  t|                  |  h                                         n1
%   |__________________|______________________________________________
%                                                                   ns
%---------------------------------------------------------------------%
%
% Anna Archetti - EPFL 2017
%
%---------------------------------------------------------------------%
% Step 1
%---------------------------------------------------------------------%
% T thickness
corDimT = 0.200;    % 2a [um] 2a 1 0.1
n1 = 2.038;       % core 1.5
n0 = 1.38;         % cladding 1
ns = 1.45;      % substrate 1.45
lambda = 0.647;  % um 1.55
%modeNUM = 3;

rootBInt = [0, 1];
k = 2*pi/lambda;
gamma = (ns^2 - n0^2)/(n1^2 - ns^2);

vT = k*corDimT/2*sqrt(n1^2 - ns^2);
% DEBUG
bx = (0:0.01:1);
by = atan(sqrt(bx./(1 - bx))) +  ...
                    atan(sqrt((bx + gamma)./(1 - bx))) - ...
                    2.*vT.*sqrt(1 - bx);               
figure, plot(bx,by), hold on % END debug

% Define the dispersion v-b function f(v, b, m, gamma)
myFunc = @(br) atan(sqrt(br./(1 - br))) +  ...
                    atan(sqrt((br + gamma)./(1 - br))) - ...
                    2.*vT.*sqrt(1 - br);

% Get the b number for the overdefined v
[bT, byRoot] = findCFRootBisection(myFunc, rootBInt, 0.001, 10^8);   

% Effective index
neT1 = sqrt(bT*(n1^2 - ns^2) + ns^2);

betaT = neT1*k;
u = vT*sqrt(1 - bT);
w =  vT*sqrt(bT);
wP = vT*sqrt(bT + gamma);

% transverse propagation constants
ki = u/(corDimT/2);        %sqrt(k^2*n1^2 - betak^2);
sigma = w/(corDimT/2);     %sqrt(betak^2 - k^2*n0^2);
xi = wP/(corDimT/2);       %sqrt(betak^2 - k^2*ns^2);

% phase
phi = 0.5*atan(w/u) - 0.5*atan(wP/u);

% power confinment 
powerConfT = (1 + (sin(u + phi))^2/(2*w) + (sin(u - phi))^2/(2*wP)) / ...
                (1 + 1/(2*w) + 1/(2*wP) );
% effective width
heT = corDimT + 1/sigma + 1/xi;

% numerical aperture
NA0Y = n1*sqrt(2*(n1 - n0));
NAsY = n1*sqrt(2*(n1 - ns));

% cut-off lambda
lambdaC0T = 2*pi/vT*corDimT/2*NA0Y;
lambdaCsT = 2*pi/vT*corDimT/2*NAsY;





%%
%------------------------------------------------------------%
% H thickness
%------------------------------------------------------------%
corDimH = 0.250;  % [um] 2a 1.2
rootBInt = [0, 1];
k = 2*pi/lambda;
gamma = (ns^2 - n0^2)/(n1^2 - ns^2);

vH = k*corDimH/2*sqrt(n1^2 - ns^2);
% DEBUG
bx = (0:0.01:1);
by = atan(sqrt(bx./(1 - bx))) +  ...
                    atan(sqrt((bx + gamma)./(1 - bx))) - ...
                    2.*vH.*sqrt(1 - bx);             
figure, plot(bx,by), hold on % END debug

myFunc = @(br) atan(sqrt(br./(1 - br))) +  ...
                    atan(sqrt((br + gamma)./(1 - br))) - ...
                    2.*vH.*sqrt(1 - br);
                
[bH, byRoot] = findCFRootBisection(myFunc, rootBInt, 0.01, 10^8);   


% effective index
neH1 = sqrt(bH*(n1^2 - ns^2) + ns^2);

betaH = neH1*k;
u = vH*sqrt(1 - bH);
w =  vH*sqrt(bH);
wP = vH*sqrt(bH + gamma);

% transverse propagation constants
ki = u/(corDimH/2);        %sqrt(k^2*n1^2 - betak^2);
sigma = w/(corDimH/2);     %sqrt(betak^2 - k^2*n0^2);
xi = wP/(corDimH/2);       %sqrt(betak^2 - k^2*ns^2);

% phase
phi = 0.5*atan(w/u) - 0.5*atan(wP/u);

% power confinment 
powerConfH = (1 + (sin(u + phi))^2/(2*w) + (sin(u - phi))^2/(2*wP)) / ...
                (1 + 1/(2*w) + 1/(2*wP) );

% effective width
heH = corDimH + 1/sigma + 1/xi;
dHEW = 1/sigma;

% penetrationd depth
dH = lambda/(2*pi)/sqrt(neH1^2 - n0^2); 

% cut-off lambda
lambdaC0H = 2*pi/vH*corDimH/2*NA0Y;
lambdaCsH = 2*pi/vH*corDimH/2*NAsY;


% Electric field fundamental mode E(x)
% y = (-corDim/2 - 4: 0.1 :corDim/2 + 4);
y = (- 5 : 0.01 : 5);
for yIdx = 1: length(y)   
    ex(yIdx) = Efield(y(yIdx), corDimH, ki, sigma, phi, xi);
end


figure,
plot(y, ex, 'k')
hold on
plot( repmat(-corDimH/2,1, length(ex)), ex, 'r--')
plot( repmat(+corDimH/2,1, length(ex)), ex, 'r--')
xlabel('y[um]')
ylabel('Efield/Emax')
title('TE0 mode profile along Y')
hold off


%%
%-----------------------------------------------------------%
% Step 2
%-----------------------------------------------------------%
corDimW = 1.5;    % 2a [um] 2a 3
n1 = neH1;      % core
% neT1 = 0;
n0 = neT1;      % cladding
ns = neT1;      % substrate
lambda = 0.647;  % um 1.55

rootBInt = [0, 1];
k = 2*pi/lambda;
gamma = (ns^2 - n0^2)/(n1^2 - ns^2);

v2 = k*corDimW/2*sqrt(n1^2 - ns^2);

% numerical aperture
NA0X = n1*sqrt(2*(n1 - n0));
NAsX = n1*sqrt(2*(n1 - ns));

% DEBUG
bx = (0:0.01:1);
by = atan(sqrt(bx./(1 - bx))) +  ...
                    atan(sqrt((bx + gamma)./(1 - bx))) - ...
                    2.*v2.*sqrt(1 - bx);               
figure, plot(bx,by), hold on % END debug 

% Grab the b number
myFunc = @(br) atan(sqrt(br./(1 - br))) +  ...
                    atan(sqrt((br + gamma)./(1 - br))) - ...
                    2.*v2.*sqrt(1 - br);
                
[b2, byRoot] = findCFRootBisection(myFunc, rootBInt, 0.01, 10^8);   


% effective index
ne2 = sqrt(b2*(n1^2 - ns^2) + ns^2);

beta = ne2*k;
u = v2*sqrt(1 - b2);
w =  v2*sqrt(b2);
wP = v2*sqrt(b2 + gamma);

% transverse propagation constants
ki = u/(corDimW/2);        %sqrt(k^2*n1^2 - betak^2);
sigma = w/(corDimW/2);     %sqrt(betak^2 - k^2*n0^2);
xi = wP/(corDimW/2);       %sqrt(betak^2 - k^2*ns^2);

% phase
phi = 0.5*atan(w/u) - 0.5*atan(wP/u);

% power confinment 
powerConf = (1 + (sin(u + phi))^2/(2*w) + (sin(u - phi))^2/(2*wP)) / ...
                (1 + 1/(2*w) + 1/(2*wP) );

% effective width
he = corDimW + 1/sigma + 1/xi;
dEW = 1/sigma;

% penetrationd depth
depth = lambda/(2*pi)/sqrt(ne2^2 - n0^2); 

% cut-off lambda
lambdaC02 = 2*pi/v2*corDimW/2*NA0X;
lambdaCs2 = 2*pi/v2*corDimW/2*NAsX;


% Electric field fundamental mode E(x)
%x = (-corDimW/2 - 2: 0.1 :corDimW/2 + 2);
x = (- 5 : 0.01 : 5);
for xIdx = 1: length(x)

    % Electric field distribution: TE mode profile in the plane perpendicular
    % to the propagation direction
    ey(xIdx) = Efield(x(xIdx), corDimW, ki, sigma, phi, xi);
    
end


figure,
plot(x, ey, 'k')
hold on
plot( repmat(-corDimW/2,1, length(ey)), ey, 'r--')
plot( repmat(+corDimW/2,1, length(ey)), ey, 'r--')
xlabel('x[um]')
ylabel('Efield/Emax')
title('TE0 mode profile along X')

% Rectangular grid in 2-D 
% [X, Y] = meshgrid (-4 : +4 , -2 : +2);
% Z = Efield(X)*Efield(Y);
pixSize = ( x(min(length(ex), length(ey))) - x(1) )/( min(length(ex), length(ey)) - 1); % um
Exy = ex( 1:min(length(ex), length(ey)) )'*ey( 1:min(length(ex), length(ey)) );

% figure, imagesc( flipud(Exy) );
figure, image( flipud(Exy),'CDataMapping','scaled'), colormap('jet')


xticklabels = (y(1) : 0.5 : y(min(length(ex), length(ey))));
xticks = linspace(1, min(length(ex), length(ey)), numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)

yticklabels = (x(1) : 0.5 : x(min(length(ex), length(ey))));
yticks = linspace(1, min(length(ex), length(ey)), numel(yticklabels));
set(gca, 'YTick', yticks, 'YTickLabel', flipud(yticklabels(:)))  
hold on

% ridgeT = 1; %um
% ridgeH = 1.2; %um
% ridgeWidth = 3; % um

ridgeT = corDimH/2 - (corDimH - corDimT); %um
ridgeH = corDimH; %um
ridgeWidth = corDimW; % um

xmax = get(gca, 'xLim');
ymax = get(gca, 'yLim');

orizXEdge1 =  xmax(2)/2 + y(1 : min(length(ex), length(ey)))./pixSize + 0.5;
orizYEdge1 =  repmat( (ymax(2)/2 - ridgeT/pixSize + 0.5), 1, min(length(ex), length(ey)));

% orizXEdge2 =  xmax(2)/2 + y(1 : min(length(ex), length(ey)))./pixSize + 0.5;
% orizYEdge2 =  repmat( (ymax(2)/2 + ridgeT/2/pixSize + 0.5), 1, min(length(ex), length(ey)));

orizXEdge3 =  xmax(2)/2 + y(1 : min(length(ex), length(ey)))./pixSize + 0.5;
orizYEdge3 =  repmat( (ymax(2)/2 - ridgeH/2/pixSize + 0.5), 1, min(length(ex), length(ey)));

orizXEdge4 =  xmax(2)/2 + y(1 : min(length(ex), length(ey)))./pixSize + 0.5;
orizYEdge4 =  repmat( (ymax(2)/2 + ridgeH/2/pixSize + 0.5), 1, min(length(ex), length(ey)));

vertYEdge1 =  xmax(2)/2 + y(1 : min(length(ex), length(ey)))./pixSize + 0.5;
vertXEdge1 =  repmat( (ymax(2)/2 - ridgeWidth/2/pixSize + 0.5), 1, min(length(ex), length(ey)));

vertYEdge2 =  xmax(2)/2 + y(1 : min(length(ex), length(ey)))./pixSize + 0.5;
vertXEdge2 =  repmat( (ymax(2)/2 + ridgeWidth/2/pixSize + 0.5), 1, min(length(ex), length(ey)));


plot(orizXEdge1, orizYEdge1, 'r--', 'LineWidth',2)
% plot(orizXEdge2, orizYEdge2, 'r--')
plot(orizXEdge3, orizYEdge3, 'r--', 'LineWidth',2)
plot(orizXEdge4, orizYEdge4, 'r--', 'LineWidth',2)
plot(vertXEdge1, vertYEdge1, 'r--', 'LineWidth',2)
plot(vertXEdge2, vertYEdge2, 'r--', 'LineWidth',2)

axis equal
xlabel('x [um]')
ylabel('y [um]')
title('Electric field propagation XY')


figure, mesh(flipud(Exy)),
hold on,
image( flipud(Exy),'CDataMapping','scaled'), colormap('jet')
xticklabels = (y(1) : 0.5 : y(min(length(ex), length(ey))));
xticks = linspace(1, min(length(ex), length(ey)), numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)

yticklabels = (x(1) : 0.5 : x(min(length(ex), length(ey))));
yticks = linspace(1, min(length(ex), length(ey)), numel(yticklabels));
set(gca, 'YTick', yticks, 'YTickLabel', flipud(yticklabels(:)))  
xlabel('x [um]')
ylabel('y [um]')
zlabel('Normalized Electric field profile')
title('Electric field propagation XY')

%% Ridge WG effective index method Y variable
% 1/Y d2/dy2 Y + [k^2*n^2 -k^2*ne(x)^2] = 0

corDim = 2;  % um

nc = 1.5; % core
na = 1;% cladding
%nr = 1.5; % rib
ns = 1.45; % substrate
lambda = 1.5; % um
d = 3; % um
h = 1.2;
s = h;

k = 2*pi/lambda;
nex = (1.451:0.001:1.499);
ney =  sin(k.*sqrt(nc^2 - nex.^2).*d).*exp(-2*(k.*sqrt(nex.^2 - ns^2).*s) + atanh(sqrt(nex.^2 - ns^2)./sqrt(nex.^2 - na^2))) - ...
    sin(k.*sqrt(nc^2 - nex.^2).*d - 2*atan(sqrt(nex.^2 - ns^2)./sqrt(nc^2 - nex.^2)));
                
figure, plot(nex, ney), hold on % END debug   

myFunc = @(ne) sin(k.*sqrt(nc^2 - ne.^2).*d).*exp(-2.*(k.*sqrt(ne.^2 - ns^2).*s) + atanh(sqrt(ne.^2 - ns^2)./sqrt(ne.^2 - na^2))) - ...
    sin(k.*sqrt(nc^2 - ne.^2).*d - 2.*atan(sqrt(ne.^2 - ns^2)./sqrt(nc^2 - ne.^2)));
[neRootH1a, yRoot] = findCFRootBisection(myFunc, [1.451, 1.47], 0.001, 10^8);
[neRootH1b, yRoot] = findCFRootBisection(myFunc, [1.48, 1.495], 0.001, 10^8);



d = 3;
t = 1;
s = t;

k = 2*pi/lambda;
nex = (1.451:0.001:1.499);
ney =  sin(k.*sqrt(nc^2 - nex.^2).*d).*exp(-2*(k.*sqrt(nex.^2 - ns^2).*s) + atanh(sqrt(nex.^2 - ns^2)./sqrt(nex.^2 - na^2))) - ...
    sin(k.*sqrt(nc^2 - nex.^2).*d - 2*atan(sqrt(nex.^2 - ns^2)./sqrt(nc^2 - nex.^2)));
                
figure, plot(nex, ney), hold on % END debug  

myFunc = @(ne) sin(k.*sqrt(nc^2 - ne.^2).*d).*exp(-2.*(k.*sqrt(ne.^2 - ns^2).*s) + atanh(sqrt(ne.^2 - ns^2)./sqrt(ne.^2 - na^2))) - ...
    sin(k.*sqrt(nc^2 - ne.^2).*d - 2.*atan(sqrt(ne.^2 - ns^2)./sqrt(nc^2 - ne.^2)));

[neRootT1a, yRoot] = findCFRootBisection(myFunc, [1.451, 1.47], 0.001, 10^8);
[neRootT1b, yRoot] = findCFRootBisection(myFunc, [1.48, 1.495], 0.001, 10^8);

%% Ridge WG effective index method X variable
% 1/X d2/dx2 X + [k^2*ne(x)^2 - beta^2] = 0
corDim = 2;  % um

nc = 1.5; % core
na = 1;% cladding
%nr = 1.5; % rib
ns = 1.45; % substrate
lambda = 1.5; % um
d = 3; % um
h = 1.2;
s = h;


k = 2*pi/lambda;
nex = (1.451:0.001:1.499);
ney =  sin(k.*sqrt(nc^2 - nex.^2).*d).*exp(-2*(k.*sqrt(nex.^2 - ns^2).*s) + atanh(sqrt(nex.^2 - ns^2).*na^2./(sqrt(nex.^2 - na^2).*ns^2))) - ...
    sin(k.*sqrt(nc^2 - nex.^2).*d - 2*atan(sqrt(nex.^2 - ns^2).*nc^2./(sqrt(nc^2 - nex.^2).*ns^2)));
                
figure, plot(nex, ney), hold on % END debug   

myFunc = @(ne) sin(k.*sqrt(nc^2 - ne.^2).*d).*exp(-2.*(k.*sqrt(ne.^2 - ns^2).*s) + atanh(sqrt(ne.^2 - ns^2).*na^2./(sqrt(ne.^2 - na^2).*ns^2))) - ...
    sin(k.*sqrt(nc^2 - ne.^2).*d - 2.*atan(sqrt(ne.^2 - ns^2).*nc^2./(sqrt(nc^2 - ne.^2).*ns^2)));
[neRootH2a, yRoot] = findCFRootBisection(myFunc, [1.451, 1.47], 0.001, 10^8);
[neRootH2b, yRoot] = findCFRootBisection(myFunc,  [1.48, 1.495], 0.001, 10^8);



d = 3;
t = 1;
s = t;

k = 2*pi/lambda;
nex = (1.451:0.001:1.499);
ney =  sin(k.*sqrt(nc^2 - nex.^2).*d).*exp(-2*(k.*sqrt(nex.^2 - ns^2).*s) + atanh(sqrt(nex.^2 - ns^2).*na^2./(sqrt(nex.^2 - na^2).*ns^2))) - ...
    sin(k.*sqrt(nc^2 - nex.^2).*d - 2*atan(sqrt(nex.^2 - ns^2).*nc^2./(sqrt(nc^2 - nex.^2).*ns^2)));
                
figure, plot(nex, ney), hold on % END debug  

myFunc = @(ne) sin(k.*sqrt(nc^2 - ne.^2).*d).*exp(-2.*(k.*sqrt(ne.^2 - ns^2).*s) + atanh(sqrt(ne.^2 - ns^2).*na^2./(sqrt(ne.^2 - na^2).*ns^2))) - ...
    sin(k.*sqrt(nc^2 - ne.^2).*d - 2.*atan(sqrt(ne.^2 - ns^2).*nc^2./(sqrt(nc^2 - ne.^2).*ns^2)));

[neRootT2a, yRoot] = findCFRootBisection(myFunc,[1.451, 1.47], 0.001, 10^8);
[neRootT2b, yRoot] = findCFRootBisection(myFunc,  [1.48, 1.495], 0.001, 10^8);

%% TE modes
% x = (-corDim/2: 0.1 :corDim/2);
% Ey = zeros(size(x), modeMum);
% 
% for mIdx = 1: modeNUM
% 
%     v(mIdx) = ( m(mIdx)*pi/2 + 1/2.*atan(sqrt(b./(1 - b))) ...
%                 + 1/2.*atan(sqrt((b + gamma(gammaIdx))./(1 - b))) )./sqrt(1 - b);
%             
%     u(mIdx) = v*sqrt(1-b);
%     w(mIdx) = v*sqrt(b);
%     wP(mIdx) = v*sqrt(b + gamma);
%     
%     phi(mIdx) = m(mIdx)*pi/2 + 0.5*atan(w/u) - 0.5*atan(wP/u);
%     
%     powerConf(mIdx)  = (1 + (sin(u + phi))^2/(2*w) + (sin(u - phi))^2/(2*wP)) / ...
%                 (1 + 1/(2*w) + 1/(2*wP) );
%     
%     aEff(mIdx) = corDim / powerConf(mIdx);
%     
%     
%     Ey(:, mIdx) = A*cos(ki*corDim/2 - phi)*exp(-sigma(x - corDim));
%             
% end
% 
% %% TM modes
% for mIdx = 1: modeNUM
% 
%     v(mIdx) = ( m(mIdx)*pi/2 + 1/2.*atan(n1^2/ns^2.*sqrt(b./(1 - b))) ...
%                 + 1/2.*atan(n1^2/n02(gammaIdx).*sqrt((b + gamma(gammaIdx))./(1 - b))) )./sqrt(1 - b);
%             
%     u(mIdx) = v*sqrt(1-b);
%     w(mIdx) = v*sqrt(b);
%     wP(mIdx) = v*sqrt(b + gamma);
%     
%     phi(mIdx) = m(mIdx)*pi/2 + 0.5*atan(w/u) - 0.5*atan(wP/u);
%     
%     powerConf(mIdx)  = (1 + (sin(u + phi))^2/(2*w) + (sin(u - phi))^2/(2*wP)) / ...
%                 (1 + 1/(2*w) + 1/(2*wP) );
%     
%     aEff(mIdx) = corDim / powerConf(mIdx);
%             
% end


%%
%--------------------------------------------------%
% Light Penetration Depth Script
%--------------------------------------------------%
% Parameters definition 
n2 = n0; % refrative index air
%n1 = 1.5; % refractive index glass
theta = 45*pi/180; % incidence angle [rad]
 
% Depth
d = lambda/(2*pi*sqrt(n1^2*sin(theta)^2-n2^2));
disp(['Penetration depth for lambda equal to ' num2str(lambda) ' n1 equal to ' num2str(n1)...
      'n2 equal to ' num2str(n2) ' and theta equal to ' num2str(theta) ' is ' num2str(d) ]);
 
disp('Varing wavelenght...')
lambda = (300:100:1000);
d = lambda./(2*pi*sqrt(n1^2*sin(theta)^2-n2^2));
plot(lambda, d, 'k*')
xlabel('wavelength [nm]')
ylabel('penetration depth [nm]')
title('Depth vs lambda')
 
disp('Varing the angle...') 
% clear z;
% clear lambda;
lambda = 0.6;
theta = (pi/4 : 0.05 : pi/2);
d = lambda./(2*pi*sqrt(n1^2.*sin(theta).^2-n2^2));
plot(theta, d, 'k*')
xlabel('angle [red]')
ylabel('penetration depth [nm]')
title('Depth vs theta')

% % Electric field fundamental mode E(x)
% % y = (-corDim/2 - 4: 0.1 :corDim/2 + 4);
% y = (- 5 : 0.1 : 5);
% for yIdx = 1: length(y)   
%     ex(yIdx) = Efield(y(yIdx), corDim, ki, sigma, phi, xi);
% end
% 
% 
% figure,
% plot(y, ex, 'k')
% hold on
% plot( repmat(-corDim/2,1, length(ex)), ex, 'r--')
% plot( repmat(+corDim/2,1, length(ex)), ex, 'r--')
% xlabel('y[um]')
% ylabel('Efield/Emax')
% title('TE0 mode profile along Y')
% hold off