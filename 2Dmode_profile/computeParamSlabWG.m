%-------------------------------------------------------------------%
%                Script for computing slab parameters               %
%-------------------------------------------------------------------%

%
%                       ________________corDimY__
%    __________________|                         |__________________no
%   |                  |  corDimX                |                  n1
%   |__________________|_________________________|_____________________
%                                                                   ns
%---------------------------------------------------------------------%
%
% Anna Archetti - EPFL 2017
%
%---------------------------------------------------------------------%

%% WG given parameters mode m = 0

n1 = 2.038;     % core Si3N4 3.48 2.038 (glass 1.5)
n0 = 1.47;     % cladding Sample Media Water 1.38
ns = 1.47;     % substrate SiO2 1.47
lambda = 0.647;   % um 1.5
NAObj = 0.55;
Fobj = 4;       % [mm]
beamD = 6.5; % [mm] % laser beam diameter at the back obj

% Mode
m = 1;

% RMS width sigma of the Gaussian approximation to the Airy disk 
% is about one-third the Airy disk radius 
sigmaLambda = 0.21*lambda/NAObj; % [um]

% Sparrow criterion
% sigmaLambda = 0.47*lambda/NAObj;  % (Abbe 0.5, Rayleigh 0.61)


beamR = 2*lambda*Fobj/(pi*beamD); % [um] w = spot size/focused beam radius where
                                    %      amplitude drop of 1/e 0.37 
                                    %      E(x) field z=0
% sigmaLambda = 0.42*lambda*Fobj/beamD; % [um]                       
% sigmaDiff = lambda/(2*NAObj); [um]

% Beam diameter
FWHM = 2.35*sigmaLambda;           % [um]
fourSigma = 4*sigmaLambda;         % [um]
fullBeamWaist = 1.67*FWHM*sqrt(2); % [um] 2w 1/e2  0.14 86% E field

%corDimY = (0.2: 0.1: 1.600);   
corDimY = 0.150;  % um 2a 2 0.01

% Number of modes along y
nMod = 2*corDimY/lambda;

for corDimIdx = 1: length(corDimY)

    corDimX = 0.15;     % um 2a 2    10
    modeNUM = 1;

    k = 2*pi/lambda;
    gamma = (ns^2 - n0^2)/(n1^2 - ns^2);

    v = k*corDimY(corDimIdx)/2*sqrt(n1^2 - ns^2);

    % DEBUG
    bx = (0:0.1:1);
    by = atan(sqrt(bx./(1 - bx))) +  ...
                        atan(sqrt((bx + gamma)./(1 - bx))) - ...
                        2.*v.*sqrt(1 - bx);                
    figure, plot(bx,by), hold on % END debug       

    % Define the dispersion b-v function(v,b,m,gamma)
    myFunc = @(br) atan(sqrt(br./(1 - br))) +  ...
                        atan(sqrt((br + gamma)./(1 - br))) - ...
                        2.*v.*sqrt(1 - br);
    % Get the b number from f(v,b,m,gamma) function roots               
    [b, byRoot] = findCFRootBisection(myFunc, [0.001, 1], 0.001, 10^10);   

    % Get the effective index
    ne = sqrt(b*(n1^2 - ns^2) + ns^2);
    
    neV(corDimIdx) = ne;

    beta = ne*k;
    u = v*sqrt(1 - b);
    w =  v*sqrt(b);
    wP = v*sqrt(b + gamma);

    % transverse propagation constants
    ki = u/(corDimY(corDimIdx)/2);        %sqrt(k^2*n1^2 - betak^2);
    sigma = w/(corDimY(corDimIdx)/2);     %sqrt(betak^2 - k^2*n0^2);
    xi = wP/(corDimY(corDimIdx)/2);       %sqrt(betak^2 - k^2*ns^2);

    % phase
    phi = 0.5*atan(w/u) - 0.5*atan(wP/u);

    % penetrationd depth
    d = lambda/(2*pi)/sqrt(ne^2 - n0^2); 

    % power confinment 
    powerConf = (1 + (sin(u + phi))^2/(2*w) + (sin(u - phi))^2/(2*wP)) / ...
                    (1 + 1/(2*w) + 1/(2*wP) );
    % effective width
    he(corDimIdx) = corDimY(corDimIdx) + 1/sigma + 1/xi;

    % penetration depth 
    dHe = (he - corDimY(corDimIdx))/2;

    % numerical aperture
    %NA = n1*sqrt(2*(n1 - n0));
    NA = sqrt(n1^2 - n0^2);

    % v cut off for fudamental mode (Okamoto ch. 2  pag 20)
    vC = m*pi/2 + 0.5*atan(sqrt(gamma)) ;

    % lambda cut off
    lambdaC = pi*corDimY(corDimIdx)*NA/v;

    % thickness cut off
    dC = lambda*v/(pi*NA);

    % coupling efficiency
    NAfiber = 1;
    % eta = (NAfiber/NA)^2;


    % From Hunsperger book
    % End-but coupling (eq. 7.4) 
%     nL = 1.468; % core refractive index laser
%     laserDim = 1;
%     m = 0;
%     etaB(corDimIdx) = 64/((m+1)^2*pi^2) * (nL*n1)/(nL+n1)^2 ...
%         * (cos(pi*corDimY(corDimIdx)/(2*laserDim)))^2 * 1/(1 - (corDimY(corDimIdx)/((m+1)*laserDim))^2 )^2 * ...
%         corDimY(corDimIdx) * (cos(m*pi/2))^2 / laserDim;

    y = (-corDimX/2 - 2: 0.01 :corDimX/2 + 2);
    for xIdx = 1: length(y)

        % Electric field distribution: TE mode profile 
        if y(xIdx) > corDimY(corDimIdx)/2 % top cladding
            Ey(xIdx) = cos(ki*corDimY(corDimIdx)/2 - phi)*exp(-sigma*(y(xIdx) - corDimY(corDimIdx)/2));
            
        elseif y(xIdx) >= - corDimY(corDimIdx)/2 && y(xIdx) <= corDimY(corDimIdx)/2 % core
            Ey(xIdx) = cos(ki*y(xIdx) - phi);
           
        else % lower cladding
            Ey(xIdx) = cos(ki*corDimY(corDimIdx)/2 + phi)*exp(xi*(y(xIdx) + corDimY(corDimIdx)/2));
            
        end
    %     By(xIdx) = exp( -y(xIdx)^2/(2*sigmaLambda^2) ) / sqrt(2*pi*sigmaLambda^2);
        
        By(xIdx) = exp( -y(xIdx)^2/(beamR^2) ); % Electric field input beam
        IBy(xIdx) = ( By(xIdx) )^2; % Intensity profile input beam
        IEy(xIdx) = Ey(xIdx)^2; 

    end

    close 
    % From Hunsperger book
    % Direct focusing (eq. 7.3)
%     etaD(corDimIdx) = (trapz(By.*Ey))^2 / ( trapz(By.^2)*trapz(Ey.^2) );
%     eta(corDimIdx) = (ne/n1)^2 / (nL/1.6)^2;

    %figure,
    plot( y, IEy, 'k')
    hold on
    plot( repmat(-corDimY(corDimIdx)/2,1, length(IEy)), IEy, 'r--')
    plot( repmat(+corDimY(corDimIdx)/2,1, length(IEy)), IEy, 'r--')
%     plot( y, By, 'b--')
    plot( y, IBy, 'g--')
    xlabel('y[um]')
    ylabel('I_{field}/I_{max}')
    title(['Intesnsity TE0 mode profile - core ' num2str(corDimY(corDimIdx)) 'um - ref idx core ' num2str(n1)])
    legend('Intensity TE0 mode profile','core edge','core edge','Intensity beam - Gaussian approxiamtion Airy pattern')
   
    
    %figure,
    plot( y, Ey, 'k')
    hold on
    plot( repmat(-corDimY(corDimIdx)/2,1, length(Ey)), Ey, 'r--')
    plot( repmat(+corDimY(corDimIdx)/2,1, length(Ey)), Ey, 'r--')
    plot( y, By, 'b--')
    %plot( y, BIy, 'g--')
    xlabel('y[um]')
    ylabel('E_{field}/E_{max}')
    title(['TE0 mode profile - core ' num2str(corDimY(corDimIdx)) 'um - ref idx core ' num2str(n1)])
    %legend('Intensity TE0 mode profile','core edge','core edge','Intensity beam - Gaussian approxiamtion Airy pattern')
    legend('TE0 mode profile','core edge','core edge','E beam - Gaussian approxiamtion Airy pattern')

    % Cutoff width
    cutOffWidth = lambda*vC/(2*pi*n1*sqrt(2*(n1-n0)/n1));
end

% figure,
% plot(corDimY, eta, '*k')
% xlabel('CorDim [um]')
% ylabel('Eta [%]')
% title('NA coupling efficinecy')
% 
% figure,
% plot(corDimY, -10.*log10(eta), '*k')
% xlabel('CorDim [um]')
% ylabel('Eta [dB]')
% title('NA coupling efficinecy')
% 
% figure,
% plot(corDimY, etaD, '*k')
% xlabel('CorDim [um]')
% ylabel('Eta [%]')
% title('Direct focusing coupling efficinecy')
% 
% figure,
% plot(corDimY, -10.*log10(etaD), '*k')
% xlabel('CorDim [um]')
% ylabel('Eta [dB]')
% title('Direct focusing coupling efficinecy')
% 
% figure,
% plot(corDimY, etaB, '*k')
% xlabel('CorDim [um]')
% ylabel('Eta [%]')
% title('Fiber end-but coupling efficiency')
% 
% figure,
% plot(corDimY, -10.*log10(etaB), '*k')
% xlabel('CorDim [um]')
% ylabel('Eta [dB]')
% title('Fiber end-but coupling efficiency')


%%
%----------------------------------------------------------
% X direction
%---------------------------------------------------------

n1 = 2.038;     % core Si3N4 3.48 2.038 (glass 1.5)
n0 = 1.47;     % cladding Sample Media Water 1.33
ns = 1.47;     % substrate SiO2 1.46

%corDimX = 5;     % um 2a 0.2
modeNUM = 1;

k = 2*pi/lambda;
gamma = (ns^2 - n0^2)/(n1^2 - ns^2);

v = k*corDimX/2*sqrt(n1^2 - ns^2);

% DEBUG
bx = (0:0.1:1);
by = atan(sqrt(bx./(1 - bx))) +  ...
                    atan(sqrt((bx + gamma)./(1 - bx))) - ...
                    2.*v.*sqrt(1 - bx);                
figure, plot(bx,by), hold on % END debug       

% Define the dispersion b-v function(v,b,m,gamma)
myFunc = @(br) atan(sqrt(br./(1 - br))) +  ...
                    atan(sqrt((br + gamma)./(1 - br))) - ...
                    2.*v.*sqrt(1 - br);
% Get the b number from f(v,b,m,gamma) function roots               
[b, byRoot] = findCFRootBisection(myFunc, [0.001, 1], 0.001, 10^8);   

% Get the effective index
ne = sqrt(b*(n1^2 - ns^2) + ns^2);

beta = ne*k;
u = v*sqrt(1 - b);
w =  v*sqrt(b);
wP = v*sqrt(b + gamma);

% transverse propagation constants
ki = u/(corDimX/2);        %sqrt(k^2*n1^2 - betak^2);
sigma = w/(corDimX/2);     %sqrt(betak^2 - k^2*n0^2);
xi = wP/(corDimX/2);       %sqrt(betak^2 - k^2*ns^2);

% phase
phi = 0.5*atan(w/u) - 0.5*atan(wP/u);
    
x = (-corDimX /2 - 2 : 0.01 : corDimX/2 + 2);
for xIdx = 1: length(x)

    % Electric field distribution: TE mode profile in the plane perpendicular
    % to the propagation direction
    if x(xIdx) > corDimX/2
        Ex(xIdx) = cos(ki*corDimX/2 - phi)*exp(-sigma*(x(xIdx) - corDimX/2));
    elseif x(xIdx) >= - corDimX/2 && x(xIdx) <= corDimX/2
        Ex(xIdx) = cos(ki*x(xIdx) - phi);
    else
        Ex(xIdx) = cos(ki*corDimX/2 + phi)*exp(xi*(x(xIdx) + corDimX/2));
    end
    Bx(xIdx) = exp( -y(xIdx)^2/(beamR^2) ); % Electric field input beam
end


% Rectangular grid in 2-D 
% [X, Y] = meshgrid (-4 : +4 , -2 : +2);
% Z = Efield(X)*Efield(Y);
pixSize = ( y(min(length(Ex), length(Ey))) - y(1) )/( min(length(Ex), length(Ey)) - 1); % um
% Exy = Ex( 1:min(length(Ex), length(Ey)) )'*Ey( 1:min(length(Ex), length(Ey)) );
Exy = Ey( 1:min(length(Ex), length(Ey)) )'* Ex( 1:min(length(Ex), length(Ey)) );
Bxy = By( 1:min(length(Bx), length(By)) )'* Bx( 1:min(length(Bx), length(By)) );




% ne = beta/k;
% b = (ne^2 - ns^2)/(n1^2 - ns^2);
% figure, imagesc(Exy)
figure, image( flipud(Exy),'CDataMapping','scaled'), colormap('jet')


yticklabels = (y(1) : 0.2 : y(min(length(Ex), length(Ey))));
yticks = linspace(1, min(length(Ex), length(Ey)), numel(yticklabels));
set(gca, 'YTick', yticks, 'YTickLabel',flipud(yticklabels(:)))

xticklabels = (x(1) : 0.5 : x(min(length(Ex), length(Ey))));
xticks = linspace(1, min(length(Ex), length(Ey)), numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)  
hold on

% ridgeT = 1; %um
% ridgeH = 1.2; %um
% ridgeWidth = 3; % um

% ridgeT = corDimY/2 - (corDimY - corDimX); %um
ridgeH = corDimY; %um
ridgeWidth = corDimX; % um

xmax = get(gca, 'xLim');
ymax = get(gca, 'yLim');

% orizXEdge1 =  xmax(2)/2 + y(1 : min(length(Ex), length(Ey)))./pixSize + 0.5;
%orizYEdge1 =  repmat( (ymax(2)/2 - ridgeT/pixSize + 0.5), 1, min(length(Ex), length(Ey)));

% orizXEdge2 =  xmax(2)/2 + y(1 : min(length(ex), length(ey)))./pixSize + 0.5;
% orizYEdge2 =  repmat( (ymax(2)/2 + ridgeT/2/pixSize + 0.5), 1, min(length(ex), length(ey)));


orizXEdge3 =  xmax(2)/2 + y(1 : min(length(Ex), length(Ey)))./pixSize + 0.5;
orizYEdge3 =  repmat( (ymax(2)/2 - ridgeH/2/pixSize + 0.5), 1, min(length(Ex), length(Ey)));

orizXEdge4 =  xmax(2)/2 + y(1 : min(length(Ex), length(Ey)))./pixSize + 0.5;
orizYEdge4 =  repmat( (ymax(2)/2 + ridgeH/2/pixSize + 0.5), 1, min(length(Ex), length(Ey)));

vertYEdge1 =  xmax(2)/2 + y(1 : min(length(Ex), length(Ey)))./pixSize + 0.5;
vertXEdge1 =  repmat( (ymax(2)/2 - ridgeWidth/2/pixSize + 0.5), 1, min(length(Ex), length(Ey)));

vertYEdge2 =  xmax(2)/2 + y(1 : min(length(Ex), length(Ey)))./pixSize + 0.5;
vertXEdge2 =  repmat( (ymax(2)/2 + ridgeWidth/2/pixSize + 0.5), 1, min(length(Ex), length(Ey)));



%plot(orizXEdge1, orizYEdge1, 'r--', 'LineWidth',2)
% plot(orizXEdge2, orizYEdge2, 'r--')
plot(orizXEdge3, orizYEdge3, 'r--', 'LineWidth',2)
plot(orizXEdge4, orizYEdge4, 'r--', 'LineWidth',2)
plot(vertXEdge1, vertYEdge1, 'r--', 'LineWidth',2)
plot(vertXEdge2, vertYEdge2, 'r--', 'LineWidth',2)

axis equal
xlabel('x [um]')
ylabel('y [um]')
title('WG Electric field propagation XY')

set(gca,'color','none')


%--------------------------------------------------------------------------
% Plot Beam 2D profile
%--------------------------------------------------------------------------

figure, image( flipud(Bxy),'CDataMapping','scaled'), colormap('jet')

yticklabels = (y(1) : 0.2 : y(min(length(Bx), length(By))));
yticks = linspace(1, min(length(Bx), length(By)), numel(yticklabels));
set(gca, 'YTick', yticks, 'YTickLabel',flipud(yticklabels(:)))

xticklabels = (x(1) : 0.5 : x(min(length(Bx), length(By))));
xticks = linspace(1, min(length(Bx), length(By)), numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)  
hold on

plot(orizXEdge3, orizYEdge3, 'r--', 'LineWidth',2)
plot(orizXEdge4, orizYEdge4, 'r--', 'LineWidth',2)
plot(vertXEdge1, vertYEdge1, 'r--', 'LineWidth',2)
plot(vertXEdge2, vertYEdge2, 'r--', 'LineWidth',2)

axis equal
xlabel('x [um]')
ylabel('y [um]')
title('Beam Electric field propagation XY')

set(gca,'color','none')

%%
% Propagation wavelenght of the trasveral mode
lambdaP = 2*pi/beta;
z = (1:1:200);
f = Exy.*cos(2.*pi/lambdaP.*z);
[X,Y] = meshgrid(1:0.5:10,1:20);
surf(X, Y, Z)

%
figure, mesh(flipud(Exy)),
hold on,
image( flipud(Exy),'CDataMapping','scaled'), colormap('jet')
xticklabels = (y(1) : 0.5 : y(min(length(Ex), length(Ey))));
xticks = linspace(1, min(length(ex), length(Ey)), numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)

yticklabels = (x(1) : 0.5 : x(min(length(Ex), length(Ey))));
yticks = linspace(1, min(length(Ex), length(Ey)), numel(yticklabels));
set(gca, 'YTick', yticks, 'YTickLabel', flipud(yticklabels(:)))  
xlabel('x [um]')
ylabel('y [um]')
zlabel('Normalized Electric field profile')
title('Electric field propagation XY')

% figure, mesh(flipud(Bxy)),
% hold on,
% image( flipud(Bxy),'CDataMapping','scaled'), colormap('jet')
% xticklabels = (y(1) : 0.5 : y(min(length(Bx), length(By))));
% xticks = linspace(1, min(length(ex), length(Ey)), numel(xticklabels));
% set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
% 
% yticklabels = (x(1) : 0.5 : x(min(length(Bx), length(By))));
% yticks = linspace(1, min(length(Bx), length(By)), numel(yticklabels));
% set(gca, 'YTick', yticks, 'YTickLabel', flipud(yticklabels(:)))  
% xlabel('x [um]')
% ylabel('y [um]')
% zlabel('Input Beam Normalized Electric field profile')
% title('Beam Electric field propagation XY')

