clear

%% Import data
filePath = 'DATA/';

freqRange = 1:35;
nParticle = 2:11;

timeVec = importTime([filePath 'N2/tek0001CH1.csv'], 18, 10017);
Ts = timeVec(2) - timeVec(1);
Fs = 1/Ts;

dataCell = cell(length(freqRange),1);

for i = 1:length(freqRange)
    for j = 1:length(nParticle)
        fileName = strcat('N',num2str(nParticle(j)),'/tek',num2str(sprintf('%04d',i)),'CH1.csv');
        dataCell{i}(:,j) = importCH1([filePath fileName], 18, 10017);
    end
end


%% Optional plot time-velocity field

fPlot = 25;  % freqRange = 1:35, thus dataCell{i} --> i [Hz]

figure
surf(nParticle,timeVec,dataCell{fPlot})
shading interp
colormap gray; colormap(flipud(colormap));
xlabel('Mass \#','interpreter','latex')
ylabel('t','interpreter','latex')
zlabel('V (1000mm/s/V)','interpreter','latex')
title(['$f_{dr}$ = ' num2str(freqRange(fPlot)) ' Hz'],'interpreter','latex')

fontSize = 20;
set(gca,'YDir','normal',...
    'fontsize',fontSize,...
    'linewidth',1,...
    'TickLabelInterpreter', 'latex',...
    'layer','top',...
    'clipping','on');
set(gca,'TickLength',[0.025, 0.01])

% view(0,0)  % spatial view

%% 2D FFT

a = 33.4E-3;  % lattice spacing
x = a*nParticle;  % Space vector [m]

for ii = 1:length(freqRange)
    
    % Zero padding
    ZPk = 2^(nextpow2(size(x,2))); %Length of zero-padded k vector
    kxZP = 1/a*(-ZPk/2:ZPk/2-1)/ZPk; %k vector
    ZPf=2^(nextpow2(length(timeVec))); %Length of zero-padded freq vector
    fZP=1/Ts*(-ZPf/2:ZPf/2-1)/ZPf; %freq vector
    
    % FFT
    W = (fftshift(fft2(dataCell{ii},ZPf,ZPk)));
    
    % Store FFT ii
    surfMat(ii,:,:) = abs(W);

    disp([num2str(100*ii/length(freqRange)) '%'])
end

%% Plot maximum values of surface to obtain

maxSurf = max(surfMat);
maxSurf = squeeze(maxSurf(1,:,:));

maxSurf = maxSurf / max(max(maxSurf)); % Normalize

% plot
figure
surf(-2*pi*kxZP,fZP,maxSurf);
shading interp
colormap gray; colormap(flipud(colormap));

view(2)

ylim([0,40]); xlim([-pi/a pi/a]); 


xlabel('$k_x$ (1/m)','interpreter','latex')
ylabel('f (Hz)','interpreter','latex')
box on; grid off;

caxis([0 1])

fontSize = 30;
set(gca,'YDir','normal',...
    'fontsize',fontSize,...
    'linewidth',1,...
    'xtick',[-pi/a 0 pi/a],...
    'xticklabel',{'-$\frac{\pi}{a}$','0','$\frac{\pi}{a}$'},...
    'TickLabelInterpreter', 'latex',...
    'layer','top',...
    'clipping','on');
set(gca,'TickLength',[0.025, 0.01])
 
