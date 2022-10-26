clear all
close all
% addpath('./RcvData_PerFrame_PI_bulle06-March-2022_19-26-02')
addpath('C:\Users\Verasonics\Documents\Michael\beamformer_vincent\bf_last_version')
% addpath('C:\Users\system-01\Documents\Vantage-4.6.2-2110271004\sequence gem5\clutterFilter')
% addpath('C:\Users\system-01\Documents\Vantage-4.6.2-2110271004\sequence gem5\parameters')
% load('PI_jo__parameters.mat');
PathData = 'D:\michaelphantom\phantomstatique100_April_13_2022_12_30_57\';
FileParam = 'test_April_13_2022_12_30_57_parameters';
addpath(PathData);
load(FileParam);
%% Inputs
%% Dataraw (from 200% sampling to 100% sampling to be able to demodulate)
ntransmit=3;
time_step = tic; % Start timewatch
ibuff=5;
fileID =fopen('test_April_13_2022_12_30_57_data7.bin'); %fileID =fopen(['Frame', sprintf('_%0*i', 2, ibuff) '.bin']); % Open file
RF = fread(fileID, Inf, 'int16=>int16'); % Read file
fclose(fileID); % Close file
RF = reshape(RF, Receive(1).endSample, [], Trans.numelements);
RF = permute(RF, [1 3 2]);
%% check pulse inversion on one channel and first two events
figure(1)
subplot 211
plot(RF(:,40,2))
hold on
plot(RF(:,40,3))
hold off
title('SHI ')

subplot 212
plot(RF(:,40,1))
hold on
plot(RF(:,40,2))
hold off
title('pulse inversion ')


%% Pulse inverse on the RFs

RFpnf = RF(:,:,1:ntransmit:end) + RF(:,:,2:ntransmit:end); % keep harmonique
RFpshinf = RF(:,:,2:ntransmit:end) + RF(:,:,3:ntransmit:end); % keep fondamentale
RFnot=RF(:,:,1:ntransmit:end);
%%

brnA = 3.5e6;
brnB = 5e6;
fs = 4*Trans.frequency*1e6;
% [B, A] = butter(5, 2*[brnA brnB]./fs, 'bandpass');
[B, A] = butter(5, 2*brnA./fs, 'high');

RFpshi = filtfilt(B, A, double(RFpshinf));
RFpshi = single(RFpshi);

RFp = filtfilt(B, A, double(RFpnf));
RFp = single(RFp);

Fs = 4*(Trans.frequency*1e6);            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = size(RFp,1);             % Length of signal
f = linspace(-Fs/2,Fs/2,L);

figure(2)
subplot 311
plot(f,real(fftshift(fft(RFp(:,40,3))))) %check pour la fft pour que ca soit plus clair
hold on 
plot(f,imag(fftshift(fft(RFp(:,40,3))))) %check pour la fft pour que ca soit plus clair
hold off
title('PI sum')
subplot 312
plot(f,real(fftshift(fft(RFpshi(:,40,3)))))
hold on
plot(f,imag(fftshift(fft(RFpshi(:,40,3))))) %check pour la fft pour que ca soit plus clair
hold off

title('SHI')
subplot 313
plot(f,real(fftshift(fft(RFnot(:,40,3)))))
hold on 
plot(f,imag(fftshift(fft(RFnot(:,40,3))))) %check pour la fft pour que ca soit plus clair
hold off
title('Normal')

%% Filtrage (Passe bande)
% [b,a]   = butter(4,[0.6,.9]);
% RFp =filtfilt(b, a, double(RFp));
% RFm=filtfilt(b, a, double(RFm));
% figure;
% subplot 211
% plot(f,abs(fftshift(fft(RFp(:,40,1))))) %check pour la fft pour que ca soit plus clair
% title('harmonique')
% subplot 212
% plot(f,abs(fftshift(fft(RFm(:,40,1)))))
% title('fondamentale')
%%  Demodulation 
simu=1;
% RF100=single(RF(1:4:end-1,:,:,1))-1j.*single(RF(2:4:end,:,:,1)); 100%
% RF100=single(RF(1:2:end-1,:,:,1))+1j.*single(RF(2:2:end,:,:,1));

fs = (Receive(1).decimSampleRate*1e6)/(Receive(1).quadDecim);

 fc = Trans.frequency*1e6;
% fc = 2*P.frequency*1e6;

t = (0:size(RF,1)-1)'/fs;
Wn = 2*fc/fs-eps;

%%%%%%%%%%%%%%% Demodulation Harmonique %%%%%%%%%%%%%%%%
RF100p = bsxfun(@times, double(RFp), exp(-1i*2*pi*fc*t));
[B, A] = butter(5, Wn);
RF100p = filtfilt(B, A, RF100p)*2; %factor 2: to preserve the envelope amplitude
% RF100p = bsxfun(@times, double(RF100p), exp(1i*2*pi*(fc-P.PI_frequency)*t));

%%%%%%%%%%%%%% Demodulation SHI %%%%%%%%%%%%%%%
RF100pshi = bsxfun(@times, double(RFpshi), exp(-1i*2*pi*fc*t));
[B, A] = butter(5, Wn);
RF100pshi = filtfilt(B, A, RF100pshi)*2; %factor 2: to preserve the envelope amplitude
% RF100m = bsxfun(@times, double(RF100m), exp(1i*2*pi*(fc-P.PI_frequency)*t));

%%%%%%%%%%%%%% Demodulation pas de PI %%%%%%%%%%%%%%%%%%%%%
RF100 = bsxfun(@times, double(RFnot), exp(-1i*2*pi*fc*t));
[B, A] = butter(5, Wn);
RF100 = filtfilt(B, A, RF100)*2; %factor 2: to preserve the envelope amplitude
% RF100 = bsxfun(@times, double(RF100), exp(1i*2*pi*(fc-P.PI_frequency)*t));

%% Check Demodulation
isteer = 3;

figure(3)
subplot 321
imagesc(abs(fftshift(fft2(RF(:,:,isteer)))))
title('Radiofrequence non demodulÈe')
subplot 322
imagesc(abs(fftshift(fft2(RF100(:,:,isteer)))))
title('RF demodulÈe')

subplot 323
imagesc(abs(fftshift(fft2(RFp(:,:,isteer)))))
title('PI non demodulÈe')
subplot 324
imagesc(abs(fftshift(fft2(RF100p(:,:,isteer)))))
title('PI demodulÈe')


subplot 325
imagesc(abs(fftshift(fft2(RFpshi(:,:,isteer)))))
title('PI non demodulÈe')
subplot 326
imagesc(abs(fftshift(fft2(RF100pshi(:,:,isteer)))))
title('SHI demodulÈe')
%% Selection des canaux
%%%%%%%%%%%%%%%%%%%%%%%%%%%Harmonique%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(' done in %.2f seconds! \n', toc(time_step)); % Elapsed time
RF100p = reshape(RF100p,size(RF100p,1),size(RF100p,2),length(TX)/ntransmit,[]);

RF100p = RF100p(:, :, :, :);
figure; imagesc(abs(RF100p(:, :, 3, 1, 1))) % verifie si mes RF sont vide ou pas
title('Radiofrequency Harmonique')

[Nsample,Nchannel,Nsteer,Nframe] = size(RF100p);
RF100p = reshape(RF100p,Nsample,Nchannel/2,2,Nsteer,Nframe);
RF100p = cat(3,(RF100p(:,:,2,:,:)),RF100p); % all channels are used
RF100p = RF100p(:,:,2,:,:); % central channel only used 
% RF100 = RF100(:,:,[1,3],:,:); % outer rows used
RF100p=reshape(RF100p,Nsample,Nchannel/2,1,[]) ;% juste si je prends les cannaux cetraux 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%SHI%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(' done in %.2f seconds! \n', toc(time_step)); % Elapsed time
RF100pm = reshape(RF100pshi,size(RF100pshi,1),size(RF100pshi,2),length(TX)/ntransmit,[]);

RF100pshi = RF100pshi(:, :, :, :);
figure; imagesc(abs(RF100pshi(:, :, 3, 1, 1))) % verifie si mes RF sont vide ou pas
title('Radiofrequency Fondamentale')
RF100pshi = reshape(RF100pshi,Nsample,Nchannel/2,2,Nsteer,Nframe);
RF100pshi = cat(3,(RF100pshi(:,:,2,:,:)),RF100pshi); % all channels are used
RF100pshi = RF100pshi(:,:,2,:,:); % central channel only used 
% RF100 = RF100(:,:,[1,3],:,:); % outer rows used
RF100pshi=reshape(RF100pshi,Nsample,Nchannel/2,1,[]) ;% juste si je prends les connaux centraux
%%%%%%%%%%%%%%%%%%%%%%% Pas de PI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RF100 = reshape(RF100,Nsample,Nchannel/2,2,Nsteer,Nframe);
RF100 = cat(3,(RF100(:,:,2,:,:)),RF100); % all channels are used
RF100 = RF100(:,:,2,:,:); % central channel only used 
% RF100 = RF100(:,:,[1,3],:,:); % outer rows used
RF100=reshape(RF100,Nsample,Nchannel/2,1,[]) ;% juste si je prends les connaux cetraux

%% Geometry
%% elementpos
% xpos = Trans.spacingMm*(-((Trans.numelements-2)/4):((Trans.numelements-2)/4));
% Trans.ElementPos(:,1) = [xpos, xpos]; % same xpos for both rows
% % If desiring to treat the GEMScD as a full 2D array (type 2), uncomment the lines below
% % to offset the outer row of elements in the y direction.  To provide zero delay focus
% % at elevationFocusMM,the outer rows can be offset in the z direction as well. This
% % allows computeTXDelays to place the elevation focus at the range focal point.
% Trans.ElementPos(81:160,2) = elevOffsetMm; % add the y-offset for outer row
% Trans.ElementPos(81:160,3) = Trans.elevationFocusMm - sqrt(elevOffsetMm^2 + Trans.elevationFocusMm^2);
% If desiring to treat the GEMScD as a full 2D array (type 2), uncomment the lines below
% to offset the outer row of elements in the y direction.  To provide zero delay focus
% at elevationFocusMM,the outer rows can be offset in the z direction as well. This
% allows computeTXDelays to place the elevation focus at the range focal point.
lbd  = 1540/(Trans.frequency*1e6);
elevOffsetMm = 4.875;
Trans.ElementPos = zeros(240,3);
xpos = Trans.spacingMm*(-((Trans.numelements-2)/4):((Trans.numelements-2)/4));
Trans.ElementPos(:,1) = [xpos, xpos, xpos]; % same xpos for both rows
Trans.ElementPos((81:160)+0,2) = 0; 
Trans.ElementPos((81:160)+0,3) = 0;

%% position in y!=0 (Type 2)
% Trans.ElementPos((81:160)-80,2) = -elevOffsetMm;
% Trans.ElementPos((81:160)+80,2) = elevOffsetMm;% add the y-offset for outer row
%% position in y ==0 (Type 0 )
Trans.ElementPos((81:160)+80,2)=0;
Trans.ElementPos((81:160)-80,2)=0;
%% position in z!=0 formalisme vera -- moins bon contraste qu'avec le formalisme de jonathan
% Trans.ElementPos((81:160)-80,3) = Trans.elevationFocusMm - sqrt(elevOffsetMm^2 + Trans.elevationFocusMm^2);
% Trans.ElementPos((81:160)+80,3) = Trans.elevationFocusMm - sqrt(elevOffsetMm^2 + Trans.elevationFocusMm^2);
%% position z !=0 (Type 2) formalisme jonathan -- Better
% Trans.ElementPos((81:160)-80,3) = sqrt(elevOffsetMm^2 + Trans.elevationFocusMm^2) - Trans.elevationFocusMm ;
% Trans.ElementPos((81:160)+80,3) = sqrt(elevOffsetMm^2 + Trans.elevationFocusMm^2) - Trans.elevationFocusMm ;
%% position in z==0 (Type 0)
Trans.ElementPos((81:160)-80,3) = 0 ;
Trans.ElementPos((81:160)+80,3) = 0 ;
Trans.ElementPos = Trans.ElementPos*1e-3./lbd;
%% only central 
Geometry.posEleX=Trans.ElementPos(81:160,1).*lbd; % 1:80 first row (outer row)--81:160 central row-- 161:240 third row (outer row)
Geometry.posEleY=Trans.ElementPos(81:160,2).*lbd;
Geometry.posEleZ=Trans.ElementPos(81:160,3).*lbd;
%% all elements
% Geometry.posEleX=Trans.ElementPos(:,1).*lbd;
% Geometry.posEleY=Trans.ElementPos(:,2).*lbd;
% Geometry.posEleZ=Trans.ElementPos(:,3).*lbd;
%% aperture precision 
Geometry.apertureY = 13e-3;
%% Grid
%-------------------grid cartesienne (Michael)-----------------------------
%%
 wavelength  = lbd;
% 
% grid.subx   = 2;    % number of samples per pitch along x-axis
% grid.subz   = 8;    % number of samples per wavelength along z-axis
% 
% grid.dx     = Trans.spacingMm*1e-3/grid.subx;   % x-pixel size (m)
% grid.dz     = wavelength/grid.subz;                    % z-pixel size (m)
% grid.Nx     = ceil(128*Trans.spacing*grid.subx);            % number of pixels along x-axis
% grid.Nz     = ceil((P.endDepth - P.startDepth)*grid.subz);  % number of pixels along z-axis
% grid.Zmin   = P.startDepth*wavelength;                 % Start Depth [m]
% grid.Xmin   = -(grid.Nx-1)*grid.dx/2;           % Left corner [m]
% 
% grid.x      = (0:grid.Nx-1).*grid.dx + grid.Xmin; % x values of the reconstruction grid [m]
% grid.y      = 0;                        % y values of the reconstruction grid [m]
% grid.z      = (0:grid.Nz-1).*grid.dz + grid.Zmin; % z values of the reconstruction grid [m]
% 
% [gridBfX,gridBfZ,gridBfY] = meshgrid(grid.x,grid.z,grid.y); % Careful with convention image
% Geometry.gridBfX     =  gridBfX;         % x-grid for beamforming (m, 3D meshgrid)
% Geometry.gridBfY     =  gridBfY;         % y-grid for beamforming (m, 3D meshgrid)
% Geometry.gridBfZ     =  gridBfZ;         % z-grid for beamforming (m, 3D meshgrid)
% Geometry.posEleX     =  Trans.ElementPos(1:80,1).*wavelength;  % position of the elements along x-axis (m, vector)
% Geometry.posEleY     =  Trans.ElementPos(1:80,2).*wavelength;  % position of the elements along y-axis (m, vector)
% Geometry.posEleZ     =  Trans.ElementPos(1:80,3).*wavelength;  % position of the elements along z-axis (m, vector)
% 
% Geometry.apertureX   =  max(Geometry.posEleX)-min(Geometry.posEleX);  % aperture along x-axis (m)        
% Geometry.apertureY   =  max(Geometry.posEleY)-min(Geometry.posEleY);  % aperture along x-axis (m)  
% Geometry.apertureZ   =  max(Geometry.posEleZ)-min(Geometry.posEleZ);  % aperture along x-axis (m)  
%%----------------------------------------grid secteur (Jonathan)--------------------------------------------------------
grid.subx   = 1;    % number of samples per pitch along x-axis
grid.subz   = 8;    % number of samples per wavelength along z-axis

grid.dx     = wavelength/grid.subx;   % x-pixel size (m)
grid.dz     = wavelength/grid.subz;   % z-pixel size (m)
grid.Nx     = 128;                    % number of pixels along x-axis
grid.Nz     = ceil((P.endDepth - P.startDepth-10)*grid.subz);  % number of pixels along z-axis
grid.Zmin   = P.startDepth*wavelength;                 % Start Depth [m]
grid.Xmin   = -(grid.Nx-1)*grid.dx/2;           % Left corner [m]

grid.x      = (0:grid.Nx-1).*grid.dx + grid.Xmin; % x values of the reconstruction grid [m]
grid.y      = 0;                        % y values of the reconstruction grid [m]
grid.z      = (0:grid.Nz-1).*grid.dz + grid.Zmin; % z values of the reconstruction grid [m]

[gridBfX,gridBfZ,gridBfY] = meshgrid(grid.x,grid.z,grid.y); % Careful with convention image

Geometry.gridBfX     =  gridBfX;         % x-grid for beamforming (m, 3D meshgrid)
Geometry.gridBfY     =  gridBfY;         % y-grid for beamforming (m, 3D meshgrid)
Geometry.gridBfZ     =  gridBfZ;         % z-grid for beamforming (m, 3D meshgrid)
Geometry.posEleX     =  Trans.ElementPos(1:80,1).*wavelength;  % position of the elements along x-axis (m, vector)
Geometry.posEleY     =  Trans.ElementPos(1:80,2).*wavelength;  % position of the elements along y-axis (m, vector)
Geometry.posEleZ     =  Trans.ElementPos(1:80,3).*wavelength;  % position of the elements along z-axis (m, vector)

Geometry.apertureX   =  max(Geometry.posEleX)-min(Geometry.posEleX);  % aperture along x-axis (m)
Geometry.apertureY   =  max(Geometry.posEleY)-min(Geometry.posEleY);  % aperture along x-axis (m)
Geometry.apertureZ   =  max(Geometry.posEleZ)-min(Geometry.posEleZ);  % aperture along x-axis (m)

%%-------------------------- Pseudo spherique (Jonathan)--------------------
%% Warning X^2+Y^2+Z^2 != R
% pmig.do     = 1/4*pi/180;%.15*pi/180;
% pmig.No     = 75*4+1;
% pmig.dr     = lbd/4;
% pmig.Nr     = abs(round(Receive(1).endDepth - 50)*8); % (8 x 216 lbd-- profondeur du RCV buffer)
% %  pmig.Nr=      1000;
% pmig.df     = 1/2*pi/180;
% pmig.Nf     = 30*2+1;
% Zmin        = P.startDepth*lbd;
% 
% [Opix,Rpix,Fpix] = meshgrid(...
%     (0:pmig.No-1).*pmig.do - (pmig.No-1).*pmig.do/2,...
%     (0:pmig.Nr-1).*pmig.dr + Zmin,...
%     (0:pmig.Nf-1).*pmig.df - (pmig.Nf-1).*pmig.df/2); %
% 
% Z=Rpix.*cos(Opix).*cos(Fpix);
% X=Rpix.*sin(Opix).*cos(Fpix);
% Y=Rpix.*cos(Opix).*sin(Fpix);
% 
% Geometry.gridBfX=X;
% Geometry.gridBfY=Y;
% Geometry.gridBfZ=Z;


% %--------------------grid pyramid(comme vera)-----------------
%% Warning X^2+Y^2+Z^2 != R
% pmig.do     = 1*pi/180;%.15*pi/180;
% pmig.No     = 75;
% pmig.dr     = lbd;
% pmig.Nr     = abs(round(Receive(1).endDepth - 40)); % (8 x 216 lbd-- profondeur du RCV buffer)
% pmig.df     = 1*pi/180;
% pmig.Nf     = 35;
% Zmin        = 40*lbd;
% 
% thetax=(0:pmig.No-1).*pmig.do - (pmig.No-1).*pmig.do/2;
% thetay=(0:pmig.Nf-1).*pmig.df - (pmig.Nf-1).*pmig.df/2;
% z=(0:pmig.Nr-1).*pmig.dr + Zmin;
% [ThX,Z,ThY]=meshgrid(thetax,z,thetay);
% X=Z.*tan(ThX);
% Y=Z.*tan(ThY);
% Geometry.gridBfX=X;
% Geometry.gridBfY=Y;
% Geometry.gridBfZ=Z;

%% Grid visualisation
figure;
plot3(gridBfX(:),gridBfY(:),gridBfZ(:),'ro')
xlabel('x')
ylabel('y')
zlabel('z')
axis equal
%% Param % le if √† changer...
if (strcmp(P.Sampling, 'NS200BW') && simu==1)
     Param.samplingFrequency = (Receive(1).decimSampleRate*1e6)/ ...
                                   (1); % Sampling frequency of IQ signals (Hz)
    
elseif (strcmp(P.Sampling, 'NS200BW'))
        Param.samplingFrequency = (Receive(1).decimSampleRate*1e6)/ ...
                                   (Receive(1).quadDecim); % Sampling frequency of IQ signals (Hz)
    
else
            Param.samplingFrequency = (Receive(1).decimSampleRate*1e6)/ ...
                                   (2*Receive(1).quadDecim); % S
end

Param.demodulationFrequency    = fc;%Receive(1).demodFrequency*1e6; % Demodulation frequency (Hz)
Param.speedOfSound      = Resource.Parameters.speedOfSound; % Speed-of-sound (m/s)
Param.startTime         = ((-2*(Receive(1).startDepth- ...
                                    Trans.lensCorrection))+ ...
                                    TW(1).peak)/ ...
                                    (Trans.frequency*1e6);  % Start time (s)
Param.fnumber               = 0; %1; % F-number, YOU CAN TUNE THIS VALUE 

%% TX Pour le phantome de flux
posx=zeros(1,length(TX));
posz=zeros(1,length(TX));
posy=zeros(1,length(TX));
for n=1:length(TX)
posx(1,n)=TX(n).FocalPt(1).*lbd;
posz(1,n)=TX(n).FocalPt(3).*lbd;
end
TXb.virtualSourceX=posx(1:ntransmit:end);
TXb.virtualSourceY=posy(1:ntransmit:end);
TXb.virtualSourceZ=posz(1:ntransmit:end);
%% Beamform
% addpath('./BF-VirtualSourceGeneric-v2')
% addpath('./BF-VirtualSourceGeneric-v3')
% addpath('./BF-VirtualSourceGeneric-v4')
addpath('./bf_last_version')

Options.flagCompounding = true;
Options.flagInterpolationLanczos = false;
Options.decreaseFactorGPU = 4;
%Param.demodulationFrequency = Param.demodFrequency;

tic
[dataBf, mask] = bfVirtualSourceGeneric(RF100(:,:,:,:),Geometry,TXb,Param, Options);
[dataBfp, mask] = bfVirtualSourceGeneric(RF100p(:,:,:,:),Geometry,TXb,Param, Options);
[dataBfshi, mask] = bfVirtualSourceGeneric(RF100pshi(:,:,:,:),Geometry,TXb,Param, Options);
toc
dataBf  = squeeze(dataBf); 
dataBfp = squeeze(dataBfp);
dataBfshi = squeeze(dataBfshi);
MaxFObf = max(max(sum(abs(dataBfshi),3)));
MaxPIbf = max(max(sum(abs(dataBfp),3)));
%%
% RF100PI = 
% RF100PI = reshape(RF100,size(RF100,1),size(RF100,2),size(RF100,3),2,[]);
% RF100PIp = RF100PI(:,:,:,1,:) + RF100PI(:,:,:,2,:);
% RF100PIp = permute(RF100PIp,[1 2 3 5 4]);
% RF100PIm = RF100PI(:,:,:,1,:) - RF100PI(:,:,:,2,:);
% RF100PIm = permute(RF100PIm,[1 2 3 5 4]);
% RFPI100    = permute(RFPI100(:,81:160,:),[1 2 4 3]);
% tic
% [dataBfPI, mask] = bfVirtualSourceGeneric(RFPI100,Geometry,TXb(1:2:end),Param, Options);
% [dataBfPIm, mask] = bfVirtualSourceGeneric(RF100PIm,Geometry,TXb(1:2:end),Param, Options);
% toc
%%
%% display Bmode image
%  figure(9); 
% colormap(gray(512));
% subplot 121
% set(pcolor (X(:,:,31),...
%             Z(:,:,31),...
%             20*log10(rescale(abs(dataBf(:,:,31))))),'Edgecolor','none');
% axis equal tight ij
%  title('Plan Y = 0')
% caxis([-60 0])
% xlabel('mm')
% ylabel('mm')
% 
% subplot 122
% set(pcolor (squeeze(Y(:,151,:)),...
%             squeeze(Z(:,151,:)),...
%             squeeze(20*log10(rescale(abs(dataBf(:,38,:)))))),'Edgecolor','none');
% axis equal tight ij
% title('Plan X = 0')
% caxis([-60 0])
% xlabel('mm')
% ylabel('mm')
% dataBf=squeeze(dataBf);
% dataBfPI=squeeze(dataBfPI);

%% Filtrage Fondamental
% dataBf2 = real(dataBf.*exp(1j*pi*Param.demodulationFrequency*gridBfZ));
% [b,a]   = butter(4,[0.6,.9]);
% dataBf2 = single(filtfilt(b,a,double(dataBf2)));

%% Demodulation
% dataBf2 = dataBf2.*exp(-1j*pi*Param.demodulationFrequency*gridBfZ);
% [b,a]   = butter(4,.5);
% dataBf2 = single(filtfilt(b,a,double(dataBf2)));

%%
figure
subplot 221
plot(abs(fftshift(fft(dataBfshi(:,128,1)))))
title('spectre donnÈe beamform shi')

subplot 222
plot(abs(fftshift(fft(dataBfp(:,128,1)))))
title('spectre donnÈe beamform harmonique')

%% Fondamental Supression in dataBfp and dataBfshi
% mask     = zeros(size(dataBfp));
% mask(size(mask,1)/2+1:end,:,:) = 1;
% mask = convn(mask,hanning(11),'same');
% 
% dataBfpf = fftshift(fft2(dataBfp)).*(mask); %masque en fourier
% dataBfpf = ifft2(ifftshift(dataBfpf));
% 
% dataBfshif = fftshift(fft2(dataBfshi)).*(mask); %masque en fourier
% dataBfshif = ifft2(ifftshift(dataBfshif));
%% recentrer le spectre
% dataBfpf2 = dataBfpf.*exp(-2j*pi*(Param.demodulationFrequency-P.PI_frequency*1e6)*gridBfZ);
% dataBfpf=dataBfpf2;  
%%
% figure
% subplot 141
% imagesc(abs(fftshift(fft2(dataBf(:,:,1)))))
% 
% subplot 142
% imagesc(abs(fftshift(fft2(dataBfp(:,:,1)))))
% 
% subplot 143
% imagesc(abs(fftshift(fft2(dataBfpf(:,:,1)))))
% 
% subplot 144
% imagesc(abs(fftshift(fft2(dataBfpf2(:,:,1)))))
%%
% figure
% subplot 141
% plot(abs(fftshift(fft(dataBfshi(:,128,1,1)))))
% title('spectre donnÈe beamform SHI')
% 
% subplot 142
% plot(abs(fftshift(fft(dataBfshif(:,128,1,1)))))
% title('spectre donnÈe beamform SHI')
% 
% subplot 143
% plot(abs(fftshift(fft(dataBfp(:,128,2,1)))))
% title('spectre donnÈe beamform PI')
% 
% subplot 144
% plot(abs(fftshift(fft(dataBfpf(:,128,2,1)))))
% title('spectre donnÈe beamform PI filtre')


%%
% figure
% subplot 211
% imagesc(squeeze(abs(fftshift(fft(dataBfPIm(:,64,:),[],3),3))))
% caxis([0 1]*1e4)
% subplot 212
% imagesc(squeeze(abs(fftshift(fft(dataBfPIp(:,64,:),[],3),3))))
% caxis([0 1]*1e4)
% 
% %%
% figure
% subplot 221
% plot(abs(fftshift(fft(RF100(:,40,1,1)))))
% subplot 222
% plot(abs(fftshift(fft(RF100(:,40,1,2)))))
% 
% subplot 223
% plot(abs(fftshift(fft(sum(RF100(:,40,1,1:2),4)))))
% subplot 224
% plot(abs(fftshift(fft(diff(RF100(:,40,1,1:2),[],4)))))
% 
% %%
% figure
% subplot 121
% imagesc(abs(fftshift(fft2(RF100(:,:,1,6)))))
% 
% subplot 122
% imagesc(abs(fftshift(fft2(dataBf(:,:,1,6)))))

%%
BmodeF  = 20.*log10(abs((dataBf)));
BmodeF  = BmodeF - max(BmodeF(:));

BmodeSHI = 20.*log10(abs((dataBfshi)));
BmodeSHI = BmodeSHI - max(BmodeSHI(:));

BmodePI  = 20.*log10(abs((dataBfp)));
BmodePI  = BmodePI - max(BmodePI(:));


% maxA = max(abs(BmodeF(:)));
% maxB = max(abs(BmodeSHI(:)));
% maxC = max(abs(BmodePI(:)));
% 
% % maxBmode = max([maxA, maxB, maxC]);
% 
% % BmodeF  = BmodeF - maxBmode;
% % BmodeSHI = BmodeSHI - maxBmode;
% % BmodePI  = BmodePI - maxBmode;
% 
% BmodeF  = BmodeF - max(BmodeF(:));
% BmodeSHI = BmodeSHI - max(BmodeSHI(:));
% BmodePI  = BmodePI - max(BmodePI(:));


% BmodeF  = BmodeF - max(BmodeF(:));
% BmodeSHI = BmodeSHI - max(BmodeSHI(:));
% BmodeSHIf  = BmodeSHIf - max(BmodeSHIf(:));
% BmodeHf  = BmodeHf - max(BmodeHf(:));


figure(133)
% for i=1:Nframe
for i=1:size(BmodeF, 3)
colormap(gray(512));

% set(pcolor (X(:,:,17),...
%             Z(:,:,17),...
%             20*log10(rescale(abs(dataBf(:,:,30,i))))),'Edgecolor','none');

subplot 131
pcolor(Geometry.gridBfX.*1e3,...
       Geometry.gridBfZ.*1e3,...
        BmodeF(:,:,i))
shading interp    
axis equal tight ij
 title('Plan Y = 0')
caxis([-60 0])
colorbar
xlabel('mm')
ylabel('mm')
title('Normal')

subplot 132
pcolor(Geometry.gridBfX.*1e3,...
       Geometry.gridBfZ.*1e3,...
        BmodeSHI(:,:,i))
shading interp    
axis equal tight ij
 title('Plan Y = 0')
caxis([-60 0])
colorbar
xlabel('mm')
ylabel('mm')
title('SHI')

subplot 133
pcolor(Geometry.gridBfX.*1e3,...
       Geometry.gridBfZ.*1e3,...
       BmodePI(:,:,i))
shading interp    
axis equal tight ij
 title('Plan Y = 0')
caxis([-60 0])
colorbar
xlabel('mm')
ylabel('mm')
title('PI ')
drawnow
end

% B = abs(dataBf(:, :, 31));
% thrB = quantile(B(:), 0.99);
% B(B>=thrB) = thrB; 
% B = rescale(B);
% 
% figure(10); 
% colormap(gray(512));
% set(pcolor (X(:,:,31),...
%             Z(:,:,31),...
%             20*log10(B)),'Edgecolor','none');
% axis equal tight ij
%  title('Plan Y = 0')
% caxis([-40 0])

%%
%% SVD clutter A revoir
disp('SVD Clutter')
Ncut = [2 size(dataBfshi,3)];
wcut = Ncut./size(dataBfshi,3);

dataBf   = dataBf./max(abs(dataBf(:))); 
[dataBf_f, eig_val] = SVDClutterFilter(dataBf,[],wcut);
dataBf_f   = dataBf_f./max(abs(dataBf_f(:))); 

dataBfshi   = dataBfshi./max(abs(dataBfshi(:))); 
[dataBfm_f, eig_val] = SVDClutterFilter(dataBfshi,[],wcut);
dataBfm_f   = dataBfm_f./max(abs(dataBfm_f(:))); 

dataBfp   = dataBfp./max(abs(dataBfp(:))); 
[dataBfp_f, eig_val] = SVDClutterFilter(dataBfp,[],wcut);
dataBfp_f   = dataBfp_f./max(abs(dataBfp_f(:))); 

dataBfpf   = dataBfpf./max(abs(dataBfpf(:))); 
[dataBfpf_f, eig_val] = SVDClutterFilter(dataBfpf,[],wcut);
dataBfpf_f   = dataBfpf_f./max(abs(dataBfpf_f(:))); 

% IQBfPI_f            = IQBfPI_f./sqrt(mask);
% IQBfPI              = squeeze(IQBfPI./mask);

%%
comp = .25

figure(10)
for it = 1:size(dataBfm_f,3)
    subplot 241
    imagesc((abs(dataBf(:,:,it))),[0 1]*comp)
    title('Fondamental')
    ylabel('wo SVD')
    
    subplot 242
    imagesc((abs(dataBfshi(:,:,it))),[0 1]*comp)
    title('Fondamental PI')
    
    subplot 243
    imagesc((abs(dataBfp(:,:,it))),[0 1]*comp)
    title('Harmonic PI')
    
    subplot 244
    imagesc((abs(dataBfpf(:,:,it))),[0 1]*comp)
    title('Harmonic PI filter')
    
    subplot 245
    imagesc((abs(dataBf_f(:,:,it))),[0 1]*comp)
    ylabel('w SVD')
    
    subplot 246
    imagesc((abs(dataBfm_f(:,:,it))),[0 1]*comp)
    
    subplot 247
    imagesc((abs(dataBfp_f(:,:,it))),[0 1]*comp)
    
    subplot 248
    imagesc((abs(dataBfpf_f(:,:,it))),[0 1]*comp)
    colormap gray
    pause(0.1)
end
%% echantillionage donn√©e brute
% figure(3);
% imagesc(abs(fftshift(fft2(RF100(10:end-10,10:end-10,3))))) % Last index indicate steer angle
% %% echantillonage donn√©e beamform
% figure(4);
% subplot 121
% imagesc(abs(fftshift(fft2(dataBf(10:end-10,10:end-10,31))))) % fixed index for plan in y (grid dep)
% subplot 122
% imagesc(abs(fftshift(fft2(squeeze(dataBf(10:end-10,151,10:end-10)))))) %fixed index for plan in x (grid dep)

%% Contrast assessment
%% Fonction contraste pour une pour une zone prendre premier click centre de la zone d'interet deuxieme click rayon troisieme click centre du backround 
Ncyst=1;
colormap(gray(512));
figure(551); 
pcolor(Geometry.gridBfX,...
       Geometry.gridBfZ,...
        BmodeF(:,:,50))
shading interp    
axis equal tight ij
 title('Plan Y = 0')
caxis([-60 0])
xlabel('mm')
ylabel('mm')
title('Fondamental')

 xInput=[];
 yInput=[];
while 0<1
    [xinput,yinput,b] = ginput(1);
    if isempty(b) 
        break; %if enter quit
    elseif b==57 % if 9 zoom out
        ax = axis; width=ax(2)-ax(1); height=ax(4)-ax(3);
        axis([xinput-width/2 xinput+width/2 yinput-height/2 yinput+height/2]);
        zoom(1/2);
    elseif b==56 % if 8 zoom in
        ax = axis; width=ax(2)-ax(1); height=ax(4)-ax(3);
        axis([xinput-width/2 xinput+width/2 yinput-height/2 yinput+height/2]);
        zoom(2);    
    else
        xInput=[xInput;xinput];
        yInput=[yInput;yinput];
    end
end
%% if coordinates are set 
% load('xInput')
% load('yInput')
%%
rin=zeros(1,Ncyst);
for i =1:Ncyst
rin(i)=sqrt((xInput(3*i-2)-xInput(3*i-1)).^2+(yInput(3*i-2)-yInput(3*i-1)).^2);
end
Nseuil = [1,2,3,4,5,6];
CNR = zeros(1,size(Nseuil,2));
CNRm = zeros(1,size(Nseuil,2));
CNRp = zeros(1,size(Nseuil,2));
CNRpf = zeros(1,size(Nseuil,2));
CNR_f = zeros(1,size(Nseuil,2));
CNRm_f = zeros(1,size(Nseuil,2));
CNRp_f = zeros(1,size(Nseuil,2));
CNRpf_f = zeros(1,size(Nseuil,2));

for i=1:size(Nseuil,2)
Ncut = [Nseuil(i) size(dataBfshi,3)];
wcut = Ncut./size(dataBfshi,3);

dataBf   = dataBf./max(abs(dataBf(:))); 
[dataBf_f, eig_val] = SVDClutterFilter(dataBf,[],wcut);
dataBf_f   = dataBf_f./max(abs(dataBf_f(:))); 

dataBfshi   = dataBfshi./max(abs(dataBfshi(:))); 
[dataBfm_f, eig_val] = SVDClutterFilter(dataBfshi,[],wcut);
dataBfm_f   = dataBfm_f./max(abs(dataBfm_f(:))); 

dataBfp   = dataBfp./max(abs(dataBfp(:))); 
[dataBfp_f, eig_val] = SVDClutterFilter(dataBfp,[],wcut);
dataBfp_f   = dataBfp_f./max(abs(dataBfp_f(:))); 

dataBfpf   = dataBfpf./max(abs(dataBfpf(:))); 
[dataBfpf_f, eig_val] = SVDClutterFilter(dataBfpf,[],wcut);

dataBfpf_f   = dataBfpf_f./max(abs(dataBfpf_f(:)));

CNR(1,i)=CNR_michael(abs(dataBfp(:,:,10)), Geometry.gridBfX, Geometry.gridBfZ, rin, Ncyst, xInput, yInput, xInput, yInput);
CNRm(1,i)=CNR_michael(abs(dataBfshi(:,:,10)), Geometry.gridBfX, Geometry.gridBfZ, rin, Ncyst, xInput, yInput, xInput, yInput);
CNRp(1,i)=CNR_michael(abs(dataBfp(:,:,10)), Geometry.gridBfX, Geometry.gridBfZ, rin, Ncyst, xInput, yInput, xInput, yInput);
CNRpf(1,i)=CNR_michael(abs(dataBfpf(:,:,10)), Geometry.gridBfX, Geometry.gridBfZ, rin, Ncyst, xInput, yInput, xInput, yInput);
CNR_f(1,i)=CNR_michael(abs(dataBf_f(:,:,10)), Geometry.gridBfX, Geometry.gridBfZ, rin, Ncyst, xInput, yInput, xInput, yInput);
CNRm_f(1,i)=CNR_michael(abs(dataBfm_f(:,:,10)), Geometry.gridBfX, Geometry.gridBfZ, rin, Ncyst, xInput, yInput, xInput, yInput);
CNRp_f(1,i)=CNR_michael(abs(dataBfp_f(:,:,10)), Geometry.gridBfX, Geometry.gridBfZ, rin, Ncyst, xInput, yInput, xInput, yInput);
CNRpf_f(1,i)=CNR_michael(abs(dataBfpf_f(:,:,10)), Geometry.gridBfX, Geometry.gridBfZ, rin, Ncyst, xInput, yInput, xInput, yInput);
end

figure;
plot(Nseuil,CNR_f)
hold on
plot(Nseuil,CNRm_f)
hold on 
plot(Nseuil,CNRp_f)
hold on
plot(Nseuil,CNRpf_f)
hold on 
legend('no PI', 'RF soustrait PI','RF somme PI','Harmonique filtree')
title('contraste vs seuil de SVD ')
%% CNR SHI
CNR=CNR_michael(abs(dataBf(:,:,10)), Geometry.gridBfX, Geometry.gridBfZ, rin, Ncyst, xInput, yInput, xInput, yInput);
CNRSHI=CNR_michael(abs(dataBfshi(:,:,10)), Geometry.gridBfX, Geometry.gridBfZ, rin, Ncyst, xInput, yInput, xInput, yInput);
CNRPI=CNR_michael(abs(dataBfp(:,:,10)), Geometry.gridBfX, Geometry.gridBfZ, rin, Ncyst, xInput, yInput, xInput, yInput);

