%% script_post_proc_PI_2022_09_26
clc
clear all

%%
% addpath(genpath('C:\Users\system-01\Desktop\MICHAELM\jopo86-cardiac_harmonic_imaging\'))
addpath(genpath('C:\Users\Verasonics\Desktop\MICHAELM\jopo86-cardiac_harmonic_imaging\'))

%% PathData
RootData    = 'D:\michael\';
DirData     = 'PI_Tissu_October_24_2022_14_41_44\'; %
% RootData    = 'F:\ProjetHarmonic\Acquisition20220926\';
% DirData     = 'A4CMICH200pc_September_26_2022_14_51_02\'; %


%% load Parameters
fileParam   = dir(fullfile(RootData,DirData,'*parameters*'));
load(fullfile(fileParam.folder,fileParam.name));

%%
idfile	= 5:6;

%% load Data
bufferRealigned = [];
TimeTagValues = [];
for ifile = 1:numel(idfile)
    ifile
    fileData    = dir([RootData,DirData,'*data',num2str(idfile(ifile)),'.bin']);
    %%
    [data, ...
        datatime] = LoadDataVera2D(fileData.folder, ...
        fileData.name, ...
        Trans, ...
        Receive, ...
        Resource);
    
    bufferRealigned = cat(3,bufferRealigned, data);
    TimeTagValues = cat(2,TimeTagValues, datatime);
end
TimeTagValues = TimeTagValues - min(TimeTagValues);

%% TimeTagValues
figure(1)
plot(TimeTagValues)
drawnow
%%
figure(2)
imagesc(squeeze(abs(bufferRealigned(:,ceil(end/4),1:P.nbAngles*P.nTXperAngle:end))))

%% General Parameters
speedOfSound    = Resource.Parameters.speedOfSound;
waveLength      = speedOfSound /(Trans.frequency * 1e6);
% centerFrequency = P.frequency*1e6;
fondamentalFrequency = P.frequency*1e6;%2/3*centerFrequency;
harmonicFrequency  = 2*fondamentalFrequency;

%% Define Transmit scheme
if strcmp(Trans.units,'wavelengths')
    for iAngle = 1:P.nbAngles
        iTX = (iAngle - 1)*P.nTXperAngle + 1;
        TXb.virtualSourceX(iAngle) = TX(iTX).FocalPt(1).*waveLength;
        TXb.virtualSourceY(iAngle) = TX(iTX).FocalPt(2).*waveLength;
        TXb.virtualSourceZ(iAngle) = TX(iTX).FocalPt(3).*waveLength;
    end
else
    for iAngle = 1:P.nbAngles
        iTX = (iAngle - 1)*P.nTXperAngle + 1;
        TXb.virtualSourceX(iAngle) = TX(iTX).FocalPt(1);
        TXb.virtualSourceY(iAngle) = TX(iTX).FocalPt(2);
        TXb.virtualSourceZ(iAngle) = TX(iTX).FocalPt(3);
    end
end

%% Beamforming grid
elevOffsetMm            = 4.875; % mm
ElementPosY             = zeros(160,1);
ElementPosY(81:160)     = elevOffsetMm.*1e-3;

%%
%-- Probe
idReceiveElements    = 1:160;
if strcmp(Trans.units,'wavelengths')
    Trans.ElementPos(:,2)= ElementPosY./waveLength;
    
    Geometry.posEleX     =  Trans.ElementPos(idReceiveElements,1).*waveLength;  % position of the elements along x-axis (m, vector)
    Geometry.posEleY     =  Trans.ElementPos(idReceiveElements,2).*waveLength;  % position of the elements along y-axis (m, vector)
    Geometry.posEleZ     =  Trans.ElementPos(idReceiveElements,3).*waveLength;  % position of the elements along z-axis (m, vector)
else
    Trans.ElementPos(:,2)=  ElementPosY;
    Geometry.posEleX     =  Trans.ElementPos(idReceiveElements,1);  % position of the elements along x-axis (m, vector)
    Geometry.posEleY     =  Trans.ElementPos(idReceiveElements,2);  % position of the elements along y-axis (m, vector)
    Geometry.posEleZ     =  Trans.ElementPos(idReceiveElements,3);  % position of the elements along z-axis (m, vector)
end

%%
Geometry.apertureX   =  max(Geometry.posEleX)-min(Geometry.posEleX);  % aperture along x-axis (m)
Geometry.apertureY   =  max(Geometry.posEleY)-min(Geometry.posEleY);  % aperture along x-axis (m)
Geometry.apertureZ   =  max(Geometry.posEleZ)-min(Geometry.posEleZ);  % aperture along x-axis (m)

%% Angular Resolution
AglRes      = waveLength/Geometry.apertureX; % [rad]
AglSampl    = AglRes/2; % For Central Frequency
AglSampl    = AglSampl/(3/2); % For Harmonic Frequency

%-- Grid
Origo   = [0 , 0, -P.apertureXMm*1e-3/tan(P.SectorScan/2)/2];
grid.do = AglSampl;   % o-pixel size (rad)
grid.dr = waveLength;   % r-pixel size (m)
grid.height = (P.endDepth - P.startDepth+20)*waveLength; % height of the sector (m)
grid.sector = P.SectorScan; %90*pi/180; sector (rad)
grid.o  = -grid.sector/2:grid.do:grid.sector/2; % (rad)
grid.r  = (0:grid.dr:grid.height) ...
    + (P.startDepth-5)*waveLength ...
    - Origo(3); % m

[gridBfO, ...
    gridBfR, ...
    gridBfY] = meshgrid(grid.o, ...
    grid.r, ...
    0); % Careful with convention image

gridBfZ = gridBfR.*cos(gridBfO) + Origo(3);
gridBfX = gridBfR.*sin(gridBfO);

Geometry.gridBfO     =  gridBfO;         % Theta-grid for beamforming (m, 3D meshgrid)
Geometry.gridBfP     =  0;               % Phi-grid for beamforming (m, 3D meshgrid)
Geometry.gridBfR     =  gridBfR;         % r-grid for beamforming (m, 3D meshgrid)

Geometry.gridBfX     =  gridBfX;         % x-grid for beamforming (m, 3D meshgrid)
Geometry.gridBfY     =  gridBfY;         % y-grid for beamforming (m, 3D meshgrid)
Geometry.gridBfZ     =  gridBfZ;         % z-grid for beamforming (m, 3D meshgrid)

Geometry.origoX      =  Origo(1);
Geometry.origoY      =  Origo(2);
Geometry.origoZ      =  Origo(3);

%% Beamformer parameters
Param.speedOfSound          = Resource.Parameters.speedOfSound; % Speed-of-sound (m/s)
Param.startTime             = ((-2*(Receive(1).startDepth- ...
    Trans.lensCorrection))+ ...
    TW(1).peak)/ ...
    (Trans.frequency*1e6);          % Start time (s)
Param.fnumber               = 0;

%-- IQ sampling frequency
if strcmp(Receive(1).sampleMode,'NS200BW')
    Param.samplingFrequency = (Receive(1).decimSampleRate*1e6);
elseif strcmp(Receive(1).sampleMode, 'BS100BW')
    Param.samplingFrequency = (Receive(1).decimSampleRate*1e6)/4;
elseif strcmp(Receive(1).sampleMode, 'BS50BW')
    Param.samplingFrequency = (Receive(1).decimSampleRate*1e6)/8;
else
    return;
end
Param.demodulationFrequency = Receive(1).demodFrequency*1e6;    % Demodulation frequency (Hz)

%-- Fondamental demodulation
ParamFondamental = Param;
ParamFondamental.demodulationFrequency = fondamentalFrequency;

%-- Harmonic demodulation
ParamHarmonic = Param;
ParamHarmonic.demodulationFrequency = harmonicFrequency;

%% Pulse Inversion
RawFondamental  = double( bufferRealigned(:,:,1:P.nTXperAngle:end) ) - ...
    double( bufferRealigned(:,:,2:P.nTXperAngle:end) );

RawHarmonic     = double( bufferRealigned(:,:,1:P.nTXperAngle:end) )+ ...
    double( bufferRealigned(:,:,2:P.nTXperAngle:end) );

clear bufferRealigned
%% Demodulation
disp('Demodulation')
tic
if strcmp(Receive(1).sampleMode,'NS200BW')
    fcut  = 1/6*Trans.frequency*1e6;
    fs    = Receive(1).decimSampleRate*1e6;
    [B,A] = butter(9, fcut/(fs/2));
        
    t = (0:size(RawFondamental, 1)-1).'/fs;
    
    IqRawFondamental = RawFondamental.*...
        exp(-1i*2*pi*ParamFondamental.demodulationFrequency*t);
    
    IqRawHarmonic  = RawHarmonic.*...
        exp(-1i*2*pi*ParamHarmonic.demodulationFrequency*t);
    
    IqRawFondamental = filtfiltGPU(B, A, IqRawFondamental)*2; % factor 2: to preserve the envelope amplitude
    IqRawHarmonic = filtfiltGPU(B, A, IqRawHarmonic)*2; % factor 2: to preserve the envelope amplitude
    
elseif strcmp(Receive(1).sampleMode, 'BS100BW')
    %-- 1/4 wavelegnth sampling
    IqRawFondamental = complex(RawFondamental(1:2:end,:,:),...
        -RawFondamental(2:2:end,:,:));
    
    IqRawHarmonic    = complex(RawHarmonic(1:2:end,:,:),...
        -RawHarmonic(2:2:end,:,:));
    
    %-- lowpass filter
    fcut  = 1/6*Trans.frequency*1e6;
    fs    = Receive(1).decimSampleRate*1e6/4;
    [B,A] = butter(9, fcut/(fs/2));
    
    %-- residual demodulation
    t = (0:size(IqRawFondamental, 1)-1).'/fs;
    
    IqRawFondamental = IqRawFondamental.*...
        exp(-1i*2*pi*ParamFondamental.demodulationFrequency*t).*...
        exp(1i*2*pi*Param.demodulationFrequency*t);
    
    IqRawHarmonic = IqRawHarmonic.*...
        exp(-1i*2*pi*ParamHarmonic.demodulationFrequency*t).*...
        exp(1i*2*pi*Param.demodulationFrequency*t);
    
    IqRawFondamental = filtfiltGPU(B, A, IqRawFondamental)*2; % factor 2: to preserve the envelope amplitude
    IqRawHarmonic = filtfiltGPU(B, A, IqRawHarmonic)*2; % factor 2: to preserve the envelope amplitude
  
elseif strcmp(Receive(1).sampleMode, 'BS50BW')
    %-- 1/4 wavelegnth sampling
    IqRawFondamental = complex(RawFondamental(1:2:end,:,:),...
        -RawFondamental(2:2:end,:,:));
    IqRawHarmonic    = complex(RawHarmonic(1:2:end,:,:),...
        -RawHarmonic(2:2:end,:,:));
    
    %-- lowpass filter
    fcut  = 1/6*Trans.frequency*1e6;
    fs    = Receive(1).decimSampleRate*1e6/4;
    [B,A] = butter(9, fcut/(fs/2));
    
    %-- residual demodulation
    t = (0:size(IqRawFondamental, 1)-1).'/fs;
    
    IqRawFondamental = IqRawFondamental.*...
        exp(-1i*2*pi*ParamFondamental.demodulationFrequency*t).*...
        exp(1i*2*pi*Param.demodulationFrequency*t);
    
    IqRawHarmonic = IqRawHarmonic.*...
        exp(-1i*2*pi*ParamHarmonic.demodulationFrequency*t).*...
        exp(1i*2*pi*Param.demodulationFrequency*t);
    
    IqRawFondamental = filtfiltGPU(B, A, IqRawFondamental)*2; % factor 2: to preserve the envelope amplitude
    IqRawHarmonic = filtfiltGPU(B, A, IqRawHarmonic)*2; % factor 2: to preserve the envelope amplitude
else
    return;
end
clear RawFondamental RawHarmonic
toc
%% Reshape data for clutter filtering and beamforming
IqRawFondamental = reshape(single(IqRawFondamental),...
    size(IqRawFondamental,1),...
    size(IqRawFondamental,2),...
    1,P.nbAngles,[]);
IqRawHarmonic = reshape(single(IqRawHarmonic),...
    size(IqRawHarmonic,1),...
    size(IqRawHarmonic,2),...
    1,P.nbAngles,[]);

%% Beamforming
disp('Beamforming')
Options.decreaseFactorGPU   = 8;
Options.flagCompounding     = false;
%-- For Bmode (without clutter filter)
tic
[IqBfFondamental, mask] = bfVirtualSourceGeneric(IqRawFondamental, ...
    Geometry, ...
    TXb, ...
    ParamFondamental,...
    Options);
IqBfFondamental = squeeze(IqBfFondamental);
mask = squeeze(mask);
IqBfFondamental(isinf(IqBfFondamental) | isnan(IqBfFondamental)) = 0;

[IqBfHarmonic, mask] = bfVirtualSourceGeneric(IqRawHarmonic, ...
    Geometry, ...
    TXb, ...
    ParamHarmonic,...
    Options);
IqBfHarmonic = squeeze(IqBfHarmonic);
mask = squeeze(mask);
IqBfHarmonic(isinf(IqBfHarmonic) | isnan(IqBfHarmonic)) = 0;
mask(isinf(mask) | isnan(mask)) =0;
toc
mask = mask+eps;

%%
BmodeFondamental = 20.*log10(squeeze(abs(mean(IqBfFondamental./mask,3))));
BmodeFondamental = BmodeFondamental - max(BmodeFondamental(:));

BmodeHarmonic   = 20.*log10(squeeze(abs(mean(IqBfHarmonic./mask,3))));
BmodeHarmonic   = BmodeHarmonic - max(BmodeHarmonic(:));

%%
MovAvg  = 4;
Dyn     = 30;
Sat     = 20;

%%
figure(4)
for it = 1:MovAvg:size(BmodeFondamental,3)-MovAvg
    subplot 211
    pcolor(Geometry.gridBfX.*1e3,...
        Geometry.gridBfZ.*1e3,...
        mean(BmodeFondamental(:,:,it+(0:MovAvg)),3))
    shading interp
    axis ij equal tight
    caxis([-Dyn 0] - Sat)
    colormap gray
    colorbar
    axis equal tight
    c = colorbar;
    title(['BmodeFondamental ',num2str(it)],'color',[1 1 1])
    set(c,'Color',[1 1 1])
    set(gca,'color',[0 0 0])
    set(gca,'XColor',[1 1 1])
    set(gca,'YColor',[1 1 1])
    set(gcf,'color',[0 0 0])
    
    subplot 212
    pcolor(Geometry.gridBfX.*1e3,...
        Geometry.gridBfZ.*1e3,...
        mean(BmodeHarmonic(:,:,it+(0:MovAvg)),3))
    shading interp
    axis ij equal tight
    caxis([-Dyn 0] - Sat)
    colormap gray
    colorbar
    axis equal tight
    title(['BmodeHarmonic ',num2str(it)], 'color',[1 1 1])
    c = colorbar;
    set(c,'Color',[1 1 1])
    set(gca,'color',[0 0 0])
    set(gca,'XColor',[1 1 1])
    set(gca,'YColor',[1 1 1])
    set(gcf,'color',[0 0 0])
    
    drawnow
end

%% 
figure(51);
  pcolor(Geometry.gridBfX.*1e3,...
        Geometry.gridBfZ.*1e3,...
        BmodeHarmonic(:,:,40))
    shading interp
    axis ij equal tight
    caxis([-50 ;-5])

    colormap gray
    colorbar
    axis equal tight
    title(['BmodeHarmonic  sans filtre vera',num2str(it)], 'color',[1 1 1])
    c = colorbar;
    set(c,'Color',[1 1 1])
    set(gca,'color',[0 0 0])
    set(gca,'XColor',[1 1 1])
    set(gca,'YColor',[1 1 1])
    set(gcf,'color',[0 0 0])
