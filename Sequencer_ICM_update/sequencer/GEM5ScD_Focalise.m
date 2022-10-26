% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use.
%
% File name: SetUpGEM5ScD_128RyLns.m - Example of phased array imaging with 128 ray lines.
%
% Description:
%   Sequence programming file for GEM5ScD phased array using virtual apex format.  128 scan lines
%   are programmed with two acquisition per scan line. The first acquisition has all 180 transmit
%   and receive channels active and the acquisition acquires to P.endDepth.  The second acquisition
%   uses only the center row of elements and overwrites the first acquisition's data down to
%   P.midDepth.
%   The natural elevation focus of the lens on the GEM5ScD transducer is 7.7cm (~142 wavelengths).
%
%   Processing is asynchronous with respect to acquisition.
%
% Last update:
%   2/5/2019 - fix bug in range change callback VTS-997

clear all
% Path to vantage version
RootPath    = cd;
%--- Path to the Vantage software
PathVantage = 'C:\Users\Verasonics\Documents\Vantage-4.6.2-2110271004';
% PathVantage = 'C:\Users\system-01\Documents\Vantage-4.6.2-2110271004';
% PathVantage = 'C:\Users\Jonathan\Desktop\Vantage-4.6.2-2110271004';

%--- DO NOT MODIFY THE FEW LINES BELOW
%--- Move to and activate software
clearvars -except PathVantage; close all; clc; % Clean up
cd(PathVantage); % Move to path
activate; % Activate soft/hard ware
addpath(genpath('./')); % Add folder to search path

% Path to code using virtual source beamforming
% PathBF = 'C:\Users\Verasonics\Documents\Vantage-4.6.2-2110271004\michael_sequenceGEM5\BF-VirtualSourceGeneric-main';
PathLib = 'C:\Users\Verasonics\Desktop\MICHAELM\jopo86-cardiac_harmonic_imaging';
addpath(genpath(PathLib));
% Path to Data Save
RootSave    = 'D:\michael\';
DirSave     = ['FOCALISE'];
DataName    = 'FOCALISE';
DataNameBF    = 'sequence_focaliseBF';

%%
P.nbFrames  = 4 ;  %number of superframes for each jump (1 superframe = 1 buffer)
P.nbJumps   = 2 ; %number of buffers 
P.startDepth = 0;
P.midDepth = 138;
P.endDepth = 276;     % Acquisition depth in wavelengths
P.theta = 75*pi/180;  % Angle of full scan
P.numRays = 128;      % no. of raylines to program
rayDelta = P.theta/(P.numRays-1);
PathVantage='C:\Users\Verasonics\Documents\Vantage-4.6.2-2110271004';

P.PingPong =1;

% Define system parameters.
Resource.Parameters.numTransmit = 256;  % number of transmit channels.
Resource.Parameters.numRcvChannels = 256;  % number of receive channels.
Resource.Parameters.speedOfSound = 1540;
Resource.Parameters.speedCorrectionFactor = 1.0;
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode = 0;
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.

% Specify Trans structure array.
Trans.name = 'GEM5ScD';
Trans.units = 'wavelengths'; % Explicit declaration avoids warning message when selected by default
Trans = computeTrans(Trans);
mm2wl = Trans.frequency/(Resource.Parameters.speedOfSound/1000);
P.txFocus = 130*mm2wl; % specify focus in wavelengths.
P.txMidFocus = 70*mm2wl; % specify focus in wavelengths.
P.frequency = Trans.frequency;
aperture = 80*Trans.spacing; % aperture in wavelengths
dapex = (aperture/4)/tan(P.theta/2); % dist. to virt. apex, with half aperture at array.

% Specify PData structure array.
PData(1).Coord = 'polar';
PData(1).PDelta = [P.theta/127,0.5,0]; % [dtheta,dr,dz]
PData(1).Origin = [0,0,-dapex];
PData(1).Size(1) = ceil((P.endDepth - PData(1).Origin(3))/PData(1).PDelta(2)) + 10; % rows
PData(1).Size(2) = ceil(P.theta/PData(1).PDelta(1) + 1); % cols
PData(1).Size(3) = 1;   % pages
% - specify 128 Region structures.
PData(1).Region = repmat(struct('Shape',struct( ...
                    'Name','SectorFT',...
                    'Position',[0,0,-dapex],...
                    'z',P.startDepth,...
                    'r',dapex+P.endDepth,...
                    'angle',rayDelta,...
                    'steer',0)),1,P.numRays);
% - set position of regions to correspond to beam spacing.
Angle = (-P.theta/2):rayDelta:(P.theta/2);
for j = 1:P.numRays
    PData(1).Region(j).Shape.steer = Angle(j);
end
PData(1).Region = computeRegions(PData(1));
% Specify Media.  Use point targets in middle of PData(1).
% Set up Media points
% - Uncomment for speckle
% Media.MP = rand(10000,4);
% Media.MP(:,2) = 0;
% Media.MP(:,4) = 0.04*Media.MP(:,4) + 0.04;  % Random amplitude
% Media.MP(:,3) = P.endDepth*Media.MP(:,3);
% Media.MP(:,1) = (dapex+Media.MP(:,3)).*tan(P.theta*(Media.MP(:,1)-0.5));
Media.MP(1,:) = [-45,0,30,1.0];
Media.MP(2,:) = [-15,0,30,1.0];
Media.MP(3,:) = [15,0,30,1.0];
Media.MP(4,:) = [45,0,30,1.0];
Media.MP(5,:) = [-15,0,60,1.0];
Media.MP(6,:) = [-15,0,90,1.0];
Media.MP(7,:) = [-15,0,120,1.0];
Media.MP(8,:) = [-15,0,150,1.0];
Media.MP(9,:) = [-15,0,180,1.0];
Media.MP(10,:) = [-15,0,210,1.0];
Media.MP(11,:) = [-45,0,120,1.0];
Media.MP(12,:) = [15,0,120,1.0];
Media.MP(13,:) = [45,0,120,1.0];
Media.MP(14,:) = [-10,0,69,1.0];
Media.MP(15,:) = [-5,0,75,1.0];
Media.MP(16,:) = [0,0,78,1.0];
Media.MP(17,:) = [5,0,80,1.0];
Media.MP(18,:) = [10,0,81,1.0];
Media.MP(19,:) = [-75,0,120,1.0];
Media.MP(20,:) = [75,0,120,1.0];
Media.MP(21,:) = [-15,0,180,1.0];
Media.numPoints = 21;
% Media.attenuation = -0.5;
Media.function = 'movePoints';

% Specify Resources.
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = P.numRays*4224;
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = 48;    % 24 frames used for RF cineloop.
Resource.ImageBuffer(1).numFrames = 48;
Resource.DisplayWindow(1).Title = 'GEM5ScD_128RyLnsPolar';
Resource.DisplayWindow(1).pdelta = 0.35;
ScrnSize = get(0,'ScreenSize');
pdeltaR = PData(1).PDelta(2);
pdeltaT = PData(1).PDelta(1);
DwWidth = ceil((2*PData.Size(1)*pdeltaR*sin(pdeltaT*(PData.Size(2)-1)/2))/Resource.DisplayWindow(1).pdelta);
DwHeight = ceil(P.endDepth/Resource.DisplayWindow(1).pdelta);

Resource.DisplayWindow(1).Type = 'Verasonics';
Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
                                      DwWidth, DwHeight];
Resource.DisplayWindow(1).ReferencePt = [-DwWidth/2*Resource.DisplayWindow(1).pdelta,0,0];   % 2D imaging is in the X,Z plane
Resource.DisplayWindow(1).numFrames = 20;
Resource.DisplayWindow(1).AxesUnits = 'mm';
Resource.DisplayWindow.Colormap = gray(256);

% Specify Transmit waveform structure.
TW.type = 'parametric';
TW.Parameters = [P.frequency,.67,2,1];   % A, B, C, D

% Specify TX structure array.
TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'focus', P.txFocus, ...
                   'Steer', [0.0,0.0], ...
                   'Apod', zeros(1,160), ...
                   'Delay', zeros(1,160)), 1, 2*P.numRays);

% - Set event specific TX attributes.
m = 200;
W = hann(m);
TXorgs = dapex*tan(Angle);
TApod = zeros(1,80); % temporary apod array
for n = 1:P.numRays
    TX(n).Origin = [TXorgs(n), 0.0, 0.0];
    ce = round(80*(TXorgs(n) - Trans.ElementPos(1,1))/aperture);
    TApod(ce) = W(m/2);
    for j=1:80
        lft=ce-j;
        if lft<1, lft=1; elseif lft<(ce-m/2), lft = ce-m/2; end
        rt=ce+j;
        if rt>80, rt=80; elseif rt>(ce+m/2), rt = ce+m/2; end
    end
    TApod(lft:rt) = W((m/2-(ce-lft)):(m/2+(rt-ce)));
    TX(n).Apod(1:80) = TApod;
    TX(n).Apod(81:160) = 0;
    TX(n).Steer = [Angle(n),0.0];
    TX(n).Delay = computeTXDelays(TX(n));
    TX(n).Apod(81:160) = TX(n).Apod(1:80);
    TX(n).Delay(81:160) = TX(n).Delay(1:80);
end
m = P.numRays;
for n = 1:P.numRays
    TX(n+m) = TX(n);
    TX(n+m).Apod(81:160) = 0; % turn off outer elements
    TX(n+m).focus = P.txMidFocus;
    TX(n+m).Delay = computeTXDelays(TX(n+m));
end

% Specify Receive structure arrays.
maxAcqLength = ceil(sqrt((aperture/2)^2 + P.endDepth^2 - 2*(aperture/2)*P.endDepth*cos(P.theta/2+pi/2)) - P.startDepth);
Receive = repmat(struct('Apod', ones(1,160), ...
                        'startDepth', P.startDepth, ...
                        'endDepth', maxAcqLength, ...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode', 'NS200BW', ...
                        'mode', 0, ...
                        'callMediaFunc',0),1,2*P.numRays*Resource.RcvBuffer(1).numFrames);
% - Set event specific Receive attributes.
for i = 1:Resource.RcvBuffer(1).numFrames
    k = 2*P.numRays*(i-1);
    Receive(k+1).callMediaFunc = 1;
    for j = 1:2:2*P.numRays
        Receive(k+j).framenum = i;
        Receive(k+j).acqNum = (j+1)/2;
        Receive(k+j+1).Apod(1:80) = 1.1;
        Receive(k+j+1).Apod(81:160) = 0;
        Receive(k+j+1).framenum = i;
        Receive(k+j+1).endDepth = P.midDepth;
        Receive(k+j+1).acqNum = (j+1)/2;
    end
end

% Specify TGC Waveform structure.
TGC.CntrlPts = [0,198,387,485,564,625,685,694];
TGC.rangeMax = P.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

% Specify Recon structure arrays.
Recon = struct('senscutoff', 0.4, ...
               'pdatanum', 1, ...
               'rcvBufFrame', -1, ...
               'ImgBufDest', [1,-1], ...
               'RINums', 1:P.numRays);

% Define ReconInfo structures.
ReconInfo = repmat(struct('mode', 'replaceIntensity', ...
                          'txnum', 1, ...
                          'rcvnum', 1, ...
                          'scaleFactor', 0.5, ...
                          'regionnum', 0), 1, P.numRays);
% - Set specific ReconInfo attributes.
for j = 1:P.numRays
    ReconInfo(j).txnum = j;
    ReconInfo(j).rcvnum = 2*j - 1;
    ReconInfo(j).regionnum = j;
end

% Specify Process structure array.
pers = 20;
Process(1).classname = 'Image';
Process(1).method = 'imageDisplay';
Process(1).Parameters = {'imgbufnum',1,...   % number of buffer to process.
                         'framenum',-1,...   % (-1 => lastFrame)
                         'pdatanum',1,...    % number of PData structure to use
                         'srcData','intensity2D',...
                         'pgain',0.5,...            % pgain is image processing gain
                         'reject',2,...      % reject level
                         'persistMethod','simple',...
                         'persistLevel',pers,...
                         'interpMethod','4pt',...  %method of interp. (1=4pt)
                         'grainRemoval','none',...
                         'processMethod','none',...
                         'averageMethod','none',...
                         'compressMethod','power',...
                         'compressFactor',40,...
                         'mappingMethod','full',...
                         'display',1,...      % display image after processing
                         'displayWindow',1};
                     
Process(2).classname = 'External';
Process(2).method = 'saveBuffer';

Process(2).Parameters = {'srcbuffer','receive',...
                         'srcbufnum',1,...
                         'srcframenum',0,...
                         'dstbuffer','none'};
Process(3).classname = 'External';
Process(3).method = 'saveBufferBF';
                     
Process(3).Parameters = {'srcbuffer','imageP',...
                         'srcbufnum',1,...
                         'srcframenum',0,...
                         'dstbuffer','none'};  

%  - Jump back to start.
SeqControl(1).command = 'jump';
SeqControl(1).argument = 1;
SeqControl(2).command = 'timeToNextAcq'; % time between rays
SeqControl(2).argument = 250;
SeqControl(3).command = 'timeToNextAcq'; % time between frames
SeqControl(3).argument = 80000 - (P.numRays-1)*SeqControl(2).argument;
SeqControl(4).command = 'returnToMatlab';
nsc =  length(SeqControl)+1;
n = 1; % n is count of Events
%Create a new variable to hold the index of the last ‘transferToHost’ SeqControl.
lastTTHnsc = 0;
TTHnsc = [];
% Acquire all frames defined in RcvBuffer

for i = 1:Resource.RcvBuffer(1).numFrames
    k = 2*P.numRays*(i-1);
    for j = 1:P.numRays
        Event(n).info = 'Acquire endDepth ray line';
        Event(n).tx = j;         % use 1st TX structure.
        Event(n).rcv = k+2*j-1;
        Event(n).recon = 0;      % no reconstruction.
        Event(n).process = 0;    % no processing
        Event(n).seqControl = 2;
        n = n+1;
        Event(n).info = 'Acquire midDepth ray line';
        Event(n).tx = P.numRays+j;         % use 1st TX structure.
        Event(n).rcv = k+2*j;
        Event(n).recon = 0;      % no reconstruction.
        Event(n).process = 0;    % no processing
        Event(n).seqControl = 2;
        n = n+1;
    end
    Event(n-1).seqControl = [3,nsc]; % time between frames & transferToHostuse
    if P.PingPong ==1
        SeqControl(nsc).command = 'transferToHost'; % transfer all acqs in one super frame
        SeqControl(nsc).condition = 'waitForProcessing';
        SeqControl(nsc).argument = nsc;
        lastTTHnsc = nsc;
        TTHnsc = [TTHnsc, lastTTHnsc];
        nsc = nsc+1;
        
        % Do reconstruction and processing for 1st sub frame
        
    end
    %     n=n+1;
    Event(n).info = 'Reconstruct';
    Event(n).tx = 0;         % no transmit
    Event(n).rcv = 0;        % no rcv
    Event(n).recon = 1;      % reconstruction
    Event(n).process = 1;    % process
    Event(n).seqControl = 0;
    %     if floor(i/3) == i/3     % Exit to Matlab every 3rd frame
    Event(n).seqControl = 4;
    %     end
    n = n+1;
end

EventSaveBF = n;
Event(n).info       = 'saveBufferBF';
Event(n).tx         = 0; % use next TX structure.
Event(n).rcv        = 0;
Event(n).recon      = 0; % no reconstruction.
Event(n).process    = 3;    % processing

%pause the software sequencer until a transfer specified in the argument field is complete.
Event(n).seqControl         = [nsc,nsc+1];
SeqControl(nsc).command     = 'waitForTransferComplete';  %
SeqControl(nsc).argument    = lastTTHnsc;
nsc = nsc + 1;
SeqControl(nsc).command     = 'markTransferProcessed';
SeqControl(nsc).argument    = lastTTHnsc;
% nsc = nsc+1;
n=n+1;

EventSaveRF = n;
Event(n).info       = 'saveBufferRF';
Event(n).tx         = 0; % use next TX structure.
Event(n).rcv        = 0;
Event(n).recon      = 0; % no reconstruction.
Event(n).process    = 2;    % processing

%pause the software sequencer until a transfer specified in the argument field is complete.
Event(n).seqControl         = [nsc,nsc+1];
SeqControl(nsc).command     = 'waitForTransferComplete';  %
SeqControl(nsc).argument    = lastTTHnsc;
nsc = nsc + 1;
SeqControl(nsc).command     = 'markTransferProcessed';
SeqControl(nsc).argument    = lastTTHnsc;
nsc = nsc+1;
n=n+1;



if P.PingPong
TTHshift = circshift(TTHnsc,1);  
for i = 1:length(TTHnsc)
    SeqControl(TTHnsc(i)).argument = TTHshift(i);
end
end

Event(n).info = 'Jump back to first event';
Event(n).tx = 0;        % no TX
Event(n).rcv = 0;       % no Rcv
Event(n).recon = 0;     % no Recon
Event(n).process = 0;
Event(n).seqControl = 1; % jump command

Event(n).info       = 'sync';
Event(n).tx         = 0;
Event(n).rcv        = 0;
Event(n).recon      = 0;
Event(n).process    = 0;
Event(n).seqControl = [1]; % sync command
n=n+1;


% User specified UI Control Elements
% - Sensitivity Cutoff
UI(1).Control =  {'UserB7','Style','VsSlider','Label','Sens. Cutoff',...
                  'SliderMinMaxVal',[0,1.0,Recon(1).senscutoff],...
                  'SliderStep',[0.025,0.1],'ValueFormat','%1.3f'};
UI(1).Callback = text2cell('%SensCutoffCallback');

% - Range Change
wls2mm = 1;
AxesUnit = 'wls';
if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
    if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm');
        AxesUnit = 'mm';
        wls2mm = Resource.Parameters.speedOfSound/1000/Trans.frequency;
    end
end
UI(2).Control = {'UserA1','Style','VsSlider','Label',['Range (',AxesUnit,')'],...
                 'SliderMinMaxVal',[64,300,P.endDepth]*wls2mm,'SliderStep',[0.1,0.2],'ValueFormat','%3.0f'};
UI(2).Callback = text2cell('%RangeChangeCallback');
% - Transmit focus change
UI(3).Control = {'UserB4','Style','VsSlider','Label',['TX Focus (',AxesUnit,')'],...
                 'SliderMinMaxVal',[20,320,P.txFocus]*wls2mm,'SliderStep',[0.1,0.2],'ValueFormat','%3.0f'};
UI(3).Callback = text2cell('%TxFocusCallback');
% UI(4).Control = {'UserC2','Style','VsPushButton','Label','SaveRFnBF'};
% UI(4).Callback = text2cell('%SaveButtonCallback');
%% external function
EF(1).Function = vsv.seq.function.ExFunctionDef('saveBuffer', @saveBuffer);
EF(2).Function = vsv.seq.function.ExFunctionDef('saveBufferBF', @saveBufferBF);

%%
%%
mydate       = datestr(now,'mmmm_dd_yyyy_HH_MM_SS');
saveFolder  = [RootSave,DirSave,'_', mydate,'\'];
mkdir(saveFolder)
saveFile    = [DataName,'_', mydate];
saveFileBF    = [DataNameBF,'_', mydate];

SeqName     = mfilename; % Name the sequence.mat file after the
% name of the current script
FileSeq     = [PathVantage,'\MatFiles\',SeqName];
filename    = FileSeq;
save(FileSeq);
disp('Saving params - work in progress...');
save(fullfile(saveFolder,[saveFile,'_parameters.mat']),...
    'P','PData','Resource','TX','TW','Trans','Media','Receive','Event','saveFolder','saveFile','-v7');
disp('Params saved !');
%%
% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 3;

% Save all the structures to a .mat file.
save('MatFiles/GEM5ScD_128RyLnsPolar');

%%
VSX;

%% Save Current script to Data folder
CurrentScriptName = mfilename('fullpath');
copyfile([CurrentScriptName,'.m'],saveFolder);

%%
% addtC:\Users\Verasonics\Desktop\MICHAELM\jopo86-cardiac_harmonic_imaging\sequencer
% GEM5ScD_PI_ECG_Save_2022_09_15;
GEM5ScD_Centralfreq;
GEM5ScD_PI_Tissu;
GEM5ScD_PI_ULM;


%%
return


% **** Callback routines to be converted by text2cell function. ****
%SensCutoffCallback - Sensitivity cutoff change
ReconL = evalin('base', 'Recon');
for i = 1:size(ReconL,2)
    ReconL(i).senscutoff = UIValue;
end
assignin('base','Recon',ReconL);
Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'Recon'};
assignin('base','Control', Control);
return
%SensCutoffCallback

%RangeChangeCallback - Range change
simMode = evalin('base','Resource.Parameters.simulateMode');
% No range change if in simulate mode 2.
if simMode == 2
    set(hObject,'Value',evalin('base','P.endDepth'));
    return
end
Trans = evalin('base','Trans');
aperture = 80*Trans.spacing; % aperture in wavelengths
Resource = evalin('base','Resource');
mm2wl = Trans.frequency/(Resource.Parameters.speedOfSound/1000);
P = evalin('base','P');
P.endDepth = UIValue;
if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
    if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm')
        P.endDepth = UIValue*mm2wl;
    end
end
assignin('base','P',P);
rayDelta = P.theta/(P.numRays-1);
dapex = (aperture/4)/tan(P.theta/2); % dist. to virt. apex
PData = evalin('base','PData');
PData(1).PDelta = [P.theta/127,0.5,0]; % [dtheta,dr,dz]
PData(1).Origin = [0,0,-dapex];
PData(1).Size(1) = ceil((P.endDepth - PData(1).Origin(3))/PData(1).PDelta(2)) + 10; % rows
PData(1).Size(2) = ceil(P.theta/PData(1).PDelta(1) + 1); % cols
PData(1).Size(3) = 1;   % pages
% - specify 128 Region structures.
PData(1).Region = repmat(struct('Shape',struct( ...
                    'Name','SectorFT',...
                    'Position',[0,0,-dapex],...
                    'z',P.startDepth,...
                    'r',dapex+P.endDepth,...
                    'angle',rayDelta,...
                    'steer',0)),1,P.numRays);
% - set position of regions to correspond to beam spacing.
Angle = (-P.theta/2):rayDelta:(P.theta/2);
for j = 1:P.numRays
    PData(1).Region(j).Shape.steer = Angle(j);
end
PData(1).Region = computeRegions(PData(1));
assignin('base','PData',PData);
% Re-size DisplayWindow
evalin('base','Resource.DisplayWindow(1).Position(4) = ceil(P.endDepth/Resource.DisplayWindow(1).pdelta);');
% Re-calculate Receive.endDepth
Receive = evalin('base', 'Receive');
maxAcqLength = ceil(sqrt((aperture/4)^2 + P.endDepth^2 - 2*(aperture/4)*P.endDepth*cos(P.theta/2+pi/2)) - P.startDepth);
for i = 1:size(Receive,2)
    Receive(i).endDepth = maxAcqLength;
end
assignin('base','Receive',Receive);
evalin('base','TGC.rangeMax = P.endDepth;');
evalin('base','TGC.Waveform = computeTGCWaveform(TGC);');
Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'PData','InterBuffer','ImageBuffer','DisplayWindow','Receive','TGC','Recon'};
assignin('base','Control', Control);
assignin('base', 'action', 'displayChange');
return
%RangeChangeCallback

%TxFocusCallback - TX focus changel
simMode = evalin('base','Resource.Parameters.simulateMode');
% No focus change if in simulate mode 2.
if simMode == 2
    set(hObject,'Value',evalin('base','P.txFocus'));
    return
end
Trans = evalin('base','Trans');
Resource = evalin('base','Resource');
mm2wl = Trans.frequency/(Resource.Parameters.speedOfSound/1000);

P = evalin('base','P');
P.txFocus = UIValue;
if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
    if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm');
        P.txFocus = UIValue*mm2wl;
    end
end
assignin('base','P',P);

TX = evalin('base', 'TX');
for n = 1:size(TX,2)   % for all transmit events
    TX(n).focus = P.txFocus;
    TX(n).Delay = computeTXDelays(TX(n));
end
assignin('base','TX', TX);
% Set Control command to update TX
Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'TX'};
assignin('base','Control', Control);
return

save('data.mat','ImageBF');

%TxFocusCallback


% %SaveButtonCallback
%     Event = evalin('base','Event');
%     EventSaveRF = evalin('base','EventSaveRF');
%     EventSaveBF = evalin('base','EventSaveBF');
%     Event(EventSaveRF).process=2;
%     Event(EventSaveBF).process=3;
%     assignin('base','Event',Event);
%     Control = evalin('base','Control');
%     Control.Command = 'update&Run';
%     Control.Parameters = {'Event'};
%     assignin('base','Control', Control);
%     return
% %SaveButtonCallback

%this function save the raw RFdata for each buffer
function saveBuffer(RData) 
    disp('Saving buffer');
    saveFolder = evalin('base','saveFolder');
    saveFile = evalin('base','saveFile');
    Trans = evalin('base','Trans');
    Receive = evalin('base','Receive');
    tstartsave = tic;
%     persistent counter
% 
%     if isempty(counter)
%         counter = 1;
%     else
%         counter = counter +1;
%     end
%     disp(counter);
    
    disp('Save Raw Data');
    RData           = RData(:);
    RData           = (reshape(RData,size(RData,1)/256,256));
    RData           = RData(1:Receive(end).endSample, Trans.Connector, :);
    fileRFDATA      = fopen([saveFolder,[saveFile,'_data.bin']],'w');%,num2str(counter),
    fwrite(fileRFDATA,RData(:),'int16');
    fclose(fileRFDATA);
    tendsave = toc(tstartsave);
    disp(['Time Save = ' , num2str(tendsave,4)])
    disp(' ')
end

function saveBufferBF(DataBF) 
    disp('Saving Beamformed Data');
    saveFolder = evalin('base','saveFolder');
    saveFileBF = evalin('base','saveFileBF');
   
    ImgBF=squeeze(DataBF);
    tstartsave = tic;
%     persistent counter
% 
%     if isempty(counter)
%         counter = 1;
%     else
%         counter = counter +1;
%     end
%     disp(counter);
    disp('Save Raw Data');
    fileRFDATA      = fopen([saveFolder,[saveFileBF,'_data.bin']],'w');%',num2str(counter),'
    fwrite(fileRFDATA,ImgBF(:),'double');
    fclose(fileRFDATA);
    tendsave = toc(tstartsave);
    disp(['Time Save = ' , num2str(tendsave,4)])
    disp(' ')
end