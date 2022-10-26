% Xing, Paul - Polytechnique Montreal
% Adapt to GEM5sCD with PI
% Last update:
% 05/04/2022
% 20/04/2022 : Jonathan Poree : Update sequence to handle 100% BW

%%
clear persistant;
clear all;
clc;

%% Path Definition

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
DirSave     = ['PI_ULM'];
DataName    = 'PI_ULM';
DataNameECG = 'ecgbullshit';

%% PULSE CODE LOADING
flag_use_PC = false;
% % pathPulseCode = 'C:\Users\Verasonics\Documents\Vantage-4.6.2-2110271004\michael_sequenceGEM5\PI_SHI_pulseCode\PI_SHI_pulseCode';
% pathPulseCode = 'E:\MichaelM\sequencer';
%
% piAName = 'pulseCode_piA.mat';
% piBName = 'pulseCode_piB.mat';
% shiAName = 'pulseCode_shiA.mat';
% shiBName = 'pulseCode_shiB.mat';
%
% matfile_piA = matfile([pathPulseCode, filesep, piAName]);
% matfile_piB = matfile([pathPulseCode, filesep, piBName]);
% matfile_shiA = matfile([pathPulseCode, filesep, shiAName]);
% matfile_shiB = matfile([pathPulseCode, filesep, shiBName]);
%
% PC_piA = matfile_piA.TW;
% PC_piB = matfile_piB.TW;
% PC_shiA = matfile_shiA.TW;
% PC_shiB = matfile_shiB.TW;
%
% PC_piA = PC_piA.PulseCode;
% PC_piB = PC_piB.PulseCode;
% PC_shiA = PC_shiA.PulseCode;
% PC_shiB = PC_shiB.PulseCode;

%% Sequence parameters
Trans.name  = 'GEM5ScD'; % Array to use
Trans.units = 'wavelengths';
Trans.maxHighVoltage = 50;       % set maximum high voltage limit for pulser supply.
Trans           = computeTrans(Trans);

%% Angular Resolution Do Not Change
P.speedOfSound  = 1540;
P.waveLengthMm  = P.speedOfSound/Trans.frequency*1e-6*1e3;
if strcmp(Trans.units, 'wavelengths')
    P.AglRes = 1/range(Trans.ElementPos(:,1)); % Angular maximum resolution [rad]
elseif strcmp(Trans.units, 'mm')
    P.AglRes = P.waveLengthMm /range(Trans.ElementPos(:,1)); % Angular maximum resolution [rad]
end
P.AglSpl = P.AglRes/2/(2/3);  % Angular sampling (Nyquist) [rad]

%%
P.startDepth    = 10;    % Acquisition depth in wavelengths
P.endDepth      = 276;
P.endDepth      = min(P.endDepth,236)+P.startDepth;  % 276 = 15cm  This should preferrably be a multiple of 128 samples.
P.frequency     = 2/3*Trans.frequency;%Trans.frequency; % in MHz

% For PI/SHI
%--- Tx strategy for PI (number of steering will be half of P.nbangles )
P.nbAngles      = 1;  % Steering angles for PW in Tx [deg]
P.frameRate     = 1000; % 450; % Frame rate [Hz]
P.nTXperAngle   = 2;   % 2 Transmit for Pulse Inversion/ 3 pulse inversion +shi
P.Npulse        = 1;   % Number of Tx pulses (1 pulse = 1 full period/cycle)

P.SectorScan    = 75*pi/180; % Beta in the paper
P.AngleStep     = P.AglSpl;  % Alpha in the paper
angles          = (0:P.nbAngles-1)*P.AngleStep;
angles          = angles - mean(angles);
angles          = [angles(1:2:end) , fliplr(angles(2:2:end))];
P.Angles        = angles; clear angles

%%
P.Sampling      = 'NS200BW';%'BS100BW'; % sampling mode : 'BS100BW' 'BS50BW'
P.Voltage       = 20;        % Voltage
P.TGC           = [0,198,387,485,564,625,685,694]/694*1023;
P.ApodTX        = ones(1, Trans.numelements);
%[zeros(16,1); ones(96,1);zeros(16,1)]; %ones(128,1) kaiser(Resource.Parameters.numTransmit,1)'

%% Jumps for continuous sequence
%for each jump, number of acquisiton and frames
P.nbAcq     = ceil(80*24/P.nbAngles/P.nTXperAngle);%50;  % no. of Acquisitions in a Receive frame (this is a "superframe") 
P.trig      = 0;    % 1 = trig for each jump
P.PingPong  = 1;
P.live      = 1;    % live display
P.TimeTagEna = 1 ;  % initial state of time tag function
P.nbFrames  = 4;   %number of superframes for each jump (1 superframe = 1 buffer)
P.nbJumps   = 4; %50;   %number of buffers
P.simu      = 0;
P.ECGrec    = 0;
P.PauseBetweenNbFrames = 2.5;

%% Display
P.Dyn           = 40;      % Dynamic range Display
P.Saturation    = 10;      % Saturation Display
P.NframeDisplay = 10;      % P.nbAcq;       % NbFrameForDisplay

%% Specify system parameters.
% Define system parameters.
Resource.Parameters.numTransmit     = 256;     % number of transmit channels.
Resource.Parameters.numRcvChannels  = 256;    % number of receive channels (from UTA).
Resource.Parameters.speedOfSound    = P.speedOfSound;    % set speed of sound in m/sec before calling computeTrans
Resource.Parameters.verbose         = 2;
Resource.Parameters.initializeOnly  = 0;
Resource.Parameters.simulateMode    = P.simu;
Resource.Parameters.Connector       = 1; % Connector 1 or 2
Resource.VDAS.dmaTimeout            = 10000; %

%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.

%% Specify RcvProfile structure array.
% RcvProfile.LnaZinSel = 31; %To  be checked
RcvProfile.PgaGain  = 30; % dB overall gain in dB of the preamp output
% buffer prior to the A/D.
RcvProfile.PgaHPF   = 80;  % RcvProfile.PgaHPF = 80 enables the integrator feedback path,
%resulting in a first order high-pass frequency response for the PGA with a breakpoint at
%approximately 80 KHz.
RcvProfile.LnaGain  = 24; % overall gain in dB of the fixed-gain low
% noise input amplifier stage.
RcvProfile.LnaHPF   = 200; % provide a programmable high-pass frequency response and null the DC offset that would
% otherwise be present from the DC-coupled signal path.

%%
%-- Central frequency should lie within the probe bandwidth
if P.frequency<Trans.Bandwidth(1) || P.frequency>Trans.Bandwidth(2)
    warning(['Central frequency ', num2str(P.frequency),...
        ' should lie within the probe bandwidth [', num2str(Trans.Bandwidth),']']);
    P.frequency  =  sum(Trans.Bandwidth)/2;
    disp(['Central frequency has been changed to ' ,...
        num2str(P.frequency)  ])
end

%-- Voltage must stay bellow the maxHighVoltage
if P.Voltage > Trans.maxHighVoltage
    warning(['Voltage ', num2str(P.Voltage),...
        ' must stay bellow the maxHighVoltage [', num2str(Trans.maxHighVoltage	),']']);
    P.Voltage	=  Trans.maxHighVoltage;
    disp(['Voltage has been changed to ' ,...
        num2str(P.Voltage)  ])
end


%% Specify PData structure array.
P.theta = -P.SectorScan/2;
P.rayDelta = 2*(-P.theta);
P.aperture = Trans.numelements / 2 * Trans.spacing; % P.aperture in wavelengths
P.radius = (P.aperture/2)/tan(-P.theta); % dist. to virt. apex

% Set up PData structure.
PData(1).PDelta = [0.875, 0, 0.5];
PData(1).Size(1) = 10 + ceil((P.endDepth-P.startDepth)/PData(1).PDelta(3));
PData(1).Size(2) = 10 + ceil(2*(P.endDepth + P.radius)*sin(-P.theta)/PData(1).PDelta(1));
PData(1).Size(3) = 1;
PData(1).Origin = [-(PData(1).Size(2)/2)*PData(1).PDelta(1),0,P.startDepth];
PData(1).Region = struct(...
    'Shape',struct('Name','SectorFT', ...
    'Position',[0,0,-P.radius], ...
    'z',P.startDepth, ...
    'r',P.radius+P.endDepth, ...
    'angle',P.rayDelta, ...
    'steer',0));
PData(1).Region = computeRegions(PData(1));
pt1;
Media.attenuation = -0.5;
%Media.function = 'movePoints';

%% Sampling
%P.maxAcqLength      = ceil(sqrt(P.endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
P.waveLengthMm	= Resource.Parameters.speedOfSound/Trans.frequency/1e3;         % Wavelgenth [mm]

if Trans.units == 'wavelengths'
    P.apertureXMm   = range(Trans.ElementPos(:,1)).*P.waveLengthMm;
    P.apertureYMm   = range(Trans.ElementPos(:,2)).*P.waveLengthMm;
    P.apertureZMm   = range(Trans.ElementPos(:,3)).*P.waveLengthMm;
elseif Trans.units == 'mm'
    P.apertureXMm   = range(Trans.ElementPos(:,1));
    P.apertureYMm   = range(Trans.ElementPos(:,2));
    P.apertureZMm   = range(Trans.ElementPos(:,3));
end

%% todo maxAcqLength = ceil(sqrt((aperture/2)^2 + P.endDepth^2 - 2*(aperture/2)*P.endDepth*cos(P.theta/2+pi/2)) - P.startDepth);
ApertureMm          = sqrt( (P.apertureXMm)^2 + ...
    (P.apertureYMm)^2 + ...
    (P.apertureXMm)^2 ) ;
P.endDepthMm        = P.endDepth.*P.waveLengthMm;

OrigoMm             = P.apertureXMm/2/tan(P.SectorScan/2);
aMm                 = (P.endDepthMm+OrigoMm)*sin(P.SectorScan/2);
bMm                 = (P.endDepthMm+OrigoMm)*cos(P.SectorScan/2);

P.maxAcqLengthMm    = sqrt( (P.apertureXMm/2 + aMm).^2 + ...
    (bMm - OrigoMm).^2 );

% P.maxAcqLengthMm    = sqrt( (P.endDepth.*P.waveLengthMm)^2 + ...
%     (P.apertureXMm)^2 + ...
%     (P.apertureYMm)^2 + ...
%     (P.apertureXMm)^2 );

P.maxAcqLength      = P.maxAcqLengthMm./P.waveLengthMm;
P.maxAcqTime        = P.maxAcqLengthMm*2/Resource.Parameters.speedOfSound*1e-3;

maxSample           = (P.maxAcqLength - P.startDepth)*2; % [wl]
switch P.Sampling
    case 'NS200BW'
        maxSample = ceil(maxSample*4);
    case 'BS100BW'
        maxSample = ceil(maxSample*2);
    case 'BS50BW'
        maxSample = ceil(maxSample);
end
P.maxSample   = ceil(maxSample/128)*128;

% Maximum Pulse Repetition Frequency [Hz]

%% PRF
P.PRPmin	= 2 * P.maxAcqLengthMm * 1e-3 /Resource.Parameters.speedOfSound;
P.PRFmax    = 1/P.PRPmin;     % Maximum Pulse Repetition Frequency [Hz]
P.PRF       = floor(P.PRFmax/100)*100;
P.PRP       = 1/P.PRF;

if P.frameRate*P.nbAngles*P.nTXperAngle > P.PRFmax
    error(['frameRate ', num2str(P.frameRate),...
        ' must stay bellow the maximum frame rate [', num2str(P.PRFmax/(P.nbAngles*P.nTXperAngle))	,']',...
        ' Either reduce the number of flatAngles, the maximum Depth, or frame rate ' ])
end

P.BufferDuration  = P.nbAcq/P.frameRate;
disp(['BufferDuration ', num2str(P.BufferDuration),' s'])

%% Specify Resources.
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame  = P.maxSample*P.nbAngles*P.nTXperAngle*P.nbAcq; % this size allows for maximum range
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames =  P.nbFrames;    % 30 frames stored in RcvBuffer.

%% Specify Transmit waveform structure.
if flag_use_PC
    TW(2).type = 'pulseCode';
    TW(2).PulseCode = PC_piB;
    TW(1).type = 'pulseCode';
    TW(1).PulseCode = PC_piA;
    TW(3).type = 'pulseCode';
    TW(3).PulseCode = PC_shiB;
    
    %     TW(3).type = 'pulseCode';
    %     TW(3).PulseCode = PC_shiA;
else
    TW(1).type          = 'parametric'; % Set type to parametric (periodic)
    TW(1).Parameters    = [P.frequency, 0.67, 2*P.Npulse, 1]; % Parametric parameters
    TW(2).type          = 'parametric'; % Set type to parametric (periodic)
    TW(2).Parameters    = [P.frequency, 0.67, 2*P.Npulse, -1]; % Parametric parameters
end
TPC.hv              = P.Voltage; % Force initial voltage [V]

%% Specify TX structure array.
TX = repmat(struct( ...
    'Apod', P.ApodTX(:).', ... % Tx apodization
    'FocalPt', [0 0 0], ... % Focal Point
    'Delay', zeros(1, Trans.numelements)), ... % Tx delays [s]
    1, P.nTXperAngle*P.nbAngles); % One structure per angle
n = 0;
%--- Specific TX attributes (PULSE INVERSION)

for ia = 1:P.nbAngles
    for itw = 1:P.nTXperAngle
        n = n+1;
        %-- High-frame-rate echocardiography using coherent compounding
        %-- with Doppler-based motion-compensation
        th = P.Angles(ia);
        beta = P.SectorScan + max(P.Angles);
        L  = P.apertureXMm/P.waveLengthMm;
        z0 = L./(tan(th-beta/2)-tan(th+beta/2));
        x0 = z0.*tan(beta/2-th)+L/2;
        %%
        TX(n).waveform = itw;
        TX(n).FocalPt = [x0,0,z0];
        [TX(n).Delay, TX(n).FocalPt] = computeTXDelays(TX(n));
    end
end

%% Specify TGC Waveform structure.
TGC.CntrlPts = P.TGC;
TGC.rangeMax = P.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

% Specify Receive structure arrays.
% - We need na Receives for every frame.
Receive = repmat(struct('Apod', ones(1,Trans.numelements), ...
    'startDepth', P.startDepth, ...
    'endDepth', P.maxAcqLength,...
    'TGC', 1, ...
    'bufnum', 1, ...
    'framenum', 1, ...
    'acqNum', 1, ...
    'sampleMode', P.Sampling, ...
    'mode', 0, ...
    'callMediaFunc', 0), 1,  P.nTXperAngle*P.nbAngles*P.nbAcq*P.nbFrames);

%     'InputFilter',[-0.00180 +0.00000 +0.00394 +0.00000 +0.00067 +0.00000 -0.01352 ...
%                      +0.00003 +0.01575 +0.00003 +0.01508 +0.00000 -0.05600 +0.00003 ...
%                      +0.03192 +0.00003 +0.10449 +0.00000 -0.28717 +0.00003 +0.37299], ...

% - Set event specific Receive attributes for each frame.
for i = 1:P.nbFrames
    % Receive(P.numAcqs*(i-1) + 1).callMediaFunc = 1;  % move points only once per frame
    for j = 1:P.nbAcq
        % -- Acquisitions for 'super' frame.
        for k = 1: P.nbAngles
            for itw = 1:P.nTXperAngle
                rcvNum = P.nbAcq*(i-1)* P.nbAngles*P.nTXperAngle + P.nbAngles*P.nTXperAngle*(j-1)+(k-1)*P.nTXperAngle + itw;
                Receive(rcvNum).framenum = i;
                Receive(rcvNum).acqNum = P.nbAngles*P.nTXperAngle*(j-1) + (k-1)*P.nTXperAngle +itw;
            end
        end
    end
end

%% Specify Process structure array.

% %process for saving each buffer
% Process(2).classname = 'External';
% Process(2).method = 'saveBuffer';
% Process(2).Parameters = {'srcbuffer','receive',...
%     'srcbufnum',1,...
%     'srcframenum',-1,...
%     'dstbuffer','none'};

% Save Data Process
for N = 1: P.nbFrames % no. of Receive frames (real-time images are produced 1 per frame)
    Process(N).classname = 'External';
    Process(N).method = 'saveBuffer'; % calls the 'saveBuffer' function
    Process(N).Parameters = {'srcbuffer','receive',... % name of buffer to process.
        'srcbufnum',1,...
        'srcframenum',N,... % process the most recent frame.
        'dstbuffer','none'};
end

%process for live display with in-house beamformer and SVD
nProcessDisplayLive  = N+1; 
Process(nProcessDisplayLive).classname = 'External';
Process(nProcessDisplayLive).method = 'displayLive';
Process(nProcessDisplayLive).Parameters = {'srcbuffer','receive',...
    'srcbufnum',1,...
    'srcframenum',-1,...
    'dstbuffer','none'};

ncloseVSX = nProcessDisplayLive+1;
%process for closing VSX at the end of the acquisition
Process(ncloseVSX).classname      = 'External';
Process(ncloseVSX).method         = 'closeVSX';
Process(ncloseVSX).Parameters     = {'srcbuffer','receive',...  % name of buffer to process.
    'srcbufnum',1,...
    'srcframenum',0,...
    'dstbuffer','none'};

nrecordECG = ncloseVSX+1;
%process for recording ECG at each frame
Process(nrecordECG).classname    = 'External';
Process(nrecordECG).method         = 'recordECG';
Process(nrecordECG).Parameters     = {'srcbuffer','none',...  % name of buffer to process.
    'srcbufnum',0,...
    'srcframenum',0,...
    'dstbuffer','none'};

%% Specify SeqControl structure arrays.
SeqControl(1).command = 'sync'; % jump back to start
SeqControl(2).command = 'timeToNextAcq';  % time between angles
SeqControl(2).argument  = P.PRP*1e6;   % usecs
SeqControl(3).command = 'timeToNextAcq';  % time between frames
SeqControl(3).argument  = 1/P.frameRate*1e6 - ...
    (P.nbAngles*P.nTXperAngle - 1)*SeqControl(2).argument;   %  usecs
SeqControl(4).command = 'returnToMatlab';
SeqControl(5).command = 'timeToNextAcq';  % time between superframes
SeqControl(5).argument = 10000;  % 1 msecs
SeqControl(6).command   = 'loopCnt';  %loop counter          %
SeqControl(6).argument  = P.nbJumps-1; %number of loop
SeqControl(6).condition = 'counter1';
SeqControl(7).command   = 'loopTst';            % loop command
SeqControl(7).argument  = [];
SeqControl(7).condition = 'counter1';
SeqControl(8).command   = 'triggerIn'; % trigger command
SeqControl(8).condition = 'Trigger_1_Rising'; % input BNC #1
SeqControl(8).argument  = 0;
SeqControl(9).command = 'triggerOut';
SeqControl(10).command = 'noop';  % time between superframes
SeqControl(10).argument = P.PauseBetweenNbFrames*1e9/200;

nsc =  length(SeqControl)+1; % nsc is count of SeqControl objects
%Create a new variable to hold the index of the last â€˜transferToHostâ€™ SeqControl.
lastTTHnsc = 0;
TTHnsc = [];

%% Specify Event structure arrays.
n = 1;
Event(n).info           = 'loopCnt'; % jump counter event
Event(n).tx             = 0; % use next TX structure.
Event(n).rcv            = 0;
Event(n).recon          = 0; % no reconstruction.
Event(n).process        = 0; % no processing
Event(n).seqControl     = 6;
n = n+1;
SeqControl(7).argument = n;% Sets the jump event no. n

if P.ECGrec %% record ECG
    Event(n).info           = 'acquire ecg';
    Event(n).tx             = 0; % use next TX structure.
    Event(n).rcv            = 0;
    Event(n).recon          = 0; % no reconstruction.
    Event(n).process        = nrecordECG;
    Event(n).seqControl     = 1;
    n=n+1;
end

% Acquire all frames defined in RcvBuffer
for i = 1:P.nbFrames
    
    if P.trig == 1    % Trig on 1 at the begining of each buffer
        Event(n).info       = 'Wait for external trigger';
        Event(n).tx         = 0;    % use next TX structure.
        Event(n).rcv        = 0;
        Event(n).recon      = 0;    % no reconstruction.
        Event(n).process    = 0;    % no processing
        Event(n).seqControl = [8, 1, 9];
        n = n+1;
    end
    
    for j = 1:P.nbAcq % Acquire frame
        for k = 1:P.nbAngles         % Acquire angles for each frame
            for itw = 1:P.nTXperAngle
                
                Event(n).info = 'Acquire RF';
                Event(n).tx = (k-1)*P.nTXperAngle + itw;
                Event(n).rcv = P.nbAcq*(i-1)* P.nbAngles*P.nTXperAngle + ...
                    P.nbAngles*(j-1)*P.nTXperAngle+(k-1)*P.nTXperAngle + itw;
                Event(n).recon = 0;
                Event(n).process = 0;
                Event(n).seqControl = [2];
                n = n+1;
                
            end
        end
        Event(n-1).seqControl = [3,9]; % modify last acquisition Event's seqControl
    end
    % Set last acquisitions SeqControl for transferToHost.
    Event(n-1).seqControl = [3,9,nsc];
    if P.PingPong ==1
        SeqControl(nsc).command = 'transferToHost'; % transfer all acqs in one super frame
        SeqControl(nsc).condition = 'waitForProcessing';
        SeqControl(nsc).argument = nsc;
        lastTTHnsc = nsc;
        TTHnsc = [TTHnsc, lastTTHnsc];
        nsc = nsc+1;
        
        % Do reconstruction and processing for 1st sub frame
        Event(n).info       = 'saveBuffer';
        Event(n).tx         = 0; % use next TX structure.
        Event(n).rcv        = 0;
        Event(n).recon      = 0; % no reconstruction.
        Event(n).process    = [i];%[2];    % processing
        %pause the software sequencer until a transfer specified in the argument field is complete.
        Event(n).seqControl         = [nsc,nsc+1];
        SeqControl(nsc).command     = 'waitForTransferComplete';  %
        SeqControl(nsc).argument    = lastTTHnsc;
        nsc = nsc + 1;
        SeqControl(nsc).command     = 'markTransferProcessed';
        SeqControl(nsc).argument    = lastTTHnsc;
        nsc = nsc+1;
        
    else
        SeqControl(nsc).command = 'transferToHost'; % transfer all acqs in one super frame
        TTHnsc = [TTHnsc, nsc];
        nsc = nsc+1;
        Event(n).info       = 'saveBuffer';
        Event(n).tx         = 0; % use next TX structure.
        Event(n).rcv        = 0;
        Event(n).recon      = 0; % no reconstruction.
        Event(n).process    = [i];%[2];    % processing
        Event(n).seqControl = [4];
    end
    n=n+1;
%     if P.live
%         Event(n).info       = 'displayLive';
%         Event(n).tx         = 0; % use next TX structure.
%         Event(n).rcv        = 0;
%         Event(n).recon      = 0; % no reconstruction.
%         Event(n).process    = nProcessDisplayLive;    % processing
%         Event(n).seqControl = 0;
%         n=n+1;
%     end
end
%SeqControl(11).argument = lastTTHnsc;
%% 
if P.live
    Event(n).info       = 'displayLive';
    Event(n).tx         = 0; % use next TX structure.
    Event(n).rcv        = 0;
    Event(n).recon      = 0; % no reconstruction.
    Event(n).process    = nProcessDisplayLive;    % processing
    Event(n).seqControl = 0;
    n=n+1;
end
    
Event(n).info = 'noop'; % noop between frames for frame rate control
Event(n).tx = 0; % no transmit
Event(n).rcv = 0; % no rcv
Event(n).recon = 0; % no reconstruction
Event(n).process = 0; % external processing function
Event(n).seqControl = 10; % reference for ‘noop’ command
n = n+1;

%% PING PONG
if P.PingPong
    TTHshift = circshift(TTHnsc,1);
    for i = 1:length(TTHnsc)
        SeqControl(TTHnsc(i)).argument = TTHshift(i);
    end
end


%%
Event(n).info       = 'test loop cnt';
Event(n).tx         = 0; % no TX
Event(n).rcv        = 0; % no Rcv
Event(n).recon      = 0; % no Recon
Event(n).process    = 0;
Event(n).seqControl = 7; % jump command
n = n+1;


Event(n).info       = 'sync';
Event(n).tx         = 0;
Event(n).rcv        = 0;
Event(n).recon      = 0;
Event(n).process    = 0;
Event(n).seqControl = [1]; % sync command
n=n+1;

%close VSX at the end of the sequence
Event(n).info       = 'closeVSX';
Event(n).tx         = 0; % use next TX structure.
Event(n).rcv        = 0;
Event(n).recon      = 0; % no reconstruction.
Event(n).process    = ncloseVSX;%length(Process)-1;    % processing
Event(n).seqControl = 4;


%% User specified UI Control Elements
UI(1).Statement = ['[result,hv] = setTpcProfileHighVoltage(',num2str(P.Voltage),',1);'];
UI(2).Statement = 'hv1Sldr = findobj(''Tag'',''hv1Sldr'');';
UI(3).Statement = 'set(hv1Sldr,''Value'',hv);';
UI(4).Statement = 'hv1Value = findobj(''Tag'',''hv1Value'');';
UI(5).Statement = 'set(hv1Value,''String'',num2str(hv,''%.1f''));';
disp(UI(1).Statement);


%% external function
EF(1).Function = vsv.seq.function.ExFunctionDef('displayLive', @displayLive);
EF(2).Function = vsv.seq.function.ExFunctionDef('saveBuffer', @saveBuffer);
EF(3).Function = vsv.seq.function.ExFunctionDef('closeVSX', @closeVSX);
EF(4).Function = vsv.seq.function.ExFunctionDef('recordECG', @recordECG);

%% Time Tag

import com.verasonics.hal.hardware.*

switch P.TimeTagEna
    case 0
        % %         %if VDAS % can't execute this command if HW is not present
        % %             % disable time tag
        % %             rc = Hardware.enableAcquisitionTimeTagging(false);
        % %             if ~rc
        % %                 error('Error from enableAcqTimeTagging')
        % %             end
        % %         %end
        tagstr = 'off';
    case 1
        %if VDAS
        % enable time tag
        rc = Hardware.enableAcquisitionTimeTagging(true);
        if ~rc
            error('Error from enableAcqTimeTagging')
        end
        %end
        tagstr = 'on';
    case 2
        %if VDAS
        % enable time tag and reset counter
        rc = Hardware.enableAcquisitionTimeTagging(true);
        if ~rc
            error('Error from enableAcqTimeTagging')
        end
        rc = Hardware.setTimeTaggingAttributes(false, true);
        if ~rc
            error('Error from setTimeTaggingAttributes')
        end
        %end
        tagstr = 'on, reset';
end
%% Save all the structures to a .mat file.
%
% the path to the folder where data will be saved
mydate       = datestr(now,'mmmm_dd_yyyy_HH_MM_SS');
saveFolder  = [RootSave,DirSave,'_', mydate,'\'];
mkdir(saveFolder)
saveFile    = [DataName,'_', mydate];
saveFileECG = [DataNameECG,'_', mydate];
SeqName     = mfilename; % Name the sequence.mat file after the
% name of the current script
FileSeq     = [PathVantage,'\MatFiles\',SeqName];
filename    = FileSeq;
save(FileSeq);
VSX;
disp('Saving params - work in progress...');
save(fullfile(saveFolder,[saveFile,'_parameters.mat']),...
    'P','Resource','TX','TW','Trans','Media','Receive','Event','saveFolder','saveFile','-v7');
disp('Params saved !');


%% Save Current script to Data folder
CurrentScriptName = mfilename('fullpath');
copyfile([CurrentScriptName,'.m'],saveFolder);



%% **** External function definition  ****
%% Record ECG and Respiration
function [chData]=recordECG()
% keyboard()
pause(15)
disp('begining of ECG');
saveFolder = evalin('base','saveFolder');
saveFileECG = evalin('base','saveFileECG');
P = evalin('base','P');
persistent ecgcounter;
if isempty(ecgcounter)
    ecgcounter = 1;
else
    ecgcounter = ecgcounter +  P.nbFrames;
end


filebuffer  = fopen([saveFolder,'ecgcount','.bin'],'w');
fwrite(filebuffer,ecgcounter,'int16');
fclose(filebuffer);
% % RootSave='D:\michael\';
% % DirSave='simu';
% % mydate       = datestr(now,'mmmm_dd_yyyy_HH_MM_SS');
% % saveFolder  = [RootSave,DirSave,'_', mydate,'\'];
% % fileid=fopen('ecgfolder.txt','wt');
% % fwrite(fileid, saveFolder);
% % fclose(fileid);
% fileRFDATA      = fopen([saveFolder,[saveFile,'_data',num2str(counter),'.bin']],'w');
% fwrite(fileRFDATA,RData(:),'int16');
% fclose(fileRFDATA);
% tendsave = toc(tstartsave);
% disp(['Time Save = ' , num2str(tendsave,4)])
disp(['recording ECG ', num2str(ecgcounter)]);
disp(' ')
tic
job=batch(@RecordIWORX,1,{saveFolder,saveFileECG,ecgcounter});
toc
disp('its parallele')
end

%% displayLive

function displayLive(RData) %display with SVD and in-house beamforming

tstartdisp = tic;
disp ('DisplayLive')
persistent myHandle

Trans = evalin('base','Trans');
P = evalin('base','P');
Receive = evalin('base','Receive');
Resource = evalin('base','Resource');
TW = evalin('base','TW');
TX = evalin('base','TX');

if isempty(myHandle)||~ishandle(myHandle)
    myHandle = figure;
end

%% General Parameters
speedOfSound    = Resource.Parameters.speedOfSound;
waveLength      = speedOfSound /(Trans.frequency * 1e6);
centerFrequency = Trans.frequency*1e6;
fondamentalFrequency = 2/3*centerFrequency;
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
grid.subx   = 1;    % number of samples per pitch along x-axis
grid.subz   = 1;    % number of samples per wavelength along z-axis

grid.dx     = waveLength/grid.subx;   % x-pixel size (m)
grid.dz     = waveLength/grid.subz;   % z-pixel size (m)
grid.width  = P.endDepth*waveLength*tan(P.SectorScan/2); % width of the sector (m)
grid.height = (P.endDepth - P.startDepth)*waveLength; % height of the sector (m)
grid.Zmin   = P.startDepth*waveLength; % Start Depth [m]
grid.Xmin   = -grid.width/2;

grid.x      = (0:grid.dx:grid.width) + grid.Xmin;   % x values of the reconstruction grid [m]
grid.y      = 0;                                    % y values of the reconstruction grid [m]
grid.z      = (0:grid.dz:grid.height) + grid.Zmin;   % z values of the reconstruction grid [m]

[gridBfX,gridBfZ,gridBfY] = meshgrid(grid.x,grid.z,grid.y); % Careful with convention image

clear Geometry
Geometry.gridBfX     =  gridBfX;         % x-grid for beamforming (m, 3D meshgrid)
Geometry.gridBfY     =  gridBfY;         % y-grid for beamforming (m, 3D meshgrid)
Geometry.gridBfZ     =  gridBfZ;         % z-grid for beamforming (m, 3D meshgrid)

Geometry.origoX      =  0;
Geometry.origoY      =  0;
Geometry.origoZ      =  -P.apertureXMm*1e-3/tan(P.SectorScan/2)/2;

idReceiveElements    = 1:160;
if strcmp(Trans.units,'wavelengths')
    Geometry.posEleX     =  Trans.ElementPos(idReceiveElements,1).*waveLength;  % position of the elements along x-axis (m, vector)
    Geometry.posEleY     =  Trans.ElementPos(idReceiveElements,2).*waveLength;  % position of the elements along y-axis (m, vector)
    Geometry.posEleZ     =  Trans.ElementPos(idReceiveElements,3).*waveLength;  % position of the elements along z-axis (m, vector)
else
    Geometry.posEleX     =  Trans.ElementPos(idReceiveElements,1);  % position of the elements along x-axis (m, vector)
    Geometry.posEleY     =  Trans.ElementPos(idReceiveElements,2);  % position of the elements along y-axis (m, vector)
    Geometry.posEleZ     =  Trans.ElementPos(idReceiveElements,3);  % position of the elements along z-axis (m, vector)
end

Geometry.apertureX   =  max(Geometry.posEleX)-min(Geometry.posEleX);  % aperture along x-axis (m)
Geometry.apertureY   =  max(Geometry.posEleY)-min(Geometry.posEleY);  % aperture along x-axis (m)
Geometry.apertureZ   =  max(Geometry.posEleZ)-min(Geometry.posEleZ);  % aperture along x-axis (m)

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

%% Realign Buffer
RData   = RData(1:Receive(end).endSample, Trans.Connector, :);
bufferRealigned = ((RData(1:Receive(P.nbAngles*P.nTXperAngle*P.NframeDisplay).endSample,:,:)));
bufferRealigned = reshape(bufferRealigned, Receive(1).endSample, [], Trans.numelements);
bufferRealigned = permute(bufferRealigned,[1 3 2]);
bufferRealigned = bufferRealigned(:,idReceiveElements,:,:); % A change pour plus tard , ici prends juste les elements centra

%% Pulse Inversion
RawFondamental  = double( bufferRealigned(:,:,1:P.nTXperAngle:end) - ...
    bufferRealigned(:,:,2:P.nTXperAngle:end) );

RawHarmonic     = double( bufferRealigned(:,:,1:P.nTXperAngle:end) + ...
    bufferRealigned(:,:,2:P.nTXperAngle:end) );

%% Demodulation
if strcmp(Receive(1).sampleMode,'NS200BW')
    fcut  = 1/6*Trans.frequency*1e6;
    fs    = Receive(1).decimSampleRate*1e6;
    [B,A] = butter(5, fcut/(fs/2));
    
    t = (0:size(bufferRealigned, 1)-1).'/fs;
    IqRawFondamental = gpuArray(RawFondamental).*...
        exp(-1i*2*pi*ParamFondamental.demodulationFrequency*t);
    
    IqRawHarmonic  = gpuArray(RawHarmonic).*...
        exp(-1i*2*pi*ParamHarmonic.demodulationFrequency*t);
    
    IqRawFondamental = filtfiltGPU(B, A, IqRawFondamental)*2; % factor 2: to preserve the envelope amplitude
    IqRawHarmonic = filtfiltGPU(B, A, IqRawHarmonic)*2; % factor 2: to preserve the envelope amplitude
    
elseif strcmp(Receive(1).sampleMode, 'BS100BW')
    %-- 1/4 wavelegnth sampling
    IqRawFondamental = gpuArray(complex(RawFondamental(1:2:end,:,:),...
        -RawFondamental(2:2:end,:,:)));
    IqRawHarmonic    = gpuArray(complex(RawHarmonic(1:2:end,:,:),...
        -RawHarmonic(2:2:end,:,:)));
    
    %-- lowpass filter
    fcut  = 1/6*Trans.frequency*1e6;
    fs    = Receive(1).decimSampleRate*1e6/4;
    [B,A] = butter(5, fcut/(fs/2));
    
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
    IqRawFondamental = gpuArray(complex(RawFondamental(1:2:end,:,:),...
        -RawFondamental(2:2:end,:,:)));
    IqRawHarmonic    = gpuArray(complex(RawHarmonic(1:2:end,:,:),...
        -RawHarmonic(2:2:end,:,:)));
    
    %-- lowpass filter
    fcut  = 1/6*Trans.frequency*1e6;
    fs    = Receive(1).decimSampleRate*1e6/4;
    [B,A] = butter(5, fcut/(fs/2));
    
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

%% Reshape data for clutter filtering and beamforming
IqRawFondamental = reshape(IqRawFondamental, size(IqRawFondamental,1),size(IqRawFondamental,2),1,P.nbAngles,[]);
IqRawHarmonic = reshape(IqRawHarmonic, size(IqRawHarmonic,1),size(IqRawHarmonic,2),1,P.nbAngles,[]);

% %% SVDClutterFilter
% SizIqRaw = size(IqRawFondamental);
% Ncut = [5 SizIqRaw(end)];
% wcut = Ncut./SizIqRaw(end);
%
% [IqRawFondamentalf, eig_val] = SVDClutterFilter(IqRawFondamental, [],wcut);
% [IqRawHarmonicf, eig_val] = SVDClutterFilter(IqRawHarmonic, [],wcut);

%% Beamforming
%-- For Bmode (without clutter filter)
[IqBfFondamental, mask] = bfVirtualSourceGeneric(IqRawFondamental, Geometry, TXb, ParamFondamental);
IqBfFondamental = squeeze(IqBfFondamental)./mask; %normalization
IqBfFondamental(isinf(IqBfFondamental) | isnan(IqBfFondamental)) = 0;

[IqBfHarmonic, mask] = bfVirtualSourceGeneric(IqRawHarmonic, Geometry, TXb, ParamHarmonic);
IqBfHarmonic = squeeze(IqBfHarmonic)./mask; %normalization
IqBfHarmonic(isinf(IqBfHarmonic) | isnan(IqBfHarmonic)) = 0;

% %-- For PWD (with clutter filter)
% [IqBfFondamentalf, mask] = bfVirtualSourceGeneric(IqRawFondamentalf, Geometry, TXb, ParamFondamental);
% IqBfFondamentalf = squeeze(IqBfFondamentalf)./mask; %normalization
% IqBfFondamentalf(isinf(IqBfFondamentalf) | isnan(IqBfFondamentalf)) = 0;
%
% [IqBfHarmonicf, mask] = bfVirtualSourceGeneric(IqRawHarmonicf, Geometry, TXb, ParamHarmonic);
% IqBfHarmonicf = squeeze(IqBfHarmonicf)./mask; %normalization
% IqBfHarmonicf(isinf(IqBfHarmonicf) | isnan(IqBfHarmonicf)) = 0;

%% Bmode & PWD
%-- Bmode (without clutter filter)
BmodeFondamental    = 20.*log10(mean(abs(IqBfFondamental),3));
BmodeHarmonic       = 20.*log10(mean(abs(IqBfHarmonic),3));

BmodeFondamental    = BmodeFondamental - max(BmodeFondamental(:));
BmodeHarmonic       = BmodeHarmonic - max(BmodeHarmonic(:));

% %-- PWD (with clutter filter)
% PWDFondamental    = 20.*log10(mean(abs(IqBfFondamentalf),3));
% PWDHarmonic       = 20.*log10(mean(abs(IqBfHarmonicf),3));
%
% PWDFondamental    = PWDFondamental - max(PWDFondamental(:));
% PWDHarmonic       = PWDHarmonic - max(PWDHarmonic(:));

%%
figure(myHandle)
subplot 121
imagesc(unique(Geometry.gridBfX(:)).*1e3,...
    unique(Geometry.gridBfZ(:)).*1e3,...
    BmodeFondamental + P.Saturation)
caxis([-P.Dyn 0])
colorbar
title('Bmode fundamental')
ylabel('Height (mm)')
axis image

subplot 122
imagesc(unique(Geometry.gridBfX(:)).*1e3,...
    unique(Geometry.gridBfZ(:)).*1e3,...
    BmodeHarmonic + P.Saturation)
caxis([-P.Dyn 0])
title('Bmode harmonic')
colorbar
axis image
colormap gray

% subplot 223
% imagesc(unique(Geometry.gridBfX(:)).*1e3,...
%     unique(Geometry.gridBfZ(:)).*1e3,...
%     PWDFondamental + P.Saturation)
% caxis([-P.Dyn 0])
% colorbar
% title('PWD fundamental')
% xlabel('Width (mm)')
% ylabel('Height (mm)')
% axis image
%
% subplot 224
% imagesc(unique(Geometry.gridBfX(:)).*1e3,...
%     unique(Geometry.gridBfZ(:)).*1e3,...
%     PWDHarmonic + P.Saturation)
% caxis([-P.Dyn 0])
% title('PWD harmonic')
% xlabel('Width (mm)')
% colorbar
% axis image
% colormap gray
drawnow

%%
totalbf = toc(tstartdisp);
disp(['Time Display = ' , num2str(totalbf,4)])
disp(['Buffer Duration = ' , num2str(P.nbAcq/P.frameRate,4)])
disp(' ')
end

%% saveBuffer

%this function save the raw RFdata for each buffer
function saveBuffer(RData)
disp('Saving buffer');
saveFolder = evalin('base','saveFolder');
saveFile = evalin('base','saveFile');
Trans = evalin('base','Trans');
Receive = evalin('base','Receive');
tstartsave = tic;
persistent counter

if isempty(counter)
    counter = 1;
else
    counter = counter +1;
end
disp(counter);
disp('Save Raw Data');
RData           = RData(:);
RData           = (reshape(RData,size(RData,1)/256,256));
RData           = RData(1:Receive(end).endSample, Trans.Connector, :);
fileRFDATA      = fopen([saveFolder,[saveFile,'_data',num2str(counter),'.bin']],'w');
fwrite(fileRFDATA,RData(:),'int16');
fclose(fileRFDATA);
tendsave = toc(tstartsave);
disp(['Time Save = ' , num2str(tendsave,4)])
disp(' ')
end

%% closeVSX

function closeVSX(RData)
% The following, until 'return' describe the content of the function
disp('CLOSING VSX');
close('VSX Control');
end


