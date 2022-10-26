 function ecgtrac=RecordIWORX(saveFolder,saveFile,ecgbuffer)
%% PATHS and files
RootSave = 'D:\michael\';
% DirSave     = ['simu'];
% DataName    = 'ecgtrac';
% mydate       ='test6'; %datestr(now,'mmmm_dd_yyyy_HH_MM_SS') ;;
% saveFolder  = [RootSave,DirSave,'_', mydate,'\'];
% mkdir(saveFolder)
% saveFile    = [DataName,'_', mydate];%[DataName];;
addpath(RootSave);
%% ECG PARAMETERS
disp('recording ECG');
TIME_SECONDS=2; %P.nbAcq.*P.frameRate
SETTINGS_FILE_NAME = 'test.iwxset';
BUFFER_SIZE_SECONDS = 0.05;
% CHANNEL_DATA_OFFSET_IN_PLOT = 10; % This is the offset of all the channels in plot
%% Dont change below here
IWORXDAQ = 'iwxDAQ';
STR_BUFFER_SIZE = 256;
max_grabs_per_second = int32(1/BUFFER_SIZE_SECONDS);

HEADER_PATH = 'iwxDAQ\x64\iwxDAQ.h'; % url for header
DLL_PATH = 'iwxDAQ\x64\iwxDAQ.dll'; % url for .dll

fprintf('Running the iWorxDAQ example script!\n')

%% ACTUAL CODE
loadlibrary(DLL_PATH, HEADER_PATH);
logfile = 'iworx.log';
bRet = calllib(IWORXDAQ, 'OpenIworxDevice', logfile);

% FindHardware
[bRet, model_number, model_name, serial_number] = calllib(IWORXDAQ, 'FindHardware', ...
	0, blanks(STR_BUFFER_SIZE), STR_BUFFER_SIZE, blanks(STR_BUFFER_SIZE), STR_BUFFER_SIZE); 

% Load a settings file that has been created with LabScribe
failed_LoadConf = calllib('iwxDAQ', 'LoadConfiguration', SETTINGS_FILE_NAME); 

% Get current sampling speed and num of channels
[failed_GetCurrentSamplingInfo, samplerate, num_analog_channels] = calllib(IWORXDAQ, 'GetCurrentSamplingInfo', 0, 0);
samplerate = int32(samplerate); % Matlab does not like singles (floats)
DATA_BUFFER_SIZE =  2^(nextpow2(samplerate * num_analog_channels / max_grabs_per_second)); % Lets just get a power of 2 for the buffer size
% Optional get some other parameters

[failed_GetSamplingSpeed, max_samplerate, min_samplerate, samplerate_is_shared] = calllib(IWORXDAQ, 'GetSamplingSpeed', 0, 0, 0);
[failed_GetNumChannels, analog_input_channels, analog_output_channels, digital_input_bits, digital_output_bits] = calllib(IWORXDAQ, 'GetNumChannels', 0, 0, 0, 0);


% Read Data and save it to file
chData = cell(1, num_analog_channels); % Allocate data array

fprintf('Current Model : %s\nSerial Number : %s\n', model_name, serial_number)
fprintf('Number of analog channels : %u\n', num_analog_channels)
fprintf('Samplerate : %u\n', samplerate)

%% Helper array for indecies
% if num_samples_per_ch > DATA_BUFFER_SIZE, a bigger buffer should be defined!)
% Notice: data is saved as [ch1(1), ch2(1), ch3(1)... chn(1), ch1(2), ch2(2), ch3(2)... chn(2), ...]
cntr = 0 : num_analog_channels : DATA_BUFFER_SIZE - num_analog_channels;
cur_len = 0;
% Setup a timer.
% Start Acquisition
failed_StartAcq = calllib(IWORXDAQ, 'StartAcq', DATA_BUFFER_SIZE); % setup internal buffer for 1 second worth of data.
tic
while (length(chData{1,2})<32*50)
	[iRet, num_samples_per_ch, trig_index, trig_string, data] = calllib(IWORXDAQ, 'ReadDataFromDevice', ...
		0, 0, blanks(STR_BUFFER_SIZE), STR_BUFFER_SIZE, zeros(1, DATA_BUFFER_SIZE), DATA_BUFFER_SIZE);
	% Notice: data is saved as [ch1(1), ch2(1), ch3(1)... chn(1), ch1(2), ch2(2), ch3(2)... chn(2), ...]

	if num_samples_per_ch < 1
		% If there is no data, call read data in iwxDAQ again.
		continue
	end

	% Here we have an assert that is often an issue for the first data call
	pts_recorded = num_samples_per_ch * num_analog_channels;
	assert(pts_recorded<=DATA_BUFFER_SIZE, ...
		sprintf('!!! pts_recorded<=DATA_BUFFER_SIZE. Points recorded: %u, DATA_BUFFER_SIZE=%u. Increase BUFFER_SIZE_SECONDS!!!'...
		, pts_recorded, DATA_BUFFER_SIZE))

	
	indices = 1 : num_samples_per_ch;
	buffer_indices = cntr(indices); % cntr starts from 0!

	% Save data to file
	for idx_ch = 1 : num_analog_channels
		% From current length of data,
		chData{idx_ch}(cur_len + indices) = data(buffer_indices + idx_ch);
    end
	% update current length of data for each channel
	cur_len = cur_len + num_samples_per_ch; 
end
toc

%% saving ECG
fileID =fopen('ecgcount.bin'); % Open file
ecgbuffer = fread(fileID, Inf, 'int16=>int16'); % Read file
fclose(fileID); 
save([saveFolder,[saveFile,'_data',num2str(ecgbuffer),'.mat']],'chData')
% fileRFDATA  = fopen([saveFolder,[saveFile,'_data',num2str(ecgbuffer),'.bin']],'w');
% fwrite(fileRFDATA,chData{1,2},'single','ieee-le');
% fclose(fileRFDATA);
%%%
toc(t)
%% Close down iwxDAQ
calllib(IWORXDAQ, 'StopAcq');
% Close the iWorx Device
calllib(IWORXDAQ, 'CloseIworxDevice');
% Unload the library
unloadlibrary(IWORXDAQ);
 end
