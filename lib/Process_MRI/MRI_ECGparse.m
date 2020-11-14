function [signal] = MRI_ECGparse(fname)
%[signal] = MRI_ECGparse(fname): takes in the *.ecg file from the Siemens
%1.5 T MRI in Snyder and parses the ECG data.
% [signal] = MRI_ECGparse(fname)
% INPUT VARIABLE:
%       fname- path+filename
% CALCULATION:
%       Reads in the two channels' data from the text file
% OUPUT:
%       signal- [time channel_1 channel_2]
%
% Thien-Khoi N. Phung (June 6, 2017- Happy Birthday Derek!)
% Base off of the PHLEM toolbox (https://sites.google.com/site/phlemtoolbox/Home) function siemens_ECGload.m

%% Load file
% If no file present- prompt to find file
if nargin < 1 | isempty(fname);
  [filename filepath] = uigetfile('*.ecg','Get ECG File');
  fname = fullfile(filepath,filename);
end
fclose('all');fid=fopen(fname);

%% Sampling rate
Hz = 400; % Hz (https://cfn.upenn.edu/aguirre/wiki/public:pulse-oximetry_during_fmri_scanning)

%% Parse through file to pull out signals
% Signal recordings start with 6002
ECGfile = textscan(fid,'%s');
starts = find(strcmp('6002',ECGfile{1}));

% Go through each starting value and pull the signal
for jz = 1:numel(starts)
    % Reopen file
    fclose('all');
    fid=fopen(fname);
    
    % Parse to start of signal reading (the 6002 spot)
    textscan(fid,'%s',starts(jz));
    
    data = textscan(fid,'%u16');
    
    % Remove the systems own evaluation of triggers.
    t_on  = find(data{1} == 5000);  % System uses identifier 5000 as trigger ON
    t_off = find(data{1} == 6000);  % System uses identifier 6000 as trigger OFF
    
    % Filter the trigger markers from the ECG data
    data_t=transpose(1:length(data{1}));
    indx = setdiff(data_t(:),union(t_on,t_off)); %Note: depending on when the scan ends, the last size(t_off)~=size(t_on).
    data_stream = data{1}(indx);
    
    % Split a single stream of ECG data into channel 1 and channel 2.
    channel_1 =   data_stream(1:2:max(size(data_stream))-1);
    channel_AVF = data_stream(2:2:max(size(data_stream)));
    
    % Make them the same length and get time estimates.
    nsamples = min(length(channel_1),length(channel_AVF));
    channel_1 = double(channel_1(1:nsamples));
    channel_AVF= double(channel_AVF(1:nsamples));
    time = double((1:nsamples)./Hz);
    
    signal{jz} = [time(:) channel_1(:) channel_AVF(:)];
end


