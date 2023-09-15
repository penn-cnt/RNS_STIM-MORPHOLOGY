clear; close all; clc;
addpath("../reference_data/")
rns_config = jsondecode(fileread("reference_data/config.JSON")); %/Users/wojemann/Documents/CNT/RNS_STIM-MORPHOLOGY/
rootpath = rns_config.paths.RNS_project_code_path;
rns_toolbox_path = rns_config.paths.RNS_processing_toolbox;
addpath(genpath(rootpath))
addpath(genpath(rns_toolbox_path))
datapath = rns_config.paths.RNS_PROCESSED_DATA_Folder;
rns_config_path = fullfile(rns_toolbox_path, 'config.JSON');
patient_info = struct2table(load(fullfile(rns_config.paths.RNS_metadata_path,'patients_Penn.mat')).patients_Penn);

ptList = {rns_config.patients.ID};
%%
for pt = 1:length(ptList)

% Read Patient Data
ptID = ptList{pt};
disp(['Starting analysis for ',ptID])
if ~exist([datapath ptID],'dir')
    mkdir(datapath,ptID)
end

[ecogT, ecogD, stims, histT, pdms] = loadRNSptData(ptID, rns_config);

PatientData = ecogD.AllData;

%% Data Setup
% This returns an index of the scheduled events
i_sched = find(strcmp(ecogT.ECoGTrigger, 'Scheduled'));
i_seiz = find(strcmp(ecogT.ECoGTrigger, 'Saturation') | strcmp(ecogT.ECoGTrigger, 'Long_Episode'));
% Get sampling frequency
fs = mean(ecogT.SamplingRate);
if sum(abs(ecogT.SamplingRate-fs)) > 0
    throw('Variable Sampling Rate')
end
%% Window Identification
% Steps:
% 1 Define a pre- and post-stimulation window to wait
% 2 Define a data capture time period
% 3 Use the stim indices and event start stop to exclude stim indices that
% are within the data + pre window or the post stim + data window of the start/end
% of a recorded event that was scheduled
% 4 Create recording windows around the accepted stimulation events
% 5 Run the filter Windows function on the data to make sure that there are
% no repeated stimulations within the event.
% 6 Create pre and post stim data windows for analysis!

% Step 3
event_idxs = ecogT{i_sched, {'EventStartIdx', 'EventEndIdx'}};
stim_idxs = stims.StimStartStopIndex;

%%% Filtering only scheduled events with no stims%%%
% Initialize logical array for stimulations in scheduled events
use_events = zeros(length(event_idxs),1);
% Loop over each scheduled event
for i = 1:length(event_idxs)
    % Want to find if there are any stims in each event
    use_events(i) = ~(logical(sum(ismember(stim_idxs(:,1),event_idxs(i,1):event_idxs(i,2)))...
        + sum(ismember(stim_idxs(:,2),event_idxs(i,1):event_idxs(i,2)))));
end
%%%

disp(['There are ' num2str(sum(use_events)) '/' num2str(length(event_idxs)) ' scheduled events with no stims']);
clean_events = event_idxs(logical(use_events),:);

%% Saving analysis windows
analysis_windows_idxs = clean_events;
analysis_windows = {};
for i_win = 1:length(analysis_windows_idxs)
    analysis_windows{i_win} = double(PatientData(analysis_windows_idxs(i_win,1):analysis_windows_idxs(i_win,2),:));
end
save([datapath,'/',ptID,'/event_windows_',ptID,'.mat'],'analysis_windows')

%% Calculating connectivity matrices for each 
% connections =  [1,2; % same
%                 3,4; % same
%                 1,3; % diff
%                 1,4; % diff
%                 2,3; % diff
%                 3,4];% diff
freq_array = 4:100;
freqs = [freq_array(1:end-1)',freq_array(2:end)'];
% freqs = [4,8;8,15;15,30;30,100];
num_freqs = length(freqs);
num_conns = 6;
options.orders = 8;
win_len = 250*30;
all_eeg = zeros(win_len,size(analysis_windows{1},2),length(analysis_windows));
all_plvs = zeros(length(analysis_windows),num_conns,num_freqs);
non_zero_events = ones(length(analysis_windows),1);
parfor i_win = 1:length(analysis_windows)
    if length(analysis_windows{i_win})<win_len
        non_zero_events(i_win) = 0;
        continue
    end
%     signal = analysis_windows{i_win}(1:win_len,:);
%     all_plvs(i_win,:,:) = clip_filtered_plv(signal,fs,freqs,options);
end
% all_plvs = all_plvs(logical(non_zero_events),:,:);
% save([datapath,'/',ptID,'/plvs_',ptID,'.mat'],'all_plvs')

clean_i_sched = i_sched(logical(use_events));
safe_i_sched = clean_i_sched(logical(non_zero_events));

clean_events = event_idxs(logical(non_zero_events),:);

save([datapath,'/',ptID,'/event_idxs_',ptID,'.mat'],'clean_events')
time_trace = ecogT{safe_i_sched,"RawUTCTimestamp"};
save([datapath,'/',ptID,'/UTC_time_trace_',ptID,'.mat'],'time_trace')
ptime_trace = posixtime(time_trace);
save([datapath,'/',ptID,'/posix_UTC_time_trace_',ptID,'.mat'],'ptime_trace')
try
%% Finding baseline marker
event_date = ecogT{safe_i_sched,{'RawLocalTimestamp'}};
stim_on_idx = find(strcmp("OFF",pdms.Tx1_B1),1)-1;
stim_on_date = pdms{stim_on_idx,'Programming_Date'};
is_baseline = event_date < stim_on_date;
save([datapath,'/',ptID,'/baseline_mask_',ptID,'.mat'],'is_baseline')
catch
end
end