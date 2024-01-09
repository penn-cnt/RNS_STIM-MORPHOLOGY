%% event_preprocess.m
% This script 
% 1) filters out a set of clean stimulations, not during scheduled events
% nor within 4h of epileptic recordings
% 2) extracts and saves data
% 3) extracts UTC/posix time of clean stims 
% 5) extracts the stim is after which visit 
%% settings
clear; close all; 
paths;
patient_info = struct2table(load(which('patients_Penn.mat')).patients_Penn);
ptList = {rns_config.patients.ID};
%%
for pt = 1:length(ptList)
    %% Read Patient Data
    ptID = ptList{pt};
    disp(['Starting analysis for ',ptID])
    % data folder
    if ~exist(fullfile(datapath,ptID),'dir')
        mkdir(datapath,ptID)
    end
    % load RNS data
    [ecogT, ecogD, stims, histT, pdms] = loadRNSptData(ptID, rns_config);
    PatientData = ecogD.AllData;
    %% Data Setup
    % Extract index of the scheduled/abnormal recordings
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
    PRE_STIM = 0.01; % (seconds)
    POST_STIM = 2; % (seconds)
    STIM_TIME = 1; % (seconds)
    pre_idx = PRE_STIM * fs;
    post_idx = POST_STIM * fs;
    stim_time_idx = 125;

    % 2 Define a data capture time period
    DATA_WINDOW = 2; % (seconds)
    data_window_idx = DATA_WINDOW * fs;

    % 3 Use the stim indices and event start stop to exclude stim during
    % scheduled recordings
    event_idxs = ecogT{i_sched, {'EventStartIdx', 'EventEndIdx'}};
    stim_idxs = stims.StimStartStopIndex;
    [~,excl] = filterWindows(stim_idxs,event_idxs);
    clean_stims = stim_idxs(excl,:);
    disp(['There are ' num2str(length(clean_stims)) '/' ...
        num2str(length(stim_idxs)) ' stimulations passing initial criteria']);
    % Step 6
    % hours after seizure to wait
    SEIZURE_BUFFER = 4;
    % create boundaries for each "seizure" in datetime
    seizure_idxs = idx2time(ecogT,ecogT{i_seiz, {'EventStartIdx', 'EventEndIdx'}});
    % add seizure buffer in hours
    seizure_idxs(:,2) = seizure_idxs(:,2) + hours(SEIZURE_BUFFER); 
    [~,excl] = filterWindows(idx2time(ecogT,clean_stims),seizure_idxs); 
    no_seiz_stims = excl; 
    disp(['There are ' num2str(sum(no_seiz_stims(:,1))) '/' num2str(length(clean_stims)) ' stimulations >4h after a seizure']);
    clean_idxs = logical(no_seiz_stims);
    clean_stim_windows = clean_stims(clean_idxs,:);

    save(fullfile(datapath,ptID,['stim_idxs_',ptID,'.mat']),'clean_stim_windows')
    time_trace = idx2time(ecogT,clean_stim_windows(:,1),...
        'timezone','UTC');
    save(fullfile(datapath,ptID,['stim_UTC_time_trace_',ptID,'.mat']),'time_trace')
    ptime_trace = posixtime(time_trace);
    save(fullfile(datapath,ptID,['stim_posix_UTC_time_trace_',ptID,'.mat']),'ptime_trace')
    %% Saving analysis windows
    analysis_windows_idxs = [clean_stim_windows(:,2)-data_window_idx, clean_stim_windows(:,2)+data_window_idx];
    analysis_windows = {};
    for i_win = 1:length(analysis_windows_idxs)
        analysis_windows{i_win} = double(PatientData(analysis_windows_idxs(i_win,1):analysis_windows_idxs(i_win,2),:));
    end
    save(fullfile(datapath,ptID,['stim_windows_',ptID,'.mat']),'analysis_windows')
    %% Create Programming Event Array
    event_id_table = ecogT(:,{'RawLocalTimestamp','EventStartIdx','EventEndIdx'});
    session_id_table = pdms(:,'Programming_Date'); 
    session_codes = 1:height(session_id_table);
    event_times = sort(posixtime(table2array(event_id_table(:,'RawLocalTimestamp'))));
    session_times = posixtime(table2array(session_id_table));

    % Make a logical array that says wether the index in this event is after
    % the time given for the clinical session where stimulation parameters were
    % set
    session_logical = zeros(length(event_times),length(session_times));
    for s = 1:length(session_times)
        session_logical(:,s) = event_times > session_times(s);
    end

    % Calculate the first event after each session
    first_event = diff(session_logical, [],1);
    [row,col] = find(first_event); % This is the row of events that are the first in the new session

    % session_idxs contains the indices for the "next session". So if the time
    % index is less than the index in session index it is in the session before
    % it
    session_idxs = [1; flip(event_id_table.EventStartIdx(row+1)); ecogT.EventEndIdx(end)+1];

    % bin the indices of the pre and post stimuli windows and create an index
    % for the plv function to use to classify the stims.
    visit_selection_array = [discretize(clean_stim_windows(:,1),session_idxs)];
    save(fullfile(datapath,ptID,['visit_selection_array_',ptID,'.mat']),'visit_selection_array')
end