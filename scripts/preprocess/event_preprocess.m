%% event_preprocess.m
% This script 
% 1) filters out a set of clean scheduled recordings (~90sec each)
% 2) extracts and saves data
% 3) cacluates its PLV connectivity matrix, at frequency bands from 4 ~ 100 Hz, 
% 4) extracts UTC/posix time of clean scheduled recordings 
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
    % Step 3 Use the stim indices and event start stop to exclude scheduled 
    % recordings with stim interruptions
    event_idxs = ecogT{i_sched, {'EventStartIdx', 'EventEndIdx'}};
    stim_idxs = stims.StimStartStopIndex;
    % Filtering only scheduled events with no stims
    [~,use_events] = filterWindows(event_idxs,stim_idxs); 
    disp(['There are ' num2str(sum(use_events)) '/' num2str(length(event_idxs)) ' scheduled events with no stims']);
    clean_events = event_idxs(logical(use_events),:);
    clean_i_sched = i_sched(logical(use_events));
    %% Saving analysis windows
    analysis_windows_idxs = clean_events;
    analysis_windows = {};
    for i_win = 1:length(analysis_windows_idxs)
        analysis_windows{i_win} = double(PatientData(analysis_windows_idxs(i_win,1):analysis_windows_idxs(i_win,2),:));
    end
    save(fullfile(datapath,ptID,['event_windows_',ptID,'.mat']),'analysis_windows')
    %% Calculating PLV connectivity matrices 
    % Setting parameters
    win_sec = 30;
    win_len = fs*win_sec;
    num_conns = 6;
    freq_range = [4,100];
    f = get_frequency(win_len,freq_range,fs);
    num_freqs = size(f,1);
    % initializing
    all_plvs = zeros(length(analysis_windows),num_conns,num_freqs);
    non_zero_events = ones(length(analysis_windows),1);
    for i_win = 1:length(analysis_windows)
        if length(analysis_windows{i_win})<win_len
            non_zero_events(i_win) = 0;
            continue
        end
        signal = analysis_windows{i_win}(1:win_len,:);
        plv = wavelet_plv(signal,fs,freq_range);
        all_plvs(i_win,:,:) = plv;
    end
    all_plvs = all_plvs(logical(non_zero_events),:,:);
    save(fullfile(datapath,ptID,['plvs_',ptID,'.mat']),'all_plvs','f')
    
    clean_events = event_idxs(logical(non_zero_events),:);
    safe_i_sched = clean_i_sched(logical(non_zero_events));
    
    save(fullfile(datapath,ptID,['event_idxs_',ptID,'.mat']),'clean_events')
    time_trace = ecogT{safe_i_sched,"RawUTCTimestamp"};
    save(fullfile(datapath,ptID,['UTC_time_trace_',ptID,'.mat']),'time_trace')
    ptime_trace = posixtime(time_trace);
    save(fullfile(datapath,ptID,['posix_UTC_time_trace_',ptID,'.mat']),'ptime_trace')
%     try
%         %% Finding baseline marker
%         event_date = ecogT{safe_i_sched,{'RawLocalTimestamp'}};
%         stim_on_idx = find(strcmp("OFF",pdms.Tx1_B1),1)-1;
%         stim_on_date = pdms{stim_on_idx,'Programming_Date'};
%         is_baseline = event_date < stim_on_date;
%         save(fullfile(datapath,ptID,['baseline_mask_',ptID,'.mat']),'is_baseline')
%     catch
%     end
end