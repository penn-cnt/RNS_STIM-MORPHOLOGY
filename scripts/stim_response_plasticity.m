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

    % Step 1
    PRE_STIM = 0.01; % (seconds)
    POST_STIM = 2; % (seconds)
    STIM_TIME = 1; % (seconds)

    pre_idx = PRE_STIM * fs;
    post_idx = POST_STIM * fs;
    stim_time_idx = 125;

    % Step 2
    DATA_WINDOW = 2; % (seconds)
    data_window_idx = DATA_WINDOW * fs;

    % Step 3
    event_idxs = ecogT{i_sched, {'EventStartIdx', 'EventEndIdx'}};
    stim_idxs = stims.StimStartStopIndex;

    %%% Filtering only stimulations during scheduled events%%%
    % Grab only stimulation start times
%     stimstarts = stim_idxs(:,1);
%     % Initialize logical array for stimulations in scheduled events
%     use_stims1 = zeros(size(stimstarts));
%     % Loop over each scheduled event
%     for i = 1:length(event_idxs)
%         % Might be able to optimize by checking event index of each stim and
%         % checking overall ecogT idx membership in ecogT scheduled event
% 
%         % for each scheduled event, check membership of each stimulation start and
%         % add to cumulative logical array (should be no value larger than 1)
% 
%         %     use_stims1 = use_stims1 + logical(ismember(stimstarts,event_idxs(i,1):event_idxs(i,2)));
%         use_stims1(find(logical(ismember(stimstarts,event_idxs(i,1):event_idxs(i,2))),1)) = true;
%     end
%     %%%
% 
%     disp(['There are ' num2str(sum(use_stims1)) '/' num2str(length(stim_idxs)) ' stimulations in scheduled events']); % started in
% 
%     %%% Filtering only stimulations without event interuption
%     % Iterate through the stimulations to see if the windows contain the
%     % start of a new event -- FUNCTION
%     % Receiver for stimulation windows exclusion logical array
%     % that don't overlap with event windows
%     use_stim = ones(length(stim_idxs),1);
%     for i = 1:length(stim_idxs)
% 
%         % Create boundaries for before buffer
%         start = stim_idxs(i,1) - (data_window_idx);
%         finish = stim_idxs(i,2) + (data_window_idx);
% 
%         % Check to see if any indices in stim buffer contain an event index
%         use_stim(i) = ~logical(sum(ismember(event_idxs(:,1),start:finish)...
%             + ismember(event_idxs(:,2),start:finish)));
%     end
    

    % Collect only the stimulation events that pass both criteria
%     clean_stims = stim_idxs(logical(use_stim.*use_stims1),:);
%     event_idxs(:,1) = event_idxs(:,1) - data_window_idx;
%     event_idxs(:,2) = event_idxs(:,2) + data_window_idx;
    [~,excl] = filterWindows(stim_idxs,event_idxs);
    clean_stims = stim_idxs(excl,:);
    disp(['There are ' num2str(length(clean_stims)) '/' ...
        num2str(length(stim_idxs)) ' stimulations passing initial criteria']);
    pstats.NumStims(pt) = length(stim_idxs);
    %% Window Identification Cont.
    % Step 4
    % Add buffers to stim start and stop times to get data windows for analysis
    % pre_stim_windows_dirty = [clean_stims(:,1) - (pre_idx + data_window_idx) ...
    %     clean_stims(:,1) - pre_idx];
    %
    % post_stim_windows_dirty = [clean_stims(:,2) + post_idx ...
    %     clean_stims(:,2) + post_idx + data_window_idx];
    %
    % dirty_stims = [clean_stims(:,1)+stim_time_idx, clean_stims(:,2)+post_idx];%%%
    %
    % % Step 5
    % % Check data windows for any stimulations in or surrounding window edges
    % [~,pre_clean_idxs] = filterWindows(pre_stim_windows_dirty,stim_idxs);
    % [~,post_clean_idxs] = filterWindows(post_stim_windows_dirty,stim_idxs);
    %
    % % Check stim windows for any second stimulations
    % [~,stim_clean_idxs] = filterWindows(dirty_stims,stim_idxs); %%%

    % Step 6
    % Check that stimulations are more than x hours after a seizure

    % hours after seizure to wait
    SEIZURE_BUFFER = 4;
%     seizure_buffer_idx = SEIZURE_BUFFER*60*60*fs;
    % create boundaries for each "seizure" in datetime
    seizure_idxs = idx2time(ecogT,ecogT{i_seiz, {'EventStartIdx', 'EventEndIdx'}});
    % add seizure buffer in hours
    seizure_idxs(:,2) = seizure_idxs(:,2) + hours(SEIZURE_BUFFER); 

%     seiz_stims = zeros(size(clean_stims));
%     for i = 1:length(seizure_idxs)
%         % for each seizure, check membership of each stimulation start converted to datetimeand
%         % add to cumulative logical array (should be no value larger than 1)
%         seiz_stims = seiz_stims + logical(ismember(idx2time(ecogT,clean_stims(:,1)),seizure_idxs(i,1):seizure_idxs(i,2)));
%     end
    [~,excl] = filterWindows(idx2time(ecogT,clean_stims),seizure_idxs); %H:
%     no_seiz_stims = ~seiz_stims(:,1);
    no_seiz_stims = excl; %H:
    % Step 7
    % Create window start stops using clean indices
    disp(['There are ' num2str(sum(no_seiz_stims(:,1))) '/' num2str(length(clean_stims)) ' stimulations removed from a seizure']);
    % Use only windows that have clean before and after
    % clean_idxs = logical(pre_clean_idxs.*post_clean_idxs.*stim_clean_idxs.*no_seiz_stims); %%%
    clean_idxs = logical(no_seiz_stims);
    % pre_stim_windows = pre_stim_windows_dirty(clean_idxs,:);
    % post_stim_windows = post_stim_windows_dirty(clean_idxs,:);
    clean_stim_windows = clean_stims(clean_idxs,:);

    % topo_idxs = intersect(topo_cleaning(PatientData,pre_stim_windows,data_window_idx),...
    %     topo_cleaning(PatientData,post_stim_windows,data_window_idx));
    %
    % pre_stim_windows = pre_stim_windows(topo_idxs,:);
    % post_stim_windows = post_stim_windows(topo_idxs,:);

    % extract the indices for the clean stimulation windows
    % clean_stim_windows = clean_stim_windows(topo_idxs,:);
    save([datapath,'/',ptID,'/stim_idxs_',ptID,'.mat'],'clean_stim_windows')
    time_trace = idx2time(ecogT,clean_stim_windows(:,1),...
        'timezone','UTC');
    save([datapath,'/',ptID,'/stim_UTC_time_trace_',ptID,'.mat'],'time_trace')
    ptime_trace = posixtime(time_trace);
    save([datapath,'/',ptID,'/stim_posix_UTC_time_trace_',ptID,'.mat'],'ptime_trace')
    disp(['There are ' num2str(length(clean_stim_windows)) '/' num2str(length(clean_stims)) ' stimulations with no repeats']);
%     continue %H:
    %% Saving analysis windows
    analysis_windows_idxs = [clean_stim_windows(:,2)-data_window_idx, clean_stim_windows(:,2)+data_window_idx];
    analysis_windows = {};
    for i_win = 1:length(analysis_windows_idxs)
        analysis_windows{i_win} = double(PatientData(analysis_windows_idxs(i_win,1):analysis_windows_idxs(i_win,2),:));
    end
    save([datapath,'/',ptID,'/stim_windows_',ptID,'.mat'],'analysis_windows')
    %% Create Programming Event Array
    event_id_table = ecogT(:,{'RawLocalTimestamp','EventStartIdx','EventEndIdx'});
    session_id_table = pdms(:,'Programming_Date'); %H: commented out
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
    save([datapath,'/',ptID,'/visit_selection_array_',ptID,'.mat'],'visit_selection_array')
end