%% stimtraces_preprocess.m
% store stim trace of certain window
close all;
paths;
patient_info = struct2table(load(which('patients_Penn.mat')).patients_Penn);
ptList = {rns_config.patients.ID};
localization = load(fullfile(datapath,"localization.mat")).localization;
base_days = 90;
rank = 3;
windows = [0, 1, 2, 3, 4, 7, 14, 30, 90, 180, 360, 720];
%% load data
for pt = 1:length(localization)
    % Read Patient Data
    ptID = localization(pt).ptID;
    pidx = strcmp(ptID,patient_info.ID);
    disp(['Starting analysis for ',ptID])

    if ~localization(pt).meets_criteria
        continue
    end
    %
    %     outcome = localization(pt).outcome;
    %     outcome_group = localization(pt).outcome_group;
    %     depth = localization(pt).depth;
    %
    load(fullfile(datapath,ptID,['working_data_',num2str(pt),suffix,'.mat']));

    ptime_trace = load(fullfile(datapath,ptID,['posix_UTC_time_trace_',ptID,'.mat'])).ptime_trace;

    % Stim
    [~,~,~, histT] = loadRNSptData(ptID, rns_config);
    stim_data = histT(:,["UTCStartTime","EpisodeStarts"]);
    [~,stim_idxs] = pdist2(posixtime(stim_data{:,"UTCStartTime"}),ptime_trace,"euclidean",'Smallest',1);
    stims = stim_data{stim_idxs,"EpisodeStarts"};
    stim_traces = cell(length(windows),1);
    for j = 1:length(windows)
        if windows(j) == 0
            stim_trace = stims;
            stim_trace(isnan(stim_trace)) = 0;
        else
            tmp = stim_data{:,"EpisodeStarts"};
            tmp(isnan(tmp)) = 0;
            cum_stims = movingCumulSum(tmp,windows(j))'; %number of stimulations in a 90day period
            stim_trace = cum_stims(stim_idxs);
        end
        stim_traces{j} = stim_trace;
    end
    save(fullfile(datapath,ptID,['working_data_',num2str(pt),suffix,'.mat']),'-append','stim_traces');
end