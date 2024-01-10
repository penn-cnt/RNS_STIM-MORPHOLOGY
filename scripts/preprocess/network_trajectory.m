%% network_trajectory.m
% This version replicates paper's method
% This script
% 1) do processing of calculated PLV data to generate PLV in certain
% frequency bands, PLV at baseline (90 days), PLV change from baseline along
% time
% 2) aggreagates data into certain time bins, e.g. data per 45 days
%% Settings
clear; close all;
paths;
patient_info = struct2table(load(which('patients_Penn.mat')).patients_Penn);
ptList = {rns_config.patients.ID};
localization = load(fullfile(datapath,"localization.mat")).localization;

fs = 250;
bin_size = 90; % subject to change
base_days = 90;
%% Connectivity Trajectories
for pt = 1:length(ptList)
    %% Read Patient Data
    ptID = ptList{pt};
    pidx = strcmp(ptID,patient_info.ID);
    disp(['Starting analysis for ',ptID])

    all_plvs = load(fullfile(datapath,ptID,['cwt_plvs_',ptID,'.mat'])).all_plvs;
    freqs = load(fullfile(datapath,ptID,['cwt_plvs_',ptID,'.mat'])).f;
    con_labels = localization(pt).con_labels;
    outcome = localization(pt).outcome;
    outcome_group = localization(pt).outcome_group;
    depth = localization(pt).depth;

    if ~localization(pt).meets_criteria
        continue
    end
    %% Calculate PLV percent change

    % load event data
    dday = patient_info{pidx,"implantDate"};
    time_trace = load(fullfile(datapath,ptID,['UTC_time_trace_',ptID,'.mat'])).time_trace;
    ptime_trace = load(fullfile(datapath,ptID,['posix_UTC_time_trace_',ptID,'.mat'])).ptime_trace;

    % Calculating the network trajecetories combined across a time window
    baseline_period = dday + days(base_days);
    baseline_mask = time_trace < baseline_period;
    % indicates no baseline data or no after baseline data
    if all(baseline_mask) || all(~baseline_mask)
        localization(pt).meets_criteria = 0;
        continue
    end
    implant_time = time_trace - dday; % get relative day of events after implantation
    implant_time = days(implant_time(~baseline_mask)); % convert to day
    time_trace(baseline_mask) = [];
    ptime_trace(baseline_mask) = [];

    % Define frequency bands and determine which bin plv data falls to
    freq_bands_mask = {(freqs > 0 & freqs <= 8),...
        (freqs > 8) & (freqs <= 15),...
        (freqs > 15) & (freqs <= 30),...
        (freqs > 30) & (freqs <= 100)};
    try
        resampled_dplv = load(fullfile(datapath,ptID,['cwt_plvs_',ptID,'.mat'])).resampled_dplv;
        baseline_plv = load(fullfile(datapath,ptID,['cwt_plvs_',ptID,'.mat'])).baseline_plv;
    catch
        % summarize raw plv into defined frequency bands
        resampled_plv = zeros(size(all_plvs,1),size(all_plvs,2),length(freq_bands_mask));
        for mask = 1:length(freq_bands_mask)
            resampled_plv(:,:,mask) = mean(all_plvs(:,:,freq_bands_mask{mask}),3,'omitnan');
        end

        % dplv
        baseline_plv = mean(resampled_plv(baseline_mask,:,:),1,'omitnan');
        dplvs = resampled_plv(~baseline_mask,:,:);
        dplvs = 100 * (dplvs ./ baseline_plv - 1);
        dplvs(abs(dplvs) == inf) = nan;
        
        resampled_dplv = dplvs;

        save(fullfile(datapath,ptID,['cwt_plvs_',ptID,'.mat']),'-append','resampled_plv','resampled_dplv','baseline_plv')
    end
    %% Calculate Rolling Mean
    % which bin the event belong to
    bin_indices = ceil((implant_time-base_days)/bin_size);
    n_bin = max(bin_indices);
    nan_inds = ~ismember([1:n_bin],unique(bin_indices));
    binned_counts = accumarray(bin_indices,1,[],@sum); % number of cleaned events in each bin
    binned_times = accumarray(bin_indices,implant_time,[],@max); % day of the latest event in each bin
    binned_counts(nan_inds) = nan;
    binned_times(nan_inds) = nan;
    binned_dplvs = nan * zeros(n_bin,size(resampled_dplv,2),size(resampled_dplv,3));
    
    % Iterating through each frequency and connection and calculating
    % trajectories heading
    for i = 1:size(resampled_dplv,2) % each chan pair
        for j = 1:size(resampled_dplv,3) % each frequency
            binned_dplvs(:,i,j) = accumarray(bin_indices,resampled_dplv(:,i,j),[],@(x) mean(x,'omitnan'));
        end
    end
    binned_dplvs(nan_inds,:,:) = nan;

    reorg_dplvs = cell(6,4);
    for f = 1:4
        reorg_dplvs{1,f} = binned_dplvs(:,1:2,f);
        reorg_dplvs{2,f} = binned_dplvs(:,3:6,f);
        reorg_dplvs{3,f} = binned_dplvs(:,con_labels == 1,f);
        reorg_dplvs{4,f} = binned_dplvs(:,con_labels == 2,f);
        reorg_dplvs{5,f} = binned_dplvs(:,con_labels == 3,f);
        reorg_dplvs{6,f} = binned_dplvs(:,con_labels == 4,f);
    end

    % store
    plasticity(pt).ptID = ptID;
    plasticity(pt).outcome = outcome;
    plasticity(pt).outcome_group = outcome_group;
    plasticity(pt).depth = depth;
    plasticity(pt).con_labels = con_labels;
    plasticity(pt).times = binned_times;
    plasticity(pt).dplv = binned_dplvs;
    plasticity(pt).reorg_dplv = reorg_dplvs;
%     plasticity(pt).stim_counts = binned_stims;
    plasticity(pt).baseline_plvs = squeeze(baseline_plv);

end
plasticity = plasticity(~cellfun('isempty', {plasticity.ptID}));
save(fullfile(datapath,['plasticity_',num2str(bin_size),'.mat']),"plasticity")
save(fullfile(datapath,"localization.mat"),'localization');