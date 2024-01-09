%% network_trajectory.m
% This version replicates paper's method
% This script
% 1) do processing of calculated PLV data to generate PLV in certain
% frequency bands, PLV at baseline (90 days), PLV change from baseline along
% time
% 2) calculates accumulated stim amount in the previous 90-day window for
% each event
% 3) aggreagates data into certain time bins, e.g. data per 45 days
% 4) collect all relevant information into a data frame
% 5) collapse plv traces along time based on electrode location, outcome,
% frequency, connection type
%% Settings
clear; close all;
paths;
patient_info = struct2table(load(which('patients_Penn.mat')).patients_Penn);
ptList = {rns_config.patients.ID};
localization = load(fullfile(datapath,"localization.mat")).localization;

fs = 250;
bin_size = 90; % subject to change
% ttdays = 1050;
base_days = 90;
n_bin = ceil((ttdays-base_days)/bin_size);
%% Connectivity Trajectories
for pt = 1:length(ptList)
    %% Read Patient Data
    ptID = ptList{pt};
    pidx = strcmp(ptID,patient_info.ID);
    disp(['Starting analysis for ',ptID])

    all_plvs = load(fullfile(datapath,ptID,['cwt_plvs_',ptID,'.mat'])).all_plvs;
    freqs = load(fullfile(datapath,ptID,['cwt_plvs_',ptID,'.mat'])).f;
    con_labels = localization(pt).con_labels;
    outcome = localization(pt).outcome_group;
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

        save(fullfile(datapath,ptID,['cwt_plvs_',ptID,'.mat']),'-append','resampled_plv','resampled_dplv')
    end
    %% Calculate Rolling Mean
%     [~,~,~, histT] = loadRNSptData(ptID, rns_config);
%     stim_data = histT(:,["UTCStartTime","EpisodeStarts"]);
%     stim_time = hours(stim_data{:,"UTCStartTime"} - dday);
%     % calculates number of stims in a fixed time window prior to this
%     % time point
%     cum_stims = movingCumulativeSum(stim_data{:,"EpisodeStarts"},180); %number of stimulations in a 90day period
%     % find out the closest stim time to a clean event
%     [~,stim_idxs] = pdist2(posixtime(stim_data{:,"UTCStartTime"}),ptime_trace,"euclidean",'Smallest',1);
%     % get the accumulated stim for each event
%     stim_trace = cum_stims(stim_idxs);
    % which bin the event belong to
    bin_indices = ceil((implant_time-base_days)/bin_size);
    binned_counts = nan*zeros(n_bin,1);
    binned_times = nan*zeros(n_bin,1);
%     binned_stims = nan*zeros(n_bin,1);
    binned_counts(1:max(bin_indices)) = accumarray(bin_indices,1,[],@sum); % number of cleaned events in each bin
    binned_times(1:max(bin_indices)) = accumarray(bin_indices,implant_time,[],@max); % day of the latest event in each bin
%     binned_stims(1:max(bin_indices)) = accumarray(bin_indices,stim_trace,[],@mean); % average cum stim in the previous period
    binned_dplvs = nan * zeros(n_bin,size(resampled_dplv,2),size(resampled_dplv,3));

    % Iterating through each frequency and connection and calculating
    % trajectories heading
    for i = 1:size(resampled_dplv,2) % each chan pair
        for j = 1:size(resampled_dplv,3) % each frequency
            binned_dplvs(1:max(bin_indices),i,j) = accumarray(bin_indices,resampled_dplv(:,i,j),[],@(x) mean(x,'omitnan'));
        end
    end

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
% %% Predicting Data4
% freq_labels = ["Theta", "Alpha","Beta","Gamma"];
% time_ticks = 90:bin_size:ttdays-1;
% year = 1;
% frequencies = [];
% times = [];
% ids = [];
% base = [];
% plvs = [];
% dplvs = [];
% outcomes = [];
% outcome_groups = [];
% conns = [];
% depth = [];
% for pt = 1:length(ptList)
%     ptID = plasticity(pt).ptID;
%     pidx = strcmp(ptID,patient_info.ID);
%     if ~localization(pt).meets_criteria
%         continue
%     end
%     con_labels = localization(pt).con_labels;
%     d = localization(pt).depth;
%     if isempty(d) || isempty(localization(pt).outcome_group)
%         continue
%     end
%     baseline_plvs = squeeze(plasticity(pt).baseline_plvs);
%     dplv = plasticity(pt).dplv; % time x conn type x freq
%     for i = 1:size(dplv,1)
%         for j = 1:size(dplv,2) % 1 = intra, 2 = inter, 3 = soz-soz, etc.
%             for k = 1:size(dplv,3)
%                 frequencies = [frequencies; k];
%                 ids = [ids; pt];
%                 times = [times; time_ticks(i)];
%                 base = [base; baseline_plvs(j,k)];
%                 dplvs = [dplvs; dplv(i,j,k)];
%                 outcomes = [outcomes; localization(pt).outcome(year)];
%                 outcome_groups = [outcome_groups; localization(pt).outcome_group(year)];
%                 conns = [conns; j];
%                 depth = [depth; d];
%             end
%         end
%     end
% end
% data_mat = [ids,frequencies,base,times,plvs,dplvs,outcomes,outcome_groups,conns,depth];
% data_table = array2table(data_mat,"VariableNames",["IDs","Freqs","Baseline","Time","DPLV","Outcome","Outcome_Group","Connection","Depth"]);
% data_table.Freqs = categorical(data_table.Freqs);
% data_table.Connection = categorical(data_table.Connection);
% data_table.Depth = categorical(data_table.Depth);
% save(fullfile(datapath,['data_table_',num2str(bin_size),'.mat']),"data_table")
% writetable(data_table,fullfile(datapath,['data_table_',num2str(bin_size),'.csv']))
% %% Combining all traces across patients
% plasticity = load(fullfile(datapath,['plasticity_',num2str(bin_size),'.mat'])).plasticity;
% dplv_traces = cell(2,6,2,4); % Lead location x Connection type x Outcome loc x Freq
% for pt = 1:length(ptList)
%     ptID = plasticity(pt).ptID;
%     pidx = strcmp(ptID,patient_info.ID);
%     if ~localization(pt).meets_criteria
%         continue
%     end
%     con_labels = localization(pt).con_labels;
%     d = localization(pt).depth;
%     if isempty(d) || isempty(localization(pt).outcome_group)
%         continue
%     end
%     o = localization(pt).outcome_group(year);
%     baseline_plvs = squeeze(plasticity(pt).baseline_plvs);
%     dplv = plasticity(pt).dplv;
%     for j = 1:size(dplv,2)
%         for k = 1:size(dplv,3)
%             dplv_traces{d,j,o,k} = [dplv_traces{d,j,o,k};dplv(:,j,k)'];
%         end
%     end
% end
% save(fullfile(datapath,['plv_traces_',num2str(bin_size),'.mat']),"dplv_traces")

%% Functions
function movingSum = movingCumulativeSum(inputArray, windowSize)
% Check if the window size is valid
if windowSize <= 0 || windowSize > length(inputArray)
    error('Window size must be a positive integer smaller than or equal to the length of the input array.');
end

conv_factor = exp(-[0:1/24:1e4]/windowSize);
% Initialize the output array to store the moving cumulative sum
movingSum = conv(inputArray,conv_factor,'full');
movingSum = movingSum(1:length(inputArray));
end

function scaledData = scaleDataToMinus1To1(data)
% Find the minimum and maximum values in the data
minValue = min(data);
maxValue = max(data);

% Scale the data to the range [-1, 1]
scaledData = ((data - minValue) * 2) ./ (maxValue - minValue) - 1;
end