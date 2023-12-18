%% network_trajectory.m
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

calc_plasticity = true;
fs = 250;
max_time_length = 0;
bin_size = 45; % subject to change
%% Connectivity Trajectories
if ~calc_plasticity
    max_time_length = 39;
else
    for pt = 1:length(ptList)
        %% Read Patient Data
        ptID = ptList{pt};
        pidx = strcmp(ptID,patient_info.ID);
        disp(['Starting analysis for ',ptID])
        
        % load event data
        dday = patient_info{pidx,"implantDate"};
        time_trace = load(fullfile(datapath,ptID,['UTC_time_trace_',ptID,'.mat'])).time_trace;
        ptime_trace = load(fullfile(datapath,ptID,['posix_UTC_time_trace_',ptID,'.mat'])).ptime_trace;
        all_plvs = load(fullfile(datapath,ptID,['plvs_',ptID,'.mat'])).all_plvs;
        freqs = load(fullfile(datapath,ptID,['plvs_',ptID,'.mat'])).f;
        con_labels = localization(pt).con_labels;
        if ~localization(pt).meets_criteria
            continue
        end
        %% Band Limited Network Trajectories
        % Calculating the network trajecetories combined across a time window
        baseline_period = dday + days(90);
        baseline_mask = time_trace < baseline_period;
        % indicates no baseline data or no after baseline data
        if all(baseline_mask) || all(~baseline_mask)
            continue
        end
        implant_time = time_trace - dday; % get relative day of events after implantation
        implant_time = days(implant_time(~baseline_mask)); % convert to day and only keep events after baseline

        % Define frequency bands and determine which bin plv data falls to
        freq_bands_mask = {(freqs > 0 & freqs <= 8),...
            (freqs > 8) & (freqs <= 15),...
            (freqs > 15) & (freqs <= 30),...
            (freqs > 30) & (freqs <= 100)};
        
        % summarize raw plv into defined frequency bands
        band_limited_plvs = zeros(size(all_plvs,1),size(all_plvs,2),length(freq_bands_mask));
        for mask = 1:length(freq_bands_mask)
            band_limited_plvs(:,:,mask) = mean(all_plvs(:,:,mask),3,'omitnan');
        end

        % calculates baseline PLV (6x4: number of connectivity pairs by freqs)
        baseline_values = mean(band_limited_plvs(baseline_mask,:,:),1,'omitnan');
        % calculates how PLV changes from the baseline along time
        plv = band_limited_plvs(~baseline_mask,:,:);
        dplv = (band_limited_plvs(~baseline_mask,:,:)-baseline_values)./baseline_values;
        %% Loading Stimulation information
        [~,~,~, histT] = loadRNSptData(ptID, rns_config);
        stim_data = histT(:,["UTCStartTime","EpisodeStarts"]);
        % calculates number of stims in a fixed time window prior to this
        % time point
        cum_stims = movingCumulativeSum(stim_data{:,"EpisodeStarts"},90*24); %number of stimulations in a 90day period
        % find out the closest stim time to a clean event
        [~,stim_idxs] = pdist2(posixtime(stim_data{:,"UTCStartTime"}),ptime_trace,"euclidean",'Smallest',1);
        % get the accumulated stim for each event
        stim_trace = cum_stims(stim_idxs);
        % which bin the event belong to
        bin_indices = floor(implant_time/bin_size);
        binned_counts = accumarray(bin_indices,1,[],@sum); % number of cleaned events in each bin
        binned_times = accumarray(bin_indices,implant_time,[],@max); % day of the latest event in each bin
        binned_stims = accumarray(bin_indices,stim_trace(~baseline_mask),[],@mean); % average cum stim in the previous period

        % Updating max # of x day bins across patients
        if length(binned_times) > max_time_length
            max_time_length = length(binned_times);
        end

        binned_dplvs = zeros(length(unique(bin_indices)),size(dplv,2),size(dplv,3));
        binned_plvs = zeros(length(unique(bin_indices)),size(plv,2),size(plv,3));
        
        % Iterating through each frequency and connection and calculating
        % trajectories heading
        for i = 1:size(band_limited_plvs,2) % each connection
            for j = 1:size(dplv,3) % each frequency
                test = accumarray(bin_indices,dplv(:,i,j),[],@mean);
                binned_dplvs(:,i,j) = test(logical(test));
                test = accumarray(bin_indices,plv(:,i,j),[],@mean);
                binned_plvs(:,i,j) = test(logical(test));
            end
        end
        plasticity(pt).ptID = ptID;
        plasticity(pt).con_labels = con_labels;
        plasticity(pt).times = binned_times(logical(binned_times));
        plasticity(pt).dplv = binned_dplvs;
        plasticity(pt).plv = binned_plvs;
        plasticity(pt).stim_counts = binned_stims(logical(binned_stims));
        plasticity(pt).baseline_plvs = baseline_values;
    end
    save(fullfile(datapath,"plasticity.mat"),"plasticity")
end
%% Predicting Data
freq_labels = ["Theta", "Alpha","Beta","Gamma"];
time_ticks = (90 + bin_size):bin_size:(90 + bin_size*max_time_length);
frequencies = [];
times = [];
ids = [];
base = [];
plvs = [];
dplvs = [];
outcomes = [];
outcome_groups = [];
conns = [];
depth = [];
for pt = 1:length(ptList)
    ptID = plasticity(pt).ptID;
    pidx = strcmp(ptID,patient_info.ID);
    if ~localization(pt).meets_criteria
        continue
    end
    con_labels = localization(pt).con_labels;
    d = localization(pt).depth;
    if isempty(d) || isnan(localization(pt).outcome_group(1))
        continue
    end
    baseline_plvs = squeeze(plasticity(pt).baseline_plvs);
    plv = plasticity(pt).plv;
    dplv = plasticity(pt).dplv;
%     times = plasticity(pt).times;
    for i = 1:size(dplv,1)
        for j = 1:size(dplv,2)
            for k = 1:size(dplv,3)
                frequencies = [frequencies; k];
                ids = [ids; pt];
                times = [times; time_ticks(i)];
                base = [base; baseline_plvs(j,k)];
                plvs = [plvs; plv(i,j,k)];
                dplvs = [dplvs; dplv(i,j,k)];
                outcomes = [outcomes; localization(pt).outcome(1)];
                outcome_groups = [outcome_groups; localization(pt).outcome_group(1)];
                conns = [conns; con_labels(j)];
                depth = [depth; d];
            end
        end
    end
end
data_mat = [ids,frequencies,base,times,plvs,dplvs,outcomes,outcome_groups,conns,depth];
data_table = array2table(data_mat,"VariableNames",["IDs","Freqs","Baseline","Time","PLV","DPLV","Outcome","Outcome_Group","Connection","Depth"]);
data_table.Freqs = categorical(data_table.Freqs);
data_table.Connection = categorical(data_table.Connection);
data_table.Depth = categorical(data_table.Depth);
data_table = data_table(data_table.Freqs == categorical(1),:);
save(fullfile(datapath,"data_table.mat"),"data_table")

% mdl = fitlme(data_table,"DPLV ~ Baseline + Connection + Outcome + Depth + (1|IDs)");
%% Combining all traces across patients
% TODO: Add implant depth as a covariate
plasticity = load(fullfile(datapath,"plasticity.mat")).plasticity;
plv_traces = cell(2,4,2,4); % Lead location x Connection type x Outcome loc x Freq 
dplv_traces = cell(2,4,2,4); % Lead location x Connection type x Outcome loc x Freq 
for pt = 1:length(ptList)
    ptID = plasticity(pt).ptID;
    pidx = strcmp(ptID,patient_info.ID);
    if ~localization(pt).meets_criteria
        continue
    end
    con_labels = localization(pt).con_labels;
    d = localization(pt).depth;
    o = localization(pt).outcome_group(1);
    if isempty(d) || isnan(o)
        continue
    end
    baseline_plvs = squeeze(plasticity(pt).baseline_plvs);
    plv = padarray(plasticity(pt).plv,[max_time_length-size(plasticity(pt).dplv,1),0,0],nan,'post');
    dplv = padarray(plasticity(pt).dplv,[max_time_length-size(plasticity(pt).dplv,1),0,0],nan,'post');
    for j = 1:size(dplv,2)
        conns = con_labels(j);
        for k = 1:size(dplv,3)
            plv_traces{d,conns,o,k} = [plv_traces{d,conns,o,k};plv(:,j,k)'];
            dplv_traces{d,conns,o,k} = [dplv_traces{d,conns,o,k};dplv(:,j,k)'];
        end
    end
end
save(fullfile(datapath,"plv_traces.mat"),"plv_traces","dplv_traces")

%% Functions
function movingSum = movingCumulativeSum(inputArray, windowSize)
% Check if the window size is valid
if windowSize <= 0 || windowSize > length(inputArray)
    error('Window size must be a positive integer smaller than or equal to the length of the input array.');
end

% Initialize the output array to store the moving cumulative sum
movingSum = zeros(1, length(inputArray));

% Calculate the cumulative sum for the first window
movingSum(1:windowSize) = cumsum(inputArray(1:windowSize),"omitnan");

% Iterate through the rest of the array to update the moving sum
for i = windowSize+1:length(inputArray)
    % Update the moving sum by adding the current element and subtracting the element that is now out of the window
    movingSum(i) = sum(inputArray((i-windowSize+1):i),"omitnan");
    %         movingSum(i) = movingSum(i-1) + inputArray(i) - inputArray(i - windowSize);
    if isnan(movingSum(i))
        disp('x')
    end
end
end

function scaledData = scaleDataToMinus1To1(data)
% Find the minimum and maximum values in the data
minValue = min(data);
maxValue = max(data);

% Scale the data to the range [-1, 1]
scaledData = ((data - minValue) * 2) ./ (maxValue - minValue) - 1;
end