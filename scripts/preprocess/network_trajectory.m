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
bin_size = 90; % subject to change
base_days = 90;
regress_dist = 0;
if regress_dist
    suffix = '_regdist';
else
    suffix = '';
end
%% Connectivity Trajectories
for pt = 1:length(localization)
    % Read Patient Data
    ptID = localization(pt).ptID;
    pidx = strcmp(ptID,patient_info.ID);
    disp(['Starting analysis for ',ptID])

    lead_labels = localization(pt).lead_labels;
    con_labels = localization(pt).con_labels;
    outcome = localization(pt).outcome;
    outcome_group = localization(pt).outcome_group;
    depth = localization(pt).depth;
    distance = localization(pt).ch_distances;

%     try
%         load(fullfile(datapath,ptID,['working_data_',num2str(pt),suffix,'.mat']));
%     catch
        %% Calculate PLV percent change
        all_plvs = load(fullfile(datapath,ptID,['cwt_plvs_',ptID,'.mat'])).all_plvs;
        freqs = load(fullfile(datapath,ptID,['cwt_plvs_',ptID,'.mat'])).f;
        
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

        implant_time_all = days(time_trace - dday); % get relative day of events after implantation
        implant_time = implant_time_all(~baseline_mask); % convert to day

        % Define frequency bands and determine which bin plv data falls to
        freq_bands_mask = {(freqs > 0 & freqs <= 8),...
            (freqs > 8) & (freqs <= 15),...
            (freqs > 15) & (freqs <= 30),...
            (freqs > 30) & (freqs <= 100)};

        % summarize raw plv into defined frequency bands
        resampled_plv = zeros(size(all_plvs,1),size(all_plvs,2),length(freq_bands_mask));
        for mask = 1:length(freq_bands_mask)
            resampled_plv(:,:,mask) = mean(all_plvs(:,:,freq_bands_mask{mask}),3,'omitnan');
        end

        if regress_dist
            % regress here
            if isempty(distance)
                continue
            else
                connections =  [1,2; 3,4; 1,3; 1,4; 2,3; 2,4];% diff
                dist = nan * zeros(1,length(connections));
                for conn = 1:length(connections)
                    dist(conn) = distance(connections(conn,1),connections(conn,2));
                end
                dist = 1./dist.^2;
                tmp_dist = repmat(dist,[size(resampled_plv,1),1,size(resampled_plv,3)]);
%                 tmp_plv = reshape(resampled_plv,[size(resampled_plv,1)*size(resampled_plv,3),6]);
%                 dist = repmat(dist,[size(tmp_plv,1),1]);
                for freq = 1:size(resampled_plv,3)
                    mdl = fitlm(reshape(tmp_dist(:,:,freq),[],1),...
                                reshape(resampled_plv(:,:,freq),[],1));
                    pred = predict(mdl, reshape(tmp_dist(:,:,freq),[],1));
                    residuals = reshape(resampled_plv(:,:,freq),[],1) - pred;
                    residuals = (residuals - min(residuals))/(max(residuals) - min(residuals));
                    resampled_plv(:,:,freq) = reshape(residuals,size(resampled_plv,1),size(resampled_plv,2));
                end

            end
        end
        
        plv_intra = squeeze(mean(resampled_plv(:,lead_labels == 1,:),2,'omitnan'));
        plv_inter = squeeze(mean(resampled_plv(:,lead_labels == 2,:),2,'omitnan'));
        
        % si_ratio
        if ~ismember(1,lead_labels) || ~ismember(2,lead_labels)
            si_ratio = nan * zeros(size(resampled_plv,1),size(resampled_plv,3));
        else
            si_ratio = plv_inter./plv_intra;
            si_ratio(abs(si_ratio) == inf) = nan;
        end

        baseline_siratio = mean(si_ratio(baseline_mask,:),1,'omitnan');
        si_ratio = si_ratio(~baseline_mask,:);
        
        % dplv
        baseline_plv = mean(resampled_plv(baseline_mask,:,:),1,'omitnan');
        std_plv = std(resampled_plv(baseline_mask,:,:),[],1,'omitnan');
        resampled_plv_all = resampled_plv;
        resampled_plv = resampled_plv(~baseline_mask,:,:);
        dplvs = 100 * (resampled_plv ./ baseline_plv - 1);
        dplvs(abs(dplvs) == inf) = nan;

        resampled_dplv = dplvs;

        % int-seg difference
        if ~ismember(1,lead_labels) || ~ismember(2,lead_labels)
            int_seg = nan * zeros(size(resampled_dplv,1),size(resampled_dplv,3));
        else
            dplv_intra = squeeze(mean(resampled_dplv(:,lead_labels == 1,:),2,'omitnan'));
            dplv_inter = squeeze(mean(resampled_dplv(:,lead_labels == 2,:),2,'omitnan'));
            int_seg = dplv_inter-dplv_intra;
        end

        % zscored_plv
        zscored_plv = (resampled_plv - baseline_plv)./ std_plv;

        % plv_slopes
        % Get Predicting Slope
        plv_slopes = nan * zeros(size(resampled_plv));
        n_baseline = find(baseline_mask);
        for d = 1:size(resampled_plv,2)
            for i = 1:size(resampled_plv,1)
                tmp_data = resampled_plv_all(1:i+n_baseline,:);
                tmp_time = implant_time_all(1:i+n_baseline);
                for j = 1:4 % each freq band
                    p = polyfit(tmp_time,zscore(tmp_data(:,j)),1);
                    plv_slopes(i,d,j) = p(1);
                end
            end
        end

        save(fullfile(datapath,ptID,['working_data_',num2str(pt),suffix,'.mat']), ...
            'resampled_plv_all','resampled_plv','resampled_dplv','baseline_plv','implant_time', 'implant_time_all', ...
            'si_ratio','baseline_siratio','int_seg','zscored_plv','plv_slopes');
%     end
    %% Calculate Rolling Mean
    % which bin the event belong to
    bin_indices = ceil((implant_time-base_days)/bin_size);
    n_bin = max(bin_indices);
    nan_inds = ~ismember([1:n_bin],unique(bin_indices));
    binned_counts = accumarray(bin_indices,1,[],@sum); % number of cleaned events in each bin
    binned_times = accumarray(bin_indices,implant_time,[],@max); % day of the latest event in each bin
    binned_counts(nan_inds) = nan;
    binned_times(nan_inds) = nan;

    binned_dplv = nan * zeros(n_bin,size(resampled_dplv,2),size(resampled_dplv,3));
    binned_zplv = nan * zeros(n_bin,size(resampled_dplv,2),size(resampled_dplv,3));
    binned_plvslopes = nan * zeros(n_bin,size(resampled_dplv,2),size(resampled_dplv,3));
    % Iterating through each frequency and connection and calculating
    % trajectories heading
    for i = 1:size(resampled_dplv,2) % each chan pair
        for j = 1:size(resampled_dplv,3) % each frequency
            binned_dplv(:,i,j) = accumarray(bin_indices,resampled_dplv(:,i,j),[],@(x) mean(x,'omitnan'));
            binned_zplv(:,i,j) = accumarray(bin_indices,zscored_plv(:,i,j),[],@(x) mean(x,'omitnan'));
            binned_plvslopes(:,i,j) = accumarray(bin_indices,plv_slopes(:,i,j),[],@(x) mean(x,'omitnan'));
        end
    end
    binned_dplv(nan_inds,:,:) = nan;
    binned_zplv(nan_inds,:,:) = nan;
    binned_plvslopes(nan_inds,:,:) = nan;

    binned_siratio = nan * zeros(n_bin,size(resampled_dplv,3));
    binned_intseg = nan * zeros(n_bin,size(resampled_dplv,3));
    for j = 1:size(resampled_dplv,3)
        binned_siratio(:,j) = accumarray(bin_indices,si_ratio(:,j),[],@(x) mean(x,'omitnan'));
        binned_intseg(:,j) = accumarray(bin_indices,int_seg(:,j),[],@(x) mean(x,'omitnan'));
    end
    binned_siratio(nan_inds,:) = nan;
    binned_intseg(nan_inds,:) = nan;

    reorg_dplvs = cell(6,4);
    for f = 1:4
        reorg_dplvs{1,f} = binned_dplv(:,lead_labels == 1,f);
        reorg_dplvs{2,f} = binned_dplv(:,lead_labels == 2,f);
        reorg_dplvs{3,f} = binned_dplv(:,con_labels == 1,f);
        reorg_dplvs{4,f} = binned_dplv(:,con_labels == 2,f);
        reorg_dplvs{5,f} = binned_dplv(:,con_labels == 3,f);
        reorg_dplvs{6,f} = binned_dplv(:,con_labels == 4,f);
    end

    reorg_zplvs = cell(6,4);
    for f = 1:4
        reorg_zplvs{1,f} = binned_zplv(:,lead_labels == 1,f);
        reorg_zplvs{2,f} = binned_zplv(:,lead_labels == 2,f);
        reorg_zplvs{3,f} = binned_zplv(:,con_labels == 1,f);
        reorg_zplvs{4,f} = binned_zplv(:,con_labels == 2,f);
        reorg_zplvs{5,f} = binned_zplv(:,con_labels == 3,f);
        reorg_zplvs{6,f} = binned_zplv(:,con_labels == 4,f);
    end

    reorg_slopes = cell(6,4);
    for f = 1:4
        reorg_slopes{1,f} = binned_plvslopes(:,lead_labels == 1,f);
        reorg_slopes{2,f} = binned_plvslopes(:,lead_labels == 2,f);
        reorg_slopes{3,f} = binned_plvslopes(:,con_labels == 1,f);
        reorg_slopes{4,f} = binned_plvslopes(:,con_labels == 2,f);
        reorg_slopes{5,f} = binned_plvslopes(:,con_labels == 3,f);
        reorg_slopes{6,f} = binned_plvslopes(:,con_labels == 4,f);
    end
    
    %%%
    % inserted block to calculate variance of PLVs
    % Calculating variance
    % Goal is to understand how underlying connectivity (baseline or structural) might moderate the responsiveness of a function connection to stimulation. The hypothesis is that more highly connected electrodes will have a more stable response (continue to increase/decrease when stimulated more and more) than weakly connected regions.
    % Calculating variance of plv(%dplv) within time bins
    % Calculating variance across times bins
    % Calculating variance in slope between time bins
    % Calculating variance in slope between events
    % Event level %dplv, calculate the slope from one event to the next
    % This “level” of analysis is new, we haven’t done this before
    % variance of PLV within time bin
    binned_plv_var = nan * zeros(n_bin,size(resampled_plv,2),size(resampled_plv,3));
    binned_dplv_var = nan * zeros(n_bin,size(resampled_dplv,2),size(resampled_dplv,3));
    % Iterating through each frequency and connection and calculating
    % trajectories heading
    for i = 1:size(resampled_dplv,2) % each chan pair
        for j = 1:size(resampled_dplv,3) % each frequency
            binned_plv_var(:,i,j) = accumarray(bin_indices,resampled_plv(:,i,j),[],@(x) var(x,'omitnan'));
            binned_dplv_var(:,i,j) = accumarray(bin_indices,resampled_dplv(:,i,j),[],@(x) var(x,'omitnan'));
        end
    end
    binned_plv_var(nan_inds,:,:) = nan;
    binned_dplv_var(nan_inds,:,:) = nan;
    % variance of PLV across time bin
    dplv_var = var(binned_dplv,1,'omitnan');
    zplv_var = var(binned_zplv,1,'omitnan');
    % Calculating variance in slope between time bins
    slope = diff(binned_dplv,[],1);
    dplv_slope_var = var(slope,1,'omitnan');
    slope = diff(binned_zplv,[],1);
    zplv_slope_var = var(slope,1,'omitnan');
    % Calculating variance in slope between events
    slope_plv = diff(resampled_plv,[],1);
    time_lapse = diff(implant_time/365,[],1);
    slope_plv = slope_plv./time_lapse;
    slope_plv = [nan * zeros(1,6,4);slope_plv];
    slope_dplv = diff(resampled_dplv,[],1);
    slope_dplv = slope_dplv./time_lapse;
    slope_dplv = [nan * zeros(1,6,4);slope_dplv];
    binned_plv_slope_var = nan * zeros(n_bin,size(resampled_plv,2),size(resampled_plv,3));
    binned_dplv_slope_var = nan * zeros(n_bin,size(resampled_dplv,2),size(resampled_dplv,3));
    % Iterating through each frequency and connection and calculating
    % trajectories heading
    for i = 1:size(resampled_dplv,2) % each chan pair
        for j = 1:size(resampled_dplv,3) % each frequency
            binned_plv_slope_var(:,i,j) = accumarray(bin_indices,slope_plv(:,i,j),[],@(x) var(x,'omitnan'));
            binned_dplv_slope_var(:,i,j) = accumarray(bin_indices,slope_dplv(:,i,j),[],@(x) var(x,'omitnan'));
        end
    end
    binned_plv_slope_var(nan_inds,:,:) = nan;
    binned_dplv_slope_var(nan_inds,:,:) = nan;
    %%% end of insertion
    

    % store
    plasticity(pt).ptID = ptID;
    plasticity(pt).outcome = outcome;
    plasticity(pt).outcome_group = outcome_group;
    plasticity(pt).depth = depth;
    plasticity(pt).con_labels = lead_labels;
    plasticity(pt).times = binned_times;
    plasticity(pt).dplv = binned_dplv;
    plasticity(pt).reorg_dplv = reorg_dplvs;
    plasticity(pt).baseline_plvs = squeeze(baseline_plv);
    plasticity(pt).baseline_siratio = baseline_siratio;
    plasticity(pt).si_ratio = binned_siratio;
    plasticity(pt).int_seg = binned_intseg;
    plasticity(pt).zplv = binned_zplv;
    plasticity(pt).reorg_zplv = reorg_zplvs;
    plasticity(pt).plv_slopes = binned_plvslopes;
    plasticity(pt).reorg_slope = reorg_slopes;
    % additional store
    plasticity(pt).plv_var = binned_plv_var;
    plasticity(pt).dplv_var = binned_dplv_var;
    plasticity(pt).dplv_var_acrossbin = dplv_var;
    plasticity(pt).zplv_var_acrossbin = zplv_var;
    plasticity(pt).plv_slope_var = binned_plv_slope_var;
    plasticity(pt).dplv_slope_var = binned_dplv_slope_var;
    plasticity(pt).dplv_slope_var_acrossbin = dplv_slope_var;
    plasticity(pt).zplv_slope_var_acrossbin = zplv_slope_var;
    %
end
save(fullfile(datapath,['plasticity_',num2str(bin_size),suffix,'.mat']),"plasticity")
save(fullfile(datapath,"localization.mat"),'localization');