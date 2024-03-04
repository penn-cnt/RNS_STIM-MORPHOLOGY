%% network_trajectory.m
% Calculates all required feature traces and do aggregation into N-day time
% bins.
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
freq_strings = {'Theta','Alpha','Beta','Gamma'};
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

    try
        load(fullfile(datapath,ptID,['working_data_',num2str(pt),suffix,'.mat']));
    catch
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
        baseline_std = std(resampled_plv(baseline_mask,:,:),[],1,'omitnan');
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
        zscored_plv = (resampled_plv - baseline_plv)./ baseline_std;

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
            'si_ratio','baseline_siratio','int_seg','zscored_plv','plv_slopes','baseline_mask');
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

    binned_plv = nan * zeros(n_bin,size(resampled_dplv,2),size(resampled_dplv,3));
    binned_dplv = nan * zeros(n_bin,size(resampled_dplv,2),size(resampled_dplv,3));
    binned_zplv = nan * zeros(n_bin,size(resampled_dplv,2),size(resampled_dplv,3));
    % Iterating through each frequency and connection and calculating
    % trajectories heading
    for i = 1:size(resampled_dplv,2) % each chan pair
        for j = 1:size(resampled_dplv,3) % each frequency
            binned_plv(:,i,j) = accumarray(bin_indices,resampled_plv(:,i,j),[],@(x) mean(x,'omitnan'));
            binned_dplv(:,i,j) = accumarray(bin_indices,resampled_dplv(:,i,j),[],@(x) mean(x,'omitnan'));
            binned_zplv(:,i,j) = accumarray(bin_indices,zscored_plv(:,i,j),[],@(x) mean(x,'omitnan'));
        end
    end
    binned_plv(nan_inds,:,:) = nan;
    binned_dplv(nan_inds,:,:) = nan;
    binned_zplv(nan_inds,:,:) = nan;

    %%%
    % inserted block to calculate std of PLVs

    % std of PLV/dPLV/zPLV event/time bin
    % within time bin, i.e. event
    binned_plv_std = nan * zeros(n_bin,size(resampled_dplv,2),size(resampled_dplv,3));
    binned_dplv_std = nan * zeros(n_bin,size(resampled_dplv,2),size(resampled_dplv,3));
    binned_zplv_std = nan * zeros(n_bin,size(resampled_dplv,2),size(resampled_dplv,3));
    % Iterating through each frequency and connection and calculating
    % trajectories heading
    for i = 1:size(resampled_dplv,2) % each chan pair
        for j = 1:size(resampled_dplv,3) % each frequency
            binned_plv_std(:,i,j) = accumarray(bin_indices,resampled_plv(:,i,j),[],@(x) std(x,'omitnan'));
            binned_dplv_std(:,i,j) = accumarray(bin_indices,resampled_dplv(:,i,j),[],@(x) std(x,'omitnan'));
            binned_zplv_std(:,i,j) = accumarray(bin_indices,zscored_plv(:,i,j),[],@(x) std(x,'omitnan'));
        end
    end
    binned_plv_std(nan_inds,:,:) = nan;
    binned_dplv_std(nan_inds,:,:) = nan;
    binned_zplv_std(nan_inds,:,:) = nan;
    % across time bin
    plv_std = squeeze(std(binned_plv,1,'omitnan'));
    dplv_std = squeeze(std(binned_dplv,1,'omitnan'));
    zplv_std = squeeze(std(binned_zplv,1,'omitnan'));
    % Calculating mean/stdiance in slope between event/time bins
    % within time bin
    slope_plv = diff(resampled_plv,[],1);
    time_lapse = diff(implant_time,[],1);
    slope_plv = slope_plv./time_lapse;
    slope_plv = [nan * zeros(1,6,4);slope_plv];
    slope_dplv = diff(resampled_dplv,[],1);
    slope_dplv = slope_dplv./time_lapse;
    slope_dplv = [nan * zeros(1,6,4);slope_dplv];
    slope_zplv = diff(zscored_plv,[],1);
    slope_zplv = slope_zplv./time_lapse;
    slope_zplv = [nan * zeros(1,6,4);slope_zplv];
    binned_plv_slope_mean = nan * zeros(n_bin,size(resampled_plv,2),size(resampled_plv,3));
    binned_dplv_slope_mean = nan * zeros(n_bin,size(resampled_dplv,2),size(resampled_dplv,3));
    binned_plv_slope_std = nan * zeros(n_bin,size(resampled_plv,2),size(resampled_plv,3));
    binned_dplv_slope_std = nan * zeros(n_bin,size(resampled_dplv,2),size(resampled_dplv,3));
    binned_zplv_slope_mean = nan * zeros(n_bin,size(resampled_plv,2),size(resampled_plv,3));
    binned_zplv_slope_std = nan * zeros(n_bin,size(resampled_dplv,2),size(resampled_dplv,3));
    
    % Iterating through each frequency and connection and calculating
    % trajectories heading
    for i = 1:size(resampled_dplv,2) % each chan pair
        for j = 1:size(resampled_dplv,3) % each frequency
            binned_plv_slope_mean(:,i,j) = accumarray(bin_indices,slope_plv(:,i,j),[],@(x) mean(x,'omitnan'));
            binned_dplv_slope_mean(:,i,j) = accumarray(bin_indices,slope_dplv(:,i,j),[],@(x) mean(x,'omitnan'));
            binned_plv_slope_std(:,i,j) = accumarray(bin_indices,slope_plv(:,i,j),[],@(x) std(x,'omitnan'));
            binned_dplv_slope_std(:,i,j) = accumarray(bin_indices,slope_dplv(:,i,j),[],@(x) std(x,'omitnan'));
            binned_zplv_slope_mean(:,i,j) = accumarray(bin_indices,slope_zplv(:,i,j),[],@(x) mean(x,'omitnan'));
            binned_zplv_slope_std(:,i,j) = accumarray(bin_indices,slope_zplv(:,i,j),[],@(x) std(x,'omitnan'));
        end
    end
    binned_plv_slope_mean(nan_inds,:,:) = nan;
    binned_dplv_slope_mean(nan_inds,:,:) = nan;
    binned_plv_slope_std(nan_inds,:,:) = nan;
    binned_dplv_slope_std(nan_inds,:,:) = nan;
    binned_zplv_slope_mean(nan_inds,:,:) = nan;
    binned_zplv_slope_std(nan_inds,:,:) = nan;
    % across time bin
    slope = diff(binned_plv,[],1);
    plv_slope_mean = squeeze(mean(slope,1,'omitnan'));
    plv_slope_std = squeeze(std(slope,1,'omitnan'));
    slope = diff(binned_dplv,[],1);
    dplv_slope_mean = squeeze(mean(slope,1,'omitnan'));
    dplv_slope_std = squeeze(std(slope,1,'omitnan'));
    slope = diff(binned_zplv,[],1);
    zplv_slope_mean = squeeze(mean(slope,1,'omitnan'));
    zplv_slope_std = squeeze(std(slope,1,'omitnan'));
    
    % fit line and get Rsquared and p, as other way to estimate stdiance
    % event level
    fit_plv_R2 = nan * zeros(size(resampled_dplv,2),size(resampled_dplv,3));
    fit_plv_p = nan * zeros(size(resampled_dplv,2),size(resampled_dplv,3));
    fit_plv_coeff = nan * zeros(size(resampled_dplv,2),size(resampled_dplv,3));
    fit_plv_pcoeff = nan * zeros(size(resampled_dplv,2),size(resampled_dplv,3));
    full_result = [];
    for i = 1:size(resampled_dplv,2)
        for j = 1:size(resampled_dplv,3)
            tbl = fitglmBS(implant_time,resampled_plv(:,i,j),100,true);
            fit_plv_R2(i,j) = tbl{1,1};
            fit_plv_p(i,j) = tbl{1,2};
            fit_plv_coeff(i,j) = tbl{2,1};
            fit_plv_pcoeff(i,j) = tbl{2,2};
            tmp_result = [tbl{1,:},tbl{2,:}];
            full_result = [full_result;tmp_result];
        end
    end
    full_result = array2table(full_result, 'RowNames',strcat(cellstr(num2str(reshape(repmat([1:6],4,1),[],1))),{'_'},repmat(freq_strings',6,1)), ...
            'VariableNames',{'R2','pR2','R2_lowerCI','R2_higherCI','b','pCoeff','b_lowerCI','b_higherCI'});
    writetable(full_result,fullfile(datapath,'stats',[ptID,'_PLV_time.csv']),'WriteRowNames',true)

    fit_dplv_R2 = nan * zeros(size(resampled_dplv,2),size(resampled_dplv,3));
    fit_dplv_p = nan * zeros(size(resampled_dplv,2),size(resampled_dplv,3));
    fit_dplv_coeff = nan * zeros(size(resampled_dplv,2),size(resampled_dplv,3));
    fit_dplv_pcoeff = nan * zeros(size(resampled_dplv,2),size(resampled_dplv,3));
    full_result = [];
    for i = 1:size(resampled_dplv,2)
        for j = 1:size(resampled_dplv,3)
            tbl = fitglmBS(implant_time,resampled_dplv(:,i,j),100,false);
            fit_dplv_R2(i,j) = tbl{1,1};
            fit_dplv_p(i,j) = tbl{1,2};
            fit_dplv_coeff(i,j) = tbl{2,1};
            fit_dplv_pcoeff(i,j) = tbl{2,2};
            tmp_result = [tbl{1,:},tbl{2,:}];
            full_result = [full_result;tmp_result];
        end
    end
    full_result = array2table(full_result, 'RowNames',strcat(cellstr(num2str(reshape(repmat([1:6],4,1),[],1))),{'_'},repmat(freq_strings',6,1)), ...
            'VariableNames',{'R2','pR2','R2_lowerCI','R2_higherCI','b','pCoeff','b_lowerCI','b_higherCI'});
    writetable(full_result,fullfile(datapath,'stats',[ptID,'_dPLV_time.csv']),'WriteRowNames',true)

    fit_plv_R2_bin = nan * zeros(size(resampled_dplv,2),size(resampled_dplv,3));
    fit_plv_p_bin = nan * zeros(size(resampled_dplv,2),size(resampled_dplv,3));
    fit_plv_coeff_bin = nan * zeros(size(resampled_dplv,2),size(resampled_dplv,3));
    fit_plv_pcoeff_bin = nan * zeros(size(resampled_dplv,2),size(resampled_dplv,3));
    full_result = [];
    for i = 1:size(resampled_dplv,2)
        for j = 1:size(resampled_dplv,3)
            tbl = fitglmBS([1:n_bin]',binned_plv(:,i,j),100,true);
            fit_plv_R2_bin(i,j) = tbl{1,1};
            fit_plv_p_bin(i,j) = tbl{1,2};
            fit_plv_coeff_bin(i,j) = tbl{2,1};
            fit_plv_pcoeff_bin(i,j) = tbl{2,2};
            tmp_result = [tbl{1,:},tbl{2,:}];
            full_result = [full_result;tmp_result];
        end
    end
    full_result = array2table(full_result, 'RowNames',strcat(cellstr(num2str(reshape(repmat([1:6],4,1),[],1))),{'_'},repmat(freq_strings',6,1)), ...
            'VariableNames',{'R2','pR2','R2_lowerCI','R2_higherCI','b','pCoeff','b_lowerCI','b_higherCI'});
    writetable(full_result,fullfile(datapath,'stats',[ptID,'_PLV_time_bin.csv']),'WriteRowNames',true)


    fit_dplv_R2_bin = nan * zeros(size(resampled_dplv,2),size(resampled_dplv,3));
    fit_dplv_p_bin = nan * zeros(size(resampled_dplv,2),size(resampled_dplv,3));
    fit_dplv_coeff_bin = nan * zeros(size(resampled_dplv,2),size(resampled_dplv,3));
    fit_dplv_pcoeff_bin = nan * zeros(size(resampled_dplv,2),size(resampled_dplv,3));
    full_result = [];
    for i = 1:size(resampled_dplv,2)
        for j = 1:size(resampled_dplv,3)
            tbl = fitglmBS([1:n_bin]',binned_dplv(:,i,j),100,false);
            fit_dplv_R2_bin(i,j) = tbl{1,1};
            fit_dplv_p_bin(i,j) = tbl{1,2};
            fit_dplv_coeff_bin(i,j) = tbl{2,1};
            fit_dplv_pcoeff_bin(i,j) = tbl{2,2};
            tmp_result = [tbl{1,:},tbl{2,:}];
            full_result = [full_result;tmp_result];
        end
    end
    full_result = array2table(full_result, 'RowNames',strcat(cellstr(num2str(reshape(repmat([1:6],4,1),[],1))),{'_'},repmat(freq_strings',6,1)), ...
            'VariableNames',{'R2','pR2','R2_lowerCI','R2_higherCI','b','pCoeff','b_lowerCI','b_higherCI'});
    writetable(full_result,fullfile(datapath,'stats',[ptID,'_dPLV_time_bin.csv']),'WriteRowNames',true)


    % fit line and get Rsquared and p, as other way to estimate stdiance
    % event level
    fitstim_plv_R2 = nan * zeros(size(resampled_dplv,2),size(resampled_dplv,3));
    fitstim_plv_p = nan * zeros(size(resampled_dplv,2),size(resampled_dplv,3));
    fitstim_plv_coeff = nan * zeros(size(resampled_dplv,2),size(resampled_dplv,3));
    fitstim_plv_pcoeff = nan * zeros(size(resampled_dplv,2),size(resampled_dplv,3));
    full_result = [];
    for i = 1:size(resampled_dplv,2)
        for j = 1:size(resampled_dplv,3)
            tbl = fitglmBS(stim_traces{10}(~baseline_mask),resampled_plv(:,i,j),100,true);
            fitstim_plv_R2(i,j) = tbl{1,1};
            fitstim_plv_p(i,j) = tbl{1,2};
            fitstim_plv_coeff(i,j) = tbl{2,1};
            fitstim_plv_pcoeff(i,j) = tbl{2,2};
            tmp_result = [tbl{1,:},tbl{2,:}];
            full_result = [full_result;tmp_result];
        end
    end
    full_result = array2table(full_result, 'RowNames',strcat(cellstr(num2str(reshape(repmat([1:6],4,1),[],1))),{'_'},repmat(freq_strings',6,1)), ...
            'VariableNames',{'R2','pR2','R2_lowerCI','R2_higherCI','b','pCoeff','b_lowerCI','b_higherCI'});
    writetable(full_result,fullfile(datapath,'stats',[ptID,'_PLV_stim.csv']),'WriteRowNames',true)


    fitstim_dplv_R2 = nan * zeros(size(resampled_dplv,2),size(resampled_dplv,3));
    fitstim_dplv_p = nan * zeros(size(resampled_dplv,2),size(resampled_dplv,3));
    fitstim_dplv_coeff = nan * zeros(size(resampled_dplv,2),size(resampled_dplv,3));
    fitstim_dplv_pcoeff = nan * zeros(size(resampled_dplv,2),size(resampled_dplv,3));
    full_result = [];
    for i = 1:size(resampled_dplv,2)
        for j = 1:size(resampled_dplv,3)
            tbl = fitglmBS(stim_traces{10}(~baseline_mask),resampled_dplv(:,i,j),100,false);
            fitstim_dplv_R2(i,j) = tbl{1,1};
            fitstim_dplv_p(i,j) = tbl{1,2};
            fitstim_dplv_coeff(i,j) = tbl{2,1};
            fitstim_dplv_pcoeff(i,j) = tbl{2,2};
            tmp_result = [tbl{1,:},tbl{2,:}];
            full_result = [full_result;tmp_result];
        end
    end
    full_result = array2table(full_result, 'RowNames',strcat(cellstr(num2str(reshape(repmat([1:6],4,1),[],1))),{'_'},repmat(freq_strings',6,1)), ...
            'VariableNames',{'R2','pR2','R2_lowerCI','R2_higherCI','b','pCoeff','b_lowerCI','b_higherCI'});
    writetable(full_result,fullfile(datapath,'stats',[ptID,'_dPLV_stim.csv']),'WriteRowNames',true)

    %%% end of insertion
    
    % store
    plasticity(pt).ptID = ptID;
    plasticity(pt).outcome = outcome;
    plasticity(pt).outcome_group = outcome_group;
    plasticity(pt).depth = depth;
    plasticity(pt).con_labels = lead_labels; % costdiate now
    plasticity(pt).times = binned_times;
    plasticity(pt).plv = binned_plv;
    plasticity(pt).dplv = binned_dplv;
    plasticity(pt).zplv = binned_zplv;
    plasticity(pt).baseline_plvs = squeeze(baseline_plv);

    % additional store
    plasticity(pt).plv_std = binned_plv_std;
    plasticity(pt).dplv_std = binned_dplv_std;
    plasticity(pt).zplv_std = binned_zplv_std;
    plasticity(pt).plv_std_bin = plv_std;
    plasticity(pt).dplv_std_bin = dplv_std;
    plasticity(pt).zplv_std_bin = zplv_std;
    plasticity(pt).plv_slope = binned_plv_slope_mean;
    plasticity(pt).dplv_slope = binned_dplv_slope_mean;
    plasticity(pt).zplv_slope = binned_zplv_slope_mean;
    plasticity(pt).plv_slope_std = binned_plv_slope_std;
    plasticity(pt).dplv_slope_std = binned_dplv_slope_std;
    plasticity(pt).zplv_slope_std = binned_zplv_slope_std;
    plasticity(pt).plv_slope_bin = plv_slope_mean;
    plasticity(pt).dplv_slope_bin = dplv_slope_mean;
    plasticity(pt).zplv_slope_bin = zplv_slope_mean;
    plasticity(pt).plv_slope_std_bin = plv_slope_std;
    plasticity(pt).dplv_slope_std_bin = dplv_slope_std;
    plasticity(pt).zplv_slope_std_bin = zplv_slope_std;
    plasticity(pt).R2_time = fit_plv_R2;
    plasticity(pt).p_time = fit_plv_p;
    plasticity(pt).Coeff_time = fit_plv_coeff;
    plasticity(pt).pCoeff_time = fit_plv_pcoeff;
    plasticity(pt).R2_time_bin = fit_plv_R2_bin;
    plasticity(pt).p_time_bin = fit_plv_p_bin;
    plasticity(pt).Coeff_time_bin = fit_plv_coeff_bin;
    plasticity(pt).pCoeff_time_bin = fit_plv_pcoeff_bin;
    plasticity(pt).R2_stim = fitstim_plv_R2;
    plasticity(pt).p_stim = fitstim_plv_p;
    plasticity(pt).Coeff_stim = fitstim_plv_coeff;
    plasticity(pt).pCoeff_stim= fitstim_plv_pcoeff;
    plasticity(pt).R2_time_noint = fit_dplv_R2;
    plasticity(pt).p_time_noint = fit_dplv_p;
    plasticity(pt).Coeff_time_noint = fit_dplv_coeff;
    plasticity(pt).pCoeff_time_noint = fit_dplv_pcoeff;
    plasticity(pt).R2_time_bin_noint = fit_dplv_R2_bin;
    plasticity(pt).p_time_bin_noint = fit_dplv_p_bin;
    plasticity(pt).Coeff_time_bin_noint = fit_dplv_coeff_bin;
    plasticity(pt).pCoeff_time_bin_noint = fit_dplv_pcoeff_bin;
    plasticity(pt).R2_stim_noint = fitstim_dplv_R2;
    plasticity(pt).p_stim_noint = fitstim_dplv_p;
    plasticity(pt).Coeff_stim_noint = fitstim_dplv_coeff;
    plasticity(pt).pCoeff_stim_noint = fitstim_dplv_pcoeff;
    %
end
save(fullfile(datapath,['plasticity_',num2str(bin_size),suffix,'.mat']),"plasticity")
save(fullfile(datapath,"localization.mat"),'localization');

plasticity(end) = [];
plasticity(21).con_labels = [1,1,0,0,0,0];
% table with time points
pt = []; depth = []; outcome = []; freq = []; conn_type = [];
time = []; bplv = [];
plv = []; dplv = []; zplv = [];
plv_std = [];dplv_std = []; zplv_std = [];
plv_slope = [];dplv_slope = []; zplv_slope = [];
plv_slope_std = [];dplv_slope_std = []; zplv_slope_std = [];

for i = 1:length(plasticity)
    if isempty(plasticity(i).ptID)
        continue
    end
    for t = 1:length(plasticity(i).times)
        for j = 1:6
            for f = 1:4
                pt = [pt;i];
                depth = [depth;plasticity(i).depth];
                if ~isempty(plasticity(i).outcome)
                    outcome = [outcome;plasticity(i).outcome(end)];
                else
                    outcome = [outcome;nan];
                end
                freq = [freq;f];
                conn_type = [conn_type;plasticity(i).con_labels(j)];
                time = [time;t];
                bplv = [bplv;plasticity(i).baseline_plvs(j,f)];
                plv = [plv;plasticity(i).plv(t,j,f)];
                dplv = [dplv;plasticity(i).dplv(t,j,f)];
                zplv = [zplv;plasticity(i).zplv(t,j,f)];
                plv_std = [plv_std;plasticity(i).plv_std(t,j,f)];
                dplv_std = [dplv_std;plasticity(i).dplv_std(t,j,f)];
                zplv_std = [zplv_std;plasticity(i).zplv_std(t,j,f)];
                plv_slope = [plv_slope;plasticity(i).plv_slope(t,j,f)];
                dplv_slope = [dplv_slope;plasticity(i).dplv_slope(t,j,f)];
                zplv_slope = [zplv_slope;plasticity(i).zplv_slope(t,j,f)];
                plv_slope_std = [plv_slope_std;plasticity(i).plv_slope_std(t,j,f)];
                dplv_slope_std = [dplv_slope_std;plasticity(i).dplv_slope_std(t,j,f)];
                zplv_slope_std = [zplv_slope_std;plasticity(i).zplv_slope_std(t,j,f)];
            end
        end
    end
end
tbl = array2table([pt,depth,outcome,freq,conn_type,time,bplv,plv,dplv,zplv,...
    plv_std,dplv_std,zplv_std,plv_slope,dplv_slope,zplv_slope,plv_slope_std,dplv_slope_std,zplv_slope_std],...
    'VariableNames',{'PtID','Depth','Outcome','Freq','Conn','Time','BaselinePLV','PLV','dPLV','zPLV',...
    'PLV_std','dPLV_std','zPLV_std','PLV_Slope','dPLV_Slope','zPLV_Slope',...
    'PLV_Slope_std','dPLV_Slope_std','zPLV_Slope_std'});

% table with no time points
pt = []; depth = []; outcome = []; freq = []; conn_type = [];bplv = [];
plv_std = []; dplv_std = []; zplv_std = [];
plv_slope = [];dplv_slope = []; zplv_slope = [];
plv_slope_std = [];dplv_slope_std = []; zplv_slope_std = [];
R2_time = []; p_time = []; Coeff_time = [];pCoeff_time = [];
R2_time_bin = []; p_time_bin = []; Coeff_time_bin = [];pCoeff_time_bin = [];
R2_stim = []; p_stim = []; Coeff_stim = [];pCoeff_stim = [];
R2_time_noint = []; p_time_noint = []; Coeff_time_noint = [];pCoeff_time_noint = [];
R2_time_bin_noint = []; p_time_bin_noint = []; Coeff_time_bin_noint = [];pCoeff_time_bin_noint = [];
R2_stim_noint = []; p_stim_noint = []; Coeff_stim_noint = [];pCoeff_stim_noint = [];
for i = 1:length(plasticity)
    if isempty(plasticity(i).ptID)
        continue
    end
    for j = 1:6
        for f = 1:4
            pt = [pt;i];
            depth = [depth;plasticity(i).depth];
            if ~isempty(plasticity(i).outcome)
                outcome = [outcome;plasticity(i).outcome(end)];
            else
                outcome = [outcome;nan];
            end
            freq = [freq;f];
            conn_type = [conn_type;plasticity(i).con_labels(j)];
            bplv = [bplv;plasticity(i).baseline_plvs(j,f)];
            plv_std = [plv_std;plasticity(i).plv_std(j,f)];
            dplv_std = [dplv_std;plasticity(i).dplv_std(j,f)];
            zplv_std = [zplv_std;plasticity(i).zplv_std(j,f)];
            plv_slope = [plv_slope;plasticity(i).plv_slope(j,f)];
            dplv_slope = [dplv_slope;plasticity(i).dplv_slope(j,f)];
            zplv_slope = [zplv_slope;plasticity(i).zplv_slope(j,f)];
            plv_slope_std = [plv_slope_std;plasticity(i).plv_slope_std(j,f)];
            dplv_slope_std = [dplv_slope_std;plasticity(i).dplv_slope_std(j,f)];
            zplv_slope_std = [zplv_slope_std;plasticity(i).zplv_slope_std(j,f)];
            R2_time = [R2_time;plasticity(i).R2_time(j,f)];
            p_time = [p_time;plasticity(i).p_time(j,f)];
            Coeff_time = [Coeff_time;plasticity(i).Coeff_time(j,f)];
            pCoeff_time = [pCoeff_time;plasticity(i).pCoeff_time(j,f)];
            R2_time_bin = [R2_time_bin;plasticity(i).R2_time_bin(j,f)];
            p_time_bin = [p_time_bin;plasticity(i).p_time_bin(j,f)];
            Coeff_time_bin = [Coeff_time_bin;plasticity(i).Coeff_time_bin(j,f)];
            pCoeff_time_bin = [pCoeff_time_bin;plasticity(i).pCoeff_time_bin(j,f)];
            R2_stim = [R2_stim;plasticity(i).R2_stim(j,f)];
            p_stim = [p_stim;plasticity(i).p_stim(j,f)];
            Coeff_stim = [Coeff_stim;plasticity(i).Coeff_stim(j,f)];
            pCoeff_stim = [pCoeff_stim;plasticity(i).pCoeff_stim(j,f)];
            R2_time_noint = [R2_time_noint;plasticity(i).R2_time_noint(j,f)];
            p_time_noint = [p_time_noint;plasticity(i).p_time_noint(j,f)];
            Coeff_time_noint = [Coeff_time_noint;plasticity(i).Coeff_time_noint(j,f)];
            pCoeff_time_noint = [pCoeff_time_noint;plasticity(i).pCoeff_time_noint(j,f)];
            R2_time_bin_noint = [R2_time_bin_noint;plasticity(i).R2_time_bin_noint(j,f)];
            p_time_bin_noint = [p_time_bin_noint;plasticity(i).p_time_bin_noint(j,f)];
            Coeff_time_bin_noint = [Coeff_time_bin_noint;plasticity(i).Coeff_time_bin_noint(j,f)];
            pCoeff_time_bin_noint = [pCoeff_time_bin_noint;plasticity(i).pCoeff_time_bin_noint(j,f)];
            R2_stim_noint = [R2_stim_noint;plasticity(i).R2_stim_noint(j,f)];
            p_stim_noint = [p_stim_noint;plasticity(i).p_stim_noint(j,f)];
            Coeff_stim_noint = [Coeff_stim_noint;plasticity(i).Coeff_stim_noint(j,f)];
            pCoeff_stim_noint = [pCoeff_stim_noint;plasticity(i).pCoeff_stim_noint(j,f)];
        end
    end
end
tbl2 = array2table([pt,depth,outcome,freq,conn_type,bplv,...
    plv_std,dplv_std,zplv_std,plv_slope,dplv_slope,zplv_slope,plv_slope_std,dplv_slope_std,zplv_slope_std, ...
    R2_time,p_time,Coeff_time,pCoeff_time,R2_time_bin,p_time_bin,Coeff_time_bin,pCoeff_time_bin, ...
    R2_stim,p_stim,Coeff_stim,pCoeff_stim,R2_time_noint,p_time_noint,Coeff_time_noint,pCoeff_time_noint, ...
    R2_time_bin_noint,p_time_bin_noint,Coeff_time_bin_noint,pCoeff_time_bin_noint, ...
    R2_stim_noint,p_stim_noint,Coeff_stim_noint,pCoeff_stim_noint],...
    'VariableNames',{'PtID','Depth','Outcome','Freq','Conn','BaselinePLV',...
    'PLV_std','dPLV_std','zPLV_std','PLV_Slope','dPLV_Slope','zPLV_Slope',...
    'PLV_Slope_std','dPLV_Slope_std','zPLV_Slope_std',...
    'R2_time','p_time','Coeff_time','pCoeff_time',...
    'R2_time_bin','p_time_bin','Coeff_time_bin','pCoeff_time_bin', ...
    'R2_stim','p_stim','Coeff_stim','pCoeff_stim', ...
    'R2_time_noint','p_time_noint','Coeff_time_noint','pCoeff_time_noint',...
    'R2_time_bin_noint','p_time_bin_noint','Coeff_time_bin_noint','pCoeff_time_bin_noint', ...
    'R2_stim_noint','p_stim_noint','Coeff_stim_noint','pCoeff_stim_noint'});

save(fullfile(datapath,['plasticity_',num2str(bin_size),suffix,'.mat']),'-append','tbl','tbl2');
