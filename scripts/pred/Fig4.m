%% NTF_stim.m
%% Settings
clear; close all;
paths;
patient_info = struct2table(load(which('patients_Penn.mat')).patients_Penn);
ptList = {rns_config.patients.ID};
localization = load(fullfile(datapath,"localization.mat")).localization;

base_days = 90;
rank = 3;
windows = [0, 1, 2, 3, 4, 7, 14, 30, 90, 180, 360, 720];
%% load data
try
    NTF_stim_pred = load(fullfile(datapath,'NTF_stim_pred.mat')).NTF_stim_pred;
catch
    for pt = 1:length(ptList)
        % Read Patient Data
        ptID = ptList{pt};
        pidx = strcmp(ptID,patient_info.ID);
        disp(['Starting analysis for ',ptID])
    
        if ~localization(pt).meets_criteria
            continue
        end
    
        outcome = localization(pt).outcome;
        outcome_group = localization(pt).outcome_group;
        depth = localization(pt).depth;
    
        % load plv
        all_plvs = load(fullfile(datapath,ptID,['cwt_plvs_',ptID,'.mat'])).all_plvs;
        freqs = load(fullfile(datapath,ptID,['cwt_plvs_',ptID,'.mat'])).f;
        
        % load event data
        dday = patient_info{pidx,"implantDate"};
        time_trace = load(fullfile(datapath,ptID,['UTC_time_trace_',ptID,'.mat'])).time_trace;
        ptime_trace = load(fullfile(datapath,ptID,['posix_UTC_time_trace_',ptID,'.mat'])).ptime_trace;
        implant_time = time_trace - dday; % get relative day of events after implantation
        implant_time = days(implant_time); % convert to day
        
        % Calculating the network trajecetories combined across a time window
        baseline_period = dday + days(base_days);
        baseline_mask = time_trace < baseline_period;

        % NTF
        [fac_time, ~, ~] = betaNTF(permute(all_plvs,[1,3,2]),rank);
        for i = 1:3
            fac_time(:,i) = smooth(fac_time(:,i),360);
            fac_time(:,i) = (fac_time(:,i) - min(fac_time(:,i)))/(max(fac_time(:,i))-min(fac_time(:,i)));
        end
        base_fac_time = squeeze(mean(fac_time(baseline_mask,:),1,'omitnan'));
    
        % Stim
        [~,~,~, histT] = loadRNSptData(ptID, rns_config);
        stim_data = histT(:,["UTCStartTime","EpisodeStarts"]);
        stim_time = hours(stim_data{:,"UTCStartTime"} - dday);
        [~,stim_idxs] = pdist2(posixtime(stim_data{:,"UTCStartTime"}),ptime_trace,"euclidean",'Smallest',1);
        stims = stim_data{stim_idxs,"EpisodeStarts"};
    
        % pred
        stim_trace_all = {};
        pValue_all = {};
        R2_all = [];
        pR2_all = [];
        R2_null_all = {};
        for j = 1:length(windows)
            if windows(j) == 0
                stim_trace = stims;
                stim_trace(isnan(stim_trace)) = 0;
            else
                tmp = stim_data{:,"EpisodeStarts"};
                tmp(isnan(tmp)) = 0;
                cum_stims = movingCumulativeSum(tmp,windows(j))'; %number of stimulations in a 90day period
                stim_trace = cum_stims(stim_idxs);
            end
            % prediction
            mdl = fitglm(fac_time,stim_trace,'Distribution','poisson');
            R2 = mdl.Rsquared.Deviance;
            p = mdl.Coefficients{2:4,'pValue'};
            stim_trace_all = [stim_trace_all,stim_trace];
            pValue_all = [pValue_all,p];
            R2_all = [R2_all,R2];
            R2_null = [];
            for r = 1:100
                inds = randperm(length(stim_trace));
                stim_rnd = stim_trace(inds);
                mdl_rnd = fitglm(fac_time,stim_rnd,'Distribution','poisson');
                R2_rnd = mdl_rnd.Rsquared.Deviance;
                p_rnd = mdl_rnd.Coefficients{2:4,'pValue'};
                R2_null = [R2_null,R2_rnd];
            end
            R2_null_all = [R2_null_all,R2_null];
            pR2_all = [pR2_all,mean(R2_null > R2)];
        end
        NTF_stim_pred(pt).ID = ptID;
        NTF_stim_pred(pt).depth = depth;
        NTF_stim_pred(pt).outcome = outcome_group;
        NTF_stim_pred(pt).fac_time = fac_time;
        NTF_stim_pred(pt).base_fac_time = base_fac_time;
        NTF_stim_pred(pt).stim_data = stims;
        NTF_stim_pred(pt).time = implant_time;
        NTF_stim_pred(pt).stim_trace = stim_trace_all;
        NTF_stim_pred(pt).pValue = pValue_all;
        NTF_stim_pred(pt).R2 = R2_all;
        NTF_stim_pred(pt).pR2 = pR2_all;
        NTF_stim_pred(pt).R2_null = R2_null_all;
    end
    NTF_stim_pred = NTF_stim_pred(~cellfun('isempty', {NTF_stim_pred.ID}));
    save(fullfile(datapath,['NTF_stim_pred.mat']),"NTF_stim_pred")
end

%% NTM plots
%% Fig.4B
for pt = 1:length(NTF_stim_pred)
    f = figure('Position',[100,100,1200,600]);
    subplot(2,1,1)
    plot(NTF_stim_pred(pt).time,NTF_stim_pred(pt).stim_data,'k')
    hold on
    plot(NTF_stim_pred(pt).time,NTF_stim_pred(pt).stim_trace{9},'r')
    ylabel('Hourly/Cumul. # therapies')
    subplot(2,1,2)
    plot(NTF_stim_pred(pt).time,NTF_stim_pred(pt).fac_time)
    xlabel('Days since implant')
    ylabel('Time expression')
    saveas(f,fullfile(figpath,'02_NTF_stim',[NTF_stim_pred(pt).ID,'.png']))
end
close all
%% Fig.4C
colors = ['g','b','k'];
outcome_strings = ["Poor Responder","Good Responder","Null"];
log_win = log10(windows+1);
NTF_stim_pred = NTF_stim_pred(~cellfun('isempty', {NTF_stim_pred.outcome}));
years = {1,2,3,'end'};
for y = 1:length(years)
    [all_data,outcome] = getYearOutcome(NTF_stim_pred,years{y},'outcome',{'R2','R2_null'});
    R2 = all_data{1};
    R2_null = all_data{2};
    R2_combined = vertcat(R2{:});
    R2_null_combined = cellfun(@(x) vertcat(x{:})',R2_null,'UniformOutput',false);
    R2_null_combined = vertcat(R2_null_combined{:});
    f = figure('Position',[100,100,1200,600]);
    for o = 1:2
        sub_data = R2_combined(outcome == o,:);
        shadedErrorBar(log_win,mean(sub_data,1,'omitnan'), ...
                            std(sub_data,[],1,'omitnan')/ sqrt(size(sub_data,1)), ...
                            "lineProps",['-',colors(o)]);
        hold on
    end
    shadedErrorBar(log_win,mean(R2_null_combined,1,'omitnan'), ...
                        std(R2_null_combined,[],1,'omitnan')/sqrt(size(R2_null_combined,1)), ...
                        "lineProps",['-',colors(3)]);
    legend(outcome_strings)
    xticks(log_win)
    xticklabels(cellstr(num2str(windows')))
    
    % stats
    for w = 1:length(windows)
        p = ranksum(R2_combined(:,w),R2_null_combined(:,w));
        if p < 0.05
            text(log_win(w),0.1,'*','FontSize',20)
        end
    end
    [fda_pv,~] = FDA_test(R2,outcome);
    text(1,0.8,['p=',num2str(fda_pv,'%.3f')],'FontSize',12)
    saveas(f,fullfile(figpath,'02_NTF_stim',['comp_',num2str(years{y}),'.png']))
end
close all
%%
function movingSum = movingCumulativeSum(inputArray, windowSize)
% Check if the window size is valid
if windowSize <= 0 || windowSize > length(inputArray)
    error('Window size must be a positive integer smaller than or equal to the length of the input array.');
end

conv_factor = exp(-[0:1/24:1e4]/windowSize);
conv_factor = conv_factor / sum(conv_factor);
% Initialize the output array to store the moving cumulative sum
movingSum = conv(inputArray,conv_factor,'full');
movingSum = movingSum(1:length(inputArray));
end