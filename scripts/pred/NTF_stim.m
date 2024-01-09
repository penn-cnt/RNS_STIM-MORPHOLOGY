%% NTF_stim.m
%% Settings
clear; close all;
paths;
patient_info = struct2table(load(which('patients_Penn.mat')).patients_Penn);
ptList = {rns_config.patients.ID};
localization = load(fullfile(datapath,"new_ver/localization.mat")).localization;

fs = 250;
base_days = 90;
rank = 3;
windows = [0, 1, 2, 3, 4, 7, 14, 30, 90, 180, 360, 720];
try
    NTF_stim_pred = load(fullfile(datapath,'NTF_stim_pred.mat')).NTF_stim_pred;
catch
    %% Connectivity Trajectories
    for pt = 1:length(ptList)
        %% Read Patient Data
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
        
        % NTF
        [fac_time, ~, ~] = betaNTF(permute(all_plvs,[1,3,2]),rank);
        for i = 1:3
            fac_time(:,i) = smooth(fac_time(:,i),360);
            fac_time(:,i) = (fac_time(:,i) - min(fac_time(:,i)))/(max(fac_time(:,i))-min(fac_time(:,i)));
        end
    
        % Stim
        [~,~,~, histT] = loadRNSptData(ptID, rns_config);
        stim_data = histT(:,["UTCStartTime","EpisodeStarts"]);
        stim_time = hours(stim_data{:,"UTCStartTime"} - dday);
        [~,stim_idxs] = pdist2(posixtime(stim_data{:,"UTCStartTime"}),ptime_trace,"euclidean",'Smallest',1);
        stims = stim_data{stim_idxs,"EpisodeStarts"};
    
        % pred
        stim_trace_all = {};
        params_all = {};
        R2_all = [];
        R2_null_all = {};
        params_null_all = {};
        for j = 1:length(windows)
            if windows(j) == 0
                stim_trace = stims;
            else
                tmp = stim_data{:,"EpisodeStarts"};
                tmp(isnan(tmp)) = 0;
                cum_stims = movingCumulativeSum(tmp,windows(j))'; %number of stimulations in a 90day period
                stim_trace = cum_stims(stim_idxs);
            end
            % prediction
            mdl = fitglm(fac_time,stim_trace,'Distribution','poisson');
            mdl_null = fitglm(array2table(stim_trace,"VariableNames","stim"),'stim ~ 1','Distribution','poisson');
            R2 = 1 - mdl.Deviance / mdl_null.Deviance;
            b = mdl.Coefficients{2:4,'Estimate'};
            stim_trace_all = [stim_trace_all,stim_trace];
            params_all = [params_all,b];
            R2_all = [R2_all,R2];
            R2_null = [];
            params_null = [];
            for r = 1:100
                inds = randperm(length(stim_trace));
                stim_rnd = stim_trace(inds);
                mdl_rnd = fitglm(fac_time,stim_rnd,'Distribution','poisson');
                mdl_null_rnd = fitglm(array2table(stim_rnd,"VariableNames","stim"),'stim ~ 1','Distribution','poisson');
                R2_rnd = 1 - mdl_rnd.Deviance / mdl_null_rnd.Deviance;
                b_rnd = mdl_rnd.Coefficients{2:4,'Estimate'};
                R2_null = [R2_null,R2_rnd];
                params_null = [params_null,b_rnd];
            end
            R2_null_all = [R2_null_all,R2_null];
            params_null_all = [params_null_all,params_null];
        end
        NTF_stim_pred(pt).ID = ptID;
        NTF_stim_pred(pt).depth = depth;
        NTF_stim_pred(pt).outcome = outcome_group;
        NTF_stim_pred(pt).fac_time = fac_time;
        NTF_stim_pred(pt).stim_data = stims;
        NTF_stim_pred(pt).time = implant_time;
        NTF_stim_pred(pt).stim_trace = stim_trace_all;
        NTF_stim_pred(pt).params = params_all;
        NTF_stim_pred(pt).R2 = R2_all;
        NTF_stim_pred(pt).R2_null = R2_null_all;
        NTF_stim_pred(pt).params_null = params_null_all;
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
    saveas(f,fullfile(datapath,'figs','NTF',[NTF_stim_pred(pt).ID,'.png']))
end
%% Fig.4C
year = 3;
colors = ['g','b','k'];
outcome_strings = ["Poor Responder","Good Responder","Null"];
log_win = log10(windows+1);
NTF_stim_pred = NTF_stim_pred(~cellfun('isempty', {NTF_stim_pred.outcome}));
outcome = cellfun(@(x) x(year), {NTF_stim_pred.outcome});
f = figure('Position',[100,100,1200,600]);
for o = 1:2
    sub_data = NTF_stim_pred(outcome == o);
    R2_combined = vertcat(sub_data.R2);
    shadedErrorBar(log_win,mean(R2_combined,1,'omitnan'), ...
                        std(R2_combined,[],1,'omitnan')/ sqrt(size(R2_combined,1)), ...
                        "lineProps",['-',colors(o)]);
    hold on
end
R2_null_combined = [];
for w = 1:length(windows)
    tmp = cellfun(@(x) x{w},{NTF_stim_pred.R2_null},'UniformOutput',false);
    R2_null_combined = [R2_null_combined,horzcat(tmp{:})'];
end
shadedErrorBar(log_win,mean(R2_null_combined,1,'omitnan'), ...
                    std(R2_null_combined,[],1,'omitnan')/sqrt(size(R2_null_combined,1)), ...
                    "lineProps",['-',colors(3)]);
legend(outcome_strings)
xticks(log_win)
xticklabels(cellstr(num2str(windows')))

% stats
R2_combined = vertcat(NTF_stim_pred.R2);
for w = 1:length(windows)
    p = ranksum(R2_combined(:,w),R2_null_combined(:,w));
    if p < 0.05
        text(log_win(w),0.1,'*','FontSize',20)
    end
end
[fda_pv,~] = FDA_test({NTF_stim_pred.R2},outcome);
text(1,0.8,['p=',num2str(fda_pv,'%.3f')],'FontSize',12)
saveas(f,fullfile(datapath,'figs','NTF',['comp_',num2str(year),'.png']))

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