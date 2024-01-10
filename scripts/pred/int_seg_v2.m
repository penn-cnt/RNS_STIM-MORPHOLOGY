%% outcome_prediction.m
% This version replicates Ankit's paper Fig 3's analysis
% This script
% 1)
%% Settings
clear; close all;
paths;
patient_info = struct2table(load(which('patients_Penn.mat')).patients_Penn);
ptList = {rns_config.patients.ID};
localization = load(fullfile(datapath,"localization.mat")).localization;

% years = {1,2,3,'end'};
ttdays = 365 * years;
base_days = 90;
rng('default');
try
    data_table = load(fullfile(datapath,'int_seg.mat')).data_table;
catch
    ID = []; d = []; y = []; o = []; o_group = [];
    seg_int = [];
    %% Connectivity Trajectories
    for pt = 1:length(ptList)
        %% Read Patient Data
        ptID = ptList{pt};
        pidx = strcmp(ptID,patient_info.ID);
        %         disp(['Starting analysis for ',ptID])

        outcome = localization(pt).outcome;
        outcome_group = localization(pt).outcome_group;
        depth = localization(pt).depth;

        if ~localization(pt).meets_criteria || isempty(outcome)
            continue
        end

        if length(outcome) > 3
            years = [1,2,3,length(outcome)];
        else
            years = [1:length(outcome)];
        end
        outcome_pred = arrayfun(@(x) max(0,x)/100,outcome(years));
        outcome_group = outcome_group(years);

        resampled_dplv = load(fullfile(datapath,ptID,['cwt_plvs_',ptID,'.mat'])).resampled_dplv;

        dday = patient_info{pidx,"implantDate"};
        time_trace = load(fullfile(datapath,ptID,['UTC_time_trace_',ptID,'.mat'])).time_trace;

        % Calculating the network trajecetories combined across a time window
        baseline_period = dday + days(base_days);
        baseline_mask = time_trace < baseline_period;
        implant_time = time_trace - dday; % get relative day of events after implantation
        implant_time = days(implant_time(~baseline_mask)); % convert to day

        %% Get Predicting Slope
        for i = 1:length(years)
            ydays = years(i) * 365;
            ydays2 = (years(i)-1) * 365;
            year_dplv = resampled_dplv(implant_time <= ydays & implant_time >= ydays2,:,:);
            tmp = diff(year_dplv,1,1);
            tmp_intra = squeeze(mean(mean(tmp(:,1:2,:),1,'omitnan'),2,'omitnan'));
            tmp_inter = squeeze(mean(mean(tmp(:,3:6,:),1,'omitnan'),2,'omitnan'));
            si = tmp_inter - tmp_intra;
            ID = [ID;pt];
            d = [d;depth];
            y = [y;years(i)];
            o = [o;outcome_pred(i)];
            o_group = [o_group;outcome_group(i)];
            seg_int = [seg_int;si'];
        end
        ID = [ID;pt];
        d = [d;depth];
        y = [y;99]; % append again for last year data, use 99 as indicator
        o = [o;outcome_pred(i)];
        o_group = [o_group;outcome_group(i)];
        seg_int = [seg_int;si'];
    end
    data_mat = [ID,d,y,o,o_group,seg_int];
    data_table = array2table(data_mat,"VariableNames",["IDs","Depth","Year","Outcome","OutcomeGroup", ...
        "Theta_Int","Alpha_Int","Beta_Int","Gamma_Int"]);
    %     data_table.OutcomeGroup = categorical(data_table.OutcomeGroup);
    %     data_table.Year = categorical(data_table.Year);
    %     data_table.Depth = categorical(data_table.Depth);
    save(fullfile(datapath,'int_seg.mat'),"data_table")
end
%% Settings
con_strings = ["Inside SOZ", "Between SOZ", "SOZ to Normal","Normal to Normal"];
depth_strings = {'Hippocampal','Neocortical'};
outcome_strings = ["Poor Responder","Good Responder"];
freq_strings = {'Theta','Alpha','Beta','Gamma'};
groupCenters = @(nGroups,nMembers,interGroupSpace) ...
    nGroups/2+.5 : nGroups+interGroupSpace : (nGroups+interGroupSpace)*nMembers-1;
ticks = groupCenters(2, 4, 1);
%% Figure 2 v2
years = {1,2,3,'end'};
for y = 1:length(years)
    if isnumeric(years{y})
        sub_data = table2array(data_table(data_table.Year == years{y},:));
        year = num2str(years{y});
    elseif strcmp(years{y},'end')
        sub_data = table2array(data_table(data_table.Year == 99,:));
        year = 'end';
    end
    outcome = sub_data(:,5);
    f = figure(y);
    f.Position = [100,100,1200,600];
    for d = 1:2
        subplot(2,1,d)
        hold on
        ind = sub_data(:,2) == d;
        data = sub_data(ind,6:9); 
        o = outcome(ind);
        box_data = {};
        sr = data(o == 1,:);
        pr = data(o == 2,:);
        sr = padarray(sr,[max(size(sr,1),size(pr,1))-size(sr,1),0],nan,'post');
        pr = padarray(pr,[max(size(sr,1),size(pr,1))-size(pr,1),0],nan,'post');
        box_data{1} = pr;
        box_data{2} = sr;
        b = boxplotGroup(box_data,'primaryLabels',repmat({''}, 2, 1), ...
            'Colors',[0,1,0;0,0,1],'GroupType','betweenGroups', ...
            'PlotStyle','traditional','BoxStyle','outline', ...
            'Symbol','o','Widths',0.7);
        ymax = get(gca, 'YLim');
        for j = 1:4
            p = ranksum(sr(:,j),pr(:,j));
            text(ticks(j),ymax(2)-0.5,['p=',num2str(p,'%.3f')])
        end
        ylabel(depth_strings{d})
        set(gca,'XTick',ticks,'XTickLabels',freq_strings)
        sgtitle('Int-seg Difference')
    end
    saveas(f,fullfile(figpath,'06_Int_seg',[year,'.png']))
end
close all
