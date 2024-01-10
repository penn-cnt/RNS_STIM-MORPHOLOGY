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
calc_R2 = @(o_true,o_pred) 1 - sum((o_true-o_pred).^2) / sum((o_true-mean(o_true,'omitnan')).^2);
rng('default');
try
    data_table = load(fullfile(datapath,'plv_slope_pred.mat')).data_table;
catch
    ID = []; d = []; y = []; o = []; o_group = [];
    ps = []; dps = [];
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

        resampled_plv = load(fullfile(datapath,ptID,['cwt_plvs_',ptID,'.mat'])).resampled_plv;
        resampled_dplv = load(fullfile(datapath,ptID,['cwt_plvs_',ptID,'.mat'])).resampled_dplv;

        dday = patient_info{pidx,"implantDate"};
        time_trace = load(fullfile(datapath,ptID,['UTC_time_trace_',ptID,'.mat'])).time_trace;
    
        % Calculating the network trajecetories combined across a time window
        baseline_period = dday + days(base_days);
        baseline_mask = time_trace < baseline_period;
        implant_time = time_trace - dday; % get relative day of events after implantation
        implant_time = days(implant_time(~baseline_mask)); % convert to day

        %% Get Predicting Slope
        plv = squeeze(mean(resampled_plv(:,1:2,:),2,'omitnan'));
        dplv = squeeze(mean(resampled_dplv(:,1:2,:),2,'omitnan'));
        plv_slopes = nan * zeros(length(years),4);
        dplv_slopes = nan * zeros(length(years),4);

        for i = 1:length(years)
            ydays = years(i) * 365;
            year_plv = plv(implant_time <= ydays,:);
            year_dplv = dplv(implant_time <= ydays,:);
            year_time = implant_time(implant_time <= ydays);
            for j = 1:4
                p = polyfit(year_time,zscore(year_plv(:,j)),1);
                plv_slopes(i,j) = p(1);
                p = polyfit(year_time,zscore(year_dplv(:,j)),1);
                dplv_slopes(i,j) = p(1);
            end
            ID = [ID;pt];
            d = [d;depth];
            y = [y;years(i)];
            o = [o;outcome_pred(i)];
            o_group = [o_group;outcome_group(i)];
            ps = [ps;plv_slopes(i,:)];
            dps = [dps;dplv_slopes(i,:)];
        end
        ID = [ID;pt];
        d = [d;depth];
        y = [y;99]; % append again for last year data, use 99 as indicator
        o = [o;outcome_pred(i)];
        o_group = [o_group;outcome_group(i)];
        ps = [ps;plv_slopes(i,:)];
        dps = [dps;dplv_slopes(i,:)];
    end
    data_mat = [ID,d,y,o,o_group,ps,dps];
    data_table = array2table(data_mat,"VariableNames",["IDs","Depth","Year","Outcome","OutcomeGroup", ...
        "Theta_Slope","Alpha_Slope","Beta_Slope","Gamma_Slope", ...
        "dTheta_Slope","dAlpha_Slope","dBeta_Slope","dGamma_Slope"]);
%     data_table.OutcomeGroup = categorical(data_table.OutcomeGroup);
%     data_table.Year = categorical(data_table.Year);
%     data_table.Depth = categorical(data_table.Depth);
    save(fullfile(datapath,'plv_slope_pred.mat'),"data_table")
end
%% prediction
R2_all = [];
params_all = {};
preds_all = {};
trues_all = {};
R2_all_null = {};
pR2_all = [];
for d = 1:2
% one model for each depth type
    for y = [1,2,3,99]
    % first 3 years first
        sub_data = table2array(data_table(data_table.Depth == d & data_table.Year == y,:));
        o_true = sub_data(:,4);
        mdl = fitglm(sub_data(:,6:9),o_true,'Distribution','binomial','Intercept',false);
        b = mdl.Coefficients{:,"Estimate"};
        o_pred = predict(mdl,sub_data(:,6:9));
        R2 = calc_R2(o_true,o_pred);

        % record
        R2_all = [R2_all,R2];
        preds_all = [preds_all,o_pred];
        trues_all = [trues_all,o_true];
        
        % permutated
        params_null = [];
        preds_null = [];
        R2_null = [];
        for r = 1:100
            inds = randperm(length(o_true));
            o_true_rnd = o_true(inds);
            mdl_rnd = fitglm(sub_data(:,6:9),o_true_rnd,'Distribution','binomial','Intercept',false);
            b_rnd = mdl_rnd.Coefficients{:,"Estimate"};
            params_null = [params_null,b_rnd];
            o_pred_rnd = predict(mdl_rnd,sub_data(:,6:9));
            preds_null = [preds_null,o_pred_rnd];
            R2_rnd = calc_R2(o_true,o_pred_rnd);
            R2_null = [R2_null,R2_rnd];
        end
        R2_all_null = [R2_all_null,R2_null];
        pR2_all = [pR2_all,mean(R2_null > R2)];

        % stats
        norm_params = [];
        for f = 1:4
            pv = mean(params_null(f,:) > b(f));
            pv = min(pv,1-pv);
            val = (b(f) - mean(params_null(f,:),'omitnan')) / std(params_null(f,:),[],'omitnan');
            norm_params = [norm_params; [val, pv]];
        end
        params_all = [params_all,norm_params];
    end
end
R2_all = reshape(R2_all,[4,2])';
params_all = reshape(params_all,[4,2])';
preds_all = reshape(preds_all,[4,2])';
trues_all = reshape(trues_all,[4,2])';
R2_all_null = reshape(R2_all_null,[4,2])';
pR2_all = reshape(pR2_all,[4,2])';
save(fullfile(datapath,'plv_slope_pred.mat'),'-append','R2_all','preds_all','trues_all','params_all','R2_all_null','pR2_all');
%% plotting
% Fig.3A/D
depth_strings = {'Hippocampal','Neocortical'};
f = figure('Position',[100,100,1200,600]);
for d = 1:2
    for y = 1:4
        subplot(2,4,(d-1)*4+y)
        scatter(100*trues_all{d,y},100*preds_all{d,y},[],'k','filled');
        hold on
        line([0,100],[0,100],'Color','k','LineStyle','--')
        xlim([-5,105])
        ylim([-5,105])
        text(100,90,{['R^{2} = ',num2str(R2_all(d,y),'%.2f')], ...
            ['R^{2}_{null} = ',num2str(mean(R2_all_null{d,y}),'%.2f')], ...
            ['p = ',num2str(pR2_all(d,y),'%.3f')]})
        if y == 1
            ylabel(depth_strings{d},'FontWeight','bold','FontSize',12)
        end
        if d == 1
            if y == 4
                title(['Last Year'],'FontWeight','bold','FontSize',12)
            else
                title(['Year',num2str(y)],'FontWeight','bold','FontSize',12)
            end
        end
    end
end
h = axes(f,'visible','off'); 
h.Title.Visible='on';
h.XLabel.Visible='on';
h.YLabel.Visible='on';
ylabel(h,'Predicted % Seizure Change','Position',[-0.05,0.500000476837158,0], ...
    'FontWeight','bold','FontSize',14);
xlabel(h,'Observed % Seizure Change','FontWeight','bold','FontSize',14);
saveas(f,fullfile(figpath,'01_PLV_outcome','Fig.3A.png'))
%% Fig.3B/E
freqs = {'Theta','Alpha','Beta','Gamma'};
cols = {'k','r'};
g = figure('Position',[100,100,1200,600]);
for d = 1:2
    for y = 1:4
        subplot(2,4,(d-1)*4+y)
        to_plot = params_all{d,y};
        hold on
        for f = 1:4
            col = cols{(to_plot(f,2) < 0.05)+1};
            errorbar(f,-to_plot(f,1),1,[col,'o'],'MarkerFaceColor',col);
        end
        line([0,5],[0,0],'Color','k','LineStyle','--')
        xticks([1:4])
        xticklabels(freqs)
        ylim([-4,4])
        xlim([0,5])
        if y == 1
            ylabel(depth_strings{d},'FontWeight','bold','FontSize',12)
        end
        if d == 1
            if y == 4
                title(['Last Year'],'FontWeight','bold','FontSize',12)
            else
                title(['Year',num2str(y)],'FontWeight','bold','FontSize',12)
            end
        end
    end
end
saveas(g,fullfile(figpath,'01_PLV_outcome','Fig.3B.png'))
close all
