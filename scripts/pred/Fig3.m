%% Fig3.m
% This version replicates Ankit's paper Fig 3's analysis
% Specifically using intra-lead PLV slope till certain time point to predict outcome
% at that point
%% Settings
function Fig3(regress_dist,outcome_option)
close all;
paths;
patient_info = struct2table(load(which('patients_Penn.mat')).patients_Penn);
ptList = {rns_config.patients.ID};
localization = load(fullfile(datapath,"localization.mat")).localization;
if regress_dist
    suffix = '_regdist';
    figpath = fullfile(figpath,'regdist');
else
    suffix = '';
end
plasticity = load(fullfile(datapath,['plasticity_90',suffix,'.mat'])).plasticity;
for i = 1:length(localization)
    if ~isempty(localization(i).outcome_group)
        localization(i).outcome_group = localization(i).outcome_group(outcome_option,:);
    end
    if ~isempty(plasticity(i).outcome_group)
        plasticity(i).outcome_group = plasticity(i).outcome_group(outcome_option,:);
    end
end
base_days = 90;
rng('default');
%% load data
% try
%     data_table = load(fullfile(datapath,'plv_slope_pred.mat')).data_table;
% catch
    ID = []; d = []; y = []; o = []; o_group = []; 
    ps_intra = []; dps_intra = []; ps_inter = []; dps_inter = [];
    for pt = 1:length(localization)
        % Read Patient Data
        ptID = localization(pt).ptID;
        pidx = strcmp(ptID,patient_info.ID);

        outcome = plasticity(pt).outcome;
        outcome_group = plasticity(pt).outcome_group;
        depth = plasticity(pt).depth;
        lead_label = localization(pt).lead_labels;

        if ~localization(pt).meets_criteria | isempty(outcome) | outcome_group == 3
            continue
        end

        if length(outcome) > 3
            years = [1,2,3,length(outcome)];
        else
            years = [1:length(outcome)];
        end
        outcome_pred = arrayfun(@(x) max(0,x)/100,outcome(years));
        outcome_group = outcome_group(years);

        resampled_plv = load(fullfile(datapath,ptID,['working_data_',num2str(pt),suffix,'.mat'])).resampled_plv_all;
%         resampled_dplv = load(fullfile(datapath,ptID,['working_data_',num2str(pt),'.mat'])).resampled_dplv;
        implant_time = load(fullfile(datapath,ptID,['working_data_',num2str(pt),suffix,'.mat'])).implant_time_all/365;

        % Get Predicting Slope
%         plv = [];
%         dplv = [];
%         plv(1,:,:) = squeeze(mean(resampled_plv(:,lead_label == 1,:),2,'omitnan'));
%         dplv(1,:,:) = squeeze(mean(resampled_dplv(:,lead_label == 1,:),2,'omitnan'));
%         plv(2,:,:) = squeeze(mean(resampled_plv(:,lead_label == 2,:),2,'omitnan'));
%         dplv(2,:,:) = squeeze(mean(resampled_dplv(:,lead_label == 2,:),2,'omitnan'));
        plv_slopes = nan * zeros(length(years),size(resampled_plv,2),4);
% %         dplv_slopes = nan * zeros(2,length(years),4);

        for i = 1:length(years)
            ydays = years(i);
            year_plv = resampled_plv(implant_time <= ydays,:,:);
%             year_dplv = dplv(:,implant_time <= ydays,:);
            year_time = implant_time(implant_time <= ydays);
            for l = 1:size(year_plv,2)
                for j = 1:4 % each freq band
                    p = polyfit(year_time,zscore(year_plv(:,l,j)),1);
                    plv_slopes(i,l,j) = p(1);
%                     p = polyfit(year_time,zscore(year_dplv(:,d,j)),1);
%                     dplv_slopes(l,i,j) = p(1);
                end
            end
            plvslope_intra = reshape(mean(plv_slopes(i,lead_label==1,:),2,'omitnan'),[1,4]);
            plvslope_inter = reshape(mean(plv_slopes(i,lead_label==2,:),2,'omitnan'),[1,4]);
            ID = [ID;pt];
            d = [d;depth];
            y = [y;years(i)];
            o = [o;outcome_pred(i)];
            o_group = [o_group;outcome_group(i)];
            ps_intra = [ps_intra;plvslope_intra];
%             dps_intra = [dps_intra;dplv_slopes(1,i,:)];
            ps_inter = [ps_inter;plvslope_inter];
%             dps_inter = [dps_inter;dplv_slopes(2,i,:)];
        end
        ID = [ID;pt];
        d = [d;depth];
        y = [y;99]; % append again for last year data, use 99 as indicator
        o = [o;outcome_pred(i)];
        o_group = [o_group;outcome_group(i)];
        ps_intra = [ps_intra;plvslope_intra];
%         dps_intra = [dps_intra;dplv_slopes(1,i,:)];
        ps_inter = [ps_inter;plvslope_inter];
%         dps_inter = [dps_inter;dplv_slopes(2,i,:)];
    end
    data_mat = [ID,d,y,o,o_group,squeeze(ps_intra),squeeze(ps_inter)];
    data_table = array2table(data_mat,"VariableNames",["IDs","Depth","Year","Outcome","OutcomeGroup", ...
        "Theta_Intra_PLV","Alpha_Intra_PLV","Beta_Intra_PLV","Gamma_Intra_PLV", ...
        "Theta_Inter_PLV","Alpha_Inter_PLV","Beta_Inter_PLV","Gamma_Inter_PLV"]);
%         "Theta_Inter_dPLV","Alpha_Inter_dPLV","Beta_Inter_dPLV","Gamma_Inter_dPLV"]);
    save(fullfile(datapath,['plv_slope_pred',suffix,'.mat']),"data_table")
% end

%% prediction
pred_range = reshape([6:13],4,2)';
for n = 1:size(pred_range,1)
    R2_all = [];
    params_all = {};
    preds_all = {};
    trues_all = {};
    R2_all_null = {};
    pR2_all = [];
    label = data_table.Properties.VariableNames{pred_range(n,1)};
    label(1:strfind(label,'_')) = [];
    for d = 1:2
        % one model for each depth type
        for y = [1,2,3,99]
            sub_data = table2array(data_table(data_table.Depth == d & data_table.Year == y,:));
            o_true = sub_data(:,4);
            mdl = fitglm(sub_data(:,pred_range(n,:)),o_true,'Distribution','binomial','Intercept',false);
            b = mdl.Coefficients{:,"Estimate"};
            o_pred = predict(mdl,sub_data(:,pred_range(n,:)));
            R2 = mdl.Rsquared.Ordinary;

            % record
            R2_all = [R2_all,R2];
            preds_all = [preds_all,o_pred];
            trues_all = [trues_all,o_true];

            % permutated
            params_null = [];
            preds_null = [];
            R2_null = [];
            for r = 1:500
                inds = randperm(length(o_true));
                o_true_rnd = o_true(inds);
                mdl_rnd = fitglm(sub_data(:,pred_range(n,:)),o_true_rnd,'Distribution','binomial','Intercept',false);
                b_rnd = mdl_rnd.Coefficients{:,"Estimate"};
                params_null = [params_null,b_rnd];
                o_pred_rnd = predict(mdl_rnd,sub_data(:,pred_range(n,:)));
                preds_null = [preds_null,o_pred_rnd];
                R2_rnd = mdl_rnd.Rsquared.Ordinary;
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
    save(fullfile(datapath,['Outcome_Pred_',label,suffix,'.mat']),'R2_all','preds_all','trues_all','params_all','R2_all_null','pR2_all');
end
%% plotting
% Fig.3A/D
depth_strings = {'Hippocampal','Neocortical'};
pred_range = reshape([6:13],4,2)';
for n = 1:size(pred_range,1)
    label = data_table.Properties.VariableNames{pred_range(n,1)};
    label(1:strfind(label,'_')) = [];
    if ~exist(fullfile(figpath,'01_PLV_outcome',num2str(outcome_option),label),'dir')
        mkdir(fullfile(figpath,'01_PLV_outcome',num2str(outcome_option),label))
    end
    load(fullfile(datapath,['Outcome_Pred_',label,suffix,'.mat']))
    f = figure('Position',[100,100,1200,600],'Visible','Off');
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
    saveas(f,fullfile(figpath,'01_PLV_outcome',num2str(outcome_option),label,'Fig.3A.png'))
end
close all
%% Fig.3B/E
freqs = {'Theta','Alpha','Beta','Gamma'};
cols = {'k','r'};
pred_range = reshape([8:15],4,2)';
for n = 1:size(pred_range,1)
    label = data_table.Properties.VariableNames{pred_range(n,1)};
    label(1:strfind(label,'_')) = [];
    load(fullfile(datapath,['Outcome_Pred_',label,suffix,'.mat']))
    g = figure('Position',[100,100,1200,600],'Visible','Off');
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
    saveas(g,fullfile(figpath,'01_PLV_outcome',num2str(outcome_option),label,'Fig.3B.png'))
end
close all
end
