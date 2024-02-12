%% Baseline_outcome.m
% Tests whether baseline connectivity predict outcome at one, two, and three years and last outcome
%% Settings
function baseline_outcome(regress_dist,outcome_option)
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
%% Baseline connectivity predict outcome
% try
%     data_table = load(fullfile(datapath,'baseline_outcome_pred.mat')).data_table;
% catch
    ID = []; d = []; y = []; o = []; 
    bp_inter = []; bp_intra = [];
    for pt = 1:length(localization)
        % Read Patient Data
        ptID = localization(pt).ptID;

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

        % baseline
        baseline_plv = plasticity(pt).baseline_plvs;
        baseline_plv_intra = squeeze(mean(baseline_plv(lead_label==1,:),1,'omitnan'));
        baseline_plv_inter = squeeze(mean(baseline_plv(lead_label==2,:),1,'omitnan'));

        for yy = 1:length(years)
            ID = [ID;pt];
            d = [d;depth];
            y = [y;years(yy)];
            o = [o;outcome_pred(yy)];
            bp_inter = [bp_inter;baseline_plv_inter];
            bp_intra = [bp_intra;baseline_plv_intra];
        end
        ID = [ID;pt];
        d = [d;depth];
        y = [y;99];
        o = [o;outcome_pred(yy)];
        bp_inter = [bp_inter;baseline_plv_inter];
        bp_intra = [bp_intra;baseline_plv_intra];
    end
    data_mat = [ID,d,y,o,bp_intra,bp_inter];
    data_table = array2table(data_mat,"VariableNames", ...
        ["IDs","Depth","Year","Outcome",...
        "bplv_intra_Theta","bplv_intra_Alpha","bplv_intra_Beta","bplv_intra_Gamma", ...
        "bplv_inter_Theta","bplv_inter_Alpha","bplv_inter_Beta","bplv_inter_Gamma"]);
    save(fullfile(datapath,['baseline_outcome_pred',suffix,'.mat']),"data_table")
% end
%% prediction
pred_range = reshape([5:12],4,2)';
label = {'intra','inter'};
for n = 1:size(pred_range,1)
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
            R2_null = [];
            for r = 1:500
                inds = randperm(length(o_true));
                o_true_rnd = o_true(inds);
                mdl_rnd = fitglm(sub_data(:,pred_range(n,:)),o_true_rnd,'Distribution','binomial','Intercept',false);
                b_rnd = mdl_rnd.Coefficients{:,"Estimate"};
                params_null = [params_null,b_rnd];
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
    save(fullfile(datapath,['baseline_outcome_pred_',label{n},suffix,'.mat']),'R2_all','preds_all','trues_all','params_all','R2_all_null','pR2_all');
end
%% plotting
% Fig.3A/D
depth_strings = {'Hippocampal','Neocortical'};
pred_range = reshape([5:12],4,2)';
label = {'intra','inter'};
for n = 1:size(pred_range,1)
    if ~exist(fullfile(figpath,'03_Baseline_outcome',num2str(outcome_option),label{n}),'dir')
        mkdir(fullfile(figpath,'03_Baseline_outcome',num2str(outcome_option),label{n}))
    end
    load(fullfile(datapath,['baseline_outcome_pred_',label{n},suffix,'.mat']))
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
    saveas(f,fullfile(figpath,'03_Baseline_outcome',num2str(outcome_option),label{n},'Fig.3A.png'))
end
%% Fig.3B/E
freqs = {'Theta','Alpha','Beta','Gamma'};
cols = {'k','r'};
pred_range = reshape([5:12],4,2)';
label = {'intra','inter'};
for n = 1:size(pred_range,1)
    load(fullfile(datapath,['baseline_outcome_pred_',label{n},suffix,'.mat']))
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
    saveas(g,fullfile(figpath,'03_Baseline_outcome',num2str(outcome_option),label{n},'Fig.3B.png'))
end
close all