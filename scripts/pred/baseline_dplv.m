%% exploration.m
% Figure 2
% Baseline connectivity does (not?) predict outcome at one, two, and three years and last outcome
% Showing relationship in different frequencies between baseline connectivity
% and %change in connectivity at different time points.
% Figure 3
% Repeat figure 2 with structural connectivity rather than functional
% Figure 4
% See how stimulation sensitivity depends on baseline connectivity
% Figure 5
% How does structural or change in functional connectivity change the acute synchronization effects of stimulation
%% Settings
clear; close all;
paths;
patient_info = struct2table(load(which('patients_Penn.mat')).patients_Penn);
ptList = {rns_config.patients.ID};
localization = load(fullfile(datapath,"localization.mat")).localization;
plasticity = load(fullfile(datapath,"plasticity_90.mat")).plasticity;
base_days = 90;
bin_size = 90; % subject to change
ttdays = 1050;
time = 90:bin_size:ttdays-1;
n_bin = length(time);
rng('default');
%% Baseline connectivity predict outcome
% try
%     data_table = load(fullfile(datapath,'baseline_plv_pred.mat')).data_table;
% catch
ID = []; d = []; o = []; o_group = [];
bp_inter = []; bp_intra = []; dp_inter = []; dp_intra = []; t = []; s = [];
%% Connectivity Trajectories
for p = 1:length(plasticity)
    %% Read Patient Data
    ptID = plasticity(p).ptID;
    pt = find(strcmp(ptID,{localization.ptID}));

    outcome = localization(pt).outcome;
    outcome_group = localization(pt).outcome_group;
    depth = localization(pt).depth;

    if ~localization(pt).meets_criteria || isempty(outcome)
        continue
    end

    if length(outcome) > 3
        years = [1,2,3,length(outcome)];
        outcome_pred = arrayfun(@(x) max(0,x)/100,outcome(years));
        outcome_group = outcome_group(years);
    else
        years = [1:length(outcome)];
        outcome_pred = arrayfun(@(x) max(0,x)/100,outcome(years));
        outcome_group = outcome_group(years);
        outcome_pred = padarray(outcome_pred,[0,3-length(years)],3,'post');
        outcome_group = padarray(outcome_group,[0,3-length(years)],3,'post');
        outcome_pred = [outcome_pred,outcome_pred(length(outcome))];
        outcome_group = [outcome_group,outcome_group(length(outcome))];
    end

    % baseline
    baseline_plv = plasticity(p).baseline_plvs;
    baseline_plv_intra = squeeze(mean(baseline_plv(1:2,:),1,'omitnan'));
    baseline_plv_inter = squeeze(mean(baseline_plv(3:6,:),1,'omitnan'));

    % dplv at different time points
    dplv_intra = cellfun(@(x) mean(x(1:min(n_bin,size(x,1)),:),2,'omitnan'),plasticity(p).reorg_dplv(1,:),'UniformOutput',false);
    dplv_intra = horzcat(dplv_intra{:});
    dplv_intra = padarray(dplv_intra,n_bin-size(dplv_intra,1),nan,'post');
    dplv_inter = cellfun(@(x) mean(x(1:min(n_bin,size(x,1)),:),2,'omitnan'),plasticity(p).reorg_dplv(2,:),'UniformOutput',false);
    dplv_inter = horzcat(dplv_inter{:});
    dplv_inter = padarray(dplv_inter,n_bin-size(dplv_inter,1),nan,'post');

    for tt = 1:n_bin
        ID = [ID;pt];
        d = [d;depth];
        o = [o;outcome_pred];
        o_group = [o_group;outcome_group];
        bp_inter = [bp_inter;baseline_plv_inter];
        bp_intra = [bp_intra;baseline_plv_intra];
        dp_inter = [dp_inter;dplv_inter(tt,:)];
        dp_intra = [dp_intra;dplv_intra(tt,:)];
        t = [t;tt];
    end
end
data_mat = [ID,d,o,o_group,bp_intra,bp_inter,dp_intra,dp_inter,t];
data_table = array2table(data_mat,"VariableNames", ...
    ["IDs","Depth","Outcome_Year1","Outcome_Year2","Outcome_Year3", "Outcome_Last", ...
    "OutcomeGroup_Year1", "OutcomeGroup_Year2", "OutcomeGroup_Year3", "OutcomeGroup_Last", ...
    "bplv_intra_Theta","bplv_intra_Alpha","bplv_intra_Beta","bplv_intra_Gamma", ...
    "bplv_inter_Theta","bplv_inter_Alpha","bplv_inter_Beta","bplv_inter_Gamma", ...
    "dplv_intra_Theta","dplv_intra_Alpha","dplv_intra_Beta","dplv_intra_Gamma", ...
    "dplv_inter_Theta","dplv_inter_Alpha","dplv_inter_Beta","dplv_inter_Gamma", ...
    "time"]);
%     data_table.OutcomeGroup = categorical(data_table.OutcomeGroup);
%     data_table.Year = categorical(data_table.Year);
%     data_table.Depth = categorical(data_table.Depth);
save(fullfile(datapath,'baseline_plv_pred.mat'),"data_table")
% end
%% prediction
R2_all = nan*zeros(2,n_bin,4,2);
pR2_all = nan*zeros(2,n_bin,4,2);
dplv_all = cell(2,n_bin,4,2);
bplv_all = cell(2,n_bin,4,2);
outcome_all = cell(2,n_bin,3);
freq = {'Theta','Alpha','Beta','Gamma'};
conn = {'intra','inter'};
for d = 1:2
    % one model for each depth type
    for y = 1:n_bin
        % first 3 years first
        sub_data = data_table(data_table.Depth == d & data_table.time == y,:);
        o1 = table2array(sub_data(:,'OutcomeGroup_Year1'));
        o2 = table2array(sub_data(:,'OutcomeGroup_Year2'));
        o3 = table2array(sub_data(:,'OutcomeGroup_Year3'));
        o4 = table2array(sub_data(:,'OutcomeGroup_Last'));
        for f = 1:4
            for c = 1:2
                o_true = table2array(sub_data(:,['dplv_',conn{c},'_',freq{f}]));
                pred = table2array(sub_data(:,['bplv_',conn{c},'_',freq{f}]));
                [r,p] = corr(pred,o_true,'Type','Spearman');
                % record
                R2_all(d,y,f,c) = r^2;
                pR2_all(d,y,f,c) = p;
                dplv_all{d,y,f,c} = o_true;
                bplv_all{d,y,f,c} = pred;
            end
        end
        outcome_all{d,y,1} = o1;
        outcome_all{d,y,2} = o2;
        outcome_all{d,y,3} = o3;
        outcome_all{d,y,4} = o4;
    end
end
save(fullfile(datapath,'baseline_plv_pred.mat'),'-append','R2_all','pR2_all','dplv_all','bplv_all','outcome_all');
%% plotting
% Fig.3A/D
depth_strings = {'Hippocampal','Neocortical'};
freq_strings = {'Theta','Alpha','Beta','Gamma'};
conn_strings = {'intra','inter'};
pos = [0.9,0.1];
col = [0,1,0;0,0,1;0,0,0];
for freq = 1:4
    for conn = 1:2
        for y = 1:n_bin
            f = figure('Position',[100,100,1200,600]);
            for d = 1:2
                for year = 1:4
                    subplot(2,4,(d-1)*4+year)
                    scatter(bplv_all{d,y,freq,conn},dplv_all{d,y,freq,conn},[],col(outcome_all{d,y,year},:),'filled');
                    hold on
                    currentXLim = xlim;
                    currentYLim = ylim;
                    % Convert normalized coordinates to data coordinates
                    xData = currentXLim(1) + pos(1) * diff(currentXLim);
                    yData = currentYLim(1) + pos(2) * diff(currentYLim);
                    text(xData,yData,{['R^{2} = ',num2str(R2_all(d,y,freq,conn),'%.2f')], ...
                        ['p = ',num2str(pR2_all(d,y,freq,conn),'%.3f')]})
                    if d == 1
                        ylabel(depth_strings{d},'FontWeight','bold','FontSize',12)
                    end
                    if d == 1
                        if year == 4
                            title(['Last Year'],'FontWeight','bold','FontSize',12)
                        else
                            title(['Year',num2str(year)],'FontWeight','bold','FontSize',12)
                        end
                    end
                end
            end
            h = axes(f,'visible','off');
            h.Title.Visible='on';
            h.XLabel.Visible='on';
            h.YLabel.Visible='on';
            ylabel(h,'Delta PLV','Position',[-0.05,0.500000476837158,0], ...
                'FontWeight','bold','FontSize',14);
            xlabel(h,'Baseline PLV','FontWeight','bold','FontSize',14);
            saveas(f,fullfile(figpath,'04_Baseline_PLV','',['time',num2str(y),'_',conn_strings{conn},'_',freq_strings{freq},'.png']))

        end
    end
end
close all
%% Fig.3B/E
freqs = {'Theta','Alpha','Beta','Gamma'};
cols = {'k','r'};
g = figure('Position',[100,100,1200,600]);
for d = 1:2
    for y = 1:3
        subplot(2,3,(d-1)*3+y)
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
            title(['Year',num2str(y)],'FontWeight','bold','FontSize',12)
        end
    end
end
saveas(g,fullfile(datapath,'figs','Fig.3B.png'))
close all
