%% baseline_dplv.m
% Tests whether baseline connectivity has relationship with %change in connectivity at different time points.
%% Settings
function baseline_dplv(bin_size,ttdays,regress_dist,outcome_option,option)
close all;
paths;
patient_info = struct2table(load(which('patients_Penn.mat')).patients_Penn);
ptList = {rns_config.patients.ID};
localization = load(fullfile(datapath,"localization.mat")).localization;
base_days = 90;
if regress_dist
    suffix = '_regdist';
    figpath = fullfile(figpath,'regdist');
else
    suffix = '';
end
plasticity = load(fullfile(datapath,['plasticity_',num2str(bin_size),suffix,'.mat'])).plasticity;
time = 90:bin_size:ttdays-1;
n_bin = length(time);
if ~exist(fullfile(figpath,'04_Baseline_PLV',num2str(outcome_option),option(7:end)),'dir')
    mkdir(fullfile(figpath,'04_Baseline_PLV',num2str(outcome_option),option(7:end)));
    mkdir(fullfile(figpath,'04_Baseline_PLV',num2str(outcome_option),option(7:end),'inter'));
    mkdir(fullfile(figpath,'04_Baseline_PLV',num2str(outcome_option),option(7:end),'intra'));
end
for i = 1:length(localization)
    if ~isempty(localization(i).outcome_group)
        localization(i).outcome_group = localization(i).outcome_group(outcome_option,:);
    end
    if ~isempty(plasticity(i).outcome_group)
        plasticity(i).outcome_group = plasticity(i).outcome_group(outcome_option,:);
    end
end
plasticity = paddata(plasticity,option,n_bin);
rng('default');
%% Baseline connectivity predict outcome
% try
%     data_table = load(fullfile(datapath,'baseline_plv_pred.mat')).data_table;
% catch
    ID = []; d = []; o_group = [];
    bp_inter = []; bp_intra = []; dp_inter = []; dp_intra = []; t = []; s = [];
    for pt = 1:length(localization)
        % Read Patient Data
        ptID = localization(pt).ptID;

        outcome_group = localization(pt).outcome_group;
        depth = localization(pt).depth;
        lead_label = localization(pt).lead_labels;

        if ~localization(pt).meets_criteria | isempty(plasticity(pt).ptID)
            continue
        end

        if length(outcome_group) > 3
            years = [1,2,3,length(outcome_group)];
            outcome_group = outcome_group(years);
        elseif isempty(outcome_group)
            outcome_group = [3,3,3,3];
        else
            years = [1:length(outcome_group)];
            tmp = outcome_group(years);
            tmp = padarray(tmp,[0,3-length(years)],3,'post');
            outcome_group = [tmp,outcome_group(end)];
        end

        % baseline
        baseline_plv = plasticity(pt).baseline_plvs;
        if any(lead_label == 1)
            baseline_plv_intra = squeeze(mean(baseline_plv(lead_label == 1,:),1,'omitnan'));
            dplv_intra = cellfun(@(x) mean(x,2,'omitnan'),plasticity(pt).(option)(1,:),'UniformOutput',false);
            dplv_intra = horzcat(dplv_intra{:});
        else
            baseline_plv_intra = nan * zeros(1,4);
            dplv_intra = nan * zeros(n_bin,4);
        end
        if any(lead_label == 2)
            baseline_plv_inter = squeeze(mean(baseline_plv(lead_label == 2,:),1,'omitnan'));
            dplv_inter = cellfun(@(x) mean(x,2,'omitnan'),plasticity(pt).(option)(2,:),'UniformOutput',false);
            dplv_inter = horzcat(dplv_inter{:});
        else
            baseline_plv_inter = nan * zeros(1,4);
            dplv_inter = nan * zeros(n_bin,4);
        end

        for tt = 1:n_bin
            ID = [ID;pt];
            d = [d;depth];
            o_group = [o_group;outcome_group];
            bp_inter = [bp_inter;baseline_plv_inter];
            bp_intra = [bp_intra;baseline_plv_intra];
            dp_inter = [dp_inter;dplv_inter(tt,:)];
            dp_intra = [dp_intra;dplv_intra(tt,:)];
            t = [t;tt];
        end
    end
    data_mat = [ID,d,o_group,bp_intra,bp_inter,dp_intra,dp_inter,t];
    data_table = array2table(data_mat,"VariableNames", ...
        ["IDs","Depth",...
        "outcome_group_Year1", "outcome_group_Year2", "outcome_group_Year3", "outcome_group_Last", ...
        "bplv_intra_Theta","bplv_intra_Alpha","bplv_intra_Beta","bplv_intra_Gamma", ...
        "bplv_inter_Theta","bplv_inter_Alpha","bplv_inter_Beta","bplv_inter_Gamma", ...
        "dplv_intra_Theta","dplv_intra_Alpha","dplv_intra_Beta","dplv_intra_Gamma", ...
        "dplv_inter_Theta","dplv_inter_Alpha","dplv_inter_Beta","dplv_inter_Gamma", ...
        "time"]);
    save(fullfile(datapath,['baseline_plv_pred',suffix,'.mat']),"data_table")
% end
%% prediction
R2_all = nan*zeros(2,n_bin,4,2);
pR2_all = nan*zeros(2,n_bin,4,2);
dplv_all = cell(2,n_bin,4,2);
bplv_all = cell(2,n_bin,4,2);
outcome_all = cell(2,n_bin,4);
freq = {'Theta','Alpha','Beta','Gamma'};
conn = {'intra','inter'};
for y = 1:n_bin
    for d = 1:2
    % one model for each depth type
        sub_data = data_table(data_table.Depth == d & data_table.time == y,:);
        outcome_all{d,y,1} = table2array(sub_data(:,'outcome_group_Year1'));
        outcome_all{d,y,2} = table2array(sub_data(:,'outcome_group_Year2'));
        outcome_all{d,y,3} = table2array(sub_data(:,'outcome_group_Year3'));
        outcome_all{d,y,4} = table2array(sub_data(:,'outcome_group_Last'));
        for f = 1:4
            for c = 1:2
                o_true = table2array(sub_data(:,['dplv_',conn{c},'_',freq{f}]));
                pred = table2array(sub_data(:,['bplv_',conn{c},'_',freq{f}]));
                [r,p] = corr(pred,o_true,'Type','Spearman','rows', 'complete');
                % record
                R2_all(d,y,f,c) = r^2;
                pR2_all(d,y,f,c) = p;
                dplv_all{d,y,f,c} = o_true;
                bplv_all{d,y,f,c} = pred;
            end
        end
    end
end
save(fullfile(datapath,['baseline_plv_pred',suffix,'.mat']),'-append','R2_all','pR2_all','dplv_all','bplv_all','outcome_all');
%% plotting
% Fig.3A/D
depth_strings = {'Hippocampal','Neocortical'};
freq_strings = {'Theta','Alpha','Beta','Gamma'};
conn_strings = {'intra','inter'};
pos = [0.9,0.1];
col = [0,1,0;0,0,1;0,0,0];
load(fullfile(datapath,['baseline_plv_pred',suffix,'.mat']))
years = {1,2,3,'end'};
for year = [2,4]
    % just plot the first year and last year
    for conn = 1:2
        for y = 1:n_bin
            f = figure('Position',[100,100,1200,600],'Visible','Off');
            for d = 1:2
                for freq = 1:4
                    subplot(2,4,(d-1)*4+freq)
                    scatter(bplv_all{d,y,freq,conn},dplv_all{d,y,freq,conn},[],col(outcome_all{d,y,year},:),'filled');
                    hold on
                    currentXLim = xlim;
                    currentYLim = ylim;
                    % Convert normalized coordinates to data coordinates
                    xData = currentXLim(1) + pos(1) * diff(currentXLim);
                    yData = currentYLim(1) + pos(2) * diff(currentYLim);
                    text(xData,yData,{['R^{2} = ',num2str(R2_all(d,y,freq,conn),'%.2f')], ...
                        ['p = ',num2str(pR2_all(d,y,freq,conn),'%.3f')]})
                    if freq == 1
                        ylabel(depth_strings{d},'FontWeight','bold','FontSize',12)
                    end
                    if d == 1
                        title(freq_strings{freq},'FontWeight','bold','FontSize',12)
                    end
                end
            end
            h = axes(f,'visible','off');
            h.Title.Visible='on';
            h.XLabel.Visible='on';
            h.YLabel.Visible='on';
            ylabel(h,['Delta',option(7:end)],'Position',[-0.05,0.500000476837158,0], ...
                'FontWeight','bold','FontSize',14);
            xlabel(h,'Baseline PLV','FontWeight','bold','FontSize',14);
            saveas(f,fullfile(figpath,'04_Baseline_PLV',num2str(outcome_option),option(7:end),conn_strings{conn},['time',num2str(y),'_year',num2str(years{year}),'.png']))
        end
    end
    close all
end
end
% not working so well
% %% prediction outcome separate
% outcome_all = load(fullfile(datapath,'baseline_plv_pred.mat')).outcome_all;
% years = {1,2,3,'end'};
% for year = 1:length(years)
%     R2_all = nan*zeros(2,2,n_bin,4,2);
%     pR2_all = nan*zeros(2,2,n_bin,4,2);
%     dplv_all = cell(2,2,n_bin,4,2);
%     bplv_all = cell(2,2,n_bin,4,2);
%     freq = {'Theta','Alpha','Beta','Gamma'};
%     conn = {'intra','inter'};
%     for y = 1:n_bin
%         for d = 1:2
%         % one model for each depth type
%             sub_data = data_table(data_table.Depth == d & data_table.time == y,:);
%             outcome = outcome_all{d,y,year};
%             for f = 1:4
%                 for c = 1:2
%                     o_true = table2array(sub_data(:,['dplv_',conn{c},'_',freq{f}]));
%                     pred = table2array(sub_data(:,['bplv_',conn{c},'_',freq{f}]));
%                     for o = 1:2
%                         o_true_sub = o_true(outcome == o);
%                         pred_sub = pred(outcome == o);
%                         if isempty(o_true_sub)
%                             continue
%                         end
%                         [r,p] = corr(pred_sub,o_true_sub,'Type','Spearman','rows', 'complete');
%                         % record
%                         R2_all(o,d,y,f,c) = r^2;
%                         pR2_all(o,d,y,f,c) = p;
%                         dplv_all{o,d,y,f,c} = o_true_sub;
%                         bplv_all{o,d,y,f,c} = pred_sub;
%                     end
%                 end
%             end
%         end
%     end
%     save(fullfile(datapath,['baseline_plv_pred_',num2str(years{year}),'.mat']),'R2_all','pR2_all','dplv_all','bplv_all');
% end
% %% line plotting
% ylims = [40,20];
% time = 90:bin_size:ttdays-1;
% n_bin = length(time);
% depth_strings = {'Hippocampal','Neocortical'};
% freq_strings = {'Theta','Alpha','Beta','Gamma'};
% con_strings = {'intra','inter'};
% outcome_strings = {'SR','PR'};
% years = {1,2,3,'end'};
% cols = 'gb';
% pos = [-0.4,-0.3];
% for year = [1:4]
%     load(fullfile(datapath,['baseline_plv_pred_',num2str(years{year}),'.mat']))
%     for conn = 1:2
%         g = figure('Position',[100,100,1200,600]);
%         for d = 1:2
%             for f = 1:4
%                 subplot(2,4,sub2ind([4,2],f,d))
%                 hold on
%                 for o = 1:2 % outcome
%                     data = R2_all(o,d,:,f,conn);
%                     p = pR2_all(o,d,:,f,conn);
%                     if ~isempty(data)
%                         shadedErrorBar(time,mean(data,1,'omitnan'), ...
%                             std(data,[],1,'omitnan')/sqrt(size(data,1)), ...
%                             "lineProps",['-',cols(o)],'transparent',1);
%                         xlim([90,ttdays])
%                         ylim([-0.5,1])
%                         scatter(time(p < 0.05),pos(o)*ones(length(time(p < 0.05)),1),[],cols(o),'filled');
%                     end
%                 end
%                 if d == 1
%                     title(freq_strings{f},'FontWeight','bold')
%                 end
%                 if f == 1
%                     ylabel(depth_strings{d},'FontWeight','bold')
%                 end
%             end
%         end
%         h = axes(g,'visible','off');
%         h.XLabel.Visible='on';
%         h.YLabel.Visible='on';
%         ylabel(h,'R^2','Position',[-0.05,0.500000476837158,0], ...
%                 'FontWeight','bold','FontSize',14);
%         xlabel(h,'Days','FontWeight','bold','FontSize',14);
%         saveas(g,fullfile(figpath,'04_Baseline_PLV',con_strings{conn},['year_',num2str(years{year}),'_outcome_comp.png']))
%     end
%     close all
% end