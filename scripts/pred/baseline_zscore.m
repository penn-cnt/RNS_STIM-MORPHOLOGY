%% baseline_zscore.m
% Tests whether baseline connectivity has relationship with zscored connectivity at different time points.
%% Settings
clear; close all;
paths;
patient_info = struct2table(load(which('patients_Penn.mat')).patients_Penn);
ptList = {rns_config.patients.ID};
localization = load(fullfile(datapath,"localization.mat")).localization;
base_days = 90;
bin_size = 90; % subject to change
ttdays = 1050;
plasticity = load(fullfile(datapath,['plasticity_',num2str(bin_size),'.mat'])).plasticity;
time = 90:bin_size:ttdays-1;
n_bin = length(time);
for n = 1:length(plasticity)
    try
        plasticity(n).reorg_zplv = cellfun(@(x) x(1:n_bin,:),plasticity(n).reorg_zplv,'UniformOutput',false);
    catch
        plasticity(n).reorg_zplv = cellfun(@(x) padarray(x,n_bin-size(x,1),nan,'post'),plasticity(n).reorg_zplv,'UniformOutput',false);
    end
end
rng('default');
%% Baseline connectivity predict outcome
try
    data_table = load(fullfile(datapath,'baseline_plv_zscore.mat')).data_table;
catch
    ID = []; d = []; o_group = [];
    bp_inter = []; bp_intra = []; zp_inter = []; zp_intra = []; t = [];
    for p = 1:length(plasticity)
        % Read Patient Data
        ptID = plasticity(p).ptID;
        pt = find(strcmp(ptID,{localization.ptID}));

        outcome_group = localization(pt).outcome_group;
        depth = localization(pt).depth;

        if ~localization(pt).meets_criteria
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

        baseline_plv = plasticity(p).baseline_plvs;
        baseline_plv_intra = squeeze(mean(baseline_plv(1:2,:),1,'omitnan'));
        baseline_plv_inter = squeeze(mean(baseline_plv(3:6,:),1,'omitnan'));

        % dplv at different time points
        zplv_intra = cellfun(@(x) mean(x,2,'omitnan'),plasticity(p).reorg_zplv(1,:),'UniformOutput',false);
        zplv_intra = horzcat(zplv_intra{:});
        zplv_inter = cellfun(@(x) mean(x,2,'omitnan'),plasticity(p).reorg_zplv(2,:),'UniformOutput',false);
        zplv_inter = horzcat(zplv_inter{:});

        for tt = 1:n_bin
            ID = [ID;pt];
            d = [d;depth];
            o_group = [o_group;outcome_group];
            bp_inter = [bp_inter;baseline_plv_inter];
            bp_intra = [bp_intra;baseline_plv_intra];
            zp_inter = [zp_inter;zplv_inter(tt,:)];
            zp_intra = [zp_intra;zplv_intra(tt,:)];
            t = [t;tt];
        end
    end
    data_mat = [ID,d,o_group,bp_intra,bp_inter,zp_intra,zp_inter,t];
    data_table = array2table(data_mat,"VariableNames", ...
        ["IDs","Depth",...
        "OutcomeGroup_Year1", "OutcomeGroup_Year2", "OutcomeGroup_Year3", "OutcomeGroup_Last", ...
        "bplv_intra_Theta","bplv_intra_Alpha","bplv_intra_Beta","bplv_intra_Gamma", ...
        "bplv_inter_Theta","bplv_inter_Alpha","bplv_inter_Beta","bplv_inter_Gamma", ...
        "zplv_intra_Theta","zplv_intra_Alpha","zplv_intra_Beta","zplv_intra_Gamma", ...
        "zplv_inter_Theta","zplv_inter_Alpha","zplv_inter_Beta","zplv_inter_Gamma", ...
        "time"]);
    save(fullfile(datapath,'baseline_plv_zscore.mat'),"data_table")
end
%% prediction
R2_all = nan*zeros(2,n_bin,4,2);
pR2_all = nan*zeros(2,n_bin,4,2);
zplv_all = cell(2,n_bin,4,2);
bplv_all = cell(2,n_bin,4,2);
outcome_all = cell(2,n_bin,4);
freq = {'Theta','Alpha','Beta','Gamma'};
conn = {'intra','inter'};
for y = 1:n_bin
    for d = 1:2
        % one model for each depth type
        sub_data = data_table(data_table.Depth == d & data_table.time == y,:);
        outcome_all{d,y,1} = table2array(sub_data(:,'OutcomeGroup_Year1'));
        outcome_all{d,y,2} = table2array(sub_data(:,'OutcomeGroup_Year2'));
        outcome_all{d,y,2} = table2array(sub_data(:,'OutcomeGroup_Year3'));
        outcome_all{d,y,4} = table2array(sub_data(:,'OutcomeGroup_Last'));
        for f = 1:4
            for c = 1:2
                o_true = table2array(sub_data(:,['zplv_',conn{c},'_',freq{f}]));
                pred = table2array(sub_data(:,['bplv_',conn{c},'_',freq{f}]));
                [r,p] = corr(pred,o_true,'Type','Spearman','rows', 'complete');
                % record
                R2_all(d,y,f,c) = r^2;
                pR2_all(d,y,f,c) = p;
                zplv_all{d,y,f,c} = o_true;
                bplv_all{d,y,f,c} = pred;
            end
        end
    end
end
save(fullfile(datapath,'baseline_plv_zscore.mat'),'-append','R2_all','pR2_all','zplv_all','bplv_all','outcome_all');
%% plotting
% Fig.3A/D
depth_strings = {'Hippocampal','Neocortical'};
freq_strings = {'Theta','Alpha','Beta','Gamma'};
conn_strings = {'intra','inter'};
years = {1,2,3,'end'};
pos = [0.9,0.1];
col = [0,1,0;0,0,1;0,0,0];
for year = [1,4]
    % just plot the first year and last year
    for conn = 1:2
        for y = 1:n_bin
            f = figure('Position',[100,100,1200,600]);
            for d = 1:2
                for freq = 1:4
                    subplot(2,4,(d-1)*4+freq)
                    scatter(bplv_all{d,y,freq,conn},zplv_all{d,y,freq,conn},[],col(outcome_all{d,y,year},:),'filled');
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
            ylabel(h,'Delta PLV','Position',[-0.05,0.500000476837158,0], ...
                'FontWeight','bold','FontSize',14);
            xlabel(h,'Baseline PLV','FontWeight','bold','FontSize',14);
            saveas(f,fullfile(figpath,'04_Baseline_PLV','zscore',conn_strings{conn},['time',num2str(y),'_year',num2str(years{year}),'.png']))
        end
    end
    close all
end