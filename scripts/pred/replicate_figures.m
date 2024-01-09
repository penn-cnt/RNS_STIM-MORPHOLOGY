%% replicate_figures.m
% This script replicates the figure 2 and analysis of Ankit's paper,
% including:
%   1) Figure 2
%% Settings
clear; close all;
paths;
ptList = {rns_config.patients.ID};
bin_size = 90; % subject to change
ttdays = 1050;
plasticity = load(fullfile(datapath,['new_ver/plasticity_',num2str(bin_size),'.mat'])).plasticity;
plasticity = plasticity(~cellfun('isempty', {plasticity.outcome}));
plasticity = plasticity(~cellfun('isempty', {plasticity.depth}));
% Lead location x Connection type x Outcome loc x Freq
con_strings = ["Inside SOZ", "Between SOZ", "SOZ to Normal","Normal to Normal"];
depth_strings = {'Hippocampal','Neocortical'};
outcome_strings = ["Poor Responder","Good Responder"];
colors = ['g','b'];
% max_time_length = max(cellfun(@(x) size(x,2),dplv_traces),[],'all');
time = 90:bin_size:ttdays-1;
%% Figure 2 v1
year = 1;
ylims = [40,20];
for d = 1:2
    f = figure(d);
    f.Position = [100,100,1200,600];
    for i = 1:2 % conntype
        for j = 1:4 % freq
            subplot(2,4,sub2ind([4,2],j,i))
            hold on
            for o = 1:2 % outcome
                tmp = plasticity(cellfun(@(x) x==d,{plasticity.depth}) & ...
                                 cellfun(@(x) x(year)==o,{plasticity.outcome}));
                traces = cellfun(@(x) x{i+2,j},{tmp.reorg_dplv},'UniformOutput',false);
                traces = horzcat(traces{:})';
                if ~isempty(traces)
                    traces = traces(:,1:length(time));
                    shadedErrorBar(time,mean(traces,1,'omitnan'), ...
                        std(traces,[],1,'omitnan')/ sqrt(size(traces,1)), ...
                        "lineProps",['-',colors(o)],'transparent',1);
                    ylabel(con_strings(i))
                    xlim([90,ttdays])
                end
            end
            % FDA test and plot
            fda_data = plasticity(cellfun(@(x) x==d,{plasticity.depth}));
            for n = 1:length(fda_data)
                fda_data(n).outcome = fda_data(n).outcome(year);
                fda_data(n).reorg_dplv = fda_data(n).reorg_dplv{i+2,j};
            end
            outcome = [fda_data.outcome];
            data = {fda_data.reorg_dplv};
            [fda_pv,sig_days] = FDA_test(data,outcome);
            text(ttdays-200,-ylims(i)+5,['p=',num2str(fda_pv,'%.3f')],'FontSize',12);
            if ~isempty(sig_days)
                scatter(time(sig_days(:,1)),zeros(size(sig_days,1),1),[],'r','filled');
            end
        end
    end
    figure(d)
    sgtitle(depth_strings(d))
    legend(outcome_strings,'Position',[0.8,0.87,0.1,0.05]);
    saveas(f,fullfile(datapath,'figs',[depth_strings{d},'_',num2str(bin_size),'_',num2str(year),'_soz.png']))
end
close all
%% Figure 2 v2
year = 3;
con_strings2 = ["Within Lead", "Between Lead"];
ylims = [40,20];
for d = 1:2
    f = figure(d);
    f.Position = [100,100,1200,600];
    for i = 1:2 % conntype
        for j = 1:4 % freq
            subplot(2,4,sub2ind([4,2],j,i))
            hold on
            for o = 1:2 % outcome
                tmp = plasticity(cellfun(@(x) x==d,{plasticity.depth}) & ...
                                 cellfun(@(x) x(year)==o,{plasticity.outcome}));
                traces = cellfun(@(x) x{i,j},{tmp.reorg_dplv},'UniformOutput',false);
                traces = horzcat(traces{:})';
                if ~isempty(traces)
                    traces = traces(:,1:length(time));
                    traces = traces(exclude_outlier(traces,1.5),:);
                    shadedErrorBar(time,mean(traces,1,'omitnan'), ...
                        std(traces,[],1,'omitnan')/ sqrt(size(traces,1)), ...
                        "lineProps",['-',colors(o)],'transparent',1);
                    ylabel(con_strings2(i))
                    xlim([90,ttdays])
                    ylim([-ylims(i),ylims(i)])
                end
            end
            % FDA test and plot
            fda_data = plasticity(cellfun(@(x) x==d,{plasticity.depth}));
            for n = 1:length(fda_data)
                fda_data(n).outcome = fda_data(n).outcome(year);
                fda_data(n).reorg_dplv = fda_data(n).reorg_dplv{i,j};
            end
            outcome = [fda_data.outcome];
            data = {fda_data.reorg_dplv};
            [fda_pv,sig_days] = FDA_test(data,outcome);
            text(ttdays-200,-ylims(i)+5,['p=',num2str(fda_pv,'%.3f')],'FontSize',12);
            if ~isempty(sig_days)
                scatter(time(sig_days(:,1)),zeros(size(sig_days,1),1),[],'r','filled');
            end
        end
    end
    figure(d)
    sgtitle(depth_strings{d})
    legend(outcome_strings,'Position',[0.8,0.87,0.1,0.05]);
    saveas(f,fullfile(datapath,'figs',[depth_strings{d},'_',num2str(bin_size),'_',num2str(year),'_lead.png']))
end
close all
%%
figure(3);
for i = 1:4 % conntype
    for j = 1:4 % freq
        subplot(4,4,sub2ind([4,4],j,i))
        hold on
        for o = 1:2 % outcome
            trace = [dplv_traces{1,i,o,j};dplv_traces{2,i,o,j}];
            if ~isempty(trace)
                shadedErrorBar(time,median(trace,1,'omitnan'), ...
                    std(trace,[],1,'omitnan')/sqrt(size(trace,1)), ...
                    "lineProps",['-',colors(o)],'transparent',1);
                ylabel(con_strings(i))
                xlim([90,800])
                ylim([-0.5,0.5])
            end
        end
    end
end

% if plot_figure
%     % get plv instead of dplv time trace
% 
% 
% end
% if plot_figure
%     figure(pt)
%     hold on
%     for i = 1:size(band_limited_plvs,2) % each connection
%         for j = 1:size(dplv,3) % each frequency
%             subplot(7,1,i)
%             hold on
%             plot(binned_times(logical(binned_times)),binned_dplvs(:,i,j))
%             xlim([90,max(binned_times)])
%             mdl = fitglm(test(logical(test)),binned_stims(logical(binned_stims)),"y~x1","Distribution","poisson");
%             title(num2str(mdl.Rsquared.Ordinary))
%         end
%     end
%     subplot(7,1,7)
%     plot(binned_times(logical(binned_times)),binned_counts(logical(binned_counts)))
%     xlim([90,max(binned_times)])
%     sgtitle(ptID)
%     %% Network NMF Trajectories
%     [W,H,Q] = betaNTF(dplv,3);
%     % figure
%     % plot(time_trace(~baseline_mask),W(:,3))
%     mdl = fitglm(W(:,3),stim_trace(~baseline_mask)','y ~ x1','Distribution',"poisson");
%     R2 = mdl.Rsquared.Ordinary;
% end
% 
