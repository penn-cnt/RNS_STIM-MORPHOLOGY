% Plot Variance
%% Settings
function var_conn_outcome(bin_size,ttdays,regress_dist,option)
% variance of plv(%dplv) within time bins, plv_var/dplv_var
% variance across times bins, dplv_var_acrossbin/zplv_var_acrossbin
% variance in slope between time bins, plv_slope_var/dplv_slope_var
% variance in slope between events, dplv_slope_var_acrossbin/zplv_slope_var_acrossbin
% clear; close all;
paths;
ptList = {rns_config.patients.ID};
if regress_dist
    suffix = '_regdist';
    figpath = fullfile(figpath,'regdist');
else
    suffix = '';
end
plasticity = load(fullfile(datapath,['plasticity_',num2str(bin_size),suffix,'.mat'])).plasticity;
plasticity = plasticity(~cellfun('isempty', {plasticity.ptID}));
plasticity = plasticity(~cellfun('isempty', {plasticity.outcome}));
plasticity = plasticity(~cellfun('isempty', {plasticity.depth}));

% Lead location x Connection type x Outcome loc x Freq
con_strings = ["1-2","3-4", "1-3","1-4","2-3","2-4"];
freq_strings = {'Theta','Alpha','Beta','Gamma'};
time = 90:bin_size:ttdays-1;
n_bin = length(time);
cols = parula(n_bin);
col = 'gb';
if ~exist(fullfile(figpath,'00_Scatter_plot',option),'dir')
    mkdir(fullfile(figpath,'00_Scatter_plot',option));
end
%% Figure 2 v2 within/between lead
PLV = {plasticity.(option)};
baseline = {plasticity.baseline_plvs};
if ndims(PLV{1}) == 2
    f = figure('Visible','Off');
    f.Position = [100,100,1200,800];
    for j = 1:4 % freq
        subplot(2,2,j)
        hold on
        plv_traces = cellfun(@(x) squeeze(x(:,j)),PLV,'UniformOutput',false);
        base = cellfun(@(x) x(:,j),baseline,'UniformOutput',false);
        outcome = cellfun(@(x) x(3,end),{plasticity.outcome_group});
        depth = cellfun(@(x) x(end),{plasticity.depth});
        for o = 1:2
            traces = vertcat(plv_traces{outcome == o});
            b = vertcat(base{outcome == o});
            scatter(b,traces,[],col(o),'filled');
            hold on
        end
        ylims = get(gca,'YLim');
        xlims = get(gca,'XLim');
        for o = 1:2
            traces = vertcat(plv_traces{outcome == o});
            b = vertcat(base{outcome == o});
            [r,p] = corr(b,traces,'Rows','complete');
            text(xlims(2)-0.25*(xlims(2)-xlims(1)),ylims(1)+0.2*o*(ylims(2)-ylims(1)),{['r = ',num2str(r,'%.2f')];['R^2 = ', num2str(r^2,'%.2f')];['pValue = ', num2str(p,'%.3f')]},'Color',col(o))
        end
        title(freq_strings{j},'FontWeight','bold')
    end
    h = axes(f,'visible','off');
    h.XLabel.Visible='on';
    h.YLabel.Visible='on';
    ylabel(h,option,'Position',[-0.05,0.500000476837158,0], ...
        'FontWeight','bold','FontSize',14);
    xlabel(h,'baseline','FontWeight','bold','FontSize',14);
    saveas(f,fullfile(figpath,'00_Scatter_plot',option,[option,'_outcome.png']))
elseif ndims(PLV{1}) == 3
    plasticity = paddata(plasticity,option,n_bin);
    PLV = {plasticity.(option)};
    baseline = {plasticity.baseline_plvs};
    for t = 1:n_bin
        f = figure('Visible','Off');
        f.Position = [100,100,1200,800];
        outcome = cellfun(@(x) x(3,end),{plasticity.outcome_group});
        depth = cellfun(@(x) x(end),{plasticity.depth});
        for j = 1:4 % freq
            subplot(2,2,j)
            plv_traces = cellfun(@(x) squeeze(x(t,:,j)),PLV,'UniformOutput',false);
            base = cellfun(@(x) x(:,j),baseline,'UniformOutput',false);
            for o = 1:2
                traces = horzcat(plv_traces{outcome == o})';
                b = vertcat(base{outcome == o});
                scatter(b,traces,[],col(o),'filled');
                hold on
            end
            ylims = get(gca,'YLim');
            xlims = get(gca,'XLim');
            title(freq_strings{j},'FontWeight','bold')
            for o = 1:2
                traces = horzcat(plv_traces{outcome == o})';
                b = vertcat(base{outcome == o});
                [r,p] = corr(b,traces,'Rows','complete');
                text(xlims(2)-0.25*(xlims(2)-xlims(1)),ylims(1)+0.2*o*(ylims(2)-ylims(1)),{['r = ',num2str(r,'%.2f')];['R^2 = ', num2str(r^2,'%.2f')];['pValue = ', num2str(p,'%.3f')]},'Color',col(o))
            end
        end
        h = axes(f,'visible','off');
        h.XLabel.Visible='on';
        h.YLabel.Visible='on';
        ylabel(h,option,'Position',[-0.05,0.500000476837158,0], ...
            'FontWeight','bold','FontSize',14);
        xlabel(h,'baseline','FontWeight','bold','FontSize',14);
        saveas(f,fullfile(figpath,'00_Scatter_plot',option,[num2str(t),'_outcome.png']))
    end
end
% if contains(option,'pValue')
%     groupCenters = @(nGroups,nMembers,interGroupSpace) ...
%     nGroups/2+.5 : nGroups+interGroupSpace : (nGroups+interGroupSpace)*nMembers-1;
%     ticks = groupCenters(2, 4, 1);
%     f = figure('Visible','Off');
%     f.Position = [100,100,1200,600];
%     sr = nan*zeros(156,4);
%     pr = nan*zeros(156,4);
%     p = nan*zeros(4,1);
%     for j = 1:4 
%         plv_traces = cellfun(@(x) squeeze(x(:,j)),outcome,'UniformOutput',false);
%         plv_traces = vertcat(plv_traces{:});
%         b = cellfun(@(x) x(:,j),baseline,'UniformOutput',false);
%         b = vertcat(b{:});
%         inds = plv_traces < 0.05;
%         sr(1:length(find(inds)),j) = b(inds);
%         pr(1:length(find(~inds)),j) = b(~inds);
%         p(j) = ranksum(b(inds),b(~inds));
%     end
%     box_data{1} = pr; % green, non-sig
%     box_data{2} = sr; % blue, sig
%     b = boxplotGroup(box_data,'primaryLabels',repmat({''}, 2, 1), ...
%             'Colors',[0,1,0;0,0,1],'GroupType','betweenGroups', ...
%             'PlotStyle','traditional','BoxStyle','outline', ...
%             'Symbol','o','Widths',0.7);
%     ymax = get(gca, 'YLim');
%     for j = 1:4
%         text(ticks(j),ymax(2)-0.1,['p=',num2str(p(j),'%.3f')])
%     end
%     ylabel('Baseline PLV')
%     set(gca,'XTick',ticks,'XTickLabels',freq_strings)
%     saveas(f,fullfile(figpath,'00_Scatter_plot',option,[option,'_bar.png']))
% end
