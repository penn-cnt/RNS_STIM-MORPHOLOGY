function var_conn_bar(bin_size,ttdays,option,model)
paths;
ptList = {rns_config.patients.ID};
load(fullfile(datapath,['plasticity_',num2str(bin_size),'.mat']));
plasticity = plasticity(~cellfun('isempty', {plasticity.ptID}));
tbl1 = tbl;
% Lead location x Connection type x Outcome loc x Freq
con_strings = ["1-2","3-4", "1-3","1-4","2-3","2-4"];
freq_strings = {'Theta','Alpha','Beta','Gamma'};
time = 90:bin_size:ttdays-1;
n_bin = length(time);
cols = parula(n_bin);
if ~exist(fullfile(figpath,'00_Scatter_plot',num2str(model)),'dir')
    mkdir(fullfile(figpath,'00_Scatter_plot',num2str(model)));
end
option2 = strrep(option,'_bin','');

    groupCenters = @(nGroups,nMembers,interGroupSpace) ...
    nGroups/2+.5 : nGroups+interGroupSpace : (nGroups+interGroupSpace)*nMembers-1;
    ticks = groupCenters(2, 4, 1);
    f = figure('Visible','Off');
    f.Position = [100,100,1200,600];
    sr = nan*zeros(156,4);
    pr = nan*zeros(156,4);
    p = nan*zeros(4,1);
    for j = 1:4 
        tbl_active = tbl2;
        tbl = tbl_active(tbl_active.Freq == j,:);
        pvalues = tbl{:,option2};
        b = tbl{:,'BaselinePLV'};
        inds = pvalues < 0.05;
        sr(1:length(find(inds)),j) = b(inds);
        pr(1:length(find(~inds)),j) = b(~inds);
        p(j) = ranksum(b(inds),b(~inds));
    end
    box_data{1} = pr; % green, non-sig
    box_data{2} = sr; % blue, sig
    b = boxplotGroup(box_data,'primaryLabels',repmat({''}, 2, 1), ...
            'Colors',[0,1,0;0,0,1],'GroupType','betweenGroups', ...
            'PlotStyle','traditional','BoxStyle','outline', ...
            'Symbol','o','Widths',0.7);
    ymax = get(gca, 'YLim');
    for j = 1:4
        text(ticks(j),ymax(2)-0.1,['p=',num2str(p(j),'%.3f')])
    end
    ylabel('Baseline PLV')
    set(gca,'XTick',ticks,'XTickLabels',freq_strings)
    saveas(f,fullfile(figpath,'00_Scatter_plot',num2str(model),[option,'_bar.png']))
end