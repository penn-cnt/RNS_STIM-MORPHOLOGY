% Plot Variance
%% Settings
function results = var_conn(bin_size,ttdays,option,model)
% variance of plv(%dplv) within time bins, plv_var/dplv_var
% variance across times bins, dplv_var_acrossbin/zplv_var_acrossbin
% variance in slope between time bins, plv_slope_var/dplv_slope_var
% variance in slope between events, dplv_slope_var_acrossbin/zplv_slope_var_acrossbin
% clear; close all;
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
if model == 1 % only random
    formula = strcat([option2,' ~ 1 + BaselinePLV + (1|PtID)']);
elseif model == 2 % add outcome as covariate
    formula = strcat([option2,' ~ 1 + BaselinePLV + Outcome + (1|PtID)']);
elseif model == 3 % add depth as covariate
    formula = strcat([option2,' ~ 1 + BaselinePLV + Depth + (1|PtID)']);
elseif model == 4 % add conntype as covariate
    formula = strcat([option2,' ~ 1 + BaselinePLV + Conn + (1|PtID)']);
elseif model == 5 % add both as covariate
    formula = strcat([option2,' ~ 1 + BaselinePLV + Outcome + Depth + (1|PtID)']);
elseif model == 6 % add both as covariate
    formula = strcat([option2,' ~ 1 + BaselinePLV + Outcome + Depth + Conn + (1|PtID)']);
end
full_result = [];
%% Figure 2 v2 within/between lead
try
    outcome = {plasticity.(option)};
catch
    outcome = {plasticity.(lower(option))};
end
if ndims(outcome{1}) == 2
    tbl_active = tbl2;
    tbl_active{:,'BaselinePLV'} = zscore(log(tbl_active{:,'BaselinePLV'}));
    f = figure('Visible','Off');
    f.Position = [100,100,1200,800];
    for j = 1:3 % freq
        subplot(2,2,j)
        hold on
        tbl = tbl_active(tbl_active.Freq == j,:);
        disp(['Model ',num2str(model),', Outcome ',option,', Freq ',num2str(j)])
        if startsWith(option,'Coeff')
            tbl.(option2) = abs(tbl.(option2));
        end
        [tmp_result,intercept,outliers] = fitglmeBS(tbl,formula,500);
        R2 = tmp_result{1,1}; pR2 = tmp_result{1,2};
        b = tmp_result{2,1}; pCoeff = tmp_result{2,2};
        tmp_result = [tmp_result{1,:},tmp_result{2,:}];
        full_result = [full_result;tmp_result];
        scatter(tbl{:,'BaselinePLV'},tbl{:,option2},[],'k','filled');
        hold on
        scatter(tbl{outliers,'BaselinePLV'},tbl{outliers,option2},[],'r','filled');
        plot(tbl{:,'BaselinePLV'},b*tbl{:,'BaselinePLV'}+intercept)
        if contains(option,'R2')
            ylim([-1,1])
        end
        ylims = get(gca,'YLim');
        xlims = get(gca,'XLim');
        text(xlims(2)-0.25*(xlims(2)-xlims(1)),ylims(1)+0.2*(ylims(2)-ylims(1)), ...
            {['R^2 = ', num2str(R2,'%.2f')];['b = ', num2str(b,'%.3f')];['p_b = ', num2str(pCoeff,'%.3f')]})
        title(freq_strings{j},'FontWeight','bold')
    end
    h = axes(f,'visible','off');
    h.XLabel.Visible='on';
    h.YLabel.Visible='on';
    ylabel(h,option,'Position',[-0.05,0.500000476837158,0], ...
        'FontWeight','bold','FontSize',14);
    xlabel(h,'baseline','FontWeight','bold','FontSize',14);
    saveas(f,fullfile(figpath,'00_Scatter_plot',num2str(model),[option,'.png']))
    results = array2table(full_result, 'RowNames',freq_strings, ...
            'VariableNames',{'R2','pR2','R2_lowerCI','R2_higherCI','b','pCoeff','b_lowerCI','b_higherCI'});
elseif ndims(outcome{1}) == 3
    tbl_active = tbl1;
    tbl_active{:,'BaselinePLV'} = zscore(log(tbl_active{:,'BaselinePLV'}));
    for t = 1:n_bin
        f = figure('Visible','Off');
        f.Position = [100,100,1200,800];
        for j = 1:4 % freq
            subplot(2,2,j)
            hold on
            tbl = tbl_active(tbl_active.Freq == j & tbl_active.Time == t,:);
            disp(['Model ',num2str(model),', Outcome ',option,', Freq ',num2str(j),' Time ',num2str(t)])
            [tmp_result,intercept,outliers] = fitglmeBS(tbl,formula,500);
            R2 = tmp_result{1,1}; pR2 = tmp_result{1,2};
            b = tmp_result{2,1}; pCoeff = tmp_result{2,2};
            tmp_result = [tmp_result{1,:},tmp_result{2,:}];
            full_result = [full_result;tmp_result];
            scatter(tbl{:,'BaselinePLV'},tbl{:,option},[],'k','filled');
            hold on
            scatter(tbl{outliers,'BaselinePLV'},tbl{outliers,option2},[],'r','filled');
            plot(tbl{:,'BaselinePLV'},b*tbl{:,'BaselinePLV'}+intercept)
            title(freq_strings{j},'FontWeight','bold')
            ylims = get(gca,'YLim');
            xlims = get(gca,'XLim');
            text(xlims(2)-0.25*(xlims(2)-xlims(1)),ylims(1)+0.2*(ylims(2)-ylims(1)), ...
                {['R^2 = ', num2str(R2,'%.2f')];['b = ', num2str(b,'%.3f')];['p_b = ', num2str(pCoeff,'%.3f')]})
        end
        h = axes(f,'visible','off');
        h.XLabel.Visible='on';
        h.YLabel.Visible='on';
        ylabel(h,option,'Position',[-0.05,0.500000476837158,0], ...
            'FontWeight','bold','FontSize',14);
        xlabel(h,'baseline','FontWeight','bold','FontSize',14);
        saveas(f,fullfile(figpath,'00_Scatter_plot',num2str(model),[option,'_',num2str(t),'.png']))
    end
    results = array2table(full_result, 'RowNames',strcat(cellstr(num2str(reshape(repmat([1:11],4,1),[],1))),{'_'},repmat(freq_strings',11,1)), ...
            'VariableNames',{'R2','pR2','R2_lowerCI','R2_higherCI','b','pCoeff','b_lowerCI','b_higherCI'});
end
end
