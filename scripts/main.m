%% settings
paths;
%% preprocesssing
soz_assignments;
event_preprocess;
network_trajectory;
%% bar plot
outcomes = {'p_time','p_stim','p_time_bin'};
for m = [1,5]
    for i = length(outcomes)
        var_conn_bar(90,1050,outcomes{i},m);
    end
end
%% variance part

outcomes1 = {'PLV','dPLV','zPLV','PLV_std','dPLV_std','zPLV_std',...
    'PLV_Slope','dPLV_Slope','zPLV_Slope',...
    'PLV_Slope_std','dPLV_Slope_std','zPLV_Slope_std'};
outcomes2 = {'PLV_std_bin','dPLV_std_bin','zPLV_std_bin',...
    'PLV_Slope_bin','dPLV_Slope_bin','zPLV_Slope_bin',...
    'PLV_Slope_std_bin','dPLV_Slope_std_bin','zPLV_Slope_std_bin',...
    'R2_time','p_time','Coeff_time',...
    'R2_stim','p_stim','Coeff_stim',...
    'R2_time_bin','p_time_bin','Coeff_time_bin'};
for m = [1,5]
    for i = 3:length(outcomes1)
        results = var_conn(90,1050,outcomes1{i},m);
        names = results.Properties.VariableNames;
        results.Properties.VariableNames = strcat(repmat({['M',num2str(m)]},1,8),'_',names);
        writetable(results,['user_data/stats/',outcomes1{i},'_M',num2str(m),'_log.csv']);
    end
end
% stats_tbl = array2table(array,'RowNames',outcomes,'VariableNames', ...
%     {'M1_R2','M1_pR2','M1_Coeff','M1_pCoeff',...
%     'M5_R2','M5_pR2','M5_Coeff','M5_pCoeff', ...
%     'M6_R2','M6_pR2','M6_Coeff','M6_pCoeff'});

%     'M2_R2','M2_pR2','M2_Coeff','M2_pCoeff', ...
%     'M3_R2','M3_pR2','M3_Coeff','M3_pCoeff', ...
%     'M4_R2','M4_pR2','M4_Coeff','M4_pCoeff', ...

%%
for i = 1:length(outcomes)
    var_conn_outcome(90,1050,0,outcomes{i});
end
%%
outcomes = {'plv_var','dplv_var','plv_slope_var','dplv_slope_var'};
for i = 1:length(outcomes)
    var_Fig2(90,1050,0,4,outcomes{i},{'end'});
end
%% variance part
outcomes = {'dplv_var_acrossbin','zplv_var_acrossbin','dplv_slope_var_acrossbin','zplv_slope_var_acrossbin'};
for i = 1:length(outcomes)
    var_Barplot(90,1050,0,4,outcomes{i},1,{'end'});
end
%% replicating figures

%% Fig2
bin_size = 90; % subject to change
ttdays = 1050; % total days to plot traces
regress_dist = 0; 
% 1. 50%
% 2. 50%+exclude
% 3. 80%
% 4. 50%+assign exclude to 1
% 5. 50%+exclude
% 6. 80% -2,50-80% exclude,others-1
% 7. 80%-2,<50% 1, others exclude
outcome_option = 5; 
option = {'reorg_dplv','reorg_zplv'};
years = {'end'};
for outcome_option = 6%3:5
    for regress_dist = 0:1
        for o = 1:2
            Fig2(bin_size,ttdays,regress_dist,outcome_option,option{o},years);
        end
    end
end
%% Fig.3
for regress_dist = 0:1
    for outcome_option = 4:5
        Fig3(regress_dist,outcome_option);
    end
end
%% Fig.4
for regress_dist = 0
    for outcome_option = 6%3:5
        Fig4(regress_dist,outcome_option);
    end
end
%% baseline outcome
for regress_dist = 0:1
    for outcome_option = 4:5
        baseline_outcome(regress_dist,outcome_option);
    end
end
%% Int seg
bin_size = 90; % subject to change
ttdays = 1050;
option = {'si_ratio','int_seg'};
years = {'end'};
for outcome_option = 6%3:5
    for regress_dist = 0:1
        for o = 1:2
            int_seg(bin_size,ttdays,regress_dist,outcome_option,option{o},years);
        end
    end
end
%% baseline plv
bin_size = 90; % subject to change
ttdays = 1050; % total days to plot traces
option = {'reorg_dplv','reorg_zplv'};
for outcome_option = 4:5
    for regress_dist = 0:1
        for o = 1:2
            baseline_dplv(bin_size,ttdays,regress_dist,outcome_option,option{o});
        end
    end
end
%% baseline stim
for regress_dist = 0
    for outcome_option = 3:5
        baseline_stim(regress_dist,outcome_option);
    end
end

