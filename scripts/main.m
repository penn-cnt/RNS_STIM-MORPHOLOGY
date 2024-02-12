%% settings
paths;
%% preprocesssing
soz_assignments;
event_preprocess;
network_trajectory;
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

%% variance part
outcomes = {'plv_var','dplv_var','plv_slope_var','dplv_slope_var'};
for i = 1:length(outcomes)
    var_Fig2(90,1050,0,4,outcomes{i},{'end'});
end
%% variance part
outcomes = {'dplv_var_acrossbin','zplv_var_acrossbin','dplv_slope_var_acrossbin','zplv_slope_var_acrossbin'};
for i = 1:length(outcomes)
    var_Barplot(90,1050,0,4,outcomes{i},1,{'end'});
end