function [fda_pv,sig_days] = FDA_test(data,outcome)
rng('default')
if ~ismember(1,outcome) || ~ismember(2,outcome)
    fda_pv = nan;
    sig_days = [];
else
[SR_mean,NR_mean] = extract_mean(data,outcome);
true_diff = mean(SR_mean - NR_mean,'omitnan');
null_diff = [];
n_person = length(outcome);
for i = 1:10000
    inds = randperm(n_person);
    outcome_rnd = outcome(inds);
    [SR_mean_rnd,NR_mean_rnd] = extract_mean(data,outcome_rnd);
    null_diff = [null_diff, mean(SR_mean_rnd - NR_mean_rnd,'omitnan')];
end
if true_diff > 0
    fda_pv = mean(null_diff > true_diff);
else
    fda_pv = mean(null_diff < true_diff);
end

sig_days = [];
n_day = size(data{1},1);
SR_data = data(outcome == 2);
NR_data = data(outcome == 1);
for i = 1:n_day
    try
        pv = ranksum(cell2mat(cellfun(@(x) x(i,:),SR_data,'UniformOutput',false)), ...
                    cell2mat(cellfun(@(x) x(i,:),NR_data,'UniformOutput',false)));
        if pv < 0.05
            sig_days = [sig_days;[i,pv]];
        end
    catch
        pv = nan;
        sig_days = [];
    end
end
end
end
function [SR_mean,NR_mean] = extract_mean(data,outcome)
traces = data(outcome == 1);
NR_mean = mean(horzcat(traces{:})',1,'omitnan');
traces = data(outcome == 2);
SR_mean = mean(horzcat(traces{:})',1,'omitnan');
end