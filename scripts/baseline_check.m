plasticity = load(fullfile(datapath,['plasticity_90.mat'])).plasticity;
plasticity = plasticity(~cellfun('isempty', {plasticity.outcome_group}));
plasticity = plasticity(~cellfun('isempty', {plasticity.depth}));
for i = 1:length(plasticity)
    plasticity(i).outcome_group = plasticity(i).outcome_group(4,:);
end

% compare baseline PLV
p_All = cell(1,2);
for d = 1:2 % depth
    data = plasticity([plasticity.depth]==d);
    [bplv,outcome] = getYearOutcome(data,'end','outcome_group','baseline_plvs');
    bplv_intra = cellfun(@(x,conn) mean(x(conn == 1,:),1,'omitnan'),bplv,{data.con_labels},'UniformOutput',false);
    bplv_inter = cellfun(@(x,conn) mean(x(conn == 2,:),1,'omitnan'),bplv,{data.con_labels},'UniformOutput',false);
    bplv_intra = vertcat(bplv_intra{:});
    bplv_inter = vertcat(bplv_inter{:});
    p_all = nan * zeros(2,4);
    for freq = 1:4
        p = ranksum(bplv_intra(outcome == 1,freq),bplv_intra(outcome == 2,freq));
        p_all(1,freq) = p;
        p = ranksum(bplv_inter(outcome == 1,freq),bplv_inter(outcome == 2,freq));
        p_all(2,freq) = p;
    end
    p_All{d} = p_all;
end

plasticity = load(fullfile(datapath,['plasticity_90_regdist.mat'])).plasticity;
plasticity = plasticity(~cellfun('isempty', {plasticity.outcome_group}));
plasticity = plasticity(~cellfun('isempty', {plasticity.depth}));
for i = 1:length(plasticity)
    plasticity(i).outcome_group = plasticity(i).outcome_group(3,:);
end

for i = 1:length(plasticity)
    [~,idx] = ismember(plasticity(i).ptID,{localization.ptID});
    plasticity(i).dist = localization(idx).ch_distances;
end
% compare distance
conns =  [1,2; % same
    3,4; % same
    1,3; % diff
    1,4; % diff
    2,3; % diff
    2,4]; % diff
p_dist_all = nan * zeros(2,2);
for d = 1:2 % depth
    data = plasticity([plasticity.depth]==d);
    [dist,outcome] = getYearOutcome(data,'end','outcome_group','dist');
    dist_intra = cellfun(@(x,conn) mean([x(1,2),x(3,4)],'omitnan'),dist,{data.con_labels},'UniformOutput',false);
    dist_inter = cellfun(@(x,conn) mean([x(1,3),x(1,4),x(2,3),x(2,4)],'omitnan'),dist,{data.con_labels},'UniformOutput',false);
    dist_intra = vertcat(dist_intra{:});
    dist_inter = vertcat(dist_inter{:});
    p = ranksum(dist_intra(outcome == 1),dist_intra(outcome == 2));
    p_dist_all(d,1) = p;
    p = ranksum(dist_inter(outcome == 1),dist_inter(outcome == 2));
    p_dist_all(d,2) = p;
end

% compare baseline PLV
p_All = cell(1,2);
for d = 1:2 % depth
    data = plasticity([plasticity.depth]==d);
    [bplv,outcome] = getYearOutcome(data,'end','outcome_group','baseline_plvs');
    bplv_intra = cellfun(@(x,conn) mean(x(conn == 1,:),1,'omitnan'),bplv,{data.con_labels},'UniformOutput',false);
    bplv_inter = cellfun(@(x,conn) mean(x(conn == 2,:),1,'omitnan'),bplv,{data.con_labels},'UniformOutput',false);
    bplv_intra = vertcat(bplv_intra{:});
    bplv_inter = vertcat(bplv_inter{:});
    p_all = nan * zeros(2,4);
    for freq = 1:4
        p = ranksum(bplv_intra(outcome == 1,freq),bplv_intra(outcome == 2,freq));
        p_all(1,freq) = p;
        p = ranksum(bplv_inter(outcome == 1,freq),bplv_inter(outcome == 2,freq));
        p_all(2,freq) = p;
    end
    p_All{d} = p_all;
end