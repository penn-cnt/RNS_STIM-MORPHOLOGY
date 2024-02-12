for i = 1:length(plasticity)
    plasticity(i).distance = localization(i).ch_distances;
    if ~isempty(plasticity(i).distance)
        connections =  [1,2; 3,4; 1,3; 1,4; 2,3; 2,4];% diff
        dist = nan * zeros(1,length(connections));
        for conn = 1:length(connections)
            dist(conn) = plasticity(i).distance(connections(conn,1),connections(conn,2));
        end
        plasticity(i).distance = dist;
    end
end

plasticity = plasticity(~cellfun('isempty', {plasticity.outcome_group}));
plasticity = plasticity(~cellfun('isempty', {plasticity.depth}));
for i = 1:length(plasticity)
    plasticity(i).outcome_group = plasticity(i).outcome_group(4,:);
end
for i = 1:length(plasticity)
    plasticity(i).intra_dist = plasticity(i).distance(plasticity(i).con_labels == 1);
    plasticity(i).inter_dist = plasticity(i).distance(plasticity(i).con_labels == 2);
end
for depth = 1:2
    figure;
    data = plasticity([plasticity.depth] == depth);
    [dist,outcome] = getYearOutcome(data,'end','outcome_group','intra_dist');
    pr = horzcat(dist{outcome == 1})';
    sr = horzcat(dist{outcome == 2})';
    sr = padarray(sr,[max(size(sr,1),size(pr,1))-size(sr,1),0],nan,'post');
    pr = padarray(pr,[max(size(sr,1),size(pr,1))-size(pr,1),0],nan,'post');
    box_data{1} = pr;
    box_data{2} = sr;
    b = boxplotGroup(box_data,'primaryLabels',repmat({''}, 2, 1), ...
            'Colors',[0,1,0;0,0,1],'GroupType','betweenGroups', ...
            'PlotStyle','traditional','BoxStyle','outline', ...
            'Symbol','o','Widths',0.7);
    ymax = get(gca, 'YLim');
    if any(outcome == 1) && any(outcome == 2)
        p = ranksum(sr,pr);
        text(1,ymax(2)-0.1,['p=',num2str(p,'%.3f')])
    end
end