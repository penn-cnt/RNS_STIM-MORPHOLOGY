%% soz_assignments.m
% This script extracts:
% 1) 3 year outcome, including raw and binary group for good/poor responder; 
% 2) SOZ region assignment for:
%       - SOZ electrodes
%       - RNS electrodes
%       - bipolar referenced RNS channels
% 3) RNS channel coords and distances,
% 4) Connectivity type between RNS channel pairs:
%    1 = Inside SOZ, 2 = Between SOZ, 3 = SOZ to normal, 4 = normal to normal
% 5) Depth location: 1 = mesial temporal lobe, 2 = neocortical
% Script also generates 3D plot of SOZ/RNS electrodes and SOZ centroids
%% settings
clear; close all; 
paths;
patient_info = struct2table(load(which('patients_Penn.mat')).patients_Penn);
ptList = {rns_config.patients.ID};
radius = 6; % in mm
plot_figure = 0;
save_figure = 0;
[x,y,z] = sphere;
%%
for pt = 1:length(ptList)
    % Read Patient Data
    ptID = ptList{pt};
    pidx = strcmp(ptID,patient_info.ID);
    % Extracts outcome data
    outcome = getOutcomeFromTimepoint([patient_info{pidx,"implantDate"}],patient_info{pidx,"outcome"},[1,2,3]);
    outcome_group = (outcome > 50) + 1;
    outcome_group(isnan(outcome)) = nan;
    localization(pt).ptID = ptID;
    localization(pt).meets_criteria = true;
    localization(pt).outcome = outcome;
    localization(pt).outcome_group = outcome_group;
    % Extract Coordinates 
    all_electrodes = patient_info{pidx,"IEEG_coords"}{1};
    rns_electrodes = patient_info{pidx,"RNS_coords"}{1};
    rns_connected = patient_info{pidx,"RNS_connected"}{1};
    rns_channels = patient_info{pidx,"RNS_channels"}{1};
    rns_coords = rns_electrodes(ismember(upper(rns_channels),upper(rns_connected)),:);
    soz_coords = patient_info{pidx,"SOZ_coords"}{1};
    soz_centroids = patient_info{pidx,"SOZ_centroid"}{1};
    n_soz = size(soz_centroids,1);

    if isempty(rns_coords) || isempty(soz_coords)
        localization(pt).meets_criteria=false;
        continue
    end

    % Assigning SOZ electrodes to centroids
    soz_distance = pdist2(soz_coords,soz_centroids);
    [~,soz_labels] = min(soz_distance,[],2);
    
    % if RNS and SOZ electrodes within distance threshold
    rns_soz_distance = pdist2(rns_coords,soz_coords);
    is_soz = rns_soz_distance < radius;
    
    % assign RNS electrodes to SOZ electrodes
    rns_labels = zeros(size(rns_coords,1),1);
    for r = 1:size(is_soz,1)
        all_assignments = [];
        for c = 1:size(is_soz,2)
            if is_soz(r,c)
                all_assignments = [all_assignments, soz_labels(c)];
            end
        end
        if isempty(all_assignments)
            rns_labels(r) = 0;
        else
            rns_labels(r) = mode(all_assignments);
        end
    end
    localization(pt).soz_labels = soz_labels;
    localization(pt).rns_labels = rns_labels;

    % Combining contacts to channels
    ch_labels = zeros(4,1);
    for i = 1:2:size(rns_coords,1)
        sub = rns_labels(i:i+1);
        label = sub(logical(sub));
        if isempty(label)
            label = 0;
        elseif length(label) > 1
            label = label(1);
        end
        ch_labels(ceil(i/2)) = label;
    end
    localization(pt).ch_labels = ch_labels;

    % median coordinates and distances
    ch_coords = zeros(size(ch_labels,1),3);
    for i = 1:2:size(rns_coords,1)
        ch_coords(ceil(i/2),:) = mean(rns_coords(i:i+1,:),1,'omitnan');
    end
    ch_distances = pdist2(ch_coords,ch_coords);
    localization(pt).ch_coords = ch_coords;
    localization(pt).ch_distances = ch_distances;

    % Assigning connections labels
    % 1 - inside soz, 2 - between soz, 3 - soz to normal, 4 - normal to normal
    conns =  [1,2; % same
        3,4; % same
        1,3; % diff
        1,4; % diff
        2,3; % diff
        2,4]; % diff
    con_labels = zeros(length(conns),1);
    for con = 1:size(conns,1)
        l1 = ch_labels(conns(con,1));
        l2 = ch_labels(conns(con,2));
        if l1 == 0 && l2 == 0
            con_labels(con) = 4;
        elseif l1 == l2
            con_labels(con) = 1;
        elseif l1 > 0 && l2 > 0
            con_labels(con) = 2;
        else
            con_labels(con) = 3;
        end
    end
    localization(pt).con_labels = con_labels;

    if all(~con_labels)
        localization(pt).meets_criteria = false; 
    end

    % Plotting
    if plot_figure
        % Setting colors for each cluster
        col_vals  = brewermap(length(unique(soz_labels)),'set2');
        cols = col_vals(soz_labels,:);
        figure(pt)
        hold on

        % Plotting electrodes and rns electrodes
        scatter3(all_electrodes(:,1),all_electrodes(:,2),all_electrodes(:,3),'blue');
        scatter3(rns_coords(:,1),rns_coords(:,2),rns_coords(:,3),50,'r','filled')

        % Plot SOZ coordinates
        for i_cord = 1:size(soz_coords,1)
            coords = soz_coords(i_cord,:);
            xs_coord = (x*radius)+coords(1);
            ys_coord = (y*radius)+coords(2);
            zs_coord = (z*radius)+coords(3);
            surf(xs_coord,ys_coord,zs_coord,'FaceColor',cols(i_cord,:),'FaceAlpha',.5);
        end
        if save_figure
            if ~exist(fullfile(datapath,ptID),'dir')
                mkdir(datapath,ptID)
            end
            saveas(gcf,fullfile(datapath,ptID,'Electrodes.fig'))
            close all
        end
    end
    % Extracts depth type
    lead_location = patient_info{pidx,"leadLocations"};
    if strcmp(lead_location,"M")
        localization(pt).depth = 1;
    elseif strcmp(lead_location,"N")
        localization(pt).depth = 2;
    end
end
% localization2 = localization;
save(fullfile(datapath,"localization.mat"),'localization')

function [L,C] = kmeansplus(X,k)
%KMEANS Cluster multivariate data using the k-means++ algorithm.
%   [L,C] = kmeans(X,k) produces a 1-by-size(X,2) vector L with one class
%   label per column in X and a size(X,1)-by-k matrix C containing the
%   centers corresponding to each class.
%   Version: 2013-02-08
%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%
%   References:
%   [1] J. B. MacQueen, "Some Methods for Classification and Analysis of
%       MultiVariate Observations", in Proc. of the fifth Berkeley
%       Symposium on Mathematical Statistics and Probability, L. M. L. Cam
%       and J. Neyman, eds., vol. 1, UC Press, 1967, pp. 281-297.
%   [2] D. Arthur and S. Vassilvitskii, "k-means++: The Advantages of
%       Careful Seeding", Technical Report 2006-13, Stanford InfoLab, 2006.
L = [];
L1 = 0;
while length(unique(L)) ~= k

    % The k-means++ initialization.
    C = X(:,1+round(rand*(size(X,2)-1)));
    L = ones(1,size(X,2));
    for i = 2:k
        D = X-C(:,L);
        D = cumsum(sqrt(dot(D,D,1)));
        if D(end) == 0, C(:,i:k) = X(:,ones(1,k-i+1)); return; end
        C(:,i) = X(:,find(rand < D/D(end),1));
        [~,L] = max(bsxfun(@minus,2*real(C'*X),dot(C,C,1).'));
    end

    % The k-means algorithm.
    while any(L ~= L1)
        L1 = L;
        for i = 1:k, l = L==i; C(:,i) = sum(X(:,l),2)/sum(l); end
        [~,L] = max(bsxfun(@minus,2*real(C'*X),dot(C,C,1).'),[],1);
    end

end
end

function [IDX,C,K]=best_kmeans(X)
dim=size(X);
% default number of test to get minimun under differnent random centriods
test_num=10;
distortion=zeros(dim(1),1);
for k_temp=1:dim(1)
    [~,~,sumd]=kmeans(X,k_temp,'emptyaction','drop');
    destortion_temp=sum(sumd);
    % try differnet tests to find minimun disortion under k_temp clusters
    for test_count=2:test_num
        [~,~,sumd]=kmeans(X,k_temp,'emptyaction','drop');
        destortion_temp=min(destortion_temp,sum(sumd));
    end
    distortion(k_temp,1)=destortion_temp;
end
variance=distortion(1:end-1)-distortion(2:end);
distortion_percent=cumsum(variance)/(distortion(1)-distortion(end));
% plot(distortion_percent,'b*--');
[r,~]=find(distortion_percent>0.9);
if size(r,1) < 2
    K = r(1)+1;
else
    K=r(1,1)+1;
end
[IDX,C]=kmeansplus(X',K);
end