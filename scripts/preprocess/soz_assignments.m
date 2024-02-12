%% soz_assignments.m
% This script extracts:
% 1) 3 year outcome, including raw and binary group for good/poor responder;
% 2) SOZ region assignment for:
%       - SOZ electrodes
%       - RNS electrodes
%       - bipolar referenced RNS channels
% 3) RNS channel coords and distances
% 4) Connectivity type between RNS channel pairs:
%    1 = Inside SOZ, 2 = Between SOZ, 3 = SOZ to normal, 4 = normal to normal
% 5) Depth location: 1 = mesial temporal lobe, 2 = neocortical
% Script also generates 3D plot of SOZ/RNS electrodes and SOZ centroids
%% settings
clear; close all;
paths;
patient_info = struct2table(load(which('patients_Penn.mat')).patients_Penn);
ptList = {rns_config.patients.ID};
radius = 10;
plot_figure = 0;
save_figure = 0;
[x,y,z] = sphere;
fluc_list = {'HUP121','HUP127','HUP136','HUP143','HUP147','HUP156','HUP108'};
%%
for pt = 1:length(ptList)
    % Read Patient Data
    ptID = ptList{pt};
    pidx = strcmp(ptID,patient_info.ID);

    % Extracts outcome data
    outcome_data = patient_info.outcome{pidx,1};
    if isempty(outcome_data) || isempty(outcome_data{1})
        outcome = [];
        outcome_group = [];
    else
        outcome_group = [];
        last_visit_year = ceil(years(outcome_data{1}(end)-patient_info.implantDate(pidx)));
        outcome = getOutcomeFromTimepoint([patient_info{pidx,"implantDate"}],{outcome_data},[1:last_visit_year]);
        outcome_group(1,:) = (outcome > 50) + 1;
        outcome_group(2,:) = (outcome > 50) + 1;
        outcome_group(4,:) = (outcome > 50) + 1;
        outcome_group(3,:) = (outcome > 80) + 1;
        outcome_group(5,:) = (outcome > 80) + 1;
        outcome_group(6,:) = (outcome > 80) + 1;
        if ismember(ptID,fluc_list)
            outcome_group(2,:) = 3*ones(1,length(outcome));
            outcome_group(4,:) = ones(1,length(outcome));
            outcome_group(5,:) = 3*ones(1,length(outcome));
        end
        outcome_group(6,logical(outcome > 50 & outcome <=80)) = 3;
    end
    localization(pt).ptID = ptID;
    localization(pt).meets_criteria = true;
    localization(pt).outcome = outcome;
    localization(pt).outcome_group = outcome_group;
    localization(pt).lead_labels = [1,1,2,2,2,2];
    %% depth
    lead_location = patient_info.leadLocations{pidx,1};
    laterality = patient_info.laterality{pidx,1};
    if strcmp(lead_location,"M")
        localization(pt).depth = 1;
        if ~strcmp(laterality,'B')
            localization(pt).lead_labels = [1,1,1,1,1,1];
        else
            localization(pt).lead_labels = [1,1,2,2,2,2];
        end
    elseif strcmp(lead_location,"B")
        switch ptID
            case 'HUP129'
                localization(pt).depth = 1;
                localization(pt).lead_labels = [1,1,1,1,1,1];
            case 'HUP131'
                localization(pt).depth = 2;
            case 'HUP136'
                localization(pt).depth = 2;
            case 'HUP182'
                localization(pt).depth = 1;% manual fix at last
            case 'HUP199'
                localization(pt).depth = 1;
                localization(pt).lead_labels = [1,1,2,2,2,2];
        end
    elseif strcmp(lead_location,"N")
        localization(pt).depth = 2;
    end
    localization(pt).con_labels = localization(pt).lead_labels;

    % Extract Coordinates
    % RNS
    rns_channels = patient_info.RNS_MNI_channels{pidx,1};

    % Stop if no RNS data available
    if isempty(rns_channels)
        continue
    end

    rns_coords = patient_info.RNS_MNI_coords{pidx,1};
%     rns_locations = patient_info.RNS_locations{pidx,1};
    rns_connected = patient_info.RNS_connected{pidx,1};
    inds = ismember(rns_channels,rns_connected);
    rns_channels = rns_channels(inds,:);
    rns_coords = rns_coords(inds,:);
%     rns_locations = rns_locations(inds,:);

    % all chan info
    all_channels = patient_info.IEEG_channels{pidx,1};
    all_coords = patient_info.IEEG_coords{pidx,1};
    all_locations = patient_info.IEEG_locations{pidx,1};

    % SOZ chan info
    soz_channels = patient_info.SOZ_channels{pidx,1};
    soz_coords = patient_info.SOZ_coords{pidx,1};
    
    %% coords & distances
    if ~isempty(rns_coords)
        % median coordinates and distances
        ch_coords = zeros(4,3);
        for i = 1:2:size(rns_coords,1)
            ch_coords(ceil(i/2),:) = mean(rns_coords(i:i+1,:),1,'omitnan');
        end
        ch_distances = pdist2(ch_coords,ch_coords);
        localization(pt).ch_coords = ch_coords;
        localization(pt).ch_distances = ch_distances;
    end
    %% connectivity type for N
    % If not coords info from both rns and soz available, directly assign
    % inter-lead as between soz and intra-lead as within soz
    if localization(pt).depth == 2
        if isempty(rns_coords) || isempty(soz_coords)
            continue
        else
            % close to soz channels
            soz_labels = [soz_channels{:,2}];
            rns_soz_distance = pdist2(rns_coords,soz_coords);
            is_soz = rns_soz_distance < 2*radius;
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
        end
    end


    % Plotting
    if plot_figure
        % Setting colors for each cluster
        col_vals  = brewermap(length(unique(soz_labels)),'set2');
        cols = col_vals(soz_labels,:);
        figure(pt)
        hold on

        % Plotting electrodes and rns electrodes
        scatter3(all_coords(:,1),all_coords(:,2),all_coords(:,3),'blue');
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
end
%% manual fix for HUP182
ptID = 'HUP182';
pt = strcmp(ptID,{localization.ptID});
% Read Patient Data
pidx = strcmp(ptID,patient_info.ID);
localization(end+1) = localization(pt);
% M instance
localization(pt).depth = 1;
localization(pt).lead_labels = [0,1,0,0,0,0];
% N instance
localization(end).depth = 2;
localization(end).lead_labels = [1,0,0,0,0,0];
% localization2 = localization;
save(fullfile(datapath,"localization.mat"),'localization')
