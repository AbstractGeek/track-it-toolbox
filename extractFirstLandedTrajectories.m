function [trackit_sorted_data_matfile] = extractFirstLandedTrajectories(matfile, force_rewrite)
% function [trackit_sorted_data_matfile] = extractFirstLandedTrajectories(matfile, force_rewrite)
%
% Sift through all the trajectories captured by the trackit system and only
% extract the ones where the flies land on the object. Among the flies that
% land, pick only the first fly.
% 
% Inputs:
%   matfile: 
%       File location of the matfile containing all extracted trajectories 
%       of experiments (output of extractTrajectories function).
%   force_rewrite: 
%       Flag to rewrite the data even if the analysis has been done before.
% 
% Outputs:
%   trackit_sorted_data_matfile: 
%       a mat file containing the trajectories of flies that land first
%       on the odor/visual object (sorted based on treatments and 
%       experiment days).
% 
% Dinesh Natesan 
% Last modified: 12th Oct 2017

if (nargin<1)
    error('extractFirstLandedTrajectories needs a rootdir input to load matfiles');
elseif nargin == 1
    force_rewrite = 0;
end

% Out file
trackit_sorted_data_mat = 'trackit_sorted_data.mat';
trackit_dist_data_mat = 'trackit_dist_data.mat';

rootdir = fileparts(matfile);
trackit_traj = load(matfile);

trackit_sorted_data_matfile = fullfile(rootdir, trackit_sorted_data_mat);
trackit_dist_data_matfile = fullfile(rootdir, trackit_dist_data_mat);

if (exist(trackit_sorted_data_matfile,'file') == 2) && (force_rewrite == 0)
    trackit_sorted_data = load(trackit_sorted_data_matfile);
    trackit_dist_data = load(trackit_dist_data_matfile);
else
    trackit_sorted_data = struct;
    trackit_dist_data = struct;
end

%% Extract First Landed Trajectory

treatments = fieldnames(trackit_traj);

for i=1:length(treatments)
    
    curr_treatment = treatments{i};
    object_num = sum(isletter(curr_treatment))/2;   
    
    days = fieldnames(trackit_traj.(treatments{i}));
    days(ismember(days, {'name'})) = [];    
    
    csvfolder = fullfile(rootdir,trackit_traj.(treatments{i}).name);
    trajfolder = fullfile(rootdir,trackit_traj.(treatments{i}).name,'Traj-plots (Daywise)');
    
    if ~isdir(trajfolder)
        mkdir(trajfolder);
    end
    
    for j=1:length(days)
        
        currDir = fullfile(rootdir,trackit_traj.(treatments{i}).name,...
            'Sorted-Data',days{j});
        
        % Check if log_file exists
        log_file = fullfile(currDir, 'trackit_first_landed_extraction.log');
        
        if (exist(log_file, 'file')==2) && (force_rewrite==0)
            % Skip extraction (already done)
            continue;
        end
        
        % Append events to log file
        logid = fopen(log_file, 'a');
        fprintf(logid,'\n#################### %s ####################\n\n',...
            datestr(now));
        
        % Perform object extraction
        [cs,trials] = extractObjectLocations(...
            trackit_traj.(treatments{i}).(days{j}),object_num);        
        
        if isempty(trials)
            fprintf('Found no object trials in day %s of treatment %s. Skipping day\n',...
                days{j}(3:end), treatments{i});
            fprintf(logid, 'Found no object trials in day %s of treatment %s. Skipping day\n',...
                cs, days{j}(3:end), treatments{i});
            fclose(logid);
            continue;
        elseif size(cs,1)> sum(isletter(curr_treatment))/2
            fprintf('Found too many visual/odor objects in day %s of treatment %s. Skipping day\n',...
                days{j}(3:end), treatments{i});
            fprintf(logid,'Found too many visual/odor objects in day %s of treatment %s. Skipping day\n',...
                days{j}(3:end), treatments{i});
            continue;
        end
        
        % Add to log
        fprintf(logid,'Object Extraction successful!. Beginning first landed trajectory extraction\n\n');
        
        for k=1:length(trials)
            
            [firstLanded, selected_fly,...
                trackit_dist_data.(treatments{i}).(days{j}).(trials{k})] = ...
                getFirstLanded(trackit_traj.(treatments{i}).(days{j}).(trials{k}), cs);
            
            if selected_fly                
                % Save the xyz details
                trackit_sorted_data.(treatments{i}).(days{j}).(trials{k}) = ...
                    firstLanded;
                
                % Save data as csv and copy to main folder
                writetable(trackit_sorted_data.(treatments{i}).(days{j}).(trials{k}),...
                    fullfile(currDir, sprintf('%s_%s_%s.csv',days{j}(3:end),...
                    trackit_traj.(treatments{i}).name,trials{k})));
                copyfile(fullfile(currDir, sprintf('%s_%s_%s.csv',days{j}(3:end),...
                    trackit_traj.(treatments{i}).name,trials{k})), csvfolder);
                
                % Add to log
                fprintf(logid, '%s: %s selected as the first landed fly.\nTrajectory successfully saved as csv\n\n',...
                    sprintf('%s_%s_%s',days{j}(3:end),...
                    trackit_traj.(treatments{i}).name,trials{k}), selected_fly);                
            else
                % Add to log
                fprintf(logid, '%s: No object selected as the first landed fly as none of the qualified all criteria\n\n',...
                    sprintf('%s_%s_%s',days{j}(3:end),...
                    trackit_traj.(treatments{i}).name,trials{k}));
            end
            
        end
        
        % Plot trajectories
        if (isfield(trackit_sorted_data.(treatments{i}),days{j}))
            TrajPlotData(cs, treatments{i},...
                trackit_sorted_data.(treatments{i}).(days{j}),...
                fullfile(trajfolder, sprintf('%s_%s.fig',days{j}(3:end),...
                trackit_traj.(treatments{i}).name)));
        end
        
        fprintf(logid, 'FIRST LANDED TRAJECTORY EXTRACTION SUCCESSFUL\n\n');
        fclose(logid);
    end
    
    % Save name
    trackit_sorted_data.(treatments{i}).name = trackit_traj.(treatments{i}).name;
    trackit_dist_data.(treatments{i}).name = trackit_traj.(treatments{i}).name;
    
    % Save first landed trajectories of this treatment as a mat-file
    temp_struct.(treatments{i}) = trackit_sorted_data.(treatments{i}); 
    save(fullfile(rootdir,trackit_traj.(treatments{i}).name,...
        sprintf('%s_sorted_data.mat',trackit_traj.(treatments{i}).name)),...
        '-struct', 'temp_struct','-v7.3');
    clearvars temp_struct;
    
    % Save distance tables of this treatment in a mat-file
    temp_struct.(treatments{i}) = trackit_dist_data.(treatments{i}); 
    save(fullfile(rootdir,trackit_traj.(treatments{i}).name,...
        sprintf('%s_dist_data.mat',trackit_traj.(treatments{i}).name)),...
        '-struct', 'temp_struct','-v7.3');
    clearvars temp_struct;
    
end

% Save mat file
save(trackit_sorted_data_matfile,'-struct','trackit_sorted_data','-v7.3');
save(trackit_dist_data_matfile,'-struct','trackit_dist_data','-v7.3');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [firstLanded, selected_fly, dist_tables] = ...
    getFirstLanded(trialData, cs)
% [firstLanded, selected_obj] = getFirstLanded(trialData, cs)
%
% Inputs: Trial Data: Trial Data structure with multiple trajectors
%         cs: Object center (X,Y,Z)
%
% Outputs: firstLanded: Table with trajectory of the first landed fly along
%                       with the object centers used for selection.
%          selected_obj: The selected fly based on first landed criterion.
%
% Dinesh Natesan
% 10th Feb, 2017

objects = fieldnames(trialData);

start_dist = nan(size(objects,1),size(cs,1));
stop_dist = nan(size(objects,1),size(cs,1));
t_stop_dist = nan(size(objects,1),size(cs,1));  % Thresholded stop distance
land_time = nan(size(objects));

dist_threshold = 0.08;      % 8 cms
dist_roundoff = 2;          % 1 cm round off is good enough
minimum_dist_threshold = 0.01;   % if greater than 1 cm, do not consider.
land_dist_threshold = 0.006;     % if less than 6 mm, the fly is considered to have landed on the object

% find minimum distance for each object
for i = 1:length(objects)    
    start_dist(i,:) = (sqrt(sum([trialData.(objects{i}).X(1) - cs(:,1),...
        trialData.(objects{i}).Y(1) - cs(:,2),...
        trialData.(objects{i}).Z(1) - cs(:,3)].^2,2)))';
    stop_dist(i,:) = sqrt(sum([trialData.(objects{i}).X(end) - cs(:,1),...
        trialData.(objects{i}).Y(end) - cs(:,2),...
        trialData.(objects{i}).Z(end) - cs(:,3)].^2,2));
    t_stop_dist(i,:) = floorn(stop_dist(i,:),dist_roundoff)';
    land_time(i) = trialData.(objects{i}).time(end);
end

% Make a table of objects
obj_table = table(objects, start_dist, stop_dist, t_stop_dist, land_time, ...
    'VariableNames',{'objects','start_dist', 'stop_dist',...
                     't_stop_dist', 'land_time'});

% Find the minimum for each object and use that to create a sorted table
[min_dist, min_ind] = min(t_stop_dist,[],2);
min_ind = sub2ind(size(t_stop_dist), (1:size(t_stop_dist,1))', min_ind);

% Find landed object
[~,landed_obj_num] = ind2sub(size(t_stop_dist),min_ind);

sorted_table = sortrows(table(objects, landed_obj_num,...
    start_dist(min_ind), stop_dist(min_ind), min_dist, land_time,...
    'VariableNames',{'objects','landed_obj_num','start_dist','stop_dist',...
                     't_stop_dist', 'land_time'}),...
    {'t_stop_dist','land_time'},{'ascend'});

% Check if it meets the minimum distance criteria
selected_fly = false;
if(boolean(size(sorted_table,1)) && ...
    (sorted_table.stop_dist(1)<=minimum_dist_threshold) &&  ...
    (sorted_table.start_dist(1)>=dist_threshold))
    selected_fly = sorted_table.objects{1};
    landed_object = sorted_table.landed_obj_num(1);
end

% Add distance table struct and return that too
dist_tables = struct('obj_table', obj_table, 'sorted_table', sorted_table);


if selected_fly
    % Create table and reset time
    firstLanded =  trialData.(selected_fly);
    firstLanded.time = firstLanded.time - firstLanded.time(1);

    % Obtain cs headers (handles multiple cs)
    xyz_tag = repmat({'_X', '_Y', '_Z'},size(cs,1),1);
    cs_strings = repmat(arrayfun(@(x) sprintf('cs%02d',x),...
        (1:size(cs,1)),'UniformOutput', false)',1,3);
    cs_headers = reshape(strcat(reshape(cs_strings,3*size(cs,1),1),...
        reshape(xyz_tag,3*size(cs,1),1)),size(cs,1),3);
    
    % reshape and format cs and its headers
    reshaped_cs_headers = reshape(cs_headers',size(cs_headers,1)*size(cs_headers,2),1)';
    reshaped_cs = reshape(cs',size(cs,1)*size(cs,2),1)';
    
    % Find trailing points based on landing distance cutoff
    source_dist = (sqrt(sum([firstLanded.X - cs(landed_object,1),...
        firstLanded.Y - cs(landed_object,2),...
        firstLanded.Z - cs(landed_object,3)].^2,2)));
    t_source_dist = (source_dist < land_dist_threshold); 
    transitions = find(diff(t_source_dist) == 1);        
    
    if isempty(transitions)    
        selected_fly = false;
        firstLanded = [];
        
        % save object transitions table
        transition_table = table(source_dist, t_source_dist,...            
            'VariableNames', {'source_dist','t_source_dist'});
        dist_tables.('transition_table') = transition_table;

    else
        if length(transitions) > 1        
            warning('More than one transitions into the landing distance cutoff sphere.');        
        end
        
        % Object transitions table
        transition_table = table(source_dist, t_source_dist,...
            [ones(transitions(1)+1,1);...
            zeros(length(source_dist) - (transitions(1)+1),1)],...
            'VariableNames', {'source_dist','t_source_dist','selected_inds'});
        % save the transition table too
        dist_tables.('transition_table') = transition_table;
        
        % Remove trailing points
        firstLanded = firstLanded(1:transitions(1)+1,:);
        
        % Create a temp table with cs
        table_rows = length(firstLanded.time);
        temp_table = array2table(nan(table_rows,3*size(cs,1)),'VariableNames',reshaped_cs_headers);
        temp_table{1,:} =  reshaped_cs;
        
        % add it to the main table
        firstLanded = [firstLanded, temp_table];
    end
    
else
    firstLanded = [];
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Xr] = floorn(X,n)
% function [X] = floorn(X,n)
%
% Floors the input at the nth decimal place (similar to round at nth
% decimal point)
%
% Dinesh Natesan
% 10th Feb 2015

Xr = floor(X*10^n)*10^-n;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [out] = iff(cond,a,b)
% function iff(cond,a,b)
% A custom written function that mimic the traditional C+ conditional
% expression: out = cond?true:false
%
% Dinesh Natesam, 6th Mar 2014

if cond
    out = a;
else
    out = b;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [cs, trials] = extractObjectLocations(trialData, object_num)
% function [cs, trials] = extractObjectLocations(trialData)
% 
% A very simple function that extracts object xyz points and
% returns it as cs. It uses object number to cluster points using k-means
% clustering. The object number is extracted from the treatment name and
% has to be an integer. If it is not, a default object number based
% clustering would be performed which will be checked for validity later
% on. 
% 
% Dinesh Natesan
% 10 Feb, 2017

% Defaults
objposition_names = {'objpositions'};

trials = fieldnames(trialData);
objtrials = trials(ismember(trials, objposition_names));

if (length(objtrials) ~= 1)
    cs = length(objtrials);
    trials = [];
    return;
end

if (floor(object_num) == object_num)
    % k-means clustering based object extraction
    object_data = cell2mat(cellfun(@(x) table2array(x),...
        struct2cell(trialData.objpositions),'UniformOutput',false));
    [~,cs,sumd,~] = kmeans(object_data(:,2:4),object_num);
    
    if any(round(sumd,3) ~= 0)
       clustering_error = true; 
    else
       clustering_error = false;
    end

end
    
if ((floor(object_num) ~= object_num) || clustering_error)
    % Standard mean based object extraction
    objects = fieldnames(trialData.(objtrials{1}));
    cs = NaN(length(objects), 3);
    
    % Find the mean X, Y, Z for each object
    for i=1:length(objects)
        cs(i,1) = nanmean(trialData.(objtrials{1}).(objects{i}).X);
        cs(i,2) = nanmean(trialData.(objtrials{1}).(objects{i}).Y);
        cs(i,3) = nanmean(trialData.(objtrials{1}).(objects{i}).Z);
    end
    
end

% Sort cs based on X
cs = sortrows(cs);
trials(ismember(trials, objposition_names)) = [];

end