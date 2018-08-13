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
land_dist_threshold = 0.010;     % if less than 6 mm, the fly is considered to have landed on the object

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