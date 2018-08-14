function [firstLanded, selected_fly, dist_tables] = ...
    getFirstLanded(trialData, cs, treatment, filename_prefix)
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

% Defaults
dist_threshold = 0.08;      % 8 cms
minimum_dist_threshold = 0.01;   % if greater than 1 cm, do not consider.
object_radius_std = 3;  % Number of standard deviations away from the mean

% Find the radius of the landed object
objectRadius = findObjectRadius(trialData, cs,...
    minimum_dist_threshold, object_radius_std, treatment, filename_prefix);

% Initalize data structures
objects = fieldnames(trialData);
start_dist = nan(size(objects,1),size(cs,1));
stop_dist = nan(size(objects,1),size(cs,1));
object_radius = repmat(objectRadius',size(objects,1),1);
t_stop_dist = nan(size(objects,1),size(cs,1));  % Thresholded stop distance
land_time = nan(size(objects));

% find minimum distance for each object
for i = 1:length(objects)    
    start_dist(i,:) = (sqrt(sum([trialData.(objects{i}).X(1) - cs(:,1),...
        trialData.(objects{i}).Y(1) - cs(:,2),...
        trialData.(objects{i}).Z(1) - cs(:,3)].^2,2)))';
    stop_dist(i,:) = sqrt(sum([trialData.(objects{i}).X(end) - cs(:,1),...
        trialData.(objects{i}).Y(end) - cs(:,2),...
        trialData.(objects{i}).Z(end) - cs(:,3)].^2,2));
    t_stop_dist(i,:) = stop_dist(i,:) > objectRadius';
    land_time(i) = trialData.(objects{i}).time(end);
end

% Make a table of objects
obj_table = table(objects, start_dist, stop_dist, object_radius,...
    t_stop_dist, land_time, ...
    'VariableNames',{'objects','start_dist', 'stop_dist',...
                     'object_radius', 't_stop_dist', 'land_time'});

% Find the minimum for each object and use that to create a sorted table
[min_dist, min_ind] = min(t_stop_dist,[],2);
min_ind = sub2ind(size(t_stop_dist), (1:size(t_stop_dist,1))', min_ind);

% Find landed object
[~,landed_obj_num] = ind2sub(size(t_stop_dist),min_ind);

% Sort the objects based on stop_distance and land_time
sorted_table = sortrows(table(objects, landed_obj_num,...
    start_dist(min_ind), stop_dist(min_ind), object_radius(min_ind),...
    min_dist, land_time,...
    'VariableNames',{'objects','landed_obj_num','start_dist','stop_dist',...
                     'object_radius', 't_stop_dist', 'land_time'}),...
    {'t_stop_dist','land_time'},{'ascend'});

% Add distance table struct and return that too
dist_tables = struct('obj_table', obj_table, 'sorted_table', sorted_table);

% Check if it meets the minimum distance criteria
selected_fly = false;
if(boolean(size(sorted_table,1)) && ...
    (sorted_table.start_dist(1)>=dist_threshold))
    selected_fly = sorted_table.objects{1};
    landed_object = sorted_table.landed_obj_num(1);
elseif (boolean(size(sorted_table,1)))
    fprintf('\t\t\tFound a landing that started inside the 8cm sphere (not considering it)\n');
else
    fprintf('\t\t\tFound no landings\n');
end


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
    t_source_dist = (source_dist < objectRadius(landed_object)); 
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
            fprintf('\t\t\tMore than one transitions into the landing cutoff sphere.\n\t\t\tConsidering the first entry as cutoff.\n');        
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

function objectRadius = findObjectRadius(trialData, cs, tdist, rstd, treatment, filename_prefix)
% function objectRadius = findObjectRadius(trialData, cs, tdist, rstd,
% filename)
% 
% 
% 

% Pool dataset (to check for clusters)
objects = fieldnames(trialData);
pooled_objects = table([],[],[],[],'VariableNames',{'time','X','Y','Z'});
for i=1:length(objects)
    % pool all points
    pooled_objects = [pooled_objects; trialData.(objects{i})];
end
pooled_objects(:,1) = [];

% Get source distance
source_distance = nan(size(pooled_objects,1),size(cs,1));
objectRadius = nan(size(cs,1),1);
idx_all = false(size(source_distance,1),1);
for j=1:size(cs,1)
    source_distance(:,j) = sqrt(sum((pooled_objects.Variables - cs(j,:)).^2,2));
    idx = source_distance(:,j) <= tdist;   % Generally 1 cm
    objectRadius(j) = nanmedian(source_distance(idx,j)) + ...
        rstd * nanstd(source_distance(idx,j));
    idx_all = idx_all | idx;
end

% Fill NaNs with the maximum object size
objectRadius(isnan(objectRadius)) = max(objectRadius);

% Plot the data
plot3(pooled_objects.X(~idx_all),pooled_objects.Y(~idx_all),pooled_objects.Z(~idx_all),'.b');
hold on, plot3(pooled_objects.X(idx_all),pooled_objects.Y(idx_all),pooled_objects.Z(idx_all),'.r');
decoratePlot(cs, treatment);
decoratePlot(cs, treatment, objectRadius);
% Save data
savefig(sprintf('%s_raw.fig',filename_prefix));
% Clear figure
clf(gcf);

% Plot individual data
% Get colors
colors = lines(length(objects));
% Plot all data
for i = 1:length(objects)    
    plot3(trialData.(objects{i}).X, trialData.(objects{i}).Y, trialData.(objects{i}).Z,...
        'Color', colors(i,:), 'DisplayName', objects{i});
    if i==1
       hold on;
       grid on;
    end
end
decoratePlot(cs, treatment);
decoratePlot(cs, treatment, objectRadius);
% Save data
savefig(sprintf('%s_obj.fig',filename_prefix));
% Clear figure
clf(gcf);

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

function [] = decoratePlot(cs, treatment, radius)
%
%

% Defaults
sphere_radius = 0.003;
cylinder_radius = 0.00025;
cylinder_height = 0.005;
redcolor = [0.8980 0 0];
palepinkcolor = [1.0000 0.8118 0.8627];

% Initialize shape
[xs,ys,zs]=sphere;
[xc,yc,zc]=cylinder(cylinder_radius);

% Initialize color
black_sphere = zeros(21,21,3);
red_sphere = reshape(repmat(redcolor,21*21,1),21,21,3);
coral_sphere = reshape(repmat(palepinkcolor,21*21,1),21,21,3);
black_cyl= zeros(2,21,3);
red_cyl = reshape(repmat(redcolor,2*21,1),2,21,3);
palepink_cyl = reshape(repmat(palepinkcolor,2*21,1),2,21,3);

% Assign table for object shape and color
obj_det = table(cell(size(cs,1),1),cell(size(cs,1),1),cell(size(cs,1),1),...
    'VariableNames',{'Shape', 'Color', 'EdgeColor'});

% Cleanup treatment name
treatment = treatment(isletter(treatment)); 
if rem(length(treatment),2) ~= 0
    % Unknown experimental treament
    axis on;    
    axis tight;   
    axis vis3d;
    return;
end
obj_num = length(treatment)/2;

% Assign shape based on treatment
low_contrast = treatment(1:obj_num)=='c';
high_contrast = treatment(1:obj_num)=='V';
obj_det.Shape(low_contrast) = {@(r) [xc*r,yc*r,-zc*cylinder_height]};    % Low contrast cylinder
obj_det.Shape(high_contrast) = {@(r) [xs,ys,zs].*r};   % High contrast sphere

% Create radius variable if it does not exist
if ~exist('radius','var')
    radius = nan(size(cs,1),1);
    radius(low_contrast) = cylinder_radius;
    radius(high_contrast) = sphere_radius;    
end

% Obtain color characteristics of the treatment
no_odor = ismember((1:obj_num),(find(treatment=='n')-obj_num));
low_odor = ismember((1:obj_num),(find(treatment=='L')-obj_num));
high_odor = ismember((1:obj_num),(find(treatment=='H')-obj_num));

% Look for all combinations of vision and odor
% vision + no odor
obj_det.Color(low_contrast & no_odor) = {black_cyl};
obj_det.Color(high_contrast & no_odor) = {black_sphere};
obj_det.EdgeColor(no_odor) = {[0,0,0]};
% vision + low odor
obj_det.Color(low_contrast & low_odor) = {palepink_cyl};
obj_det.Color(high_contrast & low_odor) = {coral_sphere};
obj_det.EdgeColor(low_odor) = {palepinkcolor};
% vision + high odor
obj_det.Color(low_contrast & high_odor) = {red_cyl};
obj_det.Color(high_contrast & high_odor) = {red_sphere};
obj_det.EdgeColor(high_odor) = {redcolor};

% Plot objects
for i = 1:size(cs,1)
    shapefn = obj_det.Shape{i};
    shape = shapefn(radius(i));
    a= shape(:,1:21)+cs(i,1);
    b= shape(:,22:42)+cs(i,2);
    c= shape(:,43:63)+cs(i,3);
    surf(a,b,c,obj_det.Color{i},'EdgeColor',obj_det.EdgeColor{i},...
        'EdgeAlpha',0.4,'FaceAlpha',0.4);
end
    
% Make axis equal and square
axis tight;
axis equal;
axis vis3d;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%