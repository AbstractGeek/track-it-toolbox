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

% Defaults
dist_threshold = 0.08;      % 8 cms

% Handle inputs
if (nargin<1)
    error('extractFirstLandedTrajectories needs a rootdir input to load matfiles');
elseif nargin == 1
    force_rewrite = 0;
end

% Load matfile
rootdir = fileparts(matfile);
trackit_traj = load(matfile);

% Out file
trackit_sorted_data_mat = 'trackit_sorted_data.mat';
trackit_dist_data_mat = 'trackit_dist_data.mat';
trackit_short_data_mat = 'trackit_short_data.mat';
log_file_global = fullfile(rootdir, 'trackit_first_landed_extraction.log');

trackit_sorted_data_matfile = fullfile(rootdir, trackit_sorted_data_mat);
trackit_dist_data_matfile = fullfile(rootdir, trackit_dist_data_mat);
trackit_short_data_matfile = fullfile(rootdir, trackit_short_data_mat);

if (exist(trackit_sorted_data_matfile,'file') == 2) && (force_rewrite == 0)
    trackit_sorted_data = load(trackit_sorted_data_matfile);
    trackit_dist_data = load(trackit_dist_data_matfile);
    trackit_short_data = load(trackit_short_data_matfile);
    logid_global = fopen(log_file_global, 'a');
else
    trackit_sorted_data = struct;
    trackit_dist_data = struct;
    trackit_short_data = struct;
    logid_global = fopen(log_file_global, 'w');
end

% Initialize global log file
fprintf(logid_global,'\n#################### %s ####################\n\n',...
    datestr(now));


%% Extract First Landed Trajectory

treatments = fieldnames(trackit_traj);

for i=1:length(treatments)
    
    fprintf('Treatment: %s\n\n', treatments{i});
    fprintf(logid_global, 'Treatment: %s\n\n', treatments{i});
    
    curr_treatment = treatments{i};
    object_num = sum(isletter(curr_treatment))/2;   
    
    days = fieldnames(trackit_traj.(treatments{i}));
    days(ismember(days, {'name'})) = [];    
    
    csvfolder = fullfile(rootdir,trackit_traj.(treatments{i}).name,'First-landings');
    trajfolder = fullfile(rootdir,trackit_traj.(treatments{i}).name,'Traj-plots');
    shortfolder = fullfile(rootdir,trackit_traj.(treatments{i}).name,'Spurned-Data','Short-Data');
    
    if ~isfolder(trajfolder)
        mkdir(trajfolder);
    end    
    if ~isfolder(csvfolder)
        mkdir(csvfolder);
    end
    if ~isfolder(shortfolder)
        mkdir(shortfolder);
    end    
    
    % Begin count of trials
    treatment_trial_total = 0;
    treatment_trial_success = 0;
    
    for j=1:length(days)
        
        fprintf('\tDay: %s\n', days{j});
        
        % Check if log_file exists
        currDir = fullfile(rootdir,trackit_traj.(treatments{i}).name,...
            'Sorted-Data',days{j});
        log_file_local = fullfile(currDir, 'trackit_first_landed_extraction.log');        
        if (exist(log_file_local, 'file')==2) && (force_rewrite==0)
            % Skip extraction (already done)
            continue;
        end
        
        % Append events to log file
        if force_rewrite
            logid_local = fopen(log_file_local, 'w');
        else
            logid_local = fopen(log_file_local, 'a');
        end
        fprintf(logid_local,'\n#################### %s ####################\n\n',...
            datestr(now));
        
        % Initiate out message string.
        out_message = sprintf('\tDay: %s\n', days{j});
        
        % Initiate day count and update treatment count
        day_trial_total = length(fieldnames(trackit_traj.(treatments{i}).(days{j}))) - 1;
        day_trial_success = 0;      
        treatment_trial_total = treatment_trial_total + day_trial_total;
        
        % Perform object extraction
        [cs,trials] = extractObjectLocations(...
            trackit_traj.(treatments{i}).(days{j}),object_num);        
        
        if isempty(trials)            
            fprintf('\t\tObject Extraction: Unsuccessful. Found no object trials\n.');
            out_message = sprintf('%s\t\tObject Extraction: Unsuccessful. Found no object trials\n.', out_message);
            fclose(logid_local);
            continue;
        elseif size(cs,1)> sum(isletter(curr_treatment))/2
            fprintf('\t\tObject Extraction: Unsuccessful. Found too many objects\n.');
            out_message = sprintf('%s\t\tObject Extraction: Unsuccessful. Found too many objects\n.', out_message);
            continue;
        end
        
        % Plot object positions and add to log
        plotObjectLocations(cs, treatments{i},...
            trackit_traj.(treatments{i}).(days{j}), trials,...
            fullfile(trajfolder,days{j}));
        fprintf('\t\tObject Extraction: Successful.\n');        
        out_message = sprintf('%s\t\tObject Extraction: Successful.\n', out_message);        
        
        for k=1:length(trials)
                       
            % Save all the tracked objects
            TrajPlotData(cs, treatments{i},...
                trackit_traj.(treatments{i}).(days{j}).(trials{k}),...
                fullfile(trajfolder, days{j}, sprintf('%s_obj.fig',trials{k})));
            
            [firstLanded, selected_fly,...
                dist_tables, fn_message] = ...
                getFirstLanded(trackit_traj.(treatments{i}).(days{j}).(trials{k}), cs,...
                treatments{i}, fullfile(trajfolder, days{j}, sprintf('%s_cutoff',trials{k})));
            trackit_dist_data.(treatments{i}).(days{j}).(trials{k}) = dist_tables;
            
            if (any(selected_fly) && ...
                    (dist_tables.sorted_table.start_dist(1)>dist_threshold))
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
                fprintf('\t\tTrial %s: Successful. Object %s landed first.\n', trials{k}, selected_fly);                
                out_message = sprintf('%s\t\tTrial %s: Successful. Object %s landed first.\n', out_message, trials{k}, selected_fly);        
                out_message = sprintf('%s%s', out_message, fn_message);                
                
                % Add one to success count
                day_trial_success = day_trial_success + 1;
                
            elseif selected_fly
                % Save to short_data mat
                trackit_short_data.(treatments{i}).(days{j}).(trials{k}) = ...
                    firstLanded;

                % Save data as csv into shortfolder
                writetable(trackit_short_data.(treatments{i}).(days{j}).(trials{k}),...
                    fullfile(shortfolder, sprintf('%s_%s_%s.csv',days{j}(3:end),...
                    trackit_traj.(treatments{i}).name,trials{k})));
                
                % Add to log
                fprintf('\t\tTrial %s: Unsuccessful. First landing criteria unmet. \n', trials{k});
                fprintf('\t\t\tFirst landed fly started from %0.2f cm from the object (minimum start distance - %0.2f cm)\n',...
                    dist_tables.sorted_table.start_dist(1)*100, dist_threshold*100);
                out_message = sprintf('%s\t\tTrial %s: Unsuccessful. First landing criteria unmet.\n',...
                    out_message, trials{k});                
                out_message = sprintf('%s\t\t\tFirst landed fly started from %0.2f cm from the object (minimum start distance - %0.2f cm)\n',...
                    out_message, dist_tables.sorted_table.start_dist(1)*100, dist_threshold*100);
            else
                % Add to log
                fprintf('\t\tTrial %s: Unsuccessful. No Landings on the object.\n', trials{k});
                out_message = sprintf('%s\t\tTrial %s: Unsuccessful. No Landings on the object.\n',...
                    out_message, trials{k});                
            end
            
        end
        
        % Plot trajectories
        if (isfield(trackit_sorted_data.(treatments{i}),days{j}))
            TrajPlotData(cs, treatments{i},...
                trackit_sorted_data.(treatments{i}).(days{j}),...
                fullfile(trajfolder, sprintf('%s_%s.fig',days{j}(3:end),...
                trackit_traj.(treatments{i}).name)));
        end
        
        % Add to treatment count
        treatment_trial_success = treatment_trial_success + day_trial_success;
        
        % Add counts
        fprintf('\tDay %s: Completed. %d of %d trajectories successfully extracted\n\n',...
            days{j}, day_trial_success, day_trial_total);
        out_message = sprintf('%s\tDay %s: Completed. %d of %d trajectories successfully extracted\n\n',...
            out_message, days{j}, day_trial_success, day_trial_total);
        
        % Update log
        fprintf(logid_global,'%s', out_message);
        fprintf(logid_local,'%s', out_message); 
        fprintf(logid_local, 'FIRST LANDED TRAJECTORY EXTRACTION SUCCESSFUL\n\n');
        fclose(logid_local);
    end
       
    % Save name
    trackit_sorted_data.(treatments{i}).name = trackit_traj.(treatments{i}).name;
    trackit_dist_data.(treatments{i}).name = trackit_traj.(treatments{i}).name;
    trackit_short_data.(treatments{i}).name = trackit_traj.(treatments{i}).name;
    
   
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
    
    % Save short data of this treatment in a mat-file
    temp_struct.(treatments{i}) = trackit_short_data.(treatments{i}); 
    save(fullfile(rootdir,trackit_traj.(treatments{i}).name,...
        sprintf('%s_short_data.mat',trackit_traj.(treatments{i}).name)),...
        '-struct', 'temp_struct','-v7.3');
    clearvars temp_struct;
    
    
    % log file
    fprintf('Treatment %s: Completed. %d of %d trajectories successfully extracted\n\n',...
        treatments{i}, treatment_trial_success, treatment_trial_total);
    fprintf(logid_global,'Treatment %s: Completed. %d of %d trajectories successfully extracted\n\n',...
        treatments{i}, treatment_trial_success, treatment_trial_total);
    
end

fprintf(logid_global, 'FIRST LANDED TRAJECTORY EXTRACTION SUCCESSFUL\n\n');
fclose(logid_global);

% Save mat file
save(trackit_sorted_data_matfile,'-struct','trackit_sorted_data','-v7.3');
save(trackit_dist_data_matfile,'-struct','trackit_dist_data','-v7.3');
save(trackit_short_data_matfile,'-struct','trackit_short_data','-v7.3');

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


function [] = plotObjectLocations(cs, treatment, trialData, trials, currDir)
% function [] = plotObjectLocations(cs, treatment, sortedData, filename)
% 
% Modified version of TrajPlotData
% 
% Dinesh Natesan 
% Last updated: 12th Oct 2017

trial_names = fieldnames(trialData);
objposname = trial_names{~ismember(trial_names,trials)};
sortedData = trialData.(objposname);
filename = fullfile(currDir, sprintf('%s_detected.fig',objposname));

objects = fieldnames(sortedData);

% Get colors
colors = lines(length(objects));

% Plot all data
for k = 1:length(objects)    
    plot3(sortedData.(objects{k}).X, sortedData.(objects{k}).Y, sortedData.(objects{k}).Z,...
        '.', 'Markersize', 25, 'Color', colors(k,:), 'DisplayName', objects{k});
    hold on;
end
grid on;

legend show;

% Decorate plot
decoratePlot(cs, treatment)

% Save data
savefig(filename);

% Clear figure
clf(gcf);


function [] = decoratePlot(cs, treatment)
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
obj_det.Shape(low_contrast) = {[xc,yc,-zc.*cylinder_height]};    % Low contrast cylinder
obj_det.Shape(high_contrast) = {[xs,ys,zs].*sphere_radius};   % High contrast sphere

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
    shape = obj_det.Shape{i};
    a= shape(:,1:21)+cs(i,1);
    b= shape(:,22:42)+cs(i,2);
    c= shape(:,43:63)+cs(i,3);
    surf(a,b,c,obj_det.Color{i},'EdgeColor','none',...
        'FaceAlpha',0.2);
end
    
% Make axis equal and square
axis tight;
axis equal;
axis vis3d;

end

end