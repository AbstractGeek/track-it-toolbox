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