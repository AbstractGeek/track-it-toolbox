function [cs, trials] = extractObjectLocations(trialData)
% function [cs, trials] = extractObjectLocations(trialData)
% 
% A very simple function (for now) that extracts object xyz points and
% returns it as cs. Currently just detects only one object and it simply
% takes a mean of all the points. Might get complicated as the treatment
% types change.
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

objects = fieldnames(trialData.(objtrials{1}));
cs = NaN(length(objects), 3);

% Find the mean X, Y, Z for each object
for i=1:length(objects)
    cs(i,1) = nanmean(trialData.(objtrials{1}).(objects{i}).X);
    cs(i,2) = nanmean(trialData.(objtrials{1}).(objects{i}).Y);
    cs(i,3) = nanmean(trialData.(objtrials{1}).(objects{i}).Z);    
end

% Sort cs based on X
cs = sortrows(cs);
trials(ismember(trials, objposition_names)) = [];

end