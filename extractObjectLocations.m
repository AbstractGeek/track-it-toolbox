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
if (length(objects) > 1)
    cs = [0, length(objects)];
    trials = [];
    return;
end

meanX = nanmean(trialData.(objtrials{1}).(objects{1}).X);
meanY = nanmean(trialData.(objtrials{1}).(objects{1}).Y);
meanZ = nanmean(trialData.(objtrials{1}).(objects{1}).Z);

cs = [meanX, meanY, meanZ];
trials(ismember(trials, objposition_names)) = [];

end