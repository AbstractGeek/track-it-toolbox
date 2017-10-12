function [] = TrajPlotTreatment(sorted_matfile)
% function [] = TrajPlotTreatment(sorted_matfile)
% 
% Plot the trajectory data for each treatment in the sorted_matfile. 
% Combines experiment days in each treatment, aligns approporiately, and
% plots the data. Saves it as figures with the name of the trajectory.
% 
% Inputs:
%   sorted_matfile: 
%       a mat file containing the trajectories of flies that land first
%       on the odor/visual object (sorted based on treatments and 
%       experiment days).
% 
% Output: 
%   Saved 3d trajectory figure per treatment.
% 
% Dinesh Natesan 
% Last updated: 12th Oct 2017

% Load data etc.
data = load(sorted_matfile);
rootdir = fileparts(sorted_matfile);
treatments = fieldnames(data);

% run through each treatment
for i=1:length(treatments)
    
    treatmentData = struct;
    days = fieldnames(data.(treatments{i}));   
    days(ismember(days,{'name'})) = [];
    
    for j=1:length(days)
        % run through each day                 
        curr_day = days{j};
        flies = fieldnames(data.(treatments{i}).(curr_day));
        
        if (j==1)
            % obtain center
            cs = reshape(data.(treatments{i}).(days{1}).(flies{1}){1,5:end}, 3, [])';
            % Make it centered to the first object
            cs = cs - repmat(cs(1,:),size(cs,1),1);            
        end

        for k=1:length(flies)            
            % add fly data to treatmentData after subtracting centers
            treatmentData.(sprintf('%s_%s', curr_day, flies{k})).X = ...
                data.(treatments{i}).(days{j}).(flies{k}).X - ...
                data.(treatments{i}).(days{j}).(flies{k}).cs01_X(1);
            
            treatmentData.(sprintf('%s_%s', curr_day, flies{k})).Y = ...
                data.(treatments{i}).(days{j}).(flies{k}).Y - ...
                data.(treatments{i}).(days{j}).(flies{k}).cs01_Y(1);

            treatmentData.(sprintf('%s_%s', curr_day, flies{k})).Z = ...
                data.(treatments{i}).(days{j}).(flies{k}).Z - ...
                data.(treatments{i}).(days{j}).(flies{k}).cs01_Z(1);            
        end
    
    end
    
    % Plot data for whole of the treatment
    TrajPlotData(cs,treatments{i},treatmentData,...
        fullfile(rootdir,sprintf('%s.fig',treatments{i})));
    
    
end


end