function [] = saveSortedTrackitData(sorted_matfile, save_folder)
% function [] = saveSortedTrackitData()
% 
% 
% Dinesh Natesan 
% Last updated: 12th Oct 2017

% Load data etc.
data = load(sorted_matfile);
treatments = fieldnames(data);

% run through each treatment
for i=1:length(treatments)
    
    days = fieldnames(data.(treatments{i}));   
    days(ismember(days,{'name'})) = [];
    
    treatmentfolder = fullfile(save_folder,treatments{i});
    if ~isdir(treatmentfolder)
        mkdir(treatmentfolder);
    end

    for j=1:length(days)
        % run through each day
        flies = fieldnames(data.(treatments{i}).(days{j}));
        
        for k=1:length(flies)
            % Convert table to csv after removing the time
            writetable(data.(treatments{i}).(days{j}).(flies{k})(:,2:end),...
                fullfile(treatmentfolder,...
                sprintf('%s_%s_%s.csv',days{j},treatments{i},flies{k})));            
        end
        
    end    

end

end