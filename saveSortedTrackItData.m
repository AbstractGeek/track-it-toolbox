function [] = saveSortedTrackItData(sorted_matfile, save_folder)
% function [] = saveSortedTrackitData()
% 
% 
% Dinesh Natesan 
% Last updated: 12th Oct 2017

% Load data etc.
data = load(sorted_matfile);
sorted_data = struct;
out_matfile = fullfile(save_folder, 'trackit_sorted_data_backward_compatible.mat');
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
            curr_table = data.(treatments{i}).(days{j}).(flies{k})(:,2:end);
            
            % Obtain centers 
            cs = curr_table{1,4:end};
            bs = cs + repmat([0,0,-0.01],1,size(cs,2)/3);
            
            % Obtain new matrix
            cs_bs = nan(size(curr_table,1),(size(curr_table,2)-3)*2);
            cs = reshape(cs, 3, []);
            bs = reshape(bs, 3, []);
            
            % Obtain cs headers (handles multiple cs)
            xyz_tag = repmat({'_X', '_Y', '_Z'},size(cs,2),1);
            cs_strings = repmat(arrayfun(@(x) sprintf('cs%02d',x),...
                (1:size(cs,2)),'UniformOutput', false)',1,3);
            cs_headers = reshape(strcat(reshape(cs_strings,3*size(cs,2),1),...
                reshape(xyz_tag,3*size(cs,2),1)),size(cs,2),3)';
            
            % Obtain bs headers (handles multiple cs)
            xyz_tag = repmat({'_X', '_Y', '_Z'},size(bs,2),1);
            bs_strings = repmat(arrayfun(@(x) sprintf('bs%02d',x),...
                (1:size(bs,2)),'UniformOutput', false)',1,3);
            bs_headers = reshape(strcat(reshape(bs_strings,3*size(bs,2),1),...
                reshape(xyz_tag,3*size(bs,2),1)),size(bs,2),3)';
            
            % new center, base matrix + variable names
            headers = reshape([cs_headers;bs_headers], 1, []);
            cs_bs(1,:) = reshape([cs;bs], 1, []);           
            % create new table
            new_table = [curr_table(:,1:3), array2table(cs_bs,...
                'VariableNames', headers)];            
            
            % Convert table to csv after removing the time
            writetable(new_table,...
                fullfile(treatmentfolder,...
                sprintf('%s_%s_%s.csv',days{j},treatments{i},flies{k})),...
                'WriteVariableNames', false);            
            
            % add table to sorted data struct
            sorted_data.(treatments{i}).(days{j}).(flies{k}) = new_table;
        end
        
    end    

end

% save structure
save(out_matfile,'-struct','sorted_data','-v7.3');

end
