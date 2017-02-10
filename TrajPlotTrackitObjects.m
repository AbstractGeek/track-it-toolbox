function [] = TrajPlotTrackitObjects(sortedData, filename)


objects = fieldnames(sortedData);

% Initialize
plot3(0,0,0,'HandleVisibility','off');
hold on;
grid on;

% Get colors
colors = lines(length(objects));

% Plot all data
for i = 1:length(objects)    
    plot3(sortedData.(objects{i}).X, sortedData.(objects{i}).Y, sortedData.(objects{i}).Z,...
        'Color', colors(i,:), 'DisplayName', objects{i})    
end

legend show;
axis square;
axis equal;
axis vis3d;

% Save data
savefig(filename);

% Clear figure
clf(gcf);
end