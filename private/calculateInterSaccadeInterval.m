function [saccade_interval, turn_interval,...
    saccade_interval_hover_excluded, turn_interval_hover_excluded] = ...
    calculateInterSaccadeInterval(saccades, turns, hover, sam)
%
%
%

false_array = false(size(saccades));

% Handle saccades
if sum(saccades) < 2
    saccade_interval = NaN;
    saccade_interval_hover_excluded = NaN;
else
    time = (0:1:length(saccades)-1)'./sam;
    saccade_instances = time(logical(saccades));
    saccade_interval = diff(saccade_instances);
    % Saccades before first hover
    indd = find(hover == 1, 1, 'first');
    if isempty(indd)
        saccade_instances_hover_excluded = saccade_instances;
        saccade_interval_hover_excluded = saccade_interval;
    else
        saccade_instances_hover_excluded = ...
            time([logical(saccades(1:indd-1));false_array(indd:end)]);
        saccade_interval_hover_excluded = diff(saccade_instances_hover_excluded);
    end
end

% Handle turns
if sum(turns) < 2
    turn_interval = NaN;
    turn_interval_hover_excluded = NaN;
else
    time = (0:1:length(turns)-1)'./sam;
    turn_instances = time(logical(turns));
    turn_interval = diff(turn_instances);
    % Saccades before first hover
    indd = find(hover == 1, 1, 'first');
    if isempty(indd)
        turn_instances_hover_excluded = turn_instances;
        turn_interval_hover_excluded = turn_interval;
    else
        turn_instances_hover_excluded = ...
            time([logical(turns(1:indd-1));false_array(indd:end)]);
        turn_interval_hover_excluded = diff(turn_instances_hover_excluded);
    end
end

end