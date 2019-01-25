function [angvel,angle,saccades,turns,t_angvel,angvel_info] = ...
    calculateAngularVelocity(pts, sam, upsample, freq_cutoff, angvel_cutoff)
% function [angvel,angle,saccades,t_angvel,saccade_info] = ...
%   calculateAngularVelocity(pts, sam, upsample, freq_cutoff, angvel_cutoff)
%
% Inputs:
% 1) XYZ points of the trajectory (pts)
% 2) Sampling Rate
% 3) Up sampling rate
% 4) Freq cutoff (for butterworth filter)
% 5) Angular velocity cutoff (for saccade calculation)
%
% Outputs: 
% 1) Angular velocity
% 2) Angle of the fly between two frames
% 3) Saccades - location of saccade in the time series
% 4) t_angvel - 1 when saccading, else 0
% 5) saccade info - time and amplitude of each saccade
%
%
% Version 3: By Dinesh Natesan
% 29 March 2017

% Defaults
len = size(pts,1);
time = (-(len-1)/sam:1/sam:0)';
upsampled_time = (time(1):1/upsample:0)';
upsampled_len = size(upsampled_time,1);

% Angle calculation between points
upsampled_angle = zeros(upsampled_len-1,1);

% Spline and filter the data
upsampled_pts = nan(upsampled_len,3);
upsampled_pts(:,1) = spline(time, pts(:,1), upsampled_time);
upsampled_pts(:,2) = spline(time, pts(:,2), upsampled_time);
upsampled_pts(:,3) = spline(time, pts(:,3), upsampled_time);

%Find angle between two points and divide it by dt
for i = 2:(size(upsampled_pts,1)-1)
    u = upsampled_pts(i,:)-upsampled_pts(i-1,:);
    v = upsampled_pts(i+1,:)-upsampled_pts(i,:);
    %  DOT product by magnitude, by dt
    upsampled_angle(i) = (atan2(norm(cross(u,v)),dot(u,v))).*180/pi;    % Double checked. Right formula. Turn angle definition exactly same as Breugel et al.
end

% Obtain turn angles
turnangle = abs(upsampled_angle - (upsampled_angle>90).*180); % Bring it to 0-90 scale. Valid?
angle = ButterFilt(upsampled_angle,upsample,freq_cutoff);
angvel = (angle)*upsample;

% Save the upsampled versions to the angvel_info structure
angvel_info = struct;
angvel_info.upsampled_time = upsampled_time;
angvel_info.unfilt_upsampled_angle = upsampled_angle;
angvel_info.filt_upsampled_angle = angle;
angvel_info.unfilt_upsampled_turn_angle = turnangle;
angvel_info.filt_upsampled_turn_angle = ...
    ButterFilt(turnangle,upsample,freq_cutoff);

% Downsample it
angvel = angvel(1:upsample/sam:end);
angle = angle(1:upsample/sam:end);

% Find turns and save them in a struct (for now)
[lmval,indd] = LocalMaxima(angle,0);
turns = zeros(size(angle));
turns(indd) = 1;
% Save amplitude and time of turns in angvel info
if size(lmval,1) == length(lmval)
    angvel_info.turns = [indd,lmval];
else
    angvel_info.turns = [indd',lmval'];
end

% Threshold angvel and find maxima
t_angvel = angvel;
t_angvel(angvel<angvel_cutoff) = 0;
[lmval,indd] = LocalMaxima(t_angvel,0);
saccades = zeros(size(angvel));
saccades(indd) = 1;
% Save amplitude and time in saccade info
 if size(lmval,1) == length(lmval)
     angvel_info.saccades = [indd,lmval];
 else
     angvel_info.saccades = [indd',lmval'];
 end

% make t_angvel a zero or one saccade indicator
t_angvel(angvel>=angvel_cutoff) = 1;

end