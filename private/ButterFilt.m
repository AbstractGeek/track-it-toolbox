function [filt_signal] = ButterFilt( signal,sampling_freq,cutoff_freq)

% Filters in input signal at the given cutoff frequency. The cutoff freq should be atleast less than sam/2 (Nyquist).
%
% [Inputs]:
% Signal (an be a matrix)
% Sampling frequency
% Cutoff frequency
%
% [Outputs]:
% Filtered Signal
%
% Version 2.001
% Last Modified on 24th May 2012
% Last Checked 1st June 2012
% by Dinesh Natesan

if (size(signal,2)>1 || size(signal,3)>1)
    filt_signal=zeros(size(signal,1),size(signal,2));
    for j = 1:size(signal,3)
    for i=1:size(signal,2)
        filt_signal(:,i,j) = ButterFilt( signal(:,i,j),sampling_freq,cutoff_freq);
    end
    end
else
    Ny_freq = sampling_freq/2;
    [b,a]= butter(4,cutoff_freq/Ny_freq,'low');
    %The order of this above designed filter is '4'. This was found out by comparing the magnitude and phase response of each order of butter filter, and
    %selecting a order which minimizes ripple effect and maximizes the passband
    filt_signal= filtfilt(b,a,signal);
end

end

