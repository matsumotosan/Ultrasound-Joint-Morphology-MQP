function filtered = lowpassfilt(signal)
%LOWPASSFILT function to filter out unwanted high frequency noise
% 
% Custom-built to deal with the data acquired from the MPU6050 IMU
%
% Input: signal data (pose data)
% Output: filtered data

%design the filter
filter = designfilt('lowpassfir', ...
    'PassbandFrequency',0.03,'StopbandFrequency',0.08, ...
    'PassbandRipple',0.5,'StopbandAttenuation',65, ...
    'DesignMethod','equiripple');

%filter the signal using filtfilt
filtered = filtfilt(filter,signal);

end
