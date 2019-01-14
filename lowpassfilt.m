function filtered = lowpassfilt(signal)

filter = designfilt('lowpassfir', ...
    'PassbandFrequency',0.03,'StopbandFrequency',0.08, ...
    'PassbandRipple',0.5,'StopbandAttenuation',65, ...
    'DesignMethod','equiripple');

filtered = filtfilt(filter,signal);

end
