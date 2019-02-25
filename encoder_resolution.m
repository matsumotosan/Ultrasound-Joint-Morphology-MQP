%% rotary encoder resolution

p = input('# Pulses: ');
ppr = input('# Pulses Per Revolution: ');

angle = p*360/ppr; %Angle

fprintf('Spatial Resolution: %0.4f\n per Degree',angle);
