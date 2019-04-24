%% IMU Global Positioning Data Processing Code
% This code tests the accuracy of the IMU's global positioning skills. It
% uses the calcDisp function written to calculate displacement from the
% provided acceleration data and a custom-built lowpass filter.

close all; clear all; clc;

%% Load relevant pose files

load('do_nothing.mat'); %newest no movement file
%load('square_test.mat'); %start and end points are same, square trial
%load('line_test.mat'); %line test
%load('new_signal.mat'); %no movement file
%load('disp_trial.mat'); %had weird data

%% Extract Acc, Gyr, Time
% Nested for double integration
acc = pose(2:end,5:7);      %extract acceleration data
ypr = pose(2:end,2:4)';     %extract gyro data
t = pose(2:end,1)/1000;     %extract time data & convert to ms

%% Calculate displacement (x,y,z)
postDisp = calcDisp(acc,t); %calculate pose data from calcDisp function

%% Plot original data
% Plot displacement in x,y,z
figure; hold on
subplot(3,1,1)  % acceleration plot
plot(t, acc(:,1), t, acc(:,2), t, acc(:,3)); grid on
title('Original Acceleration Plot')
xlabel('Time (s)')
ylabel('Acceleration (mm/s^2)');
legend('a_x','a_y','a_z');

subplot(3,1,2)  % displacement plot
plot(t, postDisp(1,:)', t, postDisp(2,:)', t, postDisp(3,:)'); grid on
title('Original Displacement Plot')
xlabel('Time (s)')
ylabel('Displacement (mm)');
legend('d_x','d_y','d_z');

subplot(3,1,3)  % angular displacement plot
plot(t, ypr(1,:), t, ypr(2,:), t, ypr(3,:)); grid on
title('Original Angular Displacement Plot')
xlabel('Time (s)');
ylabel('Angular Displacement (\circ)');
legend('Yaw', 'Pitch', 'Roll');

%% Calculate Sample Frequency
fs = 1 / mean(diff(t));             % sampling frequency

%% Median Filter
acc_fil = medfilt1(acc,10);        % median filt of acc data over 5 sample window

%% Conventional Low-pass Filter
% acc_fil = lowpass(acc,1,fs);      % conventional low pass filter

%% Custom Low-pass filter
% acc_fil = lowpassfilt(acc);       % apply custom low-pass filter

%% Butterworth Band-Pass Filtering
order = 2;     %order of the filter
fcutlow=1;     %low cut frequency in Hz
fcuthigh=5;   %high cut frequency in Hz

[b,a]=butter(order,[fcutlow,fcuthigh]/(fs/2),'bandpass');
bp_acc_fil = filtfilt(b, a, acc);

%% Create expected plot for do nothing case
a0 = zeros(length(acc_fil),3); %to compare correct '0' in plot

%% Plot filtered data

figure; hold on
ax = subplot(3,1,1);
plot(t,acc(:,1),'c',t,bp_acc_fil(:,1),'b',t,acc_fil(:,1),'r',t,a0(:,1),'k'); grid on;
title('Acceleration in x');
xlabel('Time (s)');
ylabel('a_x (mm/s^2)');
legend('original','bandpass','median','expected');

ay = subplot(3,1,2);
plot(t,acc(:,2),'c',t,bp_acc_fil(:,2),'b',t,acc_fil(:,2),'r',t,a0(:,2),'k'); grid on;
title('Acceleration in y');
xlabel('Time (s)');
ylabel('a_y (mm/s^2)');
legend('original','bandpass','median','expected');

az = subplot(3,1,3);
plot(t,acc(:,3),'c',t,bp_acc_fil(:,3),'b',t,acc_fil(:,3),'r',t,a0(:,3),'k'); grid on;
title('Acceleration in z');
xlabel('Time (s)');
ylabel('a_z (mm/s^2)');
legend('original','bandpass','median','expected'); 

linkaxes([ax ay az],'xy');

%% Calculate displacement with filtered linear acceleration

disp_fil = calcDisp(bp_acc_fil,t);

% Plot displacement in x,y,z
figure; hold on
subplot(2,1,1)  % acceleration plot
plot(t, bp_acc_fil(:,1), t, bp_acc_fil(:,2), t, bp_acc_fil(:,3)); grid on
title('Filtered Acceleration Plot')
xlabel('Time (s)')
ylabel('Acceleration (mm/s^2)');
legend('a_x','a_y','a_z');

subplot(2,1,2)  % displacement plot
plot(t, disp_fil(1,:)', t, disp_fil(2,:)', t, disp_fil(3,:)'); grid on
title('Filtered Displacement Plot')
xlabel('Time (s)')
ylabel('Displacement (mm)');
legend('d_x','d_y','d_z');

%% Calc error

avg_disp_error = mean(abs(disp_fil(:,end)));
fprintf('Averaged end displacement error of %.2f mm\n',avg_disp_error);