% This script can be run to collect IMU readings through MATLAB. The
% corresponding Arduino sketch must be uploaded and running prior to
% running this script. Correct COM port must be chosen for user's system.
% This script is a modified version of code obtained from link below:
%
% https://github.com/SayanSeth/MPU-6050-MATLAB-Toolbox/blob/master/Gyroscope_Visualization.m

clear; clc; close all

%% Create serial object for Arduino
% Change the COM Port number as needed
baudrate = 115200; %74880;%9600;

%port = '/dev/tty.usbmodem14101'; %for Shion
%port = 'COM5'; %for Rosie
port = 'COM6'; %for Olivia
s = serial(port,'BaudRate',baudrate);
s.ReadAsyncMode = 'manual';
set(s,'InputBufferSize',100);

pause(0.5);

%% Connect the serial port to Arduino
try
    fopen(s);
catch err
    fclose(instrfind);
    error('Error: Select correct COM Port where Arduino is connected.');
end

%% Prepare Figures
% Fig = figure('Position',[0 40 900 700],'ToolBar','none');
% Ax(1) = axes('Position',[.05 .75 0.90 .20]);
% grid;
% hold on;
% H = zeros(1,3);
% for k = 1:3
%     H(k) = plot(0,0);
% end
% Ax(2) = axes('Position',[.15 0.05 .6 .6],'CameraPosition',[10 -10 10]);
% hold on;
% axis([-1 1 -1 1 -1 1]);

%% Read and plot the data from Arduino
Tmax = 60;
Ts = 0.02;
i = 1;
ata = 0;
t = 0;

T(i) = 0;
FLAG_CASTING = false;
CubH = [];
ypr = zeros(1,4);
Flag_Initializing = true;

% Setup
while(Flag_Initializing)
    
    while(strcmp(s.TransferStatus,'read'))
        pause(0.1);
    end
    
    readasync(s);
    sms = fscanf(s);
    
    if ~strcmp(sms(1:3),'ypr')
        fprintf(sms)
    else
        Flag_Initializing = false;
    end
end

% Initialize pose matrix [t y p r]
pose = zeros(1,4);

%% Collect measurements
% tic % start timer here
% Stall data collection until IMU output stabilized
STABLE = false;
while ~STABLE
    sms='a';
    idx = [];
    
    if ~isempty(idx)
        idx = idx(end) + 1;
        ypr = sscanf(sms(idx:end),'%f %f %f %f',[1 4]);
    end
    
    if all(abs(ypr(2:4)) < 10)
        STABLE = true;
    end
end

% Collect IMU output
while T(end) <= 3000
    T(end+1)=T(end)+1;
    sms='a';
    idx = [];
    ypr = [0];
    
    while isempty(idx) || numel(ypr)~=4
        sms = fscanf(s);
        idx = find(sms=='r');
        if ~isempty(idx)
            idx = idx(end) + 1;
            ypr = sscanf(sms(idx:end),'%f %f %f %f', [1 4]);
%             t = toc;
        end
    end
    
    % Append to matrix containing pose information
    pose = [pose; ypr];
    fprintf('%7d %8.4f %8.4f %8.4f\n', pose(end,1), pose(end,2), pose(end,3), pose(end,4));

end

fclose(s);

%% Extract data of interest
close all; clc

ypr = pose(2:end,2:4);
t = pose(2:end,1)/1000;

figure(1); 
plot(t, ypr(:,1), t, ypr(:,2), t, ypr(:,3));
title('Rotation about IMU local z-axis 0 to 180^{\circ}')
xlabel('Time (s)');
ylabel('Degrees (\circ)');
legend('Yaw', 'Pitch', 'Roll');
grid on; 

%% Apply smoothing filter

lp_ypr = smoothdata(ypr,1,'movmedian'); %better for euler
%lp_ypr = smoothdata(ypr,1,'sgolay');

figure(2);

subplot(3,1,1)
plot(t,ypr(:,1),t,lp_ypr(:,1),'--');
ylabel('Angle (degree)')
xlabel('Time (s)');
title('Yaw');
legend('Raw','Filtered')
grid on;

subplot(3,1,2)
plot(t,ypr(:,2),t,lp_ypr(:,2),'--');
ylabel('Angle (degree)')
xlabel('Time (s)');
title('Pitch');
legend('Raw','Filtered')
grid on;

subplot(3,1,3)
plot(t,ypr(:,3),t,lp_ypr(:,3),'--');
ylabel('Angle (degree)')
xlabel('Time (s)');
title('Roll');
legend('Raw','Filtered')
grid on;

%% Compare filtered and unfiltered data

figure(5)
subplot(1,2,1)
plot(t, ypr(:,1), t, ypr(:,2), t, ypr(:,3));
title('Unfiltered Rotation')
xlabel('Time (s)');
ylabel('Degrees (\circ)');
legend('Yaw', 'Pitch', 'Roll');
grid on;

subplot(1,2,2)
plot(t, lp_ypr(:,1),'--', t, lp_ypr(:,2),'--', t, lp_ypr(:,3),'--');
title('Filtered Rotation')
xlabel('Time (s)');
ylabel('Degrees (\circ)');
legend('Yaw', 'Pitch', 'Roll');
grid on; 

figure(6)
plot(t,ypr(:,1), t, ypr(:,2), t, ypr(:,3),t, lp_ypr(:,1),'--', t, lp_ypr(:,2),'--', t, lp_ypr(:,3),'--');
title('Raw vs. Filtered Rotation')
xlabel('Time (s)');
ylabel('Degrees (\circ)');
legend('Yaw', 'Pitch', 'Roll','Filtered Yaw', 'Filtered Pitch', 'Filtered Roll');
grid on; 
