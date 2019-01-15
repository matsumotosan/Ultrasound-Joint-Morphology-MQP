%% Source
% https://github.com/SayanSeth/MPU-6050-MATLAB-Toolbox/blob/master/Gyroscope_Visualization.m
% SAYAN SETH
%
% USE THIS ONE
%
clear; clc; close all

%% Create serial object for Arduino
% Change the COM Port number as needed
baudrate = 115200;
s = serial('/dev/tty.usbmodem14101','BaudRate',baudrate);
s.ReadAsyncMode = 'manual';
set(s,'InputBufferSize',100);

pause(2);

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
yprxyz = zeros(1,7);
Flag_Initializing = true;

% Setup
while(Flag_Initializing)
    
    while(strcmp(s.TransferStatus,'read'))
        pause(0.01);
    end
    
    readasync(s);
    sms = fscanf(s);
    
    if ~strcmp(sms(1:6),'yprxyz')
        fprintf(sms)
    else
        Flag_Initializing = false;
    end
end

% Initialize displacement and velocity matrix [x y z]
disp = [0 0 0];
vel = [0 0 0];
pose = zeros(1,7);

% Initialize figure to plot
figure(1); hold on

% Collect measurements
% tic % start timer here

% Stall data collection until IMU output stabilized
STABLE = false;
while ~STABLE
    sms='a';
    idx = [];
    
    if ~isempty(idx)
        idx = idx(end) + 1;
        yprxyz = sscanf(sms(idx:end),'%f %f %f %f %f %f %f',[1 7]);
    end
    
    if all(abs(yprxyz(5:7)) < 15)
        STABLE = true;
    end
end

% Collect IMU output
while T(end) <= 3000
    T(end+1)=T(end)+1;
    sms='a';
    idx = [];
    yprxyz = [0];
    
    while isempty(idx) || numel(yprxyz)~=7
        sms = fscanf(s);
        idx = find(sms=='z');
        if ~isempty(idx)
            idx = idx(end) + 1;
            yprxyz = sscanf(sms(idx:end),'%f %f %f %f %f %f %f', [1 7]);
%             t = toc;
        end
    end
    
    % Append to matrix containing pose information
%     pose = [pose; [t, yprxyz]];
    pose = [pose; yprxyz];

    % Update velocity and displacement - comment for speed (can calculate
    % displacement post-imaging based on acceleration data)
    dt = pose(end,1) - pose(end-1,1);
    vel = vel + dt * pose(end,5:7);         % current velocity
    disp = [disp; dt * vel + disp(end,:)];  % current displacement

    % Plot
%     plot(t,pose(end,5), t,pose(end,6), t,pose(end,7));
    fprintf('%7d %8.4f %8.4f %8.4f %5d %5d %5d\n', pose(end,1), pose(end,2), pose(end,3), pose(end,4), ...
        pose(end,5), pose(end,6), pose(end,7));
%     
%     % Update plot (comment for faster sampling rate)
%     k = 1;
%     vY = get(H(k),'YData');vX = get(H(k),'XData');
%     set(H(k),'YData',[vY,yprxyz(k)]);set(H(k),'XData',[vX,T(end)]);
% 
%     k = 2;
%     vY = get(H(k),'YData');vX = get(H(k),'XData');
%     set(H(k),'YData',[vY,yprxyz(k)]);set(H(k),'XData',[vX,T(end)]);
% 
%     k = 3;
%     vY = get(H(k),'YData');vX = get(H(k),'XData');
%     set(H(k),'YData',[vY,yprxyz(k)]);set(H(k),'XData',[vX,T(end)]);
% 
%     CubH = Plot_Cube(deg2rad(-Yaw),deg2rad(Pitch),deg2rad(Roll),Ax(2),CubH);
%     drawnow;
end

fclose(s);

%% Calculate displacement (x,y,z)
% Nested for double integration
acc = pose(2:end,5:7);
ypr = pose(2:end,2:4)';
t = pose(2:end,1)/1000;

postDisp = calcDisp(acc,t);

% Plot displacement in x,y,z
figure; hold on
subplot(1,2,1)  % acceleration plot
plot(t, acc(:,1), t, acc(:,2), t, acc(:,3)); grid on
title('Acceleration Plot')
xlabel('Time (s)')
ylabel('Acceleration (mm/s^2)');
legend('a_x','a_y','a_z');

subplot(1,2,2)  % displacement plot
plot(t, postDisp(1,:)', t, postDisp(2,:)', t, postDisp(3,:)'); grid on
title('Displacement Plot')
xlabel('Time (s)')
ylabel('Displacement (mm)');
legend('d_x','d_y','d_z');

subplot(3,1,3)  % angular displacement plot
plot(t, ypr(1,:), t, ypr(2,:), t, ypr(3,:)); grid on
title('Angular Displacement Plot')
xlabel('Time (s)');
ylabel('Degrees (\circ)');
legend('Yaw', 'Pitch', 'Roll');

%% Apply filter
fs = 1 / mean(diff(t));                     % sampling frequency
% acc_fil = lowpass(acc,1,fs);              % conventional low pass filter

acc_fil = lowpassfilt(acc);
accel_zero = zeros(1,length(acc_fil));

figure; hold on
subplot(3,1,1)
plot(t,acc_fil(:,1),t,acc(:,1),t,accel_zero);
title('Angular Displacement Yaw (x)')
xlabel('Time (s)');
ylabel('Acceleration (mm/s^2)');

subplot(3,1,2)
plot(t,acc_fil(:,2),t,acc(:,2),t,accel_zero);
title('Angular Displacement Pitch (y)')
xlabel('Time (s)');
ylabel('Acceleration (mm/s^2)');

subplot(3,1,3)
plot(t,acc_fil(:,3),t,acc(:,3),t,accel_zero);
title('Angular Displacement Roll (z)')
xlabel('Time (s)');
ylabel('Acceleration (mm/s^2)');

%% Calculate displacement with filtered linear acceleration
disp_fil = calcDisp(acc_fil,t);

% Plot displacement in x,y,z
figure; hold on
subplot(2,1,1)  % acceleration plot
plot(t, acc_fil(:,1), t, acc_fil(:,2), t, acc_fil(:,3)); grid on
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


