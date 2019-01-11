%% Source
% https://github.com/SayanSeth/MPU-6050-MATLAB-Toolbox/blob/master/Gyroscope_Visualization.m
% SAYAN SETH
%
% USE THIS ONE
%
clear; clc; close all

%% Create serial object for Arduino
baudrate = 115200;

 % Change the COM Port number as needed
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
Fig = figure('Position',[0 40 900 700],'ToolBar','none');
Ax(1) = axes('Position',[.05 .75 0.90 .20]);
grid;
hold on;
H = zeros(1,3);
for k = 1:3
    H(k) = plot(0,0);
end
Ax(2) = axes('Position',[.15 0.05 .6 .6],'CameraPosition',[10 -10 10]);
hold on;
axis([-1 1 -1 1 -1 1]);

ypr = [];

%% Read and plot the data from Arduino
Tmax = 60;
Ts = 0.02;
i = 1;
ata = 0;
t = 0;
% tic % Start timer
T(i)=0;
FLAG_CASTING = false;
CubH = [];
Angles = zeros(1,3);
Flag_Initializing = true;

% Setup
while(Flag_Initializing)
    while(strcmp(s.TransferStatus,'read'))
        pause(0.01);
    end
    readasync(s);
    sms = fscanf(s);
    if ~strcmp(sms(1:3),'ypr')
        fprintf(sms)
    else
        Flag_Initializing = false;
    end
end

% Initialize vector to collect displacement
disp = [];

% Initialize variable to track velocity
vel = [0 0 0]; % [x y z]

% Collect measurements
tic % start timer here
t0 = 0;

while T(end) <= 3000

    T(end+1)=T(end)+1;
    sms='a';
    idx = [];
    Angles = 0;
    while isempty(idx) || numel(Angles)~=3
        sms = fscanf(s);
        idx = find(sms=='r');
        if ~isempty(idx)
            idx = idx(end)+1;
            Angles = sscanf(sms(idx:end),'%f %f %f');
        end
    end
    
    % Assign to yaw, pitch, roll
    Yaw = Angles(3);
    Pitch = Angles(2);
    Roll = Angles(1);
    
    % Append to matrix containing time, yaw, pitch, roll
    t1 = toc;
    ypr = [ypr; t1, Yaw, Pitch, Roll];
    
    % Update displacement - comment for speed (can calculate displacement
    % data post-imaging based on acceleration data
    vel = vel + (t1 - t0) * acc;    % current velocity
    disp = [disp; (t1 - t0) * vel + disp(end,:)];   % current displacement
    
    % Assign t1 to t0 for next iteration
    t0 = t1;
    
    % Update plot (comment for faster sampling rate)
    k = 1;
    vY = get(H(k),'YData');vX = get(H(k),'XData');
    set(H(k),'YData',[vY,Angles(k)]);set(H(k),'XData',[vX,T(end)]);

    k = 2;
    vY = get(H(k),'YData');vX = get(H(k),'XData');
    set(H(k),'YData',[vY,Angles(k)]);set(H(k),'XData',[vX,T(end)]);

    k = 3;
    vY = get(H(k),'YData');vX = get(H(k),'XData');
    set(H(k),'YData',[vY,Angles(k)]);set(H(k),'XData',[vX,T(end)]);

    CubH = Plot_Cube(deg2rad(-Yaw),deg2rad(Pitch),deg2rad(Roll),Ax(2),CubH);
    drawnow;
end

%% Calculate displacement (x,y,z)
% Nested for double integration
v = cumtrapz(t,acc,2);
postDisp = cumtrapz(t,cumtrapz(t,acc,2),2);

% Plot displacement in x,y,z
figure; hold on
subplot(1,2,1)  % acceleration plot
plot(t, acc(1,:), t, acc(2,:), t, acc(3,:)); grid on
title('Acceleration Plot')
xlabel('time (s)')
ylabel('Acceleration (mm/s^2)');
legend('a_x = sin(t)','a_y = cos(t)','a_z = sin(t/2)');

subplot(1,2,2)  % displacement plot
plot(t, postDisp(1,:), t, postDisp(2,:), t, postDisp(3,:)); grid on
title('Displacement Plot')
xlabel('time (s)')
ylabel('Displacement (mm)');
legend('d_x','d_y','d_z');

