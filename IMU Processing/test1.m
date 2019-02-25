function test1()
% Create serial object for Arduino
% Change the COM Port number as needed
baudrate = 115200;%9600;%115200;

%port = '/dev/tty.usbmodem14101'; %for Shion
port = 'COM5'; %for Rosie
%port = 'COM6'; %for Olivia
s = serial(port,'BaudRate',baudrate);
s.ReadAsyncMode = 'manual';
set(s,'InputBufferSize',100);

pause(2);

% Connect the serial port to Arduino
try
    fopen(s);
catch err
    fclose(instrfind);
    error('Error: Select correct COM Port where Arduino is connected.');
end

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
%figure(1); hold on

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
while T(end) <= 1000
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

end

fclose(s);

end
