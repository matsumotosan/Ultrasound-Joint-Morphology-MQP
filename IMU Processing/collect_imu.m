function imu_data = collect_imu(baudrate,port,n)
%baudrate = 115200;%9600;%115200;

%port = '/dev/tty.usbmodem14101'; %for Shion
%port = 'COM5'; %for Rosie
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
% Tmax = 60;
% Ts = 0.02;
i = 1;
% ata = 0;
% t = 0;

T(i) = 0;
% FLAG_CASTING = false;
% CubH = [];
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
i = 1;
imu_data = zeros(n,4);

while T(end) <= n
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
        end
    end
    
    % Append to matrix containing pose information
    imu_data(i,:) = ypr;
    i = i + 1;
    
    fprintf('%7d %8.4f %8.4f %8.4f\n',ypr(end,1),ypr(end,2),ypr(end,3),ypr(end,4));
end

fclose(s);

end
