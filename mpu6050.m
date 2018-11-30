clear; clc; close all

%% Setup Arduino
hold all
a = arduino;
mpu = i2cdev(a,'0x68'); % mpu adress is normally 0x68
writeRegister(mpu,hex2dec('B6'),hex2dec('00'),'int16'); % reset
data = zeros(10000,14,'int8'); % preallocating for speed
j = 1;
a1 = animatedline('Color',[1 0 0]); 
a2 = animatedline('Color',[0 1 0]);
a3 = animatedline('Color',[0 0 1]);
a5 = animatedline('Color',[1 1 0]);
a6 = animatedline('Color',[0 1 1]);
a7 = animatedline('Color',[1 0 1]);
legend('Accel_x','Accel_y','Accel_z','Gyro_x','Gyro_y','Gyro_z')

%% MPU 6050 Specs
% https://www.invensense.com/wp-content/uploads/2015/02/MPU-6000-Register-Map1.pdf

% Register    Parameter             Format                        Units
% 59          ACCEL_XOUT[15:8]      16-bit 2's complement value   g
% 60          ACCEL_XOUT[7:0]       16-bit 2's complement value   g
% 61          ACCEL_YOUT[15:8]      16-bit 2's complement value   g
% 62          ACCEL_YOUT[7:0]       16-bit 2's complement value   g
% 63          ACCEL_ZOUT[15:8]      16-bit 2's complement value   g
% 64          ACCEL_ZOUT[7:0]       16-bit 2's complement value   g
% 65          TEMP_OUT[15:8]        16-bit signed value           Celsius
% 66          TEMP_OUT[7:0]         16-bit signed value           Celsius
% 67          GYRO_XOUT[15:8]       16-bit 2's complement value   deg/s
% 68          GYRO_XOUT[7:0]        16-bit 2's complement value   deg/s
% 69          GYRO_YOUT[15:8]       16-bit 2's complement value   deg/s
% 70          GYRO_YOUT[7:0]        16-bit 2's complement value   deg/s
% 71          GYRO_ZOUT[15:8]       16-bit 2's complement value   deg/s
% 72          GYRO_ZOUT[7:0]        16-bit 2's complement value   deg/s

samplerate_div = readRegister(mpu,25,'uint8');  % sample rate divider
accel_fs = readRegister(mpu,28,'int8');         % accelerometer full scale
gyro_fs = readRegister(mpu,27,'int8');          % gyroscope full scale
dlpf = readRegister(mpu,26,'int8');             % DLPF

% Find gyro output rate based on DLPF
if (dlpf == 0 || dlpf == 7)
    gyrorate = 8000;
else
    gyrorate = 1000;
end

samplerate = gyrorate / (1 + samplerate_div);   % gyro sample rate
% accelerometer sample rate is 1 kHz

% Initialize matrix for storing pose information
acc = zeros(10000,3);
gyro = zeros(10000,3);

%% Read MPU 6050
while(true)
    x = 1;
    
    % Read register for acceleration, temperature, and gyroscope
    for i = 59:72
        data(j,x) = readRegister(mpu,i,'int8'); % each register is 8-bit
        x = x + 1;
    end
    
    % Litte-endian system
    y = swapbytes(typecast(data(j,:),'int16'))
    
    % Big-endian system
%     y = typecast(data(j,:),'int16')
    
    % Save acceleration and orientation data
    acc(j,:) = double([y(1), y(2), y(3)]);
    gyro(j,:) = double([y(5), y(6), y(7)]);
    
    % Plot
%     addpoints(a1,j,double(y(1)));
%     addpoints(a2,j,double(y(2)));
%     addpoints(a3,j,double(y(3)));
    addpoints(a5,j,double(y(5)));
    addpoints(a6,j,double(y(6)));
    addpoints(a7,j,double(y(7)));
    
    j = j + 1;
    drawnow limitrate
end

%%


    