%% IMU Analysis

%load in the pose data and analyze it for accuracy

close all; clear all; clc

%% Extract data of interest

%load('0_180_IMU_test_1.mat');
%load('0_180_IMU_test_2.mat');
load('0_180_IMU_test_3.mat');

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

% figure(2);
% 
% subplot(3,1,1)
% plot(t,ypr(:,1),t,lp_ypr(:,1),'--');
% ylabel('Angle (degree)')
% xlabel('Time (s)');
% title('Yaw');
% legend('Raw','Filtered')
% grid on;
% 
% subplot(3,1,2)
% plot(t,ypr(:,2),t,lp_ypr(:,2),'--');
% ylabel('Angle (degree)')
% xlabel('Time (s)');
% title('Pitch');
% legend('Raw','Filtered')
% grid on;
% 
% subplot(3,1,3)
% plot(t,ypr(:,3),t,lp_ypr(:,3),'--');
% ylabel('Angle (degree)')
% xlabel('Time (s)');
% title('Roll');
% legend('Raw','Filtered')
% grid on;
% 
% %% Compare filtered and unfiltered data
% 
% figure(5)
% subplot(1,2,1)
% plot(t, ypr(:,1), t, ypr(:,2), t, ypr(:,3));
% title('Unfiltered Rotation')
% xlabel('Time (s)');
% ylabel('Degrees (\circ)');
% legend('Yaw', 'Pitch', 'Roll');
% grid on;
% 
% subplot(1,2,2)
% plot(t, lp_ypr(:,1),'--', t, lp_ypr(:,2),'--', t, lp_ypr(:,3),'--');
% title('Filtered Rotation')
% xlabel('Time (s)');
% ylabel('Degrees (\circ)');
% legend('Yaw', 'Pitch', 'Roll');
% grid on; 
% 
% figure(6)
% plot(t,ypr(:,1), t, ypr(:,2), t, ypr(:,3),t, lp_ypr(:,1),'--', t, lp_ypr(:,2),'--', t, lp_ypr(:,3),'--');
% title('Raw vs. Filtered Rotation')
% xlabel('Time (s)');
% ylabel('Degrees (\circ)');
% legend('Yaw', 'Pitch', 'Roll','Filtered Yaw', 'Filtered Pitch', 'Filtered Roll');
% grid on; 

%% extract

yaw = lp_ypr(:,1);
high_threshold = max(yaw)-1.5;
low_threshold = min(yaw)+1.5;

high = zeros(1, length(yaw));
low = high;
high(yaw>high_threshold) = true;
low(yaw<low_threshold) = true;

thresh_y_high = diff(high);
thresh_y_low = diff(low);

h_start_idx = [1, find(thresh_y_high>0)];
h_end_idx = [find(thresh_y_high<0), length(yaw)];
l_start_idx = find(thresh_y_low>0);
l_end_idx = find(thresh_y_low<0);

figure(7)
plot(t,yaw,t(h_start_idx),yaw(h_start_idx),'or',t(h_end_idx),yaw(h_end_idx),'*r',t(l_start_idx),yaw(l_start_idx),'ok',t(l_end_idx),yaw(l_end_idx),'*k'); grid on;
xlabel('time (ms)');
ylabel('angle (degrees)');
title('segmenting data');
legend('raw data','high start','high end','low start','low end');

high_yaw = vertcat(yaw(h_start_idx(1):h_end_idx(1)),yaw(h_start_idx(2):h_end_idx(2)),yaw(h_start_idx(3):h_end_idx(3)));
low_yaw = vertcat(yaw(l_start_idx(1):l_end_idx(1)),yaw(l_start_idx(2):l_end_idx(2)));

high_theoretical = zeros(length(high_yaw),1);
low_theoretical = ones(length(low_yaw),1)*-180;

t_hy = t(1:length(high_yaw));
t_ly = t(1:length(low_yaw));

figure(8)
subplot(1,2,1)
plot(t_hy, high_yaw, t_hy, high_theoretical,'r'); grid on;
xlabel('time (s)');
ylabel('angle (degrees)');
title('high angle data');
legend('experimental','theoretical');

subplot(1,2,2)
plot(t_ly, low_yaw, t_ly, low_theoretical,'r'); grid on;
xlabel('time (s)');
ylabel('angle (degrees)');
title('low angle data');
legend('experimental','theoretical');

higherror = abs(high_theoretical-high_yaw);
lowerror = abs(low_theoretical-low_yaw);

avg_higherror = mean(higherror); %percent diff
avg_lowerror = mean(lowerror); %percent diff

fprintf('High angle error was: %.3f \n', avg_higherror);
fprintf('Low angle error was: %.3f \n', avg_lowerror);
% disp(['Min high error: ' num2str(min(higherror))]);
% disp(['Min low error: ' num2str(min(lowerror))]);
% disp(['Max high error: ' num2str(max(higherror))]);
% disp(['Max low error: ' num2str(max(lowerror))]);
disp(['SD high error: ' num2str(std(higherror))]);
disp(['SD low error: ' num2str(std(lowerror))]);