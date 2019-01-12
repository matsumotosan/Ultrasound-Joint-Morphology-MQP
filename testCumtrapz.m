% This script is meant to test the accuracy of the double integration
% scheme using cumtrapz. 
clear; close all; clc;

%% Trial dataset
t = pi * [0 0.1 0.24 0.5 0.51 0.75 0.9 1 1.2 1.23 1.5 1.6 1.8 1.85 2.0];
acc = [sin(t); cos(t); sin(t/2)];

% Nested cumtrapz for double trapezoidal integration
v = cumtrapz(t,acc,2);
postDisp = cumtrapz(t,cumtrapz(t,acc,2),2);

%% Plot displacement in x,y,z
figure; hold on
subplot(1,2,1)  % acceleration plot
plot(t, acc(1,:), '-^', t, acc(2,:), '-s', t, acc(3,:), '-d'); grid on
title('Acceleration Plot')
xlabel('time (s)')
ylabel('Acceleration (mm/s^2)');
lgnd = legend('a_x = sin(t)','a_y = cos(t)','a_z = sin(t/2)');
set(lgnd,'color','none');

subplot(1,2,2)  % displacement plot
plot(t, postDisp(1,:), '-^', t, postDisp(2,:), '-s', ...
     t, postDisp(3,:), '-d'); grid on
title('Displacement Plot')
xlabel('time (s)')
ylabel('Displacement (mm)');
legend('d_x','d_y','d_z');

%% Compare numerical to analytic solution
% Calculate analytic solution
t_ana = linspace(0, 2*pi, 100);
d_x_ana = -sin(t_ana) + t_ana;
d_y_ana = -cos(t_ana) + 1;
d_z_ana = -4 * sin(t_ana / 2) + 2 * t_ana;

% Plot numerical and analytic solutions
figure; hold on
plot(t, postDisp(1,:), '-^', t, postDisp(2,:), '-s', ...
     t, postDisp(3,:), '-d'); grid on
plot(t_ana, d_x_ana, '--', t_ana, d_y_ana, '--', t_ana, d_z_ana, '--')
xlim([0 2*pi])

xlabel('Time (s)')
ylabel('Displacement (mm)');
legend('d_x (Numerical)', 'd_y (Numerical)', 'd_z (Numerical)', ...
    'd_x (Analytical)', 'd_y (Analytical)', 'd_z (Analytical)');
title('Comparison of Numerical and Analytic Results');

% Error analysis
f_x_ana = @(x) -sin(x) + x;
f_y_ana = @(x) -cos(x) + 1;
f_z_ana = @(x) -4 * sin(x ./ 2) + 2 .* x;

d_x_err = abs((f_x_ana(t) - postDisp(1,:)) ./ f_x_ana(t));
d_y_err = abs((f_y_ana(t) - postDisp(2,:)) ./ f_y_ana(t));
d_z_err = abs((f_z_ana(t) - postDisp(3,:)) ./ f_z_ana(t));

d_x_err(isnan(d_x_err)) = 0;
d_y_err(isnan(d_y_err)) = 0;
d_z_err(isnan(d_z_err)) = 0;

figure; hold on
plot(t, d_x_err, t, d_y_err, t, d_z_err); grid on
xlabel('Time (s)')
ylabel('Relative Error');
legend('x','y','z')
title('Relative Error of Numerical Results')

