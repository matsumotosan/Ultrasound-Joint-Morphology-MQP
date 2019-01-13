function [disp, vel] = calcDisp(acc, t)
%calcDisp Composite nested trapezoidal numerical integration scheme to
%calculate linear displacement from linear acceleration

%% Check dimensions of input
% Transpose if necessary (data arranged vertically)
t_dim = size(t);
acc_dim = size(acc);

if t_dim(2) == 1
    t = t';
end

if acc_dim(2) == 3
    acc = acc';
end

%% Nested cumtrapz for double integration
vel = cumtrapz(t,acc,2);    % velocity
disp = cumtrapz(t,cumtrapz(t,acc,2),2); % displacement

end

