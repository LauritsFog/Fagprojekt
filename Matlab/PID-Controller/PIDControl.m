function [u,i] = PIDControl(i,r,y,y_prev,us,Kp,Ti,Td,Ts)

% i = Integral term
% r = glucose concentration target 
% y = current blood glucose concentration
% y_prev = previous blood glucose concentration
% us is the insulin steady state
% Kp, Ti and Td are the tuning parameters og the PID controller
% Ts = 5min is sampling time

e = y-r; % Error 
p = Kp*e; % Proportional part
i = i+(Kp*Ts/Ti)*e; % Integral part
d = (Kp*Td/Ts)*(y-y_prev); % different part 

u=us+p+i+d;

end

