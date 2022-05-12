% function A = quad_traj_track()

addpath('ddp_quad\')
addpath('utils\')

% boundary conditions in state space
x0 = [-7; 2; 0; 0; 0; 0];
xf = [0; 0; 0; 0; 0; 0];
T = 10;

%%%%%%%%% TRAJECTORY GENERATION %%%%%%%%%%%%%

% norm of initial and final velocity along desired path
% it determines how much the curves will bend
% and can be freely chosen

% quadrotor parameters
S.m = 0.2;
S.J = 0.15;
S.u1 = 1;
S.u1dot = 1;
S.u1ddot = 1;
S.u2 = 1;

% ddp trajectory generation
desired = ddp_quad_obst_nl([x0(1); 0; x0(2)], T);
close all;
%%
% modify to 2D
yd_ddp = [ desired.xs(1,:);
           desired.xs(3,:)];

yd_ddp = [yd_ddp(:,1)+rand(size(yd_ddp,1),10)*0.1 yd_ddp];

traj_ts = linspace(-1,T,size(yd_ddp,2));
A = [ polyfit(traj_ts, yd_ddp(1,:), 6)
      polyfit(traj_ts, yd_ddp(2,:), 6)];


% boundary conditions in flat output space 
% y0 = uni_h(x0);
% yf = uni_h(xf);
% dy0 = x0(4:5);
% dyf = xf(4:5);
% d2y0 = (1/S.m)*Rot(x0(3))*[0; S.u1] + [0; -9.81];
% d2yf = (1/S.m)*Rot(xf(3))*[0; S.u1] + [0; -9.81];
% d3y0 = (1/S.m)*Rot(x0(3))*[-S.u1*x0(6); S.u1dot];
% d3yf = (1/S.m)*Rot(xf(3))*[-S.u1*xf(6); S.u1dot];
% d4y0 = (1/S.m)*Rot(x0(3))*[-2*S.u1dot*x0(6); S.u1*x0(6)^2] + (1/S.m)*Rot(x0(3))*[-S.u1*S.u2/S.J; S.u1ddot];
% d4yf = (1/S.m)*Rot(xf(3))*[-2*S.u1dot*xf(6); S.u1*xf(6)^2] + (1/S.m)*Rot(xf(3))*[-S.u1*S.u2/S.J; S.u1ddot];

% compute path coefficients
% A = poly3_coeff(y0, dy0, d2y0, d3y0, d4y0, yf, dyf, d2yf, d3yf, d4yf, T);
% A = poly3_coeff(y0, dy0, d2y0, yf, dyf, d2yf, T);

% A = randn(2,7);

% plot desired path
X = A*poly3(0:.01:T);
plot(X(1,:), X(2,:), '-r', 'LineWidth', 2)
hold on
%%

%%%%%%%%% TRAJECTORY TRACKING %%%%%%%%%%%%%
S.A = A;

% gains
S.k0 = 0.33; S.k1 = 1.667; S.k2 = 120; S.k3 = 3;

% perturb initial condition
x = x0 + [0.1; 0.1; 0; 0; 0; 0];

% augmented state with dynamic compensator, i.e xi=u1
xa = [x; S.u1; S.u1dot];

% simulate system
[ts, xas] = ode45(@uni_ode, [0 T], xa, [], S);

% visualize
plot(xas(:,1), xas(:,2), '-b', 'LineWidth', 2);
legend('desired', 'executed')
title('Trajectory Tracking of Quadcopter')
xlabel('x1'); ylabel('x2');
hold off

% figure(2);
% plot(ts, xas(:,1))
% hold on
% plot(ts, xas(:,2))
% plot(ts, xas(:,3))
% legend('x', 'y', 'yaw')

% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% function A = poly3_coeff(y0, dy0, d2y0, d3y0, d4y0, yf, dyf, d2yf, d3yf, d4yf, T)
% % computes cubic curve connecting (y0,dy0) and (yf, dyf) at time T
% 
% Y = [y0, dy0, d2y0, d3y0, d4y0, yf, dyf, d2yf, d3yf, d4yf];
% L = [poly3(0), dpoly3(0), d2poly3(0), d3poly3(0), d4poly3(0) ...
%      poly3(T), dpoly3(T), d2poly3(T), d3poly3(T), d4poly3(T)];
% A = Y/L;
% end

function A = poly3_coeff(y0, dy0, d2y0, yf, dyf, d2yf, T)
% computes cubic curve connecting (y0,dy0) and (yf, dyf) at time T

Y = [y0, dy0, d2y0, yf, dyf, d2yf];
L = [poly3(0), dpoly3(0), d2poly3(0) ...
     poly3(T), dpoly3(T), d2poly3(T)];
A = Y/L;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = uni_h(x)
% output function
    y = x(1:2);
end

function f = poly3(t)
    f = [t.^6; t.^5; t.^4; t.^3; t.^2; t; ones(size(t))];
end

function f = dpoly3(t)
    f = [6*t.^5; 5*t.^4; 4*t.^3; 3*t.^2; 2*t; ones(size(t)); zeros(size(t))];
end

function f = d2poly3(t)
    f = [30*t.^4; 20*t.^3; 12*t.^2; 6*t; 2; zeros(size(t)); zeros(size(t))];
end

function f = d3poly3(t)
    f = [120*t.^3; 60*t.^2; 24*t; 6; zeros(size(t)); zeros(size(t)); zeros(size(t))];
end

function f = d4poly3(t)
    f = [360*t.^2; 120*t; 24; zeros(size(t)); zeros(size(t)); zeros(size(t)); zeros(size(t))];
end

% function f = poly3(t)
%     f = [t.^10; t.^9; t.^8; t.^7; t.^6; t.^5; t.^4; t.^3; t.^2; t; ones(size(t))];
% end
% 
% function f = dpoly3(t)
%     f = [10*t.^9; 9*t.^8; 8*t.^7; 7*t.^6; 6*t.^5; 5*t.^4; 4*t.^3; 3*t.^2; 2*t; ones(size(t)); zeros(size(t))];
% end
% 
% function f = d2poly3(t)
%     f = [90*t.^8; 72*t.^7; 56*t.^6; 42*t.^5; 30*t.^4; 20*t.^3; 12*t.^2; 6*t; 2; zeros(size(t)); zeros(size(t))];
% end
% 
% function f = d3poly3(t)
%     f = [720*t.^7; 504*t.^6; 336*t.^5; 210*t.^4; 120*t.^3; 60*t.^2; 24*t; 6; zeros(size(t)); zeros(size(t)); zeros(size(t))];
% end
% 
% function f = d4poly3(t)
%     f = [5040*t.^6; 3024*t.^5; 1680*t.^4; 840*t.^3; 360*t.^2; 120*t; 24; zeros(size(t)); zeros(size(t)); zeros(size(t)); zeros(size(t))];
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Ua = uni_ctrl(t, xa, S)
% tracking control law
R = Rot(xa(3));

% get desired outputs:
yd = S.A*poly3(t);
dyd = S.A*dpoly3(t);
d2yd = S.A*d2poly3(t);
d3yd = S.A*d3poly3(t);
d4yd = S.A*d4poly3(t);

% get current output and calculate error terms
y = uni_h(xa);
dy = xa(4:5);
d2y = (1/S.m)*R*[0; xa(end-1)] + [0; -9.81];
d3y = R*[-xa(end-1)*xa(6); xa(end)];

% errors
e = y - yd;
de = dy - dyd;
d2e = d2y - d2yd;
d3e = d3y - d3yd;

% z-state
B = zeros(8,8);
B(1:2,3:4) = eye(2);
B(3:4,5:6) = eye(2);
B(5:6,7:8) = eye(2);
B(7:8,:) = [-S.k0*eye(2), -S.k1*eye(2), -S.k2*eye(2), -S.k3*eye(2)];
z = B*[e; de; d2e; d3e];

%virtual input
v = d4yd - S.k3*d3e - S.k2*d2e - S.k1*de - S.k0*e;

Ua = [-S.J/xa(end-1), 0; 0, 1]*(S.m*R'*v - [-2*xa(end)*xa(6); -xa(end-1)*xa(6)^2]);

% This Ua that this function returns is [u2; u1ddot]

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function R = Rot(a)

R = [cos(a), -sin(a);
     sin(a), cos(a)];
end

function dxa = uni_ode(t, xa, S)
% state = [x y yaw xDot yDot yawDot u1 u1Dot]
% [u1; u1Dot] is the dynamic compensator

% quadrotor ODE
Ua = uni_ctrl(t, xa, S);
R = Rot(xa(3));

u2 = Ua(1);
u1dDot = Ua(2);

%noise added below
% ua = ua + 0.8*randn(2,1);


% dxa = [cos(xa(3))*x4;
%        sin(xa(3))*x4;
%        tan(u1)*x4/S.L;
%        u2];

dxa = [xa(4:6);
       (1/S.m)*R*[0; xa(end-1)] + [0; -9.81];
       u2/S.J;
       xa(end);
       u1dDot];

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
