addpath('ddp_quad\')
addpath('utils\')

% boundary conditions in state space
x0 = [-7; 2; pi/12; 0; 0; 0];
xf = [0; 0; 0; 0; 0; 0];
T = 30;
S.T = T;

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

% ddp trajectory generation in 3D
desired = ddp_quad_obst_nl([x0(1); 0; x0(2)], T);
close all;
%%

S.xs = desired.xs;

% modify to 2D
yd_ddp = [ desired.xs(1,:);
           desired.xs(3,:)];

yd_ddp = [yd_ddp(:,1)+rand(size(yd_ddp,1),10)*0.1 yd_ddp];

traj_ts = linspace(-1,T,size(yd_ddp,2));
A = [ polyfit(traj_ts, yd_ddp(1,:), 6)
      polyfit(traj_ts, yd_ddp(2,:), 6)];

% plot desired path
X = A*poly3(0:.01:T);
plot(X(1,:), X(2,:), '-r', 'LineWidth', 2)
hold on
%%

%%%%%%%%% TRAJECTORY TRACKING %%%%%%%%%%%%%
S.A = A;

% gains
S.k0 = 20; S.k1 = 80; S.k2 = 200; S.k3 = 3;

% perturb initial condition
x = x0 + [2; 2; zeros(4,1)];

% augmented state with dynamic compensator, i.e xi=[u1; u1Dot]
xa = [x; S.u1; S.u1dot];

global Us 
Us = [];

% simulate system
[ts, xas] = ode45(@uni_ode, [0 T], xa, [], S);

% visualize
plot(xas(:,1), xas(:,2), '-b', 'LineWidth', 2);
legend('desired', 'executed')
title('Trajectory Tracking of Quadcopter')
xlabel('x1'); ylabel('x2');
hold off

figure(2);
plot(ts, xas(:,1))
hold on
plot(ts, xas(:,2))
plot(ts, xas(:,3))
legend('x', 'y', 'yaw')
title('Variation of states over time')
xlabel('Time in seconds'); ylabel('States');
hold off

figure(3); hold on
plot(Us(:,1), Us(:,2), 'b', 'LineWidth', 2)
plot(Us(:,1), Us(:,3), 'r', 'LineWidth', 2)
xlabel('t'); ylabel('u')
legend({'u1', 'u2'})

% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

global Us 
Us = [Us; t xa(end-1) u2];

dxa = [xa(4:6);
       (1/S.m)*R*[0; xa(end-1)] + [0; -9.81];
       u2/S.J;
       xa(end);
       u1dDot];

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
