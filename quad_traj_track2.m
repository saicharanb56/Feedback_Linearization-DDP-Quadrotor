function xas = quad_traj_track2()
% EN.530.678: HW#4 sample
% 1) compute a reference path using a polynomial in flat output space
% 2) track the path using backstepping
%
% M. Kobilarov, Spring 2014


% boundary conditions in state space
x0 = [-5; -5; 0; 0; 0; 0];
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

% boundary conditions in flat output space 
y0 = uni_h(x0);
yf = uni_h(xf);
dy0 = x0(4:6);
dyf = xf(4:6);
d2y0 = (1/S.m)*Rot(x0(3))*[0; S.u1] + [0; -9.81];
d2y0 = [d2y0; S.u2/S.J];
d2yf = (1/S.m)*Rot(xf(3))*[0; S.u1] + [0; -9.81];
d2yf = [d2yf; S.u2/S.J];

% compute path coefficients
A = poly3_coeff(y0, dy0, d2y0, yf, dyf, d2yf, T);

% plot desired path
X = A*poly3(0:.01:T);
plot(X(1,:), X(2,:), '-r')
hold on


%%%%%%%%% TRAJECTORY TRACKING %%%%%%%%%%%%%
S.A = A;

% gains
S.k0 = 1; S.k1 = 1; S.k2 = 1; S.k3 = 1;

% perturb initial condition
x = x0 ;

% simulate system
[ts, xas] = ode45(@uni_ode, [0 T], x, [], S);

% visualize
plot(xas(:,1), xas(:,2), '-b');
legend('desired', 'executed')
title('Trajectory of quadcopter')
xlabel('x1'); ylabel('x2');
hold off

% figure(2);
% plot(ts, xas(:,1))
% hold on
% plot(ts, xas(:,2))
% plot(ts, xas(:,3))
% legend('x', 'y', 'yaw')

end
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
    y = x(1:3);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Ua = uni_ctrl(t, xa, S)
% tracking control law
R = Rot(xa(3));

% get desired outputs:
yd = S.A*poly3(t);
dyd = S.A*dpoly3(t);
d2yd = S.A*d2poly3(t);

% get current output and calculate error terms
y = uni_h(xa);
dy = xa(4:6);
% d2y = (1/S.m)*R*[0; xa(end-1)] + [0; -9.81];

% errors
e = y - yd;
de = dy - dyd;
% d2e = d2y - d2yd;

v = d2yd - S.k0*de - S.k1*e;

u1 = S.m*v(2) + S.m*(9.81);
u2 = S.J*(-S.k2*de(3) - S.k3*e(3));

Ua = [u1; u2];

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

u1 = Ua(1);
u2 = Ua(2);

dxa = [xa(4:6);
       (1/S.m)*R*[0; xa(end-1)] + [0; -9.81];
       u2/S.J];

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
