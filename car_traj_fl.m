function f = car_traj_fl()
% EN.530.678: HW#4 sample
% 1) compute a reference path using a polynomial in flat output space
% 2) track the path using feedback linearization
%
% M. Kobilarov, Spring 2014


% boundary conditions in state space
x0 = [-5; -3; 0.5 ; 1];
xf = [0; 0; 1 ; 1];
T = 10;

%%%%%%%%% TRAJECTORY GENERATION %%%%%%%%%%%%%

% norm of initial and final velocity along desired path
% it determines how much the curves will bend
% and can be freely chosen
S.u1 = 1;
S.L = 1;

% boundary conditions in flat output space 
y0 = uni_h(x0);
yf = uni_h(xf);
dy0 = S.u1*[cos(x0(3)); sin(x0(3))]; % desired starting velocity
dyf = S.u1*[cos(xf(3)); sin(xf(3))]; % desired end velocity

% compute path coefficients
A = poly3_coeff(y0, dy0, yf, dyf, T);

% plot desired path
X = A*poly3([0:.01:T]);
plot(X(1,:), X(2,:), '-r')
hold on


%%%%%%%%% TRAJECTORY TRACKING %%%%%%%%%%%%%
S.A = A;

% gains
S.k = [1;2];

% perturb initial condition
x = x0 + [.25;.25;.1;0];

% simulate system
[ts, xas] = ode45(@uni_ode, [0 T], x, [], S);

% visualize
plot(xas(:,1), xas(:,2), '-b');

legend('desired', 'executed')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function A = poly3_coeff(y0, dy0, yf, dyf, T)
% computes cubic curve connecting (y0,dy0) and (yf, dyf) at time T

Y = [y0, dy0, yf, dyf];
L = [poly3(0), dpoly3(0), poly3(T), dpoly3(T)];
A = Y/L;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = uni_h(x)
% output function

y = x(1:2);


function f = poly3(t)
f = [t.^3; t.^2; t; ones(size(t))];

function f = dpoly3(t)
f = [3*t.^2; 2*t; ones(size(t)); zeros(size(t))];

function f = d2poly3(t)
f = [6*t; 2; zeros(size(t)); zeros(size(t))];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function u = uni_ctrl(t, xa, S)
% tracking control law

% get desired outputs:
yd = S.A*poly3(t);
dyd = S.A*dpoly3(t);
d2yd = S.A*d2poly3(t);

% get current output
y = uni_h(xa);

% compensator, i.e.  xi=u1
xi = xa(end);

% current velocity
dy = [cos(xa(3)); sin(xa(3))]*xi;

% error state
z1 = y - yd;
z2 = dy - dyd;

% virtual inputs
v = d2yd - S.k(1)*z1 -S.k(2)*z2;

% augmented inputs ua=(dxi, u2)
ua = [-sin(xa(3))/(xa(4)^2) cos(xa(3))/(xa(4)^2); 
       sin(xa(3)) cos(xa(3))]*v;
   
u = [atan(ua(1)*S.L);
     ua(2)];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dxa = uni_ode(t, xa, S)
% unicycle ODE
ua = uni_ctrl(t, xa, S);

xi = xa(end);
u1 = ua(1);
u2 = ua(2);


dxa = [cos(xa(3))*xi;
       sin(xa(3))*xi;
       tan(u1)*xi/S.L;
       u2];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
