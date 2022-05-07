clear; clc;
% Feedback Linearization Control on quadrotor
% Original system:
%   s = [x y z a b c dx dy dz p q r]'
%   ds = F + G * U
%
% Linearized system:
%   Y = [x y z c]'
%   Ua = [d2u1 u2 u3 u4]'
%   d4Y1 = F1 + G1 * Ua
%   d2Y2 = F2 + G2 * Ua
%
% Full state:
%   sa = [x y z a b c dx dy dz p q r
%         d2x d3x d2y d3y d2z d3z dc u1 du1]'
%

% quadrotor parameters
S.m = 0.1;
S.Ix = 0.1; S.Iy = 0.1; S.Iz = 0.2;
S.g = -9.8;

% simulation time
T = 10;

% initial states
xyz0 = [2 1 0]';
abc0 = [0 0 -1]';
vel0 = [0 0 0]';
pqr0 = [0 0 0]';
dxyzc = [0 0 0 0 0 0 0]';
u1du1 = [S.m*S.g 0]';   % assume u1 compensates gravatiy at t=0
sa0 = [xyz0; abc0; vel0; pqr0; dxyzc; u1du1];

% desired output
xyzcd = [0 -1 -2 1]';
dxyzcd = [0 0 0 0]';


% A matrix for virtual control input
Y = [[xyz0;abc0(3)] [vel0;pqr0(3)] zeros(4,3) xyzcd dxyzcd zeros(4,3)];
L = [lambda(0,1:5) lambda(T,1:5)];
S.A = Y / L;


% controller
S.k0 = 1; S.k1 = 1; S.k2 = 1; S.k3 = 1; S.k4= 1; S.k5 = 1; 

[ts, sas] = ode45(@(t,sa) quadrotor_ode(t,sa,S), [0 T], sa0);

figure(1); clf; box on; grid on; hold on;
plot(ts, sas(:, 1), 'LineWidth',2);
plot(ts, sas(:, 2), 'LineWidth',2);
plot(ts, sas(:, 3), 'LineWidth',2);
plot(ts, sas(:, 6), 'LineWidth',2);
legend({'x (m)', 'y (m)', 'z (m)', 'c (rad)'}, 'Location','best')


function dsa = quadrotor_ode(t, sa, S)
% Quadrotor dynamics

    [x, y, z, a, b, c, dx, dy, dz, p, q, r, ... 
        d2x, d3x, d2y, d3y, d2z, d3z, dc, u1, du1] = get_state(sa);
    Ix = S.Ix; Iy = S.Iy; Iz = S.Iz; m = S.m; g = S.g;
    
    F = [ dx; dy; dz; 
          p + q*sin(a)*tan(b) + r*cos(a)*tan(b);
          q*cos(a) - r*sin(a);
          q*sin(a)/cos(b) + r*cos(a)/cos(b);
          0; 0; g;
          (Iy-Iz)/Ix * q * r;
          (Iz-Ix)/Iy * p * r;
          (Ix-Iy)/Iz * p * q];
    
    G = zeros(12, 4);
    G(7,1) = -1/m * (sin(a)*sin(c) + cos(a)*sin(b)*cos(c));
    G(8,1) = -1/m * (sin(a)*cos(c) - cos(a)*sin(b)*sin(c));
    G(9,1) = -1/m * (cos(a) * cos(b));
    G(10,2) = 1/Ix; G(11,3) = 1/Iy; G(12,4) = 1/Iz;
    
    
    % controls
    Ua = fl_control(t, sa, S);
    U = [u1; Ua(2:4)];
    d2u1 = Ua(1);
    

    % dynamcis
    ds = F + G * U;
    
    d4Yd = S.A * lambda(t, 5);
    d2Yd = S.A * lambda(t, 3);
    dsa = [ds; d3x; d4Yd(1); d3y; d4Yd(2); d3z; d4Yd(3); d2Yd(4); du1; d2u1];

end

function Ua = fl_control(t, sa, S)

    [x, y, z, a, b, c, dx, dy, dz, p, q, r, ... 
        d2x, d3x, d2y, d3y, d2z, d3z, dc, u1, du1] = get_state(sa);
    Ix = S.Ix; Iy = S.Iy; Iz = S.Iz; m = S.m; g = S.g;
    
    F1 = [ ...
    (Ix*Iy*p^2*u1*sin(a)*sin(c) - 2*Ix*Iy*du1*p*cos(a)*sin(c) - 2*Ix*Iy*du1*q*cos(b)*cos(c) + Ix*Iy*q^2*u1*sin(a)*sin(c) + Ix^2*p*r*u1*cos(b)*cos(c) - Iy^2*q*r*u1*cos(a)*sin(c) + 2*Ix*Iy*du1*p*cos(c)*sin(a)*sin(b) + Ix*Iy*p^2*u1*cos(a)*cos(c)*sin(b) + Ix*Iy*q^2*u1*cos(a)*cos(c)*sin(b) + Iy^2*q*r*u1*cos(c)*sin(a)*sin(b) - Ix*Iy*p*r*u1*cos(b)*cos(c) - Ix*Iz*p*r*u1*cos(b)*cos(c) + Ix*Iy*q*r*u1*cos(a)*sin(c) + Iy*Iz*q*r*u1*cos(a)*sin(c) - Ix*Iy*q*r*u1*cos(c)*sin(a)*sin(b) - Iy*Iz*q*r*u1*cos(c)*sin(a)*sin(b))/(Ix*Iy*m);
    (2*Ix*Iy*du1*q*cos(b)*sin(c) - 2*Ix*Iy*du1*p*cos(a)*cos(c) + Ix*Iy*p^2*u1*cos(c)*sin(a) + Ix*Iy*q^2*u1*cos(c)*sin(a) - Iy^2*q*r*u1*cos(a)*cos(c) - Ix^2*p*r*u1*cos(b)*sin(c) - 2*Ix*Iy*du1*p*sin(a)*sin(b)*sin(c) - Ix*Iy*p^2*u1*cos(a)*sin(b)*sin(c) - Ix*Iy*q^2*u1*cos(a)*sin(b)*sin(c) + Ix*Iy*q*r*u1*cos(a)*cos(c) + Iy*Iz*q*r*u1*cos(a)*cos(c) - Iy^2*q*r*u1*sin(a)*sin(b)*sin(c) + Ix*Iy*p*r*u1*cos(b)*sin(c) + Ix*Iz*p*r*u1*cos(b)*sin(c) + Ix*Iy*q*r*u1*sin(a)*sin(b)*sin(c) + Iy*Iz*q*r*u1*sin(a)*sin(b)*sin(c))/(Ix*Iy*m);
                                                                                                                                                                                                                                                        (2*Ix*Iy*du1*q*sin(b) - Ix^2*p*r*u1*sin(b) + Ix*Iy*p*r*u1*sin(b) + Ix*Iz*p*r*u1*sin(b) + 2*Ix*Iy*du1*p*cos(b)*sin(a) + Ix*Iy*p^2*u1*cos(a)*cos(b) + Ix*Iy*q^2*u1*cos(a)*cos(b) + Iy^2*q*r*u1*cos(b)*sin(a) - Ix*Iy*q*r*u1*cos(b)*sin(a) - Iy*Iz*q*r*u1*cos(b)*sin(a))/(Ix*Iy*m);
    ];
    
    F2 = ...
    -(Iy^2*p*q*cos(a)*cos(b) - Iz^2*p*r*cos(b)*sin(a) + 2*Iy*Iz*q*r*sin(b) - 2*Iy*Iz*q^2*cos(a)*sin(a)*sin(b) + 2*Iy*Iz*r^2*cos(a)*sin(a)*sin(b) - Ix*Iy*p*q*cos(a)*cos(b) - Iy*Iz*p*q*cos(a)*cos(b) + Ix*Iz*p*r*cos(b)*sin(a) + Iy*Iz*p*r*cos(b)*sin(a) - 4*Iy*Iz*q*r*cos(a)^2*sin(b))/(Iy*Iz*cos(b)^2);
    
    G1 = [...
    [-(sin(a)*sin(c) + cos(a)*cos(c)*sin(b))/m, -(u1*(cos(a)*sin(c) - cos(c)*sin(a)*sin(b)))/(Ix*m), -(u1*cos(b)*cos(c))/(Iy*m), 0]
    [-(cos(c)*sin(a) - cos(a)*sin(b)*sin(c))/m, -(u1*(cos(a)*cos(c) + sin(a)*sin(b)*sin(c)))/(Ix*m),  (u1*cos(b)*sin(c))/(Iy*m), 0]
    [                       -(cos(a)*cos(b))/m,                           (u1*cos(b)*sin(a))/(Ix*m),         (u1*sin(b))/(Iy*m), 0]
    ];
    
    G2 = ...
    [0, 0, sin(a)/(Iy*cos(b)), cos(a)/(Iz*cos(b))];
    
    % current output
    Y1 = [x; y; z];
    dY1 = [dx; dy; dz];
    d2Y1 = [d2x; d2y; d2z];
    d3Y1 = [d3x; d3y; d3z];
    Y2 = c;
    dY2 = dc;
    
    % desired output
    Yd = S.A * lambda(t, 1);
    dYd = S.A * lambda(t, 2);
    d2Yd = S.A * lambda(t, 3);
    d3Yd = S.A * lambda(t, 4);
    d4Yd = S.A * lambda(t, 5);

    % virtural control
    V1 = d4Yd(1:3) - S.k0*(Y1-Yd(1:3)) - S.k1*(dY1-dYd(1:3)) - S.k2*(d2Y1-d2Yd(1:3)) - S.k3*(d3Y1-d3Yd(1:3));
    V2 = d2Yd(4) - S.k4*(Y2-Yd(4)) - S.k5*(dY2-dYd(4));
    
    Ua = ([G1; G2]) \ ( [V1; V2] - [F1; F2] );

end

function [x, y, z, a, b, c, dx, dy, dz, p, q, r, ... 
          d2x, d3x, d2y, d3y, d2z, d3z, dc, u1, du1] = get_state(sa)
    
    % original dynamics
    x = sa(1); y = sa(2); z = sa(3);
    a = sa(4); b = sa(5); c = sa(6);
    dx = sa(7); dy = sa(8); dz = sa(9);
    p = sa(10); q = sa(11); r = sa(12);

    % input-output linearization
    d2x = sa(13); d3x = sa(14);
    d2y = sa(15); d3y = sa(16);
    d2z = sa(17); d3z = sa(18);
    dc = sa(19);
    u1 = sa(20); du1 = sa(21);

end

function Li = lambda(t, i)
    
    L = [...
    [   1,      0,      0,       0,        0]
    [   t,      1,      0,       0,        0]
    [ t^2,    2*t,      2,       0,        0]
    [ t^3,  3*t^2,    6*t,       6,        0]
    [ t^4,  4*t^3, 12*t^2,    24*t,       24]
    [ t^5,  5*t^4, 20*t^3,  60*t^2,    120*t]
    [ t^6,  6*t^5, 30*t^4, 120*t^3,  360*t^2]
    [ t^7,  7*t^6, 42*t^5, 210*t^4,  840*t^3]
    [ t^8,  8*t^7, 56*t^6, 336*t^5, 1680*t^4]
    [ t^9,  9*t^8, 72*t^7, 504*t^6, 3024*t^5]
    [t^10, 10*t^9, 90*t^8, 720*t^7, 5040*t^6]
    ];
    Li = L(:, i);

end

