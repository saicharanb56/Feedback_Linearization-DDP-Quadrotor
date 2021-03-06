clear; clc;

T = 20;
S.m = 0.5;
S.Ix = 0.1;
S.Iy = 0.1;
S.Iz = 0.2;
S.g = -9.8;

x0 = [-0.3 -0.2 -0.1 0.1 0.2 0.3 -0.3 -0.2 -0.1 0.1 0.2 0.3];

S.k1 = .1;
S.k2 = 1;

global Record
Record = [];

[ts, xs] = ode45(@(t,x) ode(t,x,S), [0 T], x0);

figure(1); clf;
subplot(2,1,1);
plot(ts,xs(:,1), ts,xs(:,2), ts,xs(:,3), ...
     ts,xs(:,4), ts,xs(:,5), ts,xs(:,6))
legend({'x', 'y', 'z', 'a', 'b', 'c'}, 'Location', 'best')

subplot(2,1,2);
% plot(ts,xs(:,7), ts,xs(:,8), ts,xs(:,9), ...
%      ts,xs(:,10), ts,xs(:,11), ts,xs(:,12))
% legend({'dx', 'dy', 'dz', 'p', 'q', 'r'}, 'Location', 'best')


for i = 1:3
    plot(Record(i,:))
    hold on;
    plot(Record(i+6,:))
end
legend({'1d', '1', '2d', '2', '3d', '3'}, 'Location', 'best')



function [x, y, z, a, b, c, dx, dy, dz, p, q, r] = state(s)

    x=s(1); y=s(2); z=s(3); a=s(4); b=s(5); c=s(6);
    dx=s(7); dy=s(8); dz=s(9); p=s(10); q=s(11); r=s(12);

end

%%%%%%%%%%%% t
function u = control(t, s, S, f, G, fa, Ga)

    [x, y, z, a, b, c, dx, dy, dz, p, q, r] = state(s);
    xi = [dx; dy; dz; p; q; r];

    Phi = [-S.k1 * x; -S.k1 * y; -S.k1 * z; 
           -S.k1 * (a - c*sin(b));
           -S.k1 * (c*sin(a)*cos(b) + b*cos(a));
           -S.k1 * (c*cos(a)*cos(b) - b*sin(a))];
    %%%%%%%%%%%%
    global Record
    Record = [Record [Phi; xi; t]];
    %%%%%%%%%%%%

    Z = xi - Phi;

    JPhi = zeros(6, 6);
    JPhi(1:3, 1:3) = -S.k1 * eye(3);
    JPhi(4:6, 4:6) = [ -S.k1,                                 S.k1*c*cos(b),                       S.k1*sin(b);
                        S.k1*b*sin(a) - S.k1*c*cos(a)*cos(b), S.k1*c*sin(a)*sin(b) - S.k1*cos(a), -S.k1*cos(b)*sin(a);
                        S.k1*b*cos(a) + S.k1*c*cos(b)*sin(a), S.k1*sin(a) + S.k1*c*cos(a)*sin(b), -S.k1*cos(a)*cos(b)];
    
    JV0 = [x y z a b c]';

%     u = pinv(Ga) * ( JPhi' * (f + G*xi) - fa - G'*JV0 - S.k2*Z );
    u = Ga \ ( JPhi' * (f + G*xi) - fa - G'*JV0 - S.k2*Z );

    disp(Ga*pinv(Ga))

end

function ds = ode(t, s, S)

    [x, y, z, a, b, c, dx, dy, dz, p, q, r] = state(s);
    xi = [dx; dy; dz; p; q; r];

    f = zeros(6, 1);

    G = eye(6);
    G(4:6,5:6) = [sin(a)*tan(b)  cos(a)*tan(b);
                  cos(a)        -sin(a);
                  sin(a)/cos(b)  cos(a)/cos(b)];

    fa = [0; 0; S.g;
          (S.Iy-S.Iz)/S.Ix * q * r;
          (S.Iz-S.Ix)/S.Iy * p * r;
          (S.Ix-S.Iy)/S.Iz * p * q];

    Ga = zeros(6, 4);
    Ga(1,1) = -1/S.m * (sin(a)*sin(c) + cos(a)*sin(b)*cos(c));
    Ga(2,1) = -1/S.m * (sin(a)*cos(c) - cos(a)*sin(b)*sin(c));
    Ga(3,1) = -1/S.m * (cos(a) * cos(b));
    Ga(4,2) = 1/S.Ix; Ga(5,3) = 1/S.Iy; Ga(6,4) = 1/S.Iz;
    
    %%%%% test dynamics
%     G = eye(6);
%     fa = zeros(6,1);
%     Ga = [[1; 2] zeros(2,3); eye(4)];
    %%%%% test dynamics
    
    u = control(t, s, S, f, G, fa, Ga);

    %%%%% test control
%     eta = [x; y; z; a; b; c];
%     xi = [dx; dy; dz; p; q; r];
%     u = pinv(Ga) * (-S.k2*xi + S.k1*S.k2*inv(G)*eta - fa - S.k1*inv(G)*eta);
    %%%%% test control

    deta = f + G * xi;
    dxi = fa + Ga * u;
    ds = [deta; dxi];
    
end