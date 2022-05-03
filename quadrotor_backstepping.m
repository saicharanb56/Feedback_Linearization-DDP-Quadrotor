clear; clc;

T = 10;
S.m = 0.5;
S.Ix = 0.1;
S.Iy = 0.1;
S.Iz = 0.1;
S.g = -9.8;

x0 = [3 2 1 0.1 0.1 0.1 0 0 0 0 0 0];
S.k0 = 1;
S.k1 = 1;
S.k2 = 1;

[ts, xs] = ode45(@(t,x) ode(t,x,S), [0 T], x0);

plot(ts,xs(:,1), ts,xs(:,2), ts,xs(:,3))

function [x, y, z, a, b, c, dx, dy, dz, p, q, r] = state(s)

    x=s(1); y=s(2); z=s(3); a=s(4); b=s(5); c=s(6);
    dx=s(7); dy=s(8); dz=s(9); p=s(10); q=s(11); r=s(12);

end

function u = control(s, S, f, G, fa, Ga)

    [x, y, z, a, b, c, dx, dy, dz, p, q, r] = state(s);
    xi = [dx; dy; dz; p; q; r];

    Phi = -S.k1 * [x; y; z; 
                   a - c*sin(b);
                   c*sin(a)*cos(b) + b*cos(a);
                   c*cos(a)*cos(b) - b*sin(a)];

    Z = xi - Phi;

    JPhi = zeros(6, 6);
    JPhi(1:3, 1:3) = -S.k1 * eye(3);
    JPhi(4:6, 4:6) = [ -S.k1,                                 S.k1*c*cos(b),                       S.k1*sin(b);
                        S.k1*b*sin(a) - S.k1*c*cos(a)*cos(b), S.k1*c*sin(a)*sin(b) - S.k1*cos(a), -S.k1*cos(b)*sin(a);
                        S.k1*b*cos(a) + S.k1*c*cos(b)*sin(a), S.k1*sin(a) + S.k1*c*cos(a)*sin(b), -S.k1*cos(a)*cos(b)];
    
    JV0 = [x y z a b c]';

    u = pinv(Ga) * ( JPhi' * (f + G*xi) - fa - G'*JV0 - S.k2*Z );

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
    Ga(2,1) = -1/S.m * (cos(a)*sin(c) - cos(a)*sin(b)*sin(c));
    Ga(3,1) = -1/S.m * (cos(a) * cos(b));
    Ga(4,2) = 1/S.Ix; Ga(5,3) = 1/S.Iy; Ga(6,4) = 1/S.Iz;
    
    u = control(s, S, f, G, fa, Ga);

    deta = f + G * xi;
    dxi = fa + Ga * u;
    ds = [deta; dxi];
    
end