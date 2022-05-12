clear; clc;

T = 10;
S.m = 0.5;
S.Ix = 0.1;
S.Iy = 0.1;
S.Iz = 0.1;
S.g = 9.8;

x0 = [3 -2 10 pi/4 -pi/4 pi/6];
S.k0 = 1;
S.k1 = 1;
S.k2 = 1;

[ts, xs] = ode45(@(t,x) ode(t,x,S), [0 T], x0);

figure(1); clf;
plot(ts,xs(:,1), ts,xs(:,2), ts,xs(:,3), ts,xs(:,4), ts,xs(:,5), ts,xs(:,6))

function [x, y, z, a, b, c] = state(s)

    x=s(1); y=s(2); z=s(3); a=s(4); b=s(5); c=s(6);

end

function u = control(s, S)

    [x, y, z, a, b, c] = state(s);

    u = -S.k1 * [x; y; z; 
                 a - c*sin(b);
                 c*sin(a)*cos(b) + b*cos(a);
                 c*cos(a)*cos(b) - b*sin(a)];

end

function ds = ode(t, s, S)

    [x, y, z, a, b, c] = state(s);

    f = zeros(6, 1);

    G = eye(6);
    G(4:6,5:6) = [sin(a)*tan(b) cos(a)*tan(b);
                  cos(a)        -sin(a);
                  sin(a)/cos(b) cos(a)/cos(b)];

    
    u = control(s, S);

    ds = f + G * u;
    
end