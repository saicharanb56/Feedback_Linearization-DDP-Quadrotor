clear; clc;

% initial states
xyz0 = [2 1 0]';
abc0 = [0 0 -1]';
vel0 = [0 0 0]';
pqr0 = [0 0 0]';
dxyzc = [0 0 0 0 0 0 0]';

% desired output
xyzcd = [0 -1 -2 1]';
dxyzcd = [0 0 0 0]';


% A matrix for virtual control input
T = 30;
Y = [[xyz0;abc0(3)] [vel0;pqr0(3)] zeros(4,3) xyzcd dxyzcd zeros(4,3)];
L = [lambda(0,1:5) lambda(T,1:5)];
A1 = Y / L;

pv = 1;
Y = [[xyz0;abc0(3)] [vel0;pqr0(3)] zeros(4,3) xyzcd dxyzcd zeros(4,3)];
L = [lambda(0,1:5) lambda(1,1:5)];
A2 = Y / L;


ys1 = [];
ys2 = [];

ts = linspace(0, T, 200);
for t = ts

    % desired output
    Yd = A1 * lambda(t, 1);
    dYd = A1 * lambda(t, 2);
    d2Yd = A1 * lambda(t, 3);
    d3Yd = A1 * lambda(t, 4);
    d4Yd = A1 * lambda(t, 5);

    ys1 = [ys1 d2Yd];
    
    % desired output
    p = t / T;
    Yd = A2 * lambda(p, 1);
    dYd = A2 * lambda(p, 2);
    d2Yd = A2 * lambda(p, 3);
    d3Yd = A2 * lambda(p, 4);
    d4Yd = A2 * lambda(p, 5);

    ys2 = [ys2 d2Yd];

end

idx = 3;
figure(1); clf;

plot(ts, ys1(idx,:), '-')
hold on;
% plot(ts, ys2(idx,:), '-')


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