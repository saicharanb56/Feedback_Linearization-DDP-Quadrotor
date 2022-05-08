clear; clc;
xyz0 = [-5 -5 3.5]';
T  = 10;
S.T = T;
desired = ddp_quad_obst_nl(xyz0, T);
clf;

S.yd = desired.xs;
S.ydDot1 = diff(S.yd, 1, 2)*T/size(S.yd,2);
S.ydDot1 = [zeros(size(S.yd,1),1) S.ydDot1];
S.ydDot2 = diff(S.ydDot1, 1, 2)*T/size(S.ydDot1,2);
S.ydDot2 = [zeros(size(S.ydDot2,1),1) S.ydDot2];
S.ydDot3 = diff(S.ydDot2, 1, 2)*T/size(S.ydDot2,2);
S.ydDot3 = [zeros(size(S.ydDot3,1),1) S.ydDot3];
S.ydDot4 = diff(S.ydDot3, 1, 2)*T/size(S.ydDot3,2);
S.ydDot4 = [zeros(size(S.ydDot4,1),1) S.ydDot4];

Yds = [];

ts = linspace(0, S.T, 333);
% for t = ts
% 
%     traj_tstep = S.T/size(S.yd,2);
%     idx = floor(t/traj_tstep) + 1;
%     idx = min(idx, size(S.yd,2)-1);
% 
%     oi = [1 2 3 6];
% 
%     Y = [S.yd(oi,idx) S.yd(oi,idx+1)];
% 
%     L = [lambda(0.5,1) lambda(1.5,1)];
%     A = Y / L;
% 
%     if t ~= S.T
%         percent = t/traj_tstep - (idx - 1) + 0.5;
%     else
%         percent = 1 + 0.5;
%     end
% 
%     Yd = A * lambda(percent, 1);
%     dYd = A * lambda(percent, 2);
%     d2Yd = A * lambda(percent, 3);
%     d3Yd = A * lambda(percent, 4);
%     d4Yd = A * lambda(percent, 5);
% 
% 
%     Yds = [Yds Yd];
% 
% end
% ys = Yds

idx = 1;

yd = S.yd;
traj_ts = linspace(0, T, size(S.yd,2));
p = polyfit(traj_ts,yd(idx,:),10);

ts = linspace(0, T, 200);
ys = polyval(p, ts);


figure(1); clf;

plot(traj_ts, S.yd(idx,:), '*')
hold on;
plot(ts, ys(1,:), '.-')


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