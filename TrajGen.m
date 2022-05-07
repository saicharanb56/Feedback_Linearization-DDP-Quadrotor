function Traj = TrajGen(initState)
    desired = ddp_quad_obst_nl(initState);

    T = 30; %seconds

    s_desired = desired.xs;
    for i = 2:size(s_desired, 2)
        Traj(i).A = poly3_coeff(s_desired(1:6, i-1), s_desired(7:12, i-1), ...
                    s_desired(1:6, i), s_desired(7:12, i), T);
        Traj(i).t0 = (i-1)*T/size(s_desired,2);
        Traj(i).tf = i*T/size(s_desired,2);
    end
    Traj(1).A = poly3_coeff(s_desired(1:6, i-1), s_desired(7:12, i-1), ...
                    s_desired(1:6, i), s_desired(7:12, i), 0);
    Traj(1).t0 = 0;
    Traj(1).tf = T/size(s_desired,2);
end

function A = poly3_coeff(y0, dy0, yf, dyf, T)
% computes cubic curve connecting (y0,dy0) and (yf, dyf) at time T
    Y = [y0, dy0, yf, dyf];
    L = [poly3(0), dpoly3(0), poly3(T), dpoly3(T)];
    A = Y*inv(L);
end

function f = poly3(t)
    f = [t.^10; t.^9; t.^8; t.^7; t.^6; t.^5; t.^4; t.^3; t.^2; t; ones(size(t))];
end

function f = dpoly3(t)
    f = [10*t.^9; 9*t.^8; 8*t.^7; 7*t.^6; 6*t.^5; 5*t.^4; 4*t.^3; 3*t.^2; 2*t; ones(size(t)); zeros(size(t))];
end

function f = d2poly3(t)
    f = [90*t.^8; 72*t.^7; 56*t.^6; 42*t.^5; 30*t.^4; 20*t.^3; 12*t.^2; 6*t; 2; zeros(size(t)); zeros(size(t))];
end

function f = d3poly3(t)
    f = [720*t.^7; 504*t.^6; 336*t.^5; 210*t.^4; 120*t.^3; 60*t.^2; 24*t; 6; zeros(size(t));; zeros(size(t)); zeros(size(t))];
end

function f = d4poly3(t)
    f = [720*t.^7; 504*t.^6; 336*t.^5; 210*t.^4; 120*t.^3; 60*t.^2; 24*t; 6; zeros(size(t));; zeros(size(t)); zeros(size(t))];
end