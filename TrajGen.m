function Traj = TrajGen()
    desired = ddp_quad_obst_nl();

    T = 30; %seconds

    s_desired = desired.xs;
    for i = 2:size(s_desired, 2)
        Traj(i).A = poly3_coeff(s_desired(1:6, i-1), s_desired(7:12, i-1), ...
                    s_desired(1:6, i), s_desired(7:12, i), T);
    end
end

function A = poly3_coeff(y0, dy0, yf, dyf, T)
% computes cubic curve connecting (y0,dy0) and (yf, dyf) at time T
    Y = [y0, dy0, yf, dyf];
    L = [poly3(0), dpoly3(0), poly3(T), dpoly3(T)];
    A = Y*inv(L);
end

function f = poly3(t)
    f = [t.^3; t.^2; t; ones(size(t))];
end

function f = dpoly3(t)
    f = [3*t.^2; 2*t; ones(size(t)); zeros(size(t))];
end

function f = d2poly3(t)
    f = [6*t; 2; zeros(size(t)); zeros(size(t))];
end