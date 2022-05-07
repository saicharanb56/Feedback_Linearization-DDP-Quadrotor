syms t m Ix Iy Iz x(t) y(t) z(t) a(t) b(t) c(t) u1(t) u2(t) u3(t) u4(t) real

J = diag([Ix; Iy; Iz]);

R = ROTZ(c)*ROTY(b)*ROTX(a);

% Dynamics of UAV
xyz_ddot = R*[0; 0; u1]/m + [0; 0; -9.81];
% x_ddot = xyz_ddot(1,:);
% y_ddot = xyz_ddot(2,:);
% z_ddot = xyz_ddot(3,:);

abc_ddot = J\[u2; u3; u4];
% a_ddot = abc_ddot(1,:);
% b_ddot = abc_ddot(2,:);
% c_ddot = abc_ddot(3,:);

xyz_d3dot = diff(xyz_ddot, t);
abc_d3dot = diff(abc_ddot, t);





