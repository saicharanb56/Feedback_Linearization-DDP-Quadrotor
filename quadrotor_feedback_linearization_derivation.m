clear; clc;

% Feedback Linearization
% Original system:
%   s = [x y z a b c dx dy dz p q r]'
%   ds = F + G * U
%
% Linearized system:
%   Y = [x y z c]'
%   Ya = [x y z c u1 du1]'
%   Ua = [du1 u2 u3 u4]'
%   dYa = [f1;f2] + [G1;G2] * Ua


syms x y z a b c dx dy dz p q r real
syms dx dy dz da db dc d2x d2y d2z dp dq dr real
syms u1 u2 u3 u4 real
syms m Ix Iy Iz g real
S1 = [x y z a b c dx dy dz p q r];

F = [dx; dy; dz; 
     p + q*sin(a)*tan(b) + r*cos(a)*tan(b);
     q*cos(a) - r*sin(a);
     q*sin(a)/cos(b) + r*cos(a)/cos(b);
     0; 0; g;
     (Iy-Iz)/Ix * q * r;
     (Iz-Ix)/Iy * p * r;
     (Ix-Iy)/Iz * p * q];

G = sym(zeros(12, 4));
G(7,1) = -1/m * (sin(a)*sin(c) + cos(a)*sin(b)*cos(c));
G(8,1) = -1/m * (sin(a)*cos(c) - cos(a)*sin(b)*sin(c));
G(9,1) = -1/m * (cos(a) * cos(b));
G(10,2) = 1/Ix; G(11,3) = 1/Iy; G(12,4) = 1/Iz;

U = [u1; u2; u3; u4];

syms du1 d2u1 real

Y1 = [x; y; z];
Y2 = c;
Ua = [d2u1; u2; u3; u4];

dY1 = jacobian(Y1, S1) * (F+G*U);
dY1 = simplify(dY1);

d2Y1 = jacobian(dY1, S1) * (F+G*U);
d2Y1 = simplify(d2Y1);

d3Y1 = jacobian(d2Y1, [S1 u1]) * [F+G*U; du1;];
d3Y1 = simplify(d3Y1);

d4Y1 = jacobian(d3Y1, [S1 u1 du1]) * [F+G*U; du1; d2u1];
d4Y1 = simplify(d4Y1);

G1 = jacobian(d4Y1, Ua);
G1 = simplify(G1)
f1 = d4Y1 - G1*Ua;
f1 = simplify(f1)
% check linearity, should be all 0
jacobian(f1, Ua)


dY2 = jacobian(Y2, S1) * (F+G*U);
dY2 = simplify(dY2);

d2Y2 = jacobian(dY2, S1) * (F+G*U);
d2Y2 = simplify(d2Y2);

G2 = jacobian(d2Y2, Ua);
G2 = simplify(G2)
f2 = d2Y2 - G2*Ua;
f2 = simplify(f2)
% check linearity, should be all 0
jacobian(f2, Ua)

