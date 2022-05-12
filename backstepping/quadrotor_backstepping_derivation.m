syms x y z a b c p q r real
syms u1 u2 u3 u4 real
syms m Ix Iy Iz g K0 real

syms x_(t) y_(t) z_(t) a_(t) b_(t) c_(t) p_(t) q_(t) r_(t)
syms dx dy dz da db dc dp dq dr real
S1 = [x y z a b c p q r];
S2 = [x_ y_ z_ a_ b_ c_ p_ q_ r_];
S3 = [diff(x_,t) diff(y_,t) diff(z_,t) diff(a_,t) diff(b_,t) diff(c_,t) diff(p_,t) diff(q_,t) diff(r_,t)];
S4 = [dx dy dz da db dc dp dq dr];


Eta = [x y z a b c]';


G = [ 1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 1 0 0 0;
      0 0 0 1 sin(a)*tan(b) cos(a)*tan(b);
      0 0 0 0 cos(a)        -sin(a);
      0 0 0 0 sin(a)/cos(b) cos(a)/cos(b)];

syms P [6,1]

dV0 = Eta' * (G*P);
dV0 = simplify(dV0);

%%
syms K1 real

A = [1 sin(a)*tan(b)   cos(a)*tan(b);
     0 cos(a)         -sin(a);
     0 sin(a)/cos(b)   cos(a)/cos(b)];

P456 = inv(A) * [-K1*a -K1*b -K1*c]';
P456 = simplify(P456);

P = [-K1*x; -K1*y; -K1*z; P456];
JP = jacobian(P, [x, y, z, a, b, c])

