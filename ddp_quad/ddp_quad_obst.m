function f = ddp_quad_obst()

% time horizon and segments
tf = 30;
S.N = 1000;
S.h = tf/S.N;

% cost function parameters
S.Q = 0.0*diag([20, 20, 2, ones(1,9)]);
S.R = diag([5, 10, 10, 10]);
S.Qf = diag([2, 2, 2, ones(1,9)]);

S.f = @quad_f;
S.L = @quad_L;
S.Lf = @quad_Lf;
S.mu = 0;

% initial state
x0 = [5; 5; 4; zeros(9,1)];

S.os(1).p = [-1.5;0;3];
S.os(1).r = 1;
S.os(2).p = [-4; -4; 3.5];
S.os(2).r = 1;
S.ko = 10000000000;

% initial control sequence
%us = [repmat([.1;0.1], 1, S.N/2), repmat(-[.1;0.1], 1, S.N/2)]/5;
%us = [repmat([.1;0], 1, N/2), repmat(-[.1;0], 1, N/2)]/5;
us = zeros(4,S.N);

xs = ddp_traj(x0, us, S);

Jmin = ddp_cost(xs, us,  S);

subplot(1,2,1)
plot3(xs(1,:), xs(2,:), xs(3,:), 'b')
hold on

if isfield(S, 'os')
%   da = .1;
%   a = -da:da:2*pi;
%   for i=1:length(S.os)
%     % draw obstacle
%     plot(S.os(i).p(1) + cos(a)*S.os(i).r,  S.os(i).p(2) + sin(a)*S.os(i).r, ...
%          '-r','LineWidth',2);
%   end
%  axis equal
    [sphereX,sphereY,sphereZ] = sphere;
    surf(sphereX*S.os(1).r + S.os(1).p(1), sphereY*S.os(1).r + S.os(1).p(2), sphereZ*S.os(1).r + S.os(1).p(3))
    hold on
    surf(sphereX*S.os(1).r + S.os(2).p(1), sphereY*S.os(2).r + S.os(2).p(2), sphereZ*S.os(2).r + S.os(2).p(3))
    hold on
end

S.a = 1;

for i=1:50
  [dus, V, Vn, dV, a] = ddp(x0, us, S);

  % update controls
  us = us + dus;
  
  S.a = a;   % reuse step-size for efficiency
  
  % update trajectory
  xs = ddp_traj(x0, us, S);
  
  plot3(xs(1,:), xs(2,:), xs(3,:), '-b');
end

plot3(xs(1,:), xs(2,:), xs(3,:), '-g', 'LineWidth', 3);
axis equal
hold off

J = ddp_cost(xs, us, S);

xlabel('x')
ylabel('y')

subplot(1,2,2)
plot(0:S.h:tf-S.h, us(1,:));
hold on
plot(0:S.h:tf-S.h, us(2,:));
hold on
plot(0:S.h:tf-S.h, us(3,:));
hold on
plot(0:S.h:tf-S.h, us(4,:));
hold off
xlabel('sec.')
%yline(0.01); yline(-0.05);
legend('u_1','u_2', 'u_3','u_4')

f.xs = xs;
f.us = us;
end


function [L, Lx, Lxx, Lu, Luu] = quad_L(k, x, u, S)
% quad cost (just standard quadratic cost)

if (k == S.N+1)
  L = x'*S.Qf*x/2;
  Lx = S.Qf*x;
  Lxx = S.Qf;
  Lu = [];
  Luu = [];
else
  L = S.h/2*(x'*S.Q*x + u'*S.R*u);
  Lx = S.h*S.Q*x;
  Lxx = S.h*S.Q;
  Lu = S.h*S.R*u;
  Luu = S.h*S.R;
end

% quadratic penalty term
if isfield(S, 'os')
  for i=1:length(S.os)
    g = x(1:3) - S.os(i).p;
    c = S.os(i).r - norm(g);
    if c < 0
      continue
    end
    
    L = L + S.ko/2*c^2;
    v = g/norm(g);
    Lx(1:3) = Lx(1:3) - S.ko*c*v;
    Lxx(1:3,1:3) = Lxx(1:3,1:3) + S.ko*v*v';  % Gauss-Newton appox
  end
end

end

function xs = ddp_traj(x0, us, S)

N = size(us, 2);
xs(:,1) = x0;

for k=1:N,
  xs(:, k+1) = S.f(k, xs(:,k), us(:,k), S);
end
end

function [x, A, B] = quad_f(k, x, u, S)
% quad dynamics and jacobians

h = S.h;

A = [zeros(6,6), eye(6,6);
     zeros(3,12);
     [0 -9.8; 9.8 0; 0 0], zeros(3,10)];
A = h*A + eye(12);  
   
B = [zeros(6,4);
     0 10 0 0;
     0 0 10 0;
     0 0 0 6.6667;
     0 0 0 0;
     0 0 0 0;
     5 0 0 0];

B = h*B;

x = A*x + B*u;
end