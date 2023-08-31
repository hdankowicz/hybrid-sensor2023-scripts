function [data, y] = hybridsensor_delay(prob, data, u) %#ok<INUSL>
% Evaluate the hybrid sensor zero problem for an actuator with constant
% delay, given the input signal sigmatilde.
%
% u: [sig2_tilde; u1(0); T; dm; m; k; c; mu; kc; tau]

  function dydt = simulation(t,y,m,k,kc,mu,T,a0_u2,a_u2,b_u2)
    u2 = fourier_func(t,a0_u2,a_u2,b_u2,T);
    dydt = [y(2); mu/m*(1-y(1)^2)*y(2)-(k+kc)/m*y(1)+kc/m*u2];
  end

N   = data.N;
T   = u(2*N+3);
dm  = u(2*N+4);
m   = u(2*N+5);
k   = u(2*N+6);
c   = u(2*N+7);
mu  = u(2*N+8);
kc  = u(2*N+9);
tau = u(2*N+10);

omega = 2*pi/T;

% construct sigma_{2,ss} and u_{2,ss} from sigmatilde
a0_Fl = u(1);
a_Fl = zeros(N,1);
b_Fl = zeros(N,1);
for i=1:N
  a_Fl(i) = u(1+i)*cos(i*omega*tau)-u(N+1+i)*sin(i*omega*tau);
  b_Fl(i) = u(1+i)*sin(i*omega*tau)+u(N+1+i)*cos(i*omega*tau);
end

a0_u2 = a0_Fl/k;
a_u2  = zeros(N,1);
b_u2  = zeros(N,1);
for i=1:N
  a_u2(i) = (a_Fl(i)*(k-i^2*omega^2*(m+dm))-c*i*omega*b_Fl(i))/...
    (c^2*i^2*omega^2+(k-i^2*omega^2*(m+dm))^2);
  b_u2(i) = (a_Fl(i)*c*i*omega+b_Fl(i)*(k-i^2*omega^2*(m+dm)))/...
    (c^2*i^2*omega^2+(k-i^2*omega^2*(m+dm))^2);
end

% simulate for one excitation period
opts = odeset('RelTol',1e-6,'AbsTol',1e-8);
[~,u1] = ode45(@(t,x) simulation(t,x,m,k,kc,mu,T,a0_u2,a_u2,b_u2), ...
  linspace(0,T,20),[u(2*N+2);0],opts);
[a0_u1,a_u1,b_u1] = fourier_coeff(u1(1:end-1,1),N);

y = zeros(3+2*N,1);
y(1:1+2*N) = [a0_Fl-kc*(a0_u1-a0_u2); a_Fl-kc*(a_u1-a_u2); ...
  b_Fl-kc*(b_u1-b_u2)];
y(2+2*N)  = u1(end,1)-u(2*N+2);
y(3+2*N)  = u1(end,2);

end
