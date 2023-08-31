%% Analysis of hybrid sensor design
%
% Script generating Figure 12 in Y. Mao, H. Dankowicz, "On a Principle for
% Mass Sensing Using Self-Excited Template Dynamics of Coupled Oscillators
% and Root-Finding Algorithms"
%
%% 
clear
N = 1; % order of highest harmonic
a_sigtilde = zeros(N,1);
b_sigtilde = zeros(N,1);

dm = 0.01;
m  = 1;
k  = 1;
c  = 0.01;
mu = c;
kc = 0.005;
ma = 1.28;
ka = 2.86;
ba = 0.25;
kp = 2;
kd = 2;

% u0 = [a0_sig2; a_sig2; b_sig2; x(1,1); T; dm; m; k; c; mu; kc; ma; ka; ba; kp; kd];
u0 = [0; a_sigtilde; b_sigtilde; 2; 2*pi; dm; m; k; c; mu; kc; ma; ka; ba; kp; kd];

%% Template constraints
prob = coco_prob();
prob = coco_add_func(prob, 'delay', @hybridsensor_PD, struct('N', N), ...
  'zero', 'u0', u0);
prob = coco_add_pars(prob, 'pars', 2*N+4:2*N+14, ...
  {'dm', 'm', 'k', 'c', 'mu', 'kc', 'ma', 'ka', 'ba', 'kp', 'kd'});
prob = coco_set(prob, 'corr', 'corrMX', inf, 'ResTOL', 1e-5, ...
  'MaxStep', inf);
prob = add_recording(prob);
bd = coco(prob, 'run', [], 0);

data = coco_read_solution('rec','run',1);
zhist = cell2mat(data.history(1:end-1,end)')';

%% visualization

% top-left panel
figure
hold on
plotopts = {'b--','o:','*--','k','r-.','g*','m--','kd','r:'};
lgnd = {'init.','1st iter.','2nd iter.','3rd iter.','4th iter.', ...
  '5th iter.','6th iter.','7th iter.'};
for i=1:min(size(zhist,1),9)
  t = linspace(0,zhist(i,2*N+3),30);
  plot(t,fourier_func(t,zhist(i,1),zhist(i,2:N+1),zhist(i,N+2:2*N+1),...
    t(end)),plotopts{i},'linewidth',1)
end
set(gca,'FontSize',12);
xlabel('$$t$$','fontsize',16,'interpreter','latex')
ylabel('$$\tilde{\sigma}(t)$$','fontsize',16,'interpreter','latex')
box on
grid on
hold off
xlim([0 max(zhist(:,2*N+3))])
set(gcf,'position',[0,200,550,450])
legend(lgnd{1:min(size(zhist,1),9)},'FontSize',14,'Location','southeast')

ax=axes;
set(ax,'units','normalized','position',[0.45,0.68,0.2,0.2])
box(ax,'on')
hold on
for i=1:min(size(zhist,1),9)
  t = linspace(0,zhist(i,2*N+3),30);
  plot(t,fourier_func(t,zhist(i,1),zhist(i,2:N+1),zhist(i,N+2:2*N+1),...
    t(end)),plotopts{i},'linewidth',1,'parent',ax)
end
grid on
hold off
set(ax,'xlim',[6.25 6.29],'ylim',[-0.05,0.45])
xticks(6.25:0.01:6.29)

% top-right panel
figure
al = zeros(1,size(zhist,1));
for j=1:size(zhist,1)
  T = zhist(j,2*N+3);
  omega = 2*pi/T;
  a0_Fl = zhist(j,1)/(1-1/kp);
  a_Fl = zeros(N,1);
  b_Fl = zeros(N,1);
  for i=1:N
    Fl = (zhist(j,1+i)+1j*zhist(j,N+1+i))/...
      (1-1/(kp-1j*kd*i*omega)-(ma*i*omega+1j*ba)*i*...
      omega/(kp-1j*kd*i*omega)*(1/ka+1/(k-(m+dm)*i^2*omega^2-1j*c*i*omega)));
    a_Fl(i) = real(Fl);
    b_Fl(i) = imag(Fl);
  end
  a0_u2 = a0_Fl/k;
  a_u2  = zeros(N,1);
  b_u2  = zeros(N,1);
  for i=1:N
    u2comp = (a_Fl(i)+1j*b_Fl(i))/(k-(m+dm)*i^2*omega^2-1j*c*i*omega);
    a_u2(i) = real(u2comp);
    b_u2(i) = imag(u2comp);
  end
  al(j) = zhist(j,2*N+2)/sqrt(a_u2(1)^2+b_u2(1)^2);
end
plot(0:size(zhist,1)-1,al,'linewidth',1)
set(gca,'FontSize',12);
xlabel('Iteration','fontsize',16,'interpreter','latex')
ylabel('$$\alpha$$','fontsize',16,'interpreter','latex')
box on
grid on
set(gcf,'position',[0,200,550,450])
hold off
xlim([1 3])
xticks([1 2 3])

% bottom-left panel
plotopts = {'b--','o:','*--','k','r-.','g*','m--','kd','r:'};
lgnd = {'init.','1st iter.','2nd iter.','3rd iter.','4th iter.', ...
  '5th iter.','6th iter.','7th iter.'};

figure
hold on
for j=1:size(zhist,1)
  T = zhist(j,2*N+3);
  omega = 2*pi/T;
  a0_Fl = zhist(j,1)/(1-1/kp);
  a_Fl = zeros(N,1);
  b_Fl = zeros(N,1);
  for i=1:N
    Fl = (zhist(j,1+i)+1j*zhist(j,N+1+i))...
      /(1-1/(kp-1j*kd*i*omega)-(ma*i*omega+1j*ba)*i*...
      omega/(kp-1j*kd*i*omega)*(1/ka+1/(k-(m+dm)*i^2*omega^2-1j*c*i*omega)));
    a_Fl(i) = real(Fl);
    b_Fl(i) = imag(Fl);
  end
  a0_u2 = a0_Fl/k;
  a_u2  = zeros(N,1);
  b_u2  = zeros(N,1);
  for i=1:N
    u2comp = (a_Fl(i)+1j*b_Fl(i))/(k-(m+dm)*i^2*omega^2-1j*c*i*omega);
    a_u2(i) = real(u2comp);
    b_u2(i) = imag(u2comp);
  end
  t = linspace(0,zhist(i,2*N+3),30);
  opts = odeset('RelTol',1e-6,'AbsTol',1e-8);
  [~,u1] = ode45(@(t,x) simulation(t,x,m,k,kc,mu,T,a0_u2,a_u2,b_u2), ...
    t,[zhist(j,2*N+2);0],opts);
  plot(t,kc*(u1(:,1)'-fourier_func(t,a0_u2,a_u2,b_u2,T)),plotopts{j},'linewidth',1)
end
hold off
set(gca,'FontSize',12);
xlabel('$$t$$','fontsize',16,'interpreter','latex')
ylabel('$$k_c(u_1-u_2)$$','fontsize',16,'interpreter','latex')
box on
grid on
set(gcf,'position',[0,200,550,450])
xlim([0 max(zhist(:,2*N+3))])
legend(lgnd{1:min(size(zhist,1),9)},'FontSize',14,'Location','southeast')

% bottom-right panel
figure
tol = zeros(1,size(zhist,1));
for i=1:size(zhist,1)
  [~, y] = hybridsensor_PD('',struct('N',N),zhist(i,:)');
  tol(i) = norm(y);
end
semilogy(0:size(zhist,1)-1,tol,'linewidth',1)

set(gca,'FontSize',12);
xlabel('Iteration','fontsize',16,'interpreter','latex')
ylabel('Error','fontsize',16,'interpreter','latex')
box on
grid on
set(gcf,'position',[0,200,550,450])
hold off
xlim([0 3])
xticks([0 1 2 3])

%% in silico model
function dydt = simulation(t,y,m,k,kc,mu,T,a0_u2,a_u2,b_u2)
u2 = fourier_func(t,a0_u2,a_u2,b_u2,T);
dydt = [y(2); mu/m*(1-y(1)^2)*y(2)-(k+kc)/m*y(1)+kc/m*u2];
end
