%% Analysis of hybrid sensor design
%
% Script generating Figure 15a in Y. Mao, H. Dankowicz, "On a Principle for
% Mass Sensing Using Self-Excited Template Dynamics of Coupled Oscillators
% and Root-Finding Algorithms"
%
%% 
clear
N = 1; % order of highest harmonic
a_sigtilde = zeros(N,1);
b_sigtilde = zeros(N,1);

dm  = 0.01;
m   = 1;
k   = 1;
c   = 0.005;
mu  = 0.01;
kc  = 0.005;
tau = 0.2;

% u0 = [a0_sig2; a_sig2; b_sig2; x(1,1); T; dm; m; k; c; mu; kc; tau];
u0 = [0; a_sigtilde; b_sigtilde; 2; 2*pi; dm; m; k; c; mu; kc; tau];

%% continuation and Newton iterations
data = struct('N', N);
prob = coco_prob();
prob = coco_add_func(prob, 'delay', @hybridsensor_delay, ...
  data, 'zero', 'u0', u0);
prob = coco_add_pars(prob, 'pars', 2*N+4:2*N+10, ...
  {'dm', 'm', 'k', 'c', 'mu', 'kc', 'tau'});
prob = coco_add_event(prob, 'cstarts', 'c', 0.01:0.01:0.1);
prob = coco_set(prob, 'cont', 'h_min', 1e-4, 'NPR', 100);
bdrun = coco(prob, 'run', [], 1, 'c', [0.005 0.105]);

labs = coco_bd_labs(bdrun, 'cstarts');

figure
hold on
set(gca,'FontSize',12);
i = 1;
for lab = labs
  chart = coco_read_solution('delay','run',lab,'chart');
  prob = coco_prob();
  prob = coco_add_func(prob, 'delay', @hybridsensor_delay, ...
    data, 'zero', 'u0', chart.x);
  prob = coco_add_pars(prob, 'pars', 2*N+4:2*N+10, ...
    {'dm', 'm', 'k', 'c', 'mu', 'kc', 'tau'});
  prob = coco_add_slot(prob, '', @bdalpha, data, 'bddat');
  prob = coco_add_event(prob, 'uz', 'dm', 0.0025:0.0025:0.0475);
  prob = coco_set(prob, 'cont', 'h_min', 1e-4);
  bdtest = coco(prob, 'test', [], 1, 'dm', [0 0.05]);
  
  thm = struct('lspec', {{'color',[numel(labs)-i i i]/numel(labs), ...
    'LineWidth', 2}});
  coco_plot_bd(thm, 'test', 'dm', 'alpha')
  labsuz = coco_bd_labs(bdtest, 'uz');
  for labuz = labsuz
    chart = coco_read_solution('delay','test',labuz,'chart');
    u0(2*N+4:2*N+10) = chart.x(2*N+4:2*N+10);
    prob = coco_prob();
    prob = coco_add_func(prob, 'delay', @hybridsensor_delay, ...
      data, 'zero', 'u0', u0);
    prob = coco_add_pars(prob, 'pars', 2*N+4:2*N+10, ...
      {'dm', 'm', 'k', 'c', 'mu', 'kc', 'tau'});
    prob = coco_add_slot(prob, '', @bdalpha, data, 'bddat');
    prob = coco_set(prob, 'corr', 'corrMX', inf, 'ResTOL', 1e-5, ...
      'MaxStep', inf);
    bdconv = coco(prob, 'conv', [], 0);
    al = coco_bd_col(bdconv, 'alpha');
    plot(chart.x(2*N+4),al,'ko')
    drawnow
  end
  i = i+1;
end
xlabel('Mass Ratio, $$\delta$$','fontsize',16,'interpreter','latex')
ylabel('Amplitude Ratio $$\alpha$$','fontsize',16,'interpreter','latex')
xlim([0 0.05])
xticks([0 0.01 0.02 0.03 0.04 0.05])
box on
grid on
set(gcf,'position',[0,200,550,450])
hold off

%% output amplitude ratio
function [data, res] = bdalpha(prob, data, command, varargin) %#ok<INUSL>

res = {};
switch command
  
  case 'init'
    res   = { 'alpha' };
    
  case 'data'
    chart = varargin{1};
    N   = data.N;
    T   = chart.x(2*N+3);
    dm  = chart.x(2*N+4);
    m   = chart.x(2*N+5);
    k   = chart.x(2*N+6);
    c   = chart.x(2*N+7);
    tau = chart.x(2*N+10);
    omega = 2*pi/T;
    a_Fl = zeros(N,1);
    b_Fl = zeros(N,1);
    for i=1:N
      a_Fl(i) = chart.x(1+i)*cos(i*omega*tau)-chart.x(N+1+i)*sin(i*omega*tau);
      b_Fl(i) = chart.x(1+i)*sin(i*omega*tau)+chart.x(N+1+i)*cos(i*omega*tau);
    end
    a_u2  = zeros(N,1);
    b_u2  = zeros(N,1);
    for i=1:N
      a_u2(i) = (a_Fl(i)*(k-i^2*omega^2*(m+dm))-c*i*omega*b_Fl(i))/...
        (c^2*i^2*omega^2+(k-i^2*omega^2*(m+dm))^2);
      b_u2(i) = (a_Fl(i)*c*i*omega+b_Fl(i)*(k-i^2*omega^2*(m+dm)))/...
        (c^2*i^2*omega^2+(k-i^2*omega^2*(m+dm))^2);
    end
    res  = { chart.x(2*N+2)/sqrt(a_u2(1)^2+b_u2(1)^2) };
    
end

end
