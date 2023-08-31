%% Analysis of active sensor design
%
% Script generating Figure 6 in Y. Mao, H. Dankowicz, "On a Principle for
% Mass Sensing Using Self-Excited Template Dynamics of Coupled Oscillators
% and Root-Finding Algorithms"
%
%% 
clear

%% ep problem along del=0
g0 = 0.02;
k0 = 0.01;
p0 = [k0; 0; g0; -g0+k0];

prob = coco_prob();
prob = coco_add_func(prob, 'hb', @hopfeqs, [], 'zero', ...
  'u0', [p0; 1; 0.1; 0.1; 0]);
prob = coco_add_pars(prob, 'pars', 1:8, ...
  {'kappa', 'delta', 'gamma', 'beta', 'omega', 'u1c', 'u2c', 'u2s'});
prob = coco_add_event(prob, 'GP', 'kappa', 0.005:0.005:0.05);
bd1  = coco(prob, 'alg_run1', [], 1, ...
  {'kappa', 'beta', 'omega', 'u2c', 'u2s'}, [0 0.05]);

GPlab1 = coco_bd_labs(bd1, 'GP');

for i=1:length(GPlab1)
  chart = coco_read_solution('hb', 'alg_run1', GPlab1(i), 'chart');
  prob = coco_prob();
  prob = coco_add_func(prob, 'hb', @hopfeqs, [], 'zero', 'u0', chart.x);
  prob = coco_add_pars(prob, 'pars', 1:8, ...
  {'kappa', 'delta', 'gamma', 'beta', 'omega', 'u1c', 'u2c', 'u2s'});
  prob = coco_add_func(prob, 'stab', @stab, [], 'regular', ...
    'stab', 'uidx', 1:4);
  prob = coco_set(prob, 'cont', 'h_max', 0.005, 'h_min', 0.005, 'PtMX', 500);
  coco(prob, ['alg_run1',num2str(i)], [], 1, ...
    {'delta','beta','omega', 'u2c', 'u2s'}, [0 0.05]);
end

%% ep problem along del=0.05
g0 = 0.02;
k0 = 0.01;
p0 = [k0; 0.05; g0; -g0+k0];

prob = coco_prob();
prob = coco_add_func(prob, 'hb', @hopfeqs, [], 'zero', ...
  'u0', [p0; 1; 0.1; 0.1; 0]);
prob = coco_add_pars(prob, 'pars', 1:8, ...
  {'kappa', 'delta', 'gamma', 'beta', 'omega', 'u1c', 'u2c', 'u2s'});
prob = coco_add_event(prob, 'GP', 'kappa', [0.02 0.025]);
bd1  = coco(prob, 'alg_run2', [], 1, ...
  {'kappa', 'beta', 'omega', 'u2c', 'u2s'}, [0 0.03]);

GPlab2 = coco_bd_labs(bd1, 'GP');

for i=1:length(GPlab2)
  chart = coco_read_solution('hb', 'alg_run2', GPlab2(i), 'chart');
  prob = coco_prob();
  prob = coco_add_func(prob, 'hb', @hopfeqs, [], 'zero', 'u0', chart.x);
  prob = coco_add_pars(prob, 'pars', 1:8, ...
  {'kappa', 'delta', 'gamma', 'beta', 'omega', 'u1c', 'u2c', 'u2s'});
  prob = coco_add_func(prob, 'stab', @stab, [], 'regular', ...
    'stab', 'uidx', 1:4);
  prob = coco_set(prob, 'cont', 'h_max', 0.005, 'h_min', 0.005, 'PtMX', 500);
  coco(prob, ['alg_run2',num2str(i)], [], 1, ...
    {'delta','beta','omega', 'u2c', 'u2s'}, [0 0.05]);
end

%% visualization

% top-left panel
figure
hold on
N1 = numel(GPlab1);
for i=1:N1
  thm = struct();
  thm.lspec = {{'color',[N1-i i i]/N1, 'LineWidth', 1}, ...
    {'color',[N1-i i i]/N1, 'LineWidth', 1, 'LineStyle','--'}};
  thm.ustab = 'stab';
  thm.ustabfun = @(x) (x~=0)+1;
  coco_plot_bd(thm,['alg_run1',num2str(i)],'delta','beta')
end
N2 = numel(GPlab2);
for i=1:N2
  thm = struct();
  thm.lspec = {{'color',[N1-i-3 i+3 i+3]/N1, 'LineWidth', 1}, ...
    {'color',[N1-i-3 i+3 i+3]/N1, 'LineWidth', 1, 'LineStyle','--'}};
  thm.ustab = 'stab';
  thm.ustabfun = @(x) (x~=0)+1;
  coco_plot_bd(thm,['alg_run2',num2str(i)],'delta','beta')
end
xlim([0 0.05])
xticks(0:0.01:0.05)
ylim([-0.05 -0.02])
yticks(-0.05:0.01:-0.02)
set(gca,'FontSize',12);
xlabel('Mass Ratio, $$\delta$$','fontsize',16,'interpreter','latex')
ylabel('$$\beta$$','fontsize',16,'interpreter','latex')
box on
set(gcf,'position',[0,200,550,450])
hold off
grid on

% top-right panel
figure
hold on
N1 = numel(GPlab1);
for i=1:N1
  thm = struct();
  thm.lspec = {{'color',[N1-i i i]/N1, 'LineWidth', 1}, ...
    {'color',[N1-i i i]/N1, 'LineWidth', 1, 'LineStyle','--'}};
  thm.ustab = 'stab';
  thm.ustabfun = @(x) (x~=0)+1;
  coco_plot_bd(thm,['alg_run1',num2str(i)],'delta', ...
    {'beta', 'gamma', 'omega', 'kappa'}, ...
    @(be,ga,om,ka) sqrt((be+ga).^2.*om.^2+(1+ka-om.^2).^2)./ka);
end
N2 = numel(GPlab2);
for i=1:N2
  thm = struct();
  thm.lspec = {{'color',[N1-i-3 i+3 i+3]/N1, 'LineWidth', 1}, ...
    {'color',[N1-i-3 i+3 i+3]/N1, 'LineWidth', 1, 'LineStyle','--'}};
  thm.ustab = 'stab';
  thm.ustabfun = @(x) (x~=0)+1;
  coco_plot_bd(thm,['alg_run2',num2str(i)],'delta', ...
    {'beta', 'gamma', 'omega', 'kappa'}, ...
    @(be,ga,om,ka) sqrt((be+ga).^2.*om.^2+(1+ka-om.^2).^2)./ka);
end
xlim([0 0.05])
xticks(0:0.01:0.05)
ylim([0 3])
yticks(0:1:3)
set(gca,'FontSize',12);
xlabel('Mass Ratio, $$\delta$$','fontsize',16,'interpreter','latex')
ylabel('Amplitude Ratio, $$\sqrt{u_{2,c}^2+u_{2,s}^2}/u_{1,c}$$','fontsize',16,'interpreter','latex')
box on
set(gcf,'position',[0,200,550,450])
hold off
grid on

% bottom panel
figure
hold on
N1 = numel(GPlab1);
for i=1:N1
  thm = struct();
  thm.lspec = {{'color',[N1-i i i]/N1, 'LineWidth', 1}, ...
    {'color',[N1-i i i]/N1, 'LineWidth', 1, 'LineStyle','--'}};
  thm.ustab = 'stab';
  thm.ustabfun = @(x) (x~=0)+1;
  coco_plot_bd(thm,['alg_run1',num2str(i)],'delta', 'omega');
end
N2 = numel(GPlab2);
for i=1:N2
  thm = struct();
  thm.lspec = {{'color',[N1-i i i]/N1, 'LineWidth', 1}, ...
    {'color',[N1-i i i]/N1, 'LineWidth', 1, 'LineStyle','--'}};
  thm.ustab = 'stab';
  thm.ustabfun = @(x) (x~=0)+1;
  coco_plot_bd(thm,['alg_run1',num2str(i)],'delta', 'omega');
end
xlim([0 0.05])
xticks(0:0.01:0.05)
set(gca,'FontSize',12);
xlabel('Mass Ratio, $$\delta$$','fontsize',16,'interpreter','latex')
ylabel('Angular Frequency, $$\omega$$','fontsize',16,'interpreter','latex')
box on
set(gcf,'position',[0,200,550,450])
hold off
grid on

%% hopf conditions
function [data, y] = hopfeqs(prob, data, u) %#ok<INUSL>

kap = u(1);
del = u(2);
gam = u(3);
bet = u(4);
ome = u(5);
u1c = u(6);
u2c = u(7);
u2s = u(8);

y = [-u2c*kap+u1c*(1+kap-ome^2); u2s*kap+u1c*(gam+bet)*ome; ...
  -u1c*kap+u2s*gam*ome+u2c*(1+kap-(1+del)*ome^2); ...
  -u2c*gam*ome-u1c*bet*(1+del)*ome+u2s*(1+kap-(1+del)*ome^2)];

end

%% stability monitor function
function [data, y] = stab(prob, data, u) %#ok<INUSL>

kap = u(1);
del = u(2);
gam = u(3);
bet = u(4);

y = any(real(eig([0 1 0 0; -1-kap -gam-bet kap 0; 0 0 0 1; ...
  kap/(1+del) -bet -(1+kap)/(1+del) -gam/(1+del)]))>1e-5);

end
