%% Analysis of active sensor design
%
% Script generating Figure 4 in Y. Mao, H. Dankowicz, "On a Principle for
% Mass Sensing Using Self-Excited Template Dynamics of Coupled Oscillators
% and Root-Finding Algorithms"
%
%%
clear
x0 = zeros(4,1);
g0 = 0.02;
k0 = 0.0107;
pnames = {'kappa','delta','gamma', 'beta'};
p0     =  [k0; 0; g0; -g0+0.01];

%% ep problem along u=0
prob = coco_prob();
prob = ode_isol2ep(prob, '', @linsensor, x0, pnames, p0);
bd1  = coco(prob, 'ep_run', [], 1, 'beta', [-1 1]);

HBlab = coco_bd_labs(bd1, 'HB');
HBlab = max(HBlab);

%% hb2hb
prob = coco_prob();
prob = ode_HB2HB(prob, '', 'ep_run', HBlab);
bd2  = coco(prob, 'ep_run_HB', [], 1, {'delta','beta'}, [0 0.05]);

%% alg
prob = coco_prob();
prob = coco_add_func(prob, 'hb', @hopfeqs, [], 'zero', 'u0', [p0; 1]);
prob = coco_add_pars(prob, 'pars', 1:5, ...
  {'kappa', 'delta', 'gamma', 'beta', 'omega'});
prob = coco_set(prob, 'cont', 'h_max', 0.005, 'h_min', 0.005);
bd3  = coco(prob, 'alg_run', [], 1, {'delta', 'beta', 'omega'}, [0 0.05]);

%% visualization
figure
hold on
coco_plot_bd('ep_run_HB','delta','beta')
thm.lspec = {'bo', 'LineWidth', 2};
coco_plot_bd(thm,'alg_run','delta','beta')
xlim([0 0.05])
xticks(0:0.01:0.05)
set(gca,'FontSize',12);
xlabel('Mass Ratio, $$\delta$$','fontsize',16,'interpreter','latex')
ylabel('$$\beta$$','fontsize',16,'interpreter','latex')
box on
set(gcf,'position',[0,200,550,450])
hold off
grid on

%% algebraic hopf conditions
function [data, y] = hopfeqs(prob, data, u) %#ok<INUSL>

kap = u(1);
del = u(2);
gam = u(3);
bet = u(4);
ome = u(5);

y = [gam*((2+del)*ome^2-2-2*kap)-bet*(1+(2+del)*kap-(1+del)*ome^2);
  bet+gam-(1-ome^2)*(1+2*kap-ome^2)/gam/ome^2+(1+kap-ome^2)*del/gam];

end

%% linear sensor model
function y = linsensor(x, p)

x1 = x(1,:);
x2 = x(2,:);
x3 = x(3,:);
x4 = x(4,:);

kappa = p(1,:);
delta = p(2,:);
gamma = p(3,:);
beta  = p(4,:);

y(1,:) = x2;
y(2,:) = -(gamma+beta).*x2-(1+kappa).*x1+kappa.*x3;
y(3,:) = x4;
y(4,:) = -beta.*x2-gamma./(1+delta).*x4+kappa./(1+delta).*x1-(1+kappa)./(1+delta).*x3;

end
