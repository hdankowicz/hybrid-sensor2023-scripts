%% Analysis of nonlinear sensor design
%
% Script generating Figure 10 in Y. Mao, H. Dankowicz, "On a Principle for
% Mass Sensing Using Self-Excited Template Dynamics of Coupled Oscillators
% and Root-Finding Algorithms"
%
%% 
clear
eta = 0.01;
gam = 0.01;
kap = 0.001;
pnames = {'eta', 'kap', 'gam', 'del'};
p0 = [eta; kap; gam; 0];
t0 = (0:2*pi/100:2*pi)';
r1 = 2*sqrt(1-kap^2/gam/eta);
r2 = r1*kap/gam;
x0 = [r1*cos(t0) -r1*sin(t0) r2*sin(t0) r2*cos(t0)];

%% po for varying kappa
prob = coco_prob();
prob = coco_set(prob, 'coll', 'NTST', 50);
prob = ode_isol2po(prob, '', @vdp, t0, x0, pnames, p0);
prob = coco_add_event(prob, 'uz', 'kap', 0.001:0.001:0.01);
data = coco_get_func_data(prob, 'po', 'data');
coco(prob, 'test', [], 1, 'kap', [0 0.01]);

bd = coco_bd_read('test');
labs = coco_bd_labs(bd,'uz');

figure
hold on
set(gca,'FontSize',12);
i=1;
for lab=labs
  prob = coco_prob();
  prob = coco_set(prob, 'coll', 'NTST', 50);
  prob = ode_po2po(prob, '', 'test', lab);
  prob = coco_add_event(prob, 'uz', 'del', 0:0.005:0.05);
  data = coco_get_func_data(prob, 'po', 'data');
  coco(prob, 'run', [], 1, 'del', [0 0.05]);
  
  thm = struct('lspec', {{'color',[10-i i i]/10, ...
    'LineWidth', 2}});
  coco_plot_bd(thm, 'run', 'del', {'MAX(x)', 'MIN(x)'}, ...
    @(ma,mi) (ma(1,:)-mi(1,:))./(ma(3,:)-mi(3,:)))
  i = i+1;
end

xlabel('Mass Ratio, $$\delta$$','fontsize',16,'interpreter','latex')
ylabel('Amplitude Ratio $$\alpha$$','fontsize',16,'interpreter','latex')
axis([0 0.05 1 11])
xticks([0 0.01 0.02 0.03 0.04 0.05])
box on
grid on
set(gcf,'position',[0,200,550,450])
hold off

%% oscillator model
function y = vdp(x,p)

u1  = x(1,:);
u1t = x(2,:);
u2  = x(3,:);
u2t = x(4,:);

eta = p(1,:);
kap = p(2,:);
gam = p(3,:);
del = p(4,:);

y = [u1t; eta.*(1-u1.^2).*u1t-(1+kap).*u1+kap.*u2;...
  u2t; -gam./(1+del).*u2t+kap./(1+del).*u1-(1+kap)./(1+del).*u2];
end
