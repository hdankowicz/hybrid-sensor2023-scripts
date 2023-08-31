%% Analysis of passive sensor design
%
% Script generating Figure 2 in Y. Mao, H. Dankowicz, "On a Principle for
% Mass Sensing Using Self-Excited Template Dynamics of Coupled Oscillators
% and Root-Finding Algorithms"
%
%%
clear
k = 0.01:0.01:0.1;

%% left panel
figure
hold on
for i=1:10
    gamma = 0.02;
    kappa = k(i);
    delta = 0:0.001:0.05;
    ratio = abs(sqrt(1+delta+delta.^2/4/kappa^2*(1+kappa)^2)-delta/2/kappa*(1+kappa));
    plot(delta,ratio,'color',[10-i i i]/10,'linewidth',1)
end
hold off
xlim([0 0.05])
xticks(0:0.01:0.05)
set(gca,'FontSize',12);
xlabel('Mass Ratio, $$\delta$$','fontsize',16,'interpreter','latex')
ylabel('Amplitude Ratio','fontsize',16,'interpreter','latex')
box on
grid on
set(gcf,'position',[0,200,550,450])

%% right panel
figure
hold on
for i=1:10
    gamma = 0.02;
    kappa = k(i);
    delta = 0:0.001:0.05;
    ratio = abs(-sqrt(1+delta+delta.^2/4/kappa^2*(1+kappa)^2)-delta/2/kappa*(1+kappa));
    plot(delta,ratio,'color',[10-i i i]/10,'linewidth',1)
end
hold off
xlim([0 0.05])
xticks(0:0.01:0.05)
set(gca,'FontSize',12);
xlabel('Mass Ratio, $$\delta$$','fontsize',16,'interpreter','latex')
ylabel('Amplitude Ratio','fontsize',16,'interpreter','latex')
box on
grid on
set(gcf,'position',[0,200,550,450])
