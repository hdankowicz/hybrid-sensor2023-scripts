%% Analysis of nonlinear sensor design
%
% Script generating Figures 8 and 9 in Y. Mao, H. Dankowicz, "On a
% Principle for Mass Sensing Using Self-Excited Template Dynamics of
% Coupled Oscillators and Root-Finding Algorithms"
%
%% 
clear
etahat = 1;
kaphat = 1/2;
gamhat = 1;

%% po for varying epsilon and delta
for eps = [0.01 0.05 0.1 0.2]
  pnames = {'eta', 'kap', 'gam', 'del'};
  p0 = eps*[etahat; kaphat; gamhat; 0];
  t0 = (0:2*pi/100:2*pi)';
  x0 = [1.73*cos(t0) -1.73*sin(t0) 0.87*sin(t0) 0.87*cos(t0)];
  
  prob = coco_prob();
  prob = coco_set(prob, 'coll', 'NTST', 50);
  prob = ode_isol2po(prob, '', @vdp, t0, x0, pnames, p0);
  prob = coco_add_event(prob, 'uz', 'del', 0:0.005:0.05);
  data = coco_get_func_data(prob, 'po', 'data');
  prob = coco_add_slot(prob, '', @bdcoup, data, 'bddat');
  coco(prob, 'test', [], 1, 'del', [0 0.05]);
  
  % amplitude ratio against mass ratio
  figure
  hold on
  set(gca,'FontSize',12);
  thm = struct('lspec', {{'b-', 'LineWidth', 2}});
  coco_plot_bd(thm, 'test', 'del', {'MAX(x)', 'MIN(x)'}, ...
    @(ma,mi) (ma(1,:)-mi(1,:))./(ma(3,:)-mi(3,:)))
  xlabel('Mass Ratio, $$\delta$$','fontsize',16,'interpreter','latex')
  ylabel('Amplitude Ratio $$\alpha$$','fontsize',16,'interpreter','latex')
  xlim([0 0.05])
  xticks([0 0.01 0.02 0.03 0.04 0.05])
  ylims = get(gca, 'ylim');
  ylim([2 ylims(2)])
  alp = 2:(ylims(2)-2)/50:ylims(2);
  dpl = (1-1./alp.^2).*sqrt(alp.^2*kaphat^2*eps^2-gamhat^2*eps^2);
  plot((1-1./alp.^2).*sqrt(alp.^2*kaphat^2*eps^2-gamhat^2*eps^2), alp, ...
    'r--','linewidth', 2)
  box on
  grid on
  set(gcf,'position',[0,200,550,450])
  hold off
  
  % coupling amplitude against mass ratio
  figure
  hold on
  set(gca,'FontSize',12);
  coco_plot_bd(thm, 'test', 'del', {'MAXDIFF(x)', 'MINDIFF(x)'}, ...
    @(ma,mi) eps*kaphat*(ma(1,:)-mi(1,:))/2)
  set(gca,'FontSize',12);
  xlabel('Mass Ratio, $$\delta$$','fontsize',16,'interpreter','latex')
  ylabel('Amplitude of $$\kappa(u_2-u_1)$$','fontsize',16,'interpreter','latex')
  xlim([0 0.05])
  xticks([0 0.01 0.02 0.03 0.04 0.05])
  cpl = sqrt(4*eps*kaphat./alp.^2.*(1-eps*gamhat./(eps*etahat*alp.^2))...
    .*(2*(1-1./alp.^2).*sqrt(alp.^2*kaphat^2*eps^2-gamhat^2*eps^2).*alp.^2 ...
    +eps*kaphat*(alp.^4-1))./(alp.^2-1));
  plot(dpl, cpl, 'r--','linewidth', 2)
  box on
  grid on
  set(gcf,'position',[0,200,550,450])
  hold off
end

%% oscillator model
function dydt = vdp(y,p)

u1  = y(1,:);
u1t = y(2,:);
u2  = y(3,:);
u2t = y(4,:);

eta = p(1,:);
kap = p(2,:);
gam = p(3,:);
del = p(4,:);

dydt = [u1t; eta.*(1-u1.^2).*u1t-(1+kap).*u1+kap.*u2;...
  u2t; -gam./(1+del).*u2t+kap./(1+del).*u1-(1+kap)./(1+del).*u2];
end

%% output coupling difference
function [data, res] = bdcoup(prob, data, command, varargin)

res = {};
switch command
  
  case 'init'
    xid   = data.po_bddat.xid;
    maxid = sprintf('MAXDIFF(%s)', xid);
    minid = sprintf('MINDIFF(%s)', xid);
    res   = { maxid, minid };
    
  case 'data'
    chart = varargin{1};
    [fdata, uidx] = coco_get_func_data(prob, data.cid, 'data', 'uidx');
    maps = fdata.coll_seg.maps;
    xbp  = reshape(chart.x(uidx(maps.xbp_idx)), maps.xbp_shp);
    res  = { max(xbp(1,:)-xbp(3,:)) , min(xbp(1,:)-xbp(3,:)) };
    
end

end
