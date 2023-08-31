function varargout=sym_linode(action,varargin)
%% Automatically generated with matlabFunction
%#ok<*DEFNU,*INUSD,*INUSL>

switch action
  case 'nargs'
   varargout{1}=3;
   return
  case 'nout'
   varargout{1}=2;
   return
  case 'argrange'
   varargout{1}=struct('t',1:1,'x',2:3,'p',4:6);
   return
  case 'argsize'
   varargout{1}=struct('t',1,'x',2,'p',3);
   return
  case 'vector'
   varargout{1}=struct('t',0,'x',1,'p',1);
   return
  case 'extension'
   varargout{1}='rhs';
   return
  case 'maxorder'
   varargout{1}=2;
   return
end
nout=2;
order=varargin{1};
f=str2func(sprintf('sym_linode_%s_%d',action,order));
varargout=cell(nout,1);
[varargout{:}]=f(varargin{2:end});




function [out1,out2] = sym_linode_rhs_0(t,x1,x2,k,f,theta,t_dev,x1_dev,x2_dev,k_dev,f_dev,theta_dev)
%SYM_LINODE_RHS_0
%    [OUT1,OUT2] = SYM_LINODE_RHS_0(T,X1,X2,K,F,THETA,T_DEV,X1_DEV,X2_DEV,K_DEV,F_DEV,THETA_DEV)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    10-Apr-2020 22:15:37

out1 = x2;
if nargout > 1
    out2 = -x2-k.*x1+f.*cos(t+theta);
end


function [out1,out2] = sym_linode_rhs_1(t,x1,x2,k,f,theta,t_dev,x1_dev,x2_dev,k_dev,f_dev,theta_dev)
%SYM_LINODE_RHS_1
%    [OUT1,OUT2] = SYM_LINODE_RHS_1(T,X1,X2,K,F,THETA,T_DEV,X1_DEV,X2_DEV,K_DEV,F_DEV,THETA_DEV)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    10-Apr-2020 22:15:38

out1 = x2_dev;
if nargout > 1
    t2 = t+theta;
    out2 = -x2_dev-k.*x1_dev-k_dev.*x1+f_dev.*cos(t2)-f.*sin(t2).*(t_dev+theta_dev);
end


function [out1,out2] = sym_linode_rhs_2(t,x1,x2,k,f,theta,t_dev,x1_dev,x2_dev,k_dev,f_dev,theta_dev)
%SYM_LINODE_RHS_2
%    [OUT1,OUT2] = SYM_LINODE_RHS_2(T,X1,X2,K,F,THETA,T_DEV,X1_DEV,X2_DEV,K_DEV,F_DEV,THETA_DEV)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    10-Apr-2020 22:15:38

out1 = 0.0;
if nargout > 1
    t2 = t_dev+theta_dev;
    t3 = t+theta;
    out2 = k_dev.*x1_dev.*-2.0-f_dev.*t2.*sin(t3).*2.0-f.*t2.^2.*cos(t3);
end
