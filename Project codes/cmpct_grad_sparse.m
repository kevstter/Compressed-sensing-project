function[x,gx]=cmpct_grad_sparse(s,k,N,kmax,plt)
% Creates a gradient sparse signal of length N and sparsity s that 
% is supported over at most k*s<N indices. 
%
% Inputs:
%   s: sparsity level
%   k: positive integer; k*s is the length of the support of the signal
%   N: signal length
%   kmax: (optional) positive integer; max signal amplitude
%   plt,set as 1 to plot,
% Outputs:
%   x: the gradient sparse signal
%   gx: s-sparse discrete gradient of x. 
%
% Ex:
%   cmpct_grad_sparse(10,8,512);
%
% Created Nov 2015,Kevin Chow
% Last modified Nov 2016,Kevin Chow
%

%% Check inputs
  if nargin<3,error('Requires at least 3 inputs');
  elseif nargin<4, kmax=10; plt=1;
  elseif nargin<5, plt=1;
  elseif nargin>5, error('Too many inputs'); 
  end
  
%% Sensible inputs?  
  if k*s>=N
    error('Desired support length is greater than signal length')
  end

%% Generate grad sparse signal of length N
% First generate grad sparse signal of length k*s
  x=grad_sparse(s-2,k*s,kmax,0);
% Sandwich with zeros  
  x0=randi([0 N-k*s]); x=[0*(1:x0-1) x 0*(1:N-k*s-x0+1)]; x=x(:);
  gx=0*x; gx(1:end-1)=x(2:end)-x(1:end-1);
   
%% Visualize  
  if plt==1
    fignum=800; figure(fignum); clf;
  % Plot signal  
    plot(x,'x-','markersize',4); hold on;
  % Plot gradient
    plot(gx,'ro')
    
    axis([0 N min(min(x),min(gx))-2 max(max(x),max(gx))+2])
    lgd=legend('\bf Signal','\bf Discrete gradient','location','southeast');
    set(lgd,'color','none');
  end
end
    
