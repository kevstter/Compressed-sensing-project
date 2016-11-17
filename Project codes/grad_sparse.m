function[x,gx,pos]=grad_sparse(s,N,kmax,plt)
% Creates a gradient sparse signal of length N with sparsity s.
%
% Inputs:
%   s: sparsity level
%   N: signal length.
%   kmax: positive integer; max signal amplitude
%   plt: set as 1 to plot,
% Outputs:
%   x: gradient sparse signal
%   gx: s-sparse discrete gradient of x.
%
% Ex:
%   grad_sparse(10,512);
%
% Created Nov 2015, Kevin Chow
% Last modifed Nov 2016, Kevin Chow
%

%% Check inputs
  if nargin < 2, error('Require at least 2 inputs');
  elseif nargin < 3, kmax=10; plt=1;
  elseif nargin < 4, plt=1;
  elseif nargin > 4, error('Too many inputs');
  end

%% Sanity  
  if s >= N,
    error('Desired sparsity level is greater than signal length')
  end

%% Construct signal  
  x=zeros(1,N);
  pos=[1 randsample(N,s)' N]; pos=sort(pos);
  val=randi([-kmax,kmax],[1,s+1]);
    
  for k=1:s+1,
    x(pos(k):pos(k+1))=val(k); 
  end
  gx=0*x; gx(1:end-1)=x(2:end)-x(1:end-1); 

%% Visualize  
  if plt == 1
    fignum=801; figure(fignum); clf
    plot(x,'x-','markersize',4); hold on ;
    plot(gx,'ro')
    
    axis([0 N min(min(x),min(gx))-2 max(max(x),max(gx))+2])
    lgd=legend('\bf Signal','\bf Discrete gradient','location','southeast');
    set(lgd,'color','none');
  end
end
    
