function[u,outer]=split_breg_2D(R,f,mu,lambda,nInner,epsilon)
% Perform CS reconstruction on a 2D signal by the split Bregman algorithm.
%   
% Inputs:
%  R: Sampling map. R(n,k)=1, if the freq (n,k) is known, 0 otherwise
%  f: Fourier data. Should have R.*f==f
%  mu: Fidelity term parameter. Best choice depends on how the data is scaled.
%  lambda: Coefficient of the constraint term in the Split Bregman model. 
%     For most problems, it is suggested to use lambda=mu
%  nInner: Number of "inner" loops the Split Bregman method. nInner=30 is
%     suggested to be safe; nInner=5-10 may also be possible.
%  epsilon: noise level/stopping criterion 
% Outputs: 
%   u: recovered data in state space
%   outer: number of outer iterations
%
% References:
% [1] Goldstein and Osher, Split Bregman for l1-regularized problems
% [2] Tom Goldstein's webpage, http://www.cs.umd.edu/~tomg/project/

% Created Dec 2015, Kevin Chow
% Last modified Nov 2016, Kevin Chow
%

%% Preprocessing; storage
  [rows,cols]=size(f);
  f0=f;
  u=zeros(rows,cols);
  x=zeros(rows,cols); y=zeros(rows,cols);
  bx=zeros(rows,cols); by=zeros(rows,cols);

%% Build Kernels
  scale=sqrt(rows*cols);
  murf=ifftn(mu*(conj(R).*f))*scale;
    
  uker=zeros(rows,cols);
  uker(1,1)=4;  uker(1,2)=-1;
  uker(2,1)=-1;   uker(rows,1)=-1;
  uker(1,cols)=-1;
  uker=mu*(conj(R).*R) + lambda*fftn(uker);
  
%% Reconstruction
  outer=0; 
  while norm(R.*fftn(u)/scale-f0,2)/norm(f0,2)>epsilon && outer<1000,
  % Inner loops  
    for inner=1:nInner;
    % Update u   
      rhs=murf+lambda*Dxt(x-bx)+lambda*Dyt(y-by);
      u=ifftn(fftn(rhs)./uker);
    % Update x and y
      dx=Dx(u);
      dy=Dy(u);
      [x,y]=shrink2(dx+bx,dy+by,1/lambda);
    % Update bregman parameters
      bx=bx+dx-x;
      by=by+dy-y;
    end

  % Update Fourier data  
    f=f+f0-R.*fftn(u)/scale;
    murf=ifftn(mu*R.*f)*scale;
    outer=outer+1;
  end
end

function d=Dx(u)
    d=zeros(size(u));
    d(:,2:end)=u(:,2:end) - u(:,1:end-1);
    d(:,1)=u(:,1) - u(:,end);
end

function d=Dxt(u)
    d=zeros(size(u));
    d(:,1:end-1)=u(:,1:end-1) - u(:,2:end);
    d(:,end)=u(:,end) - u(:,1);
end

function d=Dy(u)
    d=zeros(size(u));
    d(2:end,:)=u(2:end,:) - u(1:end-1,:);
    d(1,:)=u(1,:) - u(end,:);
end

function d=Dyt(u)
    d=zeros(size(u));
    d(1:end-1,:)=u(1:end-1,:) - u(2:end,:);
    d(end,:)=u(end,:) - u(1,:);
end

function [xs,ys]=shrink2(x,y,lambda)
    s=sqrt(x.*conj(x)+y.*conj(y));
    ss=s-lambda; ss=ss.*(ss>0);
    s=s+(s<lambda); ss=ss./s;

    xs=ss.*x; ys=ss.*y;
end