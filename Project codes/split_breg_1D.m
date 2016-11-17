function[u,outer]=split_breg_1D(R,f,mu,lambda,nInner,epsilon)
% Perform CS reconstruction on a 1D signal by the split Bregman algorithm.
%   
% Inputs:
%  R: Sampling map. R(k)=1 if the frequency k is known,and 0 otherwise.
%  f: Fourier data. Should have R.*f==f;
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
%
% Modified from mrics.m, Dec 2015, Kevin Chow
% Last modified Nov 2016, Kevin Chow
%

%% Preprocessing; storage
  [rows,cols]=size(f);
  f0=f;
  u=zeros(rows,cols);
  x=zeros(rows,cols);
  bx=zeros(rows,cols);
    
%% Build kernels
  scale=sqrt(rows*cols);
  murf=ifft(mu*(conj(R).*f))*scale;
    
  uker=zeros(rows,cols);
  uker(1)=2; uker(2)=-1; uker(end)=-1;
  uker=mu*(conj(R).*R) + lambda*fft(uker);
  
%% Reconstruction  
  outer=0; % count while loop its
  while norm(R.*fft(u)/scale-f0,2)/norm(f0,2)>epsilon && outer<1000
  % Inner loops  
    for inner=1:nInner,
    % Update u   
      rhs=murf+lambda*Dxt(x-bx);
      u=ifft(fft(rhs)./uker);

    % Update x 
      dx=Dx(u);
      x=shrink(dx+bx,1/lambda);

    % Update Bregman parameter
      bx=bx+dx-x;
    end

  % Update Fourier data  
    f=f+f0-R.*fft(u)/scale;
    murf=ifft(mu*R.*f)*scale;
    outer=outer+1;
  end
end

function d=Dx(u)
% Gradient; periodic
  d=zeros(size(u));
  d(:,2:end)=u(:,2:end) - u(:,1:end-1);
  d(:,1)=u(:,1) - u(:,end);
end

function d=Dxt(u)
% Gradient "transposed"; periodic
  d=zeros(size(u));
  d(:,1:end-1)=u(:,1:end-1) - u(:,2:end);
  d(:,end)=u(:,end) - u(:,1);
end

function[xs]=shrink(x,lambda)
% Shrinkage (4.4)-(4.6) of Split Bregman paper
  s=sqrt(x.*conj(x));
  ss=s-lambda; ss=ss.*(ss>0); 
  s=s+(s<lambda);
  xs=ss.*x./s;
end