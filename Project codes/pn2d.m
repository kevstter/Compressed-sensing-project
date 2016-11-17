function[pnk]=pn2d(n,k,N,C)
% Computes the probability based on the density
%   p(n,k)=C*(log(N)*max(1, n^2+k^2))^{-1}, for k,n = -N/2+1, ... , N/2,
% where C = C(N) is chosen so that p is a probability distribution.
%
% p(n,k) is the probability that frequency (n,k) is to be sampled.
%
% Inputs:
%   n,k: frequency number
%   N: frequencies under consideration
%   C: normalizing constant
% Output:
%   pnk: probability that frequency indexed by (n,k) is to be sampled
% 
% References: 
% [1] Poon, On the role of TV in compressed sensing
% [2] Krahmer and Ward, Stable and robust sampling stragies for CS
%

%% Check inputs
  if nargin < 3
    error('Require minimum 3 inputs');
  elseif nargin < 4
    C = 1+10/N^2;
    for q=1:N/2-1
      C=C+16./(N^2+4*q^2) + 6./q^2 + 8*sum(1./((q+1:N/2-1).^2+q^2));
    end
    C=C/log(N); C=1/C;
  end
  
  if max(n)>N/2 || min(n)<-N/2+1 || max(k)>N/2 || min(k)<-N/2+1
    error('(n,k) outside of expected frequencies');
  end
  
%% Compute probability  
    pnk = C./(log(N).*max(1,n.^2+k.^2));
end