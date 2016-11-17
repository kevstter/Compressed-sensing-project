function[pn,n]=pn1d(n,N,C)
% Computes the probability dictated by the density function
%   p(n) = C*(log(N)*max(1,abs(n)))^{-1}, for n=-N/2+1,...,N/2,
% where C=C(N) is chosen so that p is a probability distribution.
%
% p(n) is the probability that frequency n is to be sampled.
%
% Inputs: 
%   n: frequency between -N/2+1 to N/2
%   N: Number of frequencies 
%   C: the normalizing constant
% Output: 
%   p(n): probability the frequency n is to be sampled.
%
% References: 
% [1] Poon, On the role of TV in compressed sensing
% [2] Krahmer and Ward, Stable and robust sampling stragies for CS 
%
% Created Dec 2015, Kevin Chow
% Last modified Nov 2016, Kevin Chow
%

%% Check inputs
  if nargin < 2
    error('Minimum 2 inputs');
  elseif nargin < 3
    C=log(N)/(1 + 2/N + 2*sum(1./(1:N/2-1)));
  end

%% Parameter validity  
  if max(n) > N/2 || min(n) < -N/2+1
    error('Given n is outside of expected frequencies');
  end

%% Compute prob  
  n=abs(n); pn=C./(log(N).*max(1,n));
end