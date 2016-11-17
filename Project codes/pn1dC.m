function[C] = pn1dC(N)
% Returns the normalizing constant needed for pn1d.m
% See pn1d for details.
%

%% Check inputs
  if nargin < 1, error('Supply N'); end
    
%% Compute constant  
  C=log(N)/(1 + 2/N + 2*sum(1./(1:N/2-1)));
end
  