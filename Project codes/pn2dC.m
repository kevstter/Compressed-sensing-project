function[C]=pn2dC(N)
% Returns the normalizing constant needed for pn2d.m
% See pn2d for details.
%

%% Check inputs
  if nargin < 1, error('Supply N'); end
  
%% Compute constant  
  C=1+10/N^2;
  for q=1:N/2-1
    C=C+16./(N^2+4*q^2) + 6./q^2 + 8*sum(1./((q+1:N/2-1).^2+q^2));
  end
  C=C/log(N); C=1/C;
end