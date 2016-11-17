function[samples]=sample_m(m,N,prob)
% Draw m samples from frequencies -N/2+1,...,N/2 according to the
% distribution prob.
% Inputs:
%   N: frequencies -N/2+1,..,N/2,
%   m: sampling ratio of m/N.
%   prob: a probability distribution
% Outputs:
%   samples: binary vector; approx m 1's; marks the freqs to be sampled
%
% ex. 
%   N=128;
%   abab=sample_m(10,N,@(n) pn1d(n,N)); 
%   stem(-N/2+1:N/2, fftshift(abab)); axis('tight')
%
% Last modified Nov 2016, Kevin Chow
%

%% Generate unif random
  mm=rand(1,m); mm=sort(mm); samples=zeros(1,N);
  
%% Probabilites indicates interval [a,b] to check; if [a,b]\mm~=[a,b], sample 
% DC
  pr=prob(0);
  if nnz(mm<=pr)>0, samples(1)=1; end
% Each index corresp. to pair of frequencies  
  for j=1:N/2-1,
    pj=prob(j); freq=[j+1 N-j+1];
    for i=freq,
      if nnz((mm<pr+pj).*(mm>pr))>0,
        samples(i)=1;
      end
      pr=pr+pj;
    end
    if mm(end)<pr, return; end
  end
  samples(N/2)=1;
  
%% Faster, but maybe less clear  
%   pr=prob(0);
%   idx_at = nnz(mm <= pr)+1;
%   if idx_at > 1, samples(1) = 1; end 
%   if idx_at > m, return; end
%   for jj = 1:N/2-1
%       pjj = prob(jj);
%       loc = [jj+1 N-jj+1];
%       for ii = loc;
%           pr = pr+pjj;
%           if mm(idx_at) < pr
%               samples(ii) = 1;
%               while idx_at <= m && mm(idx_at) <= pr
%                       idx_at = idx_at + 1;
%               end
%               if idx_at > m, return; end
%           end
%       end
%   end
%   samples(N/2) = 1;
end