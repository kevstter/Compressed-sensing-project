function[samples]=sample_m2d(m,N,prob)
% Draw m samples from frequencies 
%   (n,k) \in {-N/2+1,...,N/2}x{-N/2+1,...,N/2},
% according to the distribution prob.
%
% Inputs: 
%   m,N: sampling ratio of m/N; N describes range of frequencies
%   prob: density function
% Output: 
%   samples: NxN matrix of 1s and 0s
%
% ex. 
%   N=128; C=pn2dC(N);
%   abab=sample_m2d(round(N^2/10),N,@(n,k) pn2d(n,k,N,C)); 
%   figure,imagesc(fftshift(abab)); 
%
% References: 
% [1] Poon, On the role of TV in compressed sensing
%
% Created Dec 2015, Kevin Chow
% Last modified Nov 2016, Kevin Chow
%

%% Generate m unif random numbers; Create storage.
% Mult by 2.6 because of repeated sampling of the same freq
  m=round(2.6*m); m=min(m,N^2); 
  if m==N^2; samples=ones(N,N); return; end
  mm=rand(1,m); mm=sort(mm);
  samples=zeros(N,N);
    
%% Construct binary sampling matrix
% DC
  pr=prob(0,0);
  idx_at=nnz(mm<=pr)+1;
  if idx_at>1, samples(1,1)=1; end 
  if idx_at>m, return; end

% Other:
% Sub-blocks in the for loop does the following resp.:
% Sample frequencies of the type (\pm j,0),(0,\pm j) 
% Sample frequencies of the type (\pm j,\pm j)
% Sample frequencies of the type (N/2,\pm j),(\pm j,N/2)
% Sample frequencies of the type (\pm j,\pm ii),(\pm ii,\pm j)
  for j=1:N/2-1
    row=[j+1 1 N-j+1 1 j+1];
    [pr,idx_at,samples]=subblock(j,row,prob,pr,m,mm,idx_at,samples);
    if idx_at>m, return; end

    row=[j+1 j+1 N-j+1 N-j+1 j+1];
    [pr,idx_at,samples]=subblock(j,row,prob,pr,m,mm,idx_at,samples);
    if idx_at>m, return; end


    row=[j+1 N/2+1 N-j+1 N/2+1 j+1];
    [pr,idx_at,samples]=subblock(j,row,prob,pr,m,mm,idx_at,samples);
    if idx_at>m, return; end

    for ii=j+1:N/2-1
      row=[j+1 ii+1 N-j+1 N-ii+1 j+1 N-ii+1 N-j+1 ii+1 j+1];
      [pr,idx_at,samples]=subblock(j,row,prob,pr,m,mm,idx_at,samples);
      if idx_at>m, return; end
    end
  end
  samples(N/2+1,N/2+1)=1;
end

function[pr,idx_at,samples]=subblock(j,row,prob,pr,m,mm,idx_at,samples)
  pj=prob(j,row(2)-1);
  for jj=1:length(row)-1,
    pr=pr+pj; 
    if mm(idx_at)<pr,
      samples(row(jj),row(jj+1))=1; 
      idx_at=mv_idx_at(idx_at,m,mm,pr);
      if idx_at>m, return; end
    end
  end
end
function[idx_at]=mv_idx_at(idx_at,m,mm,pr)
  while idx_at<=m && mm(idx_at)<=pr,
    idx_at=idx_at+1;
  end
end