function[recovered,R,mnoise,outer]=recovery_1D(sig,sampling,m,epsilon,std)
% For given 1D signal, get measurements by a subsampled Fourier operator
% and perform CS reconstruction by TV min. using the split Bregman.
%
% Inputs:
%   sig: 1D signal
%   sampling: sampling strategy, uniform, low, power, unif power
%   m: sampling ratio
%   epsilon: stopping criteria, halts SB when norm(Az - y) < epsilon
%   std: std dev of gaussian noise added to measurements to simulate
%       corrupted data
% Outputs: 
%   recovered: signal recovered by TV minimization
%   R: sampling map.  
%   mnoise: noise added to measurements
%   outer: # of iter of the SB outer loop
%
% Modified from test_mrics.m, Dec 2015, Kevin Chow
% Last modified Nov 2016, Kevin Chow
%

%% Check inputs
  if nargin < 4
    error('First 4 arguments must be specified');
  elseif nargin < 5
    std=0;
  end
  
%% Parameters
% Split Bregman 
  mu=1.1; lambda=0.10;
% Signal length  
  N=length(sig);

%% Build sampling matrix
% 4 cases:
%   uniform random sampling             'uniform'
%   low freq sampling                   'low'
%   power law density sampling          'power'
%   unif + power law density sampling   'unif power'
%
% Adjust first parameter of sample.m to get the right sampling ratio
  switch sampling
    case 'uniform'
      R=rand(size(sig));
      R=double(R < m);
      
    case 'low'
      low=round(m*N/2);
      R=zeros(size(sig));
      R(1:low)=1;
      R(end-low+1:end)=1;
      
    case 'power'
      C=pn1dC(N);
      R=sample_m(round(1.95*m*N), N, @(n) pn1d(n,N,C));
      
    case 'unif power'
      C=pn1dC(N);
      R=rand(size(sig));
      R1=double(R < m/2);
      R2=sample_m(round(0.86*m*N), N, @(n) pn1d(n,N,C));
      R=double(R1|R2);
      
  end
% Always sample 0-frequency/flat vector
  R(1)=1; 

%% Form CS data with noise
  F=R.*fft(sig)/sqrt(numel(sig));
  mnoise=R.*randn(size(F))*std; F=F+mnoise; 

%% Recover using the split Bregman
  [recovered, outer]=split_breg_1D(R, F, mu, lambda, 10, epsilon);
end
