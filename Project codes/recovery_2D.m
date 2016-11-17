function[recovered,R,F,outer]=recovery_2D(img,sampling,m,epsilon,std)
% Given an image, get measurements by a subsampled Fourier operator
% and perform CS reconstruction by TV min. using the split Bregman 
%
% Inputs:
%   img: image
%   sampling: sampling strategy, uniform, low, power, unif power
%   m: sampling ratio=# of measurements / numel(img)
%   epsilon: tolerance, halts SB when norm(Az - y) < epsilon
%   std: standard deviation of gaussian noise added to measurements
% Outputs: 
%   recovered: signal recovered by TV minimization using the SB
%   R: sampling map.  
%   F: (Noisy) Fourier samples
%   outer: # of iters of the SB outer loop
%
% Created Dec 2015, Kevin Chow
% Last modifed Nov 2016, Kevin Chow
%

%% Check inputs
  if nargin < 4
    error('First 4 arguments must be specified');
  elseif nargin < 5
    std=0;
  end
  
%% Parameters
% Split Bregman
  mu=0.50; lambda=0.10;

%% Build sampling matrix
% 4 cases:
%   uniform random sampling             'uniform'
%   low freq sampling                   'low'
%   power law density sampling          'power'
%   unif + power law density sampling   'unif power'
%
% Adjust first parameter of sample_m2d.m to get the right sampling ratio 
  switch sampling
    case 'uniform'
      R=rand(size(img)); R=double(R < m);
    case 'low'
      R=zeros(size(img));
      low=sqrt(m);
      km=round(low*size(img, 1)/2);
      kn=round(low*size(img, 2)/2);
      R(1:km, 1:kn)=1;
      R(end-km:end, end-kn:end)=1;
      R(1:km, end-kn:end)=1;
      R(end-km:end, 1:kn)=1;
    case 'power'
      N=size(img, 1); C=pn2dC(N);
      R=sample_m2d(round(m*N^2), N, @(n,k) pn2d(n,k,N,C));
    case 'unif power'
      N=size(img, 1); C=pn2dC(N);
      R=rand(N, N);
      R1=double(R < m/2);
      R2=sample_m2d(round(m/2*N^2), N, @(n,k) pn2d(n,k,N,C));
      R=double(R1|R2);
    otherwise
      error('Specify available sampling strategy');
  end
% Always sample 0-frequency  
  R(1,1)=1; 

%% Form CS data with noise
  F=R.*fftn(img)/sqrt(numel(img));
  mnoise=R.*randn(size(F))*std; F=F+mnoise; 

%% Recover by TV minimization using the split Bregman
  [recovered,outer]=split_breg_2D(R,F,mu,lambda,10,epsilon);

