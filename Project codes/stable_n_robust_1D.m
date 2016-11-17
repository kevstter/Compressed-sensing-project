function[m,l2err,l2err_grad]=stable_n_robust_1D(N,s,m,sampling,epsilon, ...
  noise,scale,std,sig_type,cmpct)
% Tests stability and robustness of recovery of 1D signals by TV minimization;
% stable to inexact sparsity; robust to noise.
%
% Inputs:
%   N: signal length; must be even
%   s: (gradient) sparsity of the signal
%   m: sampling ratio
%   sampling strategy. See recovery_1d.m for supported strategies.
%   epsilon: recon stopping crit,ie. stop when norm(Az-y) < epsilon.  
%   noise: adding noise to the signal,'no','noise'
%   scale: scale magnitude of noise,noise <- noise/scale.
%   std: standard deviation of (gaussian) measurement noise
%   sig_type: global or local signal,'global','local'
%   cmpct: local signal support length; given as multiple of s
% Outputs:
%   Generates a collection of plots
%
% Ex. N=1024;
% stable_n_robust_1D(N,round(N/50),140/N,'low',1e-3,'noise',1*5,0.1,'global',8) 
% stable_n_robust_1D(N,round(N/50),140/N,'unif power',1e-3,'noise',1*5,0*0.1,'global',8) 
% stable_n_robust_1D(N,round(N/50),140/N,'uniform',1e-3,'no',1,0,'local',8)
%
% References:
% [1] Poon, On the role of TV in compressed sensing
% [2] Goldstein and Osher, Split Bregman for l1-regularized problems
%
% Created Dec 2015, Kevin Chow
% Last modified Nov 2016, Kevin Chow
%

%% Check parameters
  if nargin < 5, error('Must specify first 5 inputs');
  elseif nargin < 7, noise='no'; scale=1;
  elseif nargin < 8, std=0;
  elseif nargin < 10, sig_type='global'; cmpct=0;
  end
  if mod(N,2)==1, error('Even signal length only'); end

%% Construct signal; gradient s-sparse  
  switch sig_type
    case 'global'
      sig=grad_sparse(s,N,10,0);
    case 'local'
      sig=cmpct_grad_sparse(s,cmpct,N,10,0);
  end
  sig0=sig;
  sig_grad=sig(2:end)-sig(1:end-1);

%% Adding noise to the signal
  switch noise 
    case 'noise'
      noise=randn(size(sig)); noise=noise/scale; 
      sig=sig+noise;
    case 'no'
      noise=zeros(size(sig));
  end

%% Sample and reconstruct
  [recovered,R,mnoise,~]=recovery_1d(sig,sampling,m,epsilon,std);
  recovered=real(recovered);
  rec_grad=recovered(2:end)-recovered(1:end-1);

%% Display some statistics
  imgName='1D gradient sparse'; m=nnz(R);
  l2err=norm(recovered-sig,2)/norm(sig,2); 
  l2err_grad=norm(rec_grad-sig_grad,2)/norm(sig_grad,2); 
  
  fprintf(['%s:\n\t Sampling ratio=%.2g\n\t l2error=%.3g\n\t' ...
    ' l2error_grad=%.3g\n'],imgName,m/numel(sig),l2err,l2err_grad);

%% Visualize
  fignum=810; figure(fignum); % suptitle('TV recovery')

  ax1=subplot(4,2,1); plot(sig);
  axis([0 N min(sig0)-2 max(sig0)+2]); 
  title('Signal');

  ax3=subplot(4,2,3); plot(sig_grad);
  axis([0 N min(sig_grad)-2 max(sig_grad)+2]); 
  title('Discrete Gradient');

  ax5=subplot(4,2,5); plot(-N/2+1:N/2,fftshift(mnoise));
  axis([-N/2+1 N/2 min(1.1*min(mnoise),-0.5) max(1.1*max(mnoise),0.5)]);
  title('Noise (added to the measurements)');

  ax7=subplot(4,2,7); plot(noise);
  axis([0 N min(1.1*min(noise),-1) max(1.1*max(noise),1)]); 
  title('Noise (added to the sparse signal)');

  ax2=subplot(4,2,2); plot(recovered);
  axis([0 N min(sig)-2 max(sig)+2]); 
  title('Recovered Signal');

  ax4=subplot(4,2,4); plot(rec_grad);
  axis([0 N min(sig_grad)-2 max(sig_grad)+2]); 
  title('Recovered Gradient');

  ax6=subplot(4,2,6); stem(-N/2+1:N/2,fftshift(R),'d','marker','none');
  axis([-N/2+1 N/2 0 1]); set(ax6,'ytick',[]);
  title(['$P_\Omega$, sampling ratio ',num2str(m/numel(sig))]);

  ax8=subplot(4,2,8); plot(sig0);
  axis([0 N min(sig0)-2 max(sig0)+2]); 
  title('Gradient Sparse Signal');

  linkaxes([ax1,ax2,ax8],'xy');
  linkaxes([ax3,ax4],'xy');
end