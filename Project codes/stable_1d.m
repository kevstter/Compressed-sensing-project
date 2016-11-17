function[] = stable_1d(N, s, m, sampling, epsilon, noise, scale, std, ...
    sig_type, cmpct)
% FOR LATEST IMPLEMENTATION, USE stable_n_robust_1D.m INSTEAD  
%
% test stability of recovery by TV for 1D signals
%
% Inputs:
%   N: signal length, N must be even
%   s: (gradient) sparsity of the signal
%   m: sampling ratio
%   sampling strategy. See test_mrics1d.m for supported strategies.
%   epsilon: recon stopping crit, ie. stop when norm(Az-y) < epsilon.  
%   noise: adding noise to the signal, 'no', 'noise'
%   scale: scale magnitude of noise, noise <- noise/scale.
%   std: standard deviation of measurement noise
%   sig_type: global or local signal, 'global', 'local'
%   cmpct: local signal support length; given as multiple of s
% Outputs:
%   Generates a collection of plots
%
% Ex. N = 1024;
% stable_1d(N, round(N/50), 140/N, 'low', 1e-3, 'noise', 1*5, 0.1, 'global', 8) 
% stable_1d(N, round(N/50), 140/N, 'uniform', 1e-3, 'no', 1, 0, 'local', 8)
%
% Dec 2015, Kevin Chow
if nargin < 5
    error('Specify first 5 inputs');
elseif nargin < 7
    noise = 'no'; scale = 1;
elseif nargin < 8
    std = 0;
elseif nargin < 10
    sig_type = 'global'; cmpct = 0;
end
if mod(N, 2) == 1
    error('Even signal length only')
end

% % % % % % % % % % % % % % % % % % % % % % % % 1D gradient sparse signal:
switch sig_type
    case 'global'
        sig = grad_sparse(s, N, 10, 0);
    case 'local'
        sig = cmpct_grad_sparse(s, cmpct, N, 10, 0);
end

sig0 = sig;
sig_grad = sig(2:end)-sig(1:end-1);

switch noise 
    case 'noise'
        noise = randn(size(sig)); noise = noise/scale; 
        sig = sig+noise;
    case 'no'
        noise = zeros(size(sig));
end

% % % % % % % % % % % % % % % Call test_mrics1d to sample and reconstruct:
[recovered, R, mnoise, outer]=test_mrics1d(sig, sampling, m, epsilon, std);
recovered = real(recovered);
rec_grad = recovered(2:end)-recovered(1:end-1);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % Some statistics: 
imgName = '1D gradient sparse';
l2error = norm(recovered - sig, 2)/norm(sig, 2);
l2error_grad = norm(rec_grad - sig_grad, 2)/norm(sig_grad, 2);
m = nnz(R);
disp([imgName, ':', 10, 9, 'Sampling ratio = ', num2str(m/numel(sig), 3), ...
    10, 9, 'l2error = ', num2str(l2error), ...
    10, 9, 'grad_l2error = ', num2str(l2error_grad)]);

% % % % % % % % % % % % % % % % % % % % build a figure to display results:
figure; 
% suptitle('TV recovery')
ax1 = subplot(4,2,1); plot(sig);
axis([0 N min(sig0)-2 max(sig0)+2]); title('Signal');

ax3 = subplot(4,2,3); plot(sig_grad);
axis([0 N min(sig_grad)-2 max(sig_grad)+2]); title('Discrete Gradient');

ax5 = subplot(4,2,5); plot(-N/2+1:N/2, mnoise);
axis([-N/2+1 N/2 min(1.1*min(mnoise),-1) max(1.1*max(mnoise), 1)]);
title('Noise (added to the measurements)');

ax7 = subplot(4,2,7); plot(noise);
axis([0 N min(1.1*min(noise),-1) max(1.1*max(noise),1)]); 
title('Noise (added to the sparse signal)');

ax2 = subplot(4,2,2); plot(recovered);
axis([0 N min(sig)-2 max(sig)+2]); title('Recovered Signal');

ax4 = subplot(4,2,4); plot(rec_grad);
axis([0 N min(sig_grad)-2 max(sig_grad)+2]); title('Recovered Gradient');

ax6 = subplot(4,2,6); stem(-N/2+1:N/2, fftshift(R),'d','marker','none');
axis([-N/2+1 N/2 0 1]); title(['R, sampling ratio ', num2str(m/numel(sig))]);

ax8 = subplot(4,2,8); plot(sig0);
axis([0 N min(sig0)-2 max(sig0)+2]); title('Gradient Sparse Signal');

linkaxes([ax1,ax2,ax8],'xy');
linkaxes([ax3,ax4],'xy');