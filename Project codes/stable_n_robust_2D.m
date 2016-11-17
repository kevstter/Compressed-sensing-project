function[R]=stable_n_robust_2D(img,m,sampling,epsilon,noise,scale,std)
% Tests stability and robustness of recovery of images by TV minimization
%
% Inputs:
%   img: specify test image
%   m: sampling ratio
%   sampling strategy. See recovery_2d.m for supported strategies.
%   epsilon: recon stopping crit, ie. stop when norm(Az-y) < epsilon.  
%   noise: adding noise to the image, 'no', 'noise'
%   scale: scale magnitude of noise, noise <- noise/scale.
%   std: standard deviation of (gaussian) measurement noise
% Outputs:
%   Generates a collection of plots
%
% Ex. 
%   stable_n_robust_2D('pears', 0.15, 'power', 1e-3, 'noise', 1, 0);
%   stable_n_robust_2D('peppers', 0.15, 'uniform', 1e-3, 'no', 1, 0);
%
% Created Dec 2015, Kevin Chow
% Last modified Nov 2016, Kevin Chow
%

%% Check parameters
  if nargin < 3, error('Specify first 3 inputs');
  elseif nargin < 4, epsilon=1e-3; noise='no'; scale=1; std=0;  
  elseif nargin < 5, noise='no'; scale=1; std=0;
  elseif nargin < 6, scale=1; std=0;
  elseif nargin < 7, std=0;
  end
  
%% Select image  
%   spine.tif, pears.png, ... see toolbox/images/imdata
% Modify to add your own image
switch img
  case 'spine'
    name=img; img=imread('spine.tif'); img=img(1:end-1,63:end-62);
  case 'pears'
    name=img; img=imread('pears.png'); img=rgb2gray(img);
  case 'peppers'
    name=img; img=imread('peppers.png'); img=rgb2gray(img);
  case 'boats'
    name=img; img=imread('boats.png'); 
  case 'cameraman'
    name=img; img=imread('cameraman.tif');
  otherwise
    error('No such image');
end

% 'squarify' the image; limitation of sample_m2d.m; square image; NxN
  N=min(size(img)); N=N-mod(N,2); img=img(1:N,1:N);
  img=im2double(img)*255;

%% Add noise to image
  switch noise 
    case 'noise'
      noise=randn(size(img)); noise=noise/scale; 
      img=img+noise;
    case 'no'
%       noise=sparse(size(img));
  end

%% Sample and reconstruct
  [rec,R,F,~]=recovery_2D(img,sampling,m,epsilon,std);
% outer

%% Display some stats
  l2err=norm(rec-img,2)/norm(img,2); m=nnz(R); m=m/numel(img);
  name=regexprep(name,'(\<[a-z])','${upper($1)}');
  fprintf('%s:\n\t Samping ratio=%.2g\n\t l2err=%.3g\n',name,m,l2err);
  
%% Visualize
  fignum=815; figure(fignum); % suptitle('TV recovery')

  ax1=subplot(2,2,1); imagesc(real(ifftn(F))); axis('image','off'); 
  title('Image reconstruction with unknowns set to zero'); 
  
  ax2=subplot(2,2,2); imagesc(real(rec)); 
  title('Split Bregman Recovery'); axis('image','off');

  ax3=subplot(2,2,3); imagesc(abs(img)); 
  title(name); axis('image','off');

  ax4=subplot(2,2,4); imagesc(abs(fftshift(R)));
  title(['R: sampling ratio=',num2str(m,3)]); axis('image','off');

  linkaxes([ax1, ax2, ax3],'xy')
end