function[grad_h, grad_v] = im_gradient(img)
% Returns discrete gradients; assumes grayscale image
%
% Inputs: 
%   img: an image, in png, jpg, etc. 
% Outputs: 
%   grad_h, grad_v: horizontal and vertical discrete gradient of the 
%   input image.
%
% Ex: im_gradient(imread('pears.png')));
%

%% Preprocess image
  if size(img,3)~=1; 
    img=rgb2gray(img); 
  end
  
%% Compute gradients  
  grad_v=img(2:end,:) - img(1:end-1,:);
  grad_h=img(:,2:end) - img(:,1:end-1);
  
%% Visualize
  fignum=802; figure(fignum); clf
  imagesc(grad_h); axis('image','off');
  figure(fignum+1); clf; 
  imagesc(grad_v); axis('image','off'); 
end
