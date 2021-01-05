function [E,cmin,cmax,ns,nt] = geterrorimages(X, x_true, fighandle)
%
% This function constructs error images and displays them on the same
% colormap
% x_true needs to a 2D image
% X contains all of the images to compare to (in the columns of X)
% fighandle is the figure handle

if nargin < 3
  fighandle = figure(119);
end

% Size of the image
[ns, nt] = size(x_true);

% Number of images
[~, nimages] = size(X);

% Absolute errors
E = abs(bsxfun(@minus, X, x_true(:)));
cmin = min(E(:));
cmax = max(E(:));
E = cmax-E;

% Plot in current figure
for i = 1:nimages
  figure(fighandle),
  subplot(2,nimages,i+nimages), imshow(reshape(E(:,i),ns, nt),[])
  caxis([cmin cmax]), colorbar
end