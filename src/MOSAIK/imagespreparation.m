
%====================================================================== 
%
% IMAGES PREPARATION: normalize all images to global minimum and global
%                     maximum
% SYNTAX:  imagestack = imagespreparation(filestub, ext, init, final,viz,BoxFilterW,GaussFilter_lambda)
%
% INPUTS:  filestub             stub of frame image file names incl. path
%          ext                  extension of file names (without the dot)
%          init                 number of initial frame
%          final                number of final frame
%          viz                  visualization mode (range 1,0)
%          BoxFilterW           box filter window (2*BoxFilterW+1)
%          GaussFilter_lambda   gauss filter lambda
%
% The input files are expected to be called <filestub>[init..final]
% with subsequent numbering referring to the corresponding frame
% in the movie.
%
% code based on matlab version of Mosaik by Ivo Sbalzarini, 12.2.2003
% update  21.12.2018
%====================================================================== 


function images=imagespreparation(fname,init,final,BoxFilter,GaussFilter_lambda)

%=============================================================
% read images and determine global extrema
%=============================================================
info = imfinfo(fname);
num_images = numel(info);

maxint = -Inf;
minint = Inf;
numsat = 0;

images=zeros(2048,2048,final-init+1);

for img=init:final
    orig = double(imread(fname, img, 'Info', info));
    locmax = max(max(orig));
    locmin = min(min(orig));
    
    if locmax > maxint, maxint = locmax; end;
    if locmin < minint, minint = locmin; end;
    images(:,:,img-init+1) = orig;
end;

%=============================================================
% normalize all images
%=============================================================
disp('normilize files to the global max and min...')
 nimg = final-init+1;
for img=1:nimg,
    image=(images(:,:,img)-minint)./(maxint-minint);
    % drop values less than zero otherwise after filter image looks strange
    image(image<0)=0;
    images(:,:,img)=image;
end;

%====================================================================== 
% Image restoration
%====================================================================== 
w=BoxFilter;
% correlation length of camera noise (usu. set to unity)
lambdan = GaussFilter_lambda;

% some often used quantities
idx =[-w:1:w];     % index vector
dm = 2*w+1;         % diameter
im = repmat(idx',1,dm);
jm = repmat(idx,dm,1);
imjm2 = im.^2+jm.^2;
siz=size(images(:,:,1)); % size of image


% build kernel K for background extraction and noise removal
% (eq. [4])
B = sum(exp(-(idx.^2/(4*lambdan^2))));
B = B^2;
K0 = 1/B*sum(exp(-(idx.^2/(2*lambdan^2))))^2-(B/(dm^2));
K = (exp(-(imjm2/(4*lambdan^2)))/B-(1/(dm^2)))/K0;


warning('off', 'Images:initSize:adjustingMag')
viz=0;

for cimg=1:nimg,
    % apply convolution filter
    filtered = imfilter(images(:,:,cimg),K,'conv','replicate');
    % Set every value that is smaller than 0 to 0 
    filtered(filtered<0)=0;
    images(:,:,cimg)=filtered;
    % cut bright artificial spot on the left
    images(1:70,1:70,cimg)=0;
end;


disp(sprintf('%d images successfully restored',nimg))

%viz_image_stack(nimg,images)
%viz_image_stack(nimg,images)
%1
return