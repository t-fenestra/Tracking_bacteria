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


function images_restored=imagespreparation(filestub, ext, init, final,viz,BoxFilter,GaussFilter_lambda)

%=============================================================
% read images and determine global extrema
%=============================================================

if viz ==1,
    nfig = 1;
else
    nfig = -1;
end;

disp('Scanning files for minimum and maximum intensity values...')
maxint = -Inf;
minint = Inf;
numsat = 0;

for img=init:final
    file = sprintf('%s%d%s',filestub,img,ext);
    orig = double(imread(file));
    locmax = max(max(orig));
    locmin = min(min(orig));
    if locmax > 254, numsat = numsat+1; end;
    if locmax > maxint, maxint = locmax; end;
    if locmin < minint, minint = locmin; end;
    images(:,:,img-init+1) = orig;
end;
nimg = final-init+1;
disp(sprintf('%d frame images successfully read',nimg))
disp(sprintf('Minimum intensity value is: %f',minint))
disp(sprintf('Maximum intensity value is: %f',maxint))
if numsat > 0,
    disp(sprintf('WARNING: found %d saturated pixels !',numsat))
end;

%=============================================================
% normalize all images
%=============================================================
disp('normilize files to the global max and min...')
for img=1:nimg,
    images(:,:,img) = (images(:,:,img)-minint)./(maxint-minint);
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

for cimg=1:nimg,
    % apply convolution filter
    filtered = conv2(images(:,:,cimg),K,'same');
    % Set every value that is smaller than 0 to 0 
    filtered(filtered<0)=0;
    
    % get rid of the kernel-boader artifacts
    filtered_boader=zeros(siz(1),siz(2));
    filtered_boader(w+1:end-w,w+1:end-w)=filtered(w+1:end-w,w+1:end-w);
    if viz == 1,
        figure(img);imshowpair(images(:,:,cimg),filtered_boader,'montage');
     end;
    
     images_restored(:,:,cimg)=filtered_boader;
end;


disp(sprintf('%d images successfully restored',nimg))
%viz_image_stack(nimg,images)

return
