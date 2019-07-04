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


function imagesFTT=imagespreparation_FTT(filestub, ext, init, final,viz,LowFreqBand,HighFreqBand)

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
    image=(images(:,:,img)-minint)./(maxint-minint);
    % drop values less than zero otherwise after filter image looks strange
    image(image<0)=0;
    images(:,:,img)=image;
end;

%=============================================================
% filter all images
%=============================================================
imsiz=size(images(:,:,img));
H=gaussianbpf(imsiz,LowFreqBand,HighFreqBand);
% HS=real(ifft2(H));
%imshow(fftshift(HS),[])
% imshow(fftshift(H));
%imshow(fftshift(H));


% gaussian filter have an artifacts on the boader
% size boader to cut off from image
Border=100;
Xdirection=Border:(imsiz(1)-Border);
Ydirection=Border:(imsiz(2)-Border);
imagesFTT=zeros(length(Xdirection),length(Ydirection),nimg);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for img=1:nimg,
    % obtain FTT with padded imput
    imFTT=fft2(images(:,:,img),2*imsiz(1)-1,2*imsiz(2)-1);%,2*imsiz(1)-1,2*imsiz(2)-1
    % perform filtering
    imFTT=imFTT.*H;
    % make inverse FTT
    im=real(ifft2(imFTT));
    % crop the image
    %im=im(1:imsiz(1),1:imsiz(2));
    im(im<0)=0;
    imagesFTT(:,:,img)=im(Xdirection,Ydirection);
end;


% ftt_norms=imhist(abs(F))/numel(F);
% cdf=cumsum(ftt_norms);
% im=real(ifft2(F)); 
% im2=im(1:imsiz(1),1:imsiz(2));
% F2=fftshift(F);      
% figure,imtool(log(1+abs(F2)),[])

%viz_image_stack(nimg,imagesFTT)
%viz_image_stack(nimg,images)



return
