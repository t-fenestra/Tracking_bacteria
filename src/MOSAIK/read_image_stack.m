function images=read_image_stack(fname,init,final)
info = imfinfo(fname);
num_images = numel(info);
images=zeros(2048,2048,final-init+1);

maxint = -Inf;
minint = Inf;
numsat = 0;

images_adjust=zeros(2048,2048,final-init+1);

%=============================================================
% read images and determine global extrema
%===========================================================
for img=init:final
    orig = imread(fname, img, 'Info', info);
    locmax = max(max(orig));
    locmin = min(min(orig));
    
    if locmax > maxint, maxint = locmax; end;
    if locmin < minint, minint = locmin; end;
    images(:,:,img-init+1) = orig;
end;

end