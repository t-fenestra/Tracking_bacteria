function images=read_image_stack(fname,init,final)
info = imfinfo(fname);
num_images = numel(info);
images=zeros(2048,2048,final-init+1);

for img=init:final
    orig = double(imread(fname, img, 'Info', info));
    
    images(:,:,img-init+1) = orig;
end;


end