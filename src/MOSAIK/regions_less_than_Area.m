orig=imagesFTT(:,:,1);
level = graythresh(orig)
BW = imbinarize(orig,level);
AreaLevel_top=200
AreaLevel_bottom=10

CC=bwlabel(BW);
stats=regionprops(CC,'Area','PixelId','MajorAxisLength','MinorAxisLength');
orig_select=zeros(size(orig));

% histogramme 
MajAxis=[stats.MajorAxisLength];histogram(MajAxis(MajAxis<200));
MinAxis=[stats.MinorAxisLength];histogram(MinAxis(MinAxis<200));

idx=find(Area>AreaLevel_bottom & Area<AreaLevel_top);
npart=length(idx);
for ii=1:npart
    orig_select(stats(idx(ii)).PixelIdxList)=1;
end;


figure,imshowpair(orig,orig_select),title('Compare orig with selected by area');

