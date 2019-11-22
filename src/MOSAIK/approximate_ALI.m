function ALI=approximate_ALI(img)
%======================================================================
% Approximate Ali by the straight line
% INPUTS: images

% Returns:Y=coordinats
%https://de.mathworks.com/help/images/texture-segmentation-using-texture-filters.html
%======================================================================

% assume that the AlI can not be deeper than 300px
img=img(1:400,:);
% entropy filter

E = entropyfilt(img);

%calculate median entropy on the top 50pixel stripe
TopLayer=E(1:50,500:end);
Threshold=median(TopLayer(:));
BW = imbinarize(E, Threshold);
%imshow(BW);

BWao = bwareaopen(BW,1000);
%imshow(BWao)

nhood = true(9);
closeBWao = imclose(BWao,nhood);
%imshow(closeBWao)

roughMask = imfill(closeBWao,'holes');
%imshow(roughMask);

boundary = bwperim(roughMask);
segmentResults = E;
segmentResults(boundary) = 255;
%imshow(boundary)

%boundary= bwareaopen(boundary,20);

[y,x]=find(boundary>0);
ALI=zeros(2048,1);

Perimeter=[x,y];
X=unique(Perimeter(:,1));


%if ~(isempty(Perimeter))
    for i=1:1:length(X)
        ALI(X(i))=max(Perimeter(Perimeter(:,1)==X(i),2));
    end
%end

