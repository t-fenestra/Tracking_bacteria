%======================================================================
%
% DETECT_PARTICLES: detect particle-shaped features in frame images
%
% SYNTAX:  peak = detect_particles(orig,w,cutoff,pth,v)
%
% INPUTS:  orig     original image to detect features in
%          w        global size parameter, >particle radius and
%                   <interparticle spacing
%          cutoff   probability cutoff for non-particle discrimination
%          pth      percentile threshold for maxima selection
%          v        visualization parameters [viz, nfig] as:
%                   viz=1 if intermediate visualization is needed
%                   nfig: figure number for first image
%
% After detection, a 1-cell list of a matrix is returned in
% "peak". peak{1} contains the particle information
% for the present frame stored in a matrix:
%
%         peak{1}(:,1)    x (col)-positions of particles
%         peak{1}(:,2)    y (row)-positions of particles
%         peak{1}(:,3)    zero order intensity moments
%         peak{1}(:,4)    second order intensity moments
%         peak{1}(:,5)    *** unused ***
%         peak{1}(:,6)    *** empty *** to be filled by linker
%
% USES:  saturate.m
%
% Ivo Sbalzarini, 12.2.2003
% Institute of Computational Science, Swiss Federal
% Institute of Technology (ETH) Zurich.
% E-mail: sbalzarini@inf.ethz.ch
%
% based on an algorithm by Crocker & Grier:
%     Crocker, J.C. & Grier, D.G., Methods of digital video microscopy
%     for colloidal studies, J. colloid interface sci., 1996,
%     179: 298-310.
%======================================================================

function [peak,segImg] =  detect_particles(orig,w,v,AreaLevel_top,AreaLevel_bottom)
viz = v(1);
nfig = v(2);


% % some often used quantities
idx = [-w:1:w];     % index vector
dm = 2*w+1;         % diameter
im = repmat(idx',1,dm);
jm = repmat(idx,dm,1);
imjm2 = im.^2+jm.^2;
siz = size(orig);   % image size

%======================================================================
% STEP 1: Locating particles
%======================================================================
orig_grey=im2uint8(orig);

% calculate thresh for each part of the image on the grid
div_col=3; % size of the grid in col direction
div_row=3; % dix of the grid in row direction

step_row=floor(siz(1)/div_col);
step_col=floor(siz(2)/div_col);

row_index=[1:step_row];
col_index=[1:step_col];

thresh_list=[];
orig_grey_cut=orig(row_index,col_index);

for i=1:div_col
    if(i~=1) col_index=col_index+step_col;end
    row_index=[1:step_row];
    for j=1:div_row
        if(j~=1) row_index=row_index+step_row;end
        
        %[min(row_index),max(row_index),min(col_index),max(col_index)]
        orig_grey_cut=orig(row_index,col_index);
        [thresh,orig_bw_cut]=maxentropie(orig_grey_cut);
        thresh_list=[thresh_list,thresh];
    end;
end;

%check threshold proportion WhitePixel
thresh_list=unique(thresh_list);
WhitePixel_ratio=zeros(length(thresh_list),3);
for thr=1:length(thresh_list)
    % calculate pixel ration only to the top  100 layer of image (200 is
    % not enougth
    %this layer is air-liquid boarder
    WhitePixel_ratio(thr,1)=size(find(orig_grey(1:50,:)>thresh_list(thr)),1)/siz(1)/siz(2); % proportion across "air" layer
    WhitePixel_ratio(thr,2)=size(find(orig_grey(:)>thresh_list(thr)),1)/siz(1)/siz(2); % proportion across all image;
    WhitePixel_ratio(thr,3)=thresh_list(thr);
    
end;
%WhitePixel_ratio
CheckedIndexes=find(WhitePixel_ratio(:,1)<0.0001);

if isempty(CheckedIndexes)
    orig_bw=zeros(siz);
    'no proper thereshold is found'
    peak = {};
    segImg=orig_bw;
else
    thresh=min(WhitePixel_ratio(CheckedIndexes,3))
    orig_bw=orig_grey>thresh;
    thresh
    orig_label=bwlabel(orig_bw);
    stats=regionprops(orig_label,'Area','Centroid','PixelIdxList');
    Area=[stats.Area];
    Centroids = cat(1,stats.Centroid);
    idx=find(Area<AreaLevel_top & Area>AreaLevel_bottom);
    orig_select=zeros(size(orig));
    
    npart=length(idx);
    % orig_select=zeros(size(orig));
    % for ii=1:npart
    %     orig_select(stats(idx(ii)).PixelIdxList)=1;
    % end;
    % figure;imshow(orig_select);
    CentroidsNew=Centroids(idx,:);
    
    %======================================================================
    % STEP 2: Calculate zero and second order intensity moments of selected particles
    %======================================================================
    
    C=round(CentroidsNew(:,1));
    R=round(CentroidsNew(:,2));
    m0 = zeros(npart,1);
    m2 = zeros(npart,1);
    AreaSelect=Area(idx);
    
    % % generate circular mask of radius w
    mask = zeros(dm,dm);
    mask(find(imjm2 <= w*w)) = 1;
    
    % for each particle: compute zero and second order moments
    % lower and upper index bounds for all particle neighborhoods
    % in local coordinates.
    li = 1-(R-w-saturate(R-w,1,siz(1)));
    lj = 1-(C-w-saturate(C-w,1,siz(2)));
    ui = dm-(R+w-saturate(R+w,1,siz(1)));
    uj = dm-(C+w-saturate(C+w,1,siz(2)));
    
    
    
    for ipart=1:npart,
        % masked image part containing the particle
        Aij = orig(R(ipart)+li(ipart)-w-1:R(ipart)+ui(ipart)-w-1,...
            C(ipart)+lj(ipart)-w-1:C(ipart)+uj(ipart)-w-1).* ...
            mask(li(ipart):ui(ipart),lj(ipart):uj(ipart));
        % moments
        m0(ipart) = sum(sum(Aij));    % eq. [6]
        % eq. [7]
        m2(ipart) = sum(sum(imjm2(li(ipart):ui(ipart),lj(ipart):uj(ipart))...
            .*Aij))/m0(ipart);
    end;
    peak = zeros(npart,6);
    peak(:,1) = C;       % col position
    peak(:,2) = R;       % row position
    peak(:,3) = m0;      % zero order moment
    peak(:,4) = m2;      % second order moment
    peak(:,5)=AreaSelect;
    %======================================================================
    % STEP 4: Visualization
    %======================================================================
    viz=0;
    if viz == 1,
        % plot crosses at particle positions
        C = peak(:,1);
        R = peak(:,2);
        X = [[C'-2; C'+2], [C'; C']];
        Y = [[R'; R'], [R'-2; R'+2]];
        
        figure(nfig)
        imshow(imbinarize(orig,thresh))
        hold on
        hand = line(X,Y);
        set(hand(:),'Color',[1 0 0]);
        set(hand(:),'LineWidth',[3.0]);
        hold off
        nfig = nfig + 1
    end;
    
    peak = {peak};
    segImg=orig_bw;
end

return




