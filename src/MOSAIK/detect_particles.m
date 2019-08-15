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

function [peak,segImg,FirstPeak] =  detect_particles(orig,w,v,AreaLevel_top,AreaLevel_bottom,FirstPeak)

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

% determining upper pth-th percentile of intensity values
disp('determining upper pth-th percentile of intensity values');
pth=0.10
[cnts,bins] = imhist(orig);
l = length(cnts);
k = 1;
while sum(cnts(l-k:l))/sum(cnts) < pth,
    k = k + 1;
end;
%thresh= bins(l-k+1);



thresh_hist = bins(l-k+1);
% % proportion ones to zeros
thresh_hp=length(find(orig<thresh_hist))/length(find(orig>thresh_hist))
thresh_outsu = adaptthresh(orig); %graythresh(orig);
thresh_op=length(find(orig<thresh_outsu))/length(find(orig>thresh_outsu))

% kernel = [-1, -1, -1; -1, 8, -1; -1, -1,-1]/8;
% diffImage = conv2(orig, kernel, 'same');
% orig=orig+abs(diffImage);
% cpp = median(diffImage(:))

if (thresh_op>thresh_hp)
    thresh=thresh_outsu;
else
    thresh=thresh_hist;
end;

orig_bw=orig>thresh;
% ones_image=length(find(orig_bw>0));
% zeros_image=length(find(orig_bw==0));
% thresh_proportion=zeros_image/ones_image;

% if thresh_proportion<5
%     J=histeq(orig);
%     imshow(J);
% end

orig_bw=bwlabel(orig_bw);
%imshow(orig_bw);
stats=regionprops(orig_bw,'Area','Centroid','PixelIdxList');
Area=[stats.Area];
Centroids = cat(1,stats.Centroid);
idx=find(Area<AreaLevel_top & Area>AreaLevel_bottom);
orig_select=zeros(size(orig));

npart=length(idx);

%for ii=1:npart
%    orig_select(stats(idx(ii)).PixelIdxList)=1;
%end;


% %======================================================================
%% radial_distribution
CentroidsNew=Centroids(idx,:);
%[CurrentFisrtPeak,Output]=radial_distribution(orig,CentroidsNew);
%CentroidsNew=int64(CentroidsNew);
%%figure,imshow(orig_select),title('Bacteria Area')

FirstPeak=15; %[FirstPeak,CurrentFisrtPeak];
%====================================================================== 
% STEP 2: Calculate zero and second order intensity moments of selected particles
%======================================================================

C=round(CentroidsNew(:,1));
R=round(CentroidsNew(:,2));
m0 = zeros(npart,1);
m2 = zeros(npart,1);

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
%

% %====================================================================== 
% % STEP 3: Non-particle discrimination
% %====================================================================== 
% 
% %-------------------------------------------%
% % check if particles have overlapped neighbourhood
% % if yes it indicates that peacks refinment may be not valid
% d2=w^2;
% for i=1:length(R)
%     for j=(i+1):length(R)
%         dist=(R(i)-R(j))^2+(C(i)-C(j))^2;
%         if dist<d2
%             disp(sprintf('Warning in peaks refinement: peaks i=%d j=%d overlapped',i,j))
%             hold on;
%             imshow(particles)
%             plot(R(i),C(i),'ro');
%             plot(R(j),C(j),'go');
%             hold off;
%         end;
%     end;    
%             
% end
% %-------------------------------------------%
% sigx = 0.1;
% sigy = 0.1;
% prob = zeros(size(m0));
% Nm = length(m0);
% 
% 
% for i=1:Nm,
%     prob(i)=sum(exp(-((m0(i)-m0).^2./(2*sigx*sigx))-((m2(i)-m2).^2./...
%         (2*sigy*sigy)))/(2*pi*sigx*sigy*Nm));
% end;
%     
% if viz == 1,
%     figure(nfig)
%     clf;
%     nfig = nfig + 1;
%     subplot(2,2,1)
%     hold on
%     m0in = m0(find(prob >= cutoff));
%     m2in = m2(find(prob >= cutoff));
%     plot(m0in,m2in,'go')
%     m0in = m0(find(prob < cutoff));
%     m2in = m2(find(prob < cutoff));
%     plot(m0in,m2in,'ro')
%     hold off
%     xlabel('m0')
%     ylabel('m2')
%     subplot(2,2,2)
%     hist(m0,50)
%     xlabel('m0')
%     subplot(2,2,3)
%     hist(m2,50)
%     xlabel('m2')
% end;
% 
% % indices of valid particles
% tmp = find(prob>=cutoff);  
% % pack data into return value
% npart = length(tmp);
% peak = zeros(npart,6);
% peak(:,2) = R(tmp);       % row position
% peak(:,1) = C(tmp);       % col position
% peak(:,3) = m0(tmp);      % zero order moment
% peak(:,4) = m2(tmp);      % second order moment
% field 5: unused
% field 6: used by linker to store linked list indices


peak = zeros(npart,6);
peak(:,1) = C;       % col position
peak(:,2) = R;       % row position
peak(:,3) = m0;      % zero order moment
peak(:,4) = m2;      % second order moment
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
segImg=imbinarize(orig,thresh);

  
return

% % generate circular mask of radius w
% mask = zeros(dm,dm);
% mask(find(imjm2 <= w*w)) = 1;
% 
% 
% % identify individual particles as local maxima in a
% % w-neighborhood that are larger than thresh
% dil = imdilate(orig,mask);
% [Rp,Cp] = find((dil-orig)==0);
% particles = zeros(siz);
% V = find(orig(sub2ind(siz,Rp,Cp))>thresh);
% R = Rp(V); 
% C = Cp(V);
% particles(sub2ind(siz,R,C)) = 1;
% npart = length(R);
% 
% viz=0;
% if viz == 1,
%     figure(nfig)
%     nfig = nfig + 1;
%     imshow(orig,[])
%     hold on;
%     plot(C,R,'r+');
%     title('intensity maxima of particles');
% end;


