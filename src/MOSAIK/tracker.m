%====================================================================== 
%
% TRACKER: track virus particles in a time series of images
%
% SYNTAX:  peaks = tracker(images,w,pth,cutoff,L)
%
% INPUTS:  images       stack of images
%          w            radius of neighborhood: > particle radius, < interparticle distance
%          pth          upper intensity percentile for particles
%          cutoff       probability cutoff for non-particle discrimination
%          L            maximum displacement between frames
%
% After detection, a cell list of matrices is returned in
% "peaks". peaks{iframe} contains the particle information
% for frame i stored in a matrix. Two different routines
% both for feature detection and particle matching (linking
% of the trajectories) are available:
%
%    Detection: detect_particles
%    Linking:   link_trajectories
%
% The final output that is returned is:
%         peaks{iframe}(:,1)    x (col)-positions of particles
%         peaks{iframe}(:,2)    y (row)-positions of particles
%         peaks{iframe}(:,3)    zero order intensity moments
%         peaks{iframe}(:,4)    second order intensity moments
%         peaks{iframe}(:,5)    *** unused ***
%         peaks{iframe}(:,6)    link list index of same part.
%                               in next frame. -1 if none.

%         w =  % radius of neighborhood: > particle radius, 
               % < interparticle distance
%
%
% This cell list is also stored in the file "trackdata" for later use.
%
% USES:    detect_particles, link_trajectories
%
% based on matlab code by Ivo Sbalzarini, Feb. 12, 2003
% update 21.12.2018
%====================================================================== 

function [peaks,SegmentedImageStack] = tracker(images,w,LinkedDistance)

nimg=size(images,3);
siz=size(images(:,:,1));
viz = 0;

%=============================================================
% detect particles in all image frames
%=============================================================
peaks = [];
nfig=1;
SegmentedImageStack=zeros(siz(1),siz(2),nimg);



for img=1:nimg,
    disp(sprintf('\nParticle recoginition in image %d of %d',img,nimg))
    if viz == 1,
       figure(nfig)
       nfig = nfig + 1;
       imshow(images(:,:,img))
       title('Original micrograph');
    end;
     
   viz=0;
   
   [peak,segmImg] = detect_particles(images(:,:,img),w,[viz,nfig]);
   peaks = [peaks, peak];
   SegmentedImageStack(:,:,img)=segmImg;
end;


%=============================================================
% assemble paths across frames as linked list
%=============================================================
% L maximum displacement between frames determined by radial distribution
% function
disp('Linking distance')
L=LinkedDistance;

if ~isempty(peaks)
    peaks = link_trajectories(peaks, L, viz, 100);
end
% save data for later use
%save trackdata peaks;

%=============================================================
% visualize paths  
%=============================================================
%viz=0;
%orig = images(:,:,1);
%[orig] = normalize(orig);
%nframe = length(peaks);

% figure(200)
% clf
% imshow(orig)
% hold on
% 
% C = peaks{1}(:,1);
% R = peaks{1}(:,2);
% X = [[C'-2; C'+2], [C'; C']];
% Y = [[R'; R'], [R'-2; R'+2]];
% hand = line(X,Y);
% set(hand(:),'Color',[1 0 0]);
% for iframe=2:nframe,
%     oldind = find(peaks{iframe-1}(:,6)>0);
%     curind = peaks{iframe-1}(oldind,6);
%     X = [peaks{iframe-1}(oldind,1), peaks{iframe}(curind,1)];
%     Y = [peaks{iframe-1}(oldind,2), peaks{iframe}(curind,2)];
%     hand = line(X',Y');
%     set(hand(:),'Color',[1 0 0]);
%     set(hand(:),'LineWidth',[1.0]);
% end;
% hold off
