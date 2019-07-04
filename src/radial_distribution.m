function [FistPeak,Output]=radial_distribution(im,thresh)
    Lx=size(im,1);
    Ly=size(im,2);
    
    % calculate centrioids of image particles
    bw=imbinarize(im,thresh);
    CC=bwconncomp(bw);
    S = regionprops(CC,'Centroid');
    centroids = cat(1, S.Centroid);
    
    figure; imshow(bw);
    hold on;
    plot(centroids(:,1),centroids(:,2),'r*')
    hold off;
    
    
    nPart=size(centroids,1);
    rho = nPart/(Lx*Ly); % general density of the particles
    if (nPart<5)
        warning('too few particles to calculate radial distribution function');
    end;    
    
    % pair-pair distance
    count=1;
    for partA=1:size(centroids,1)-1
        for partB=(partA+1):size(centroids,1)
            % calculate particle-particle distance
            dr=centroids(partA,:) - centroids(partB,:);
            % Get the size of this distance vector
            r(count) = sqrt(sum(dr.*dr));
            count=count+1;
        end
    end

    
% we are interested at the particle-particle distance less than 100 px
r=r(r>10 & r<100);
figure()
hr=histogram(r,100);
counts=hr.Values;
Edges=hr.BinEdges;
position=Edges+hr.BinWidth/2; % center of bins
normalize=rho*pi*(Edges(2:end).^2-Edges(1:end-1).^2);
counts_normilize=counts./normalize;

% plot(position(1:end-1),counts_normilize,'r');
% xlabel('distance');
% ylabel('g(R)');

% simple median smoothing
RDD = movmedian(counts_normilize,3);
[maxv,index]=max(RDD);
FistPeak=position(index);
Output=[position(1:end-1)',RDD'];

figure();
plot(position(1:end-1),RDD,'b','LineWidth',3);
xlabel('distance');
ylabel('g(R)');
title('Radial Distribution Function');
saveas(gcf,'RDF.png');
    
end