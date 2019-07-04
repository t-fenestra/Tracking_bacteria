function filter3= gaussianbpf(imsiz,d0,d1)
%   d0 = Lower cut off frequency
%   d1 = Higher cut off frequency
%
% The function makes use of the simple principle that a bandpass filter
% can be obtained by multiplying a lowpass filter with a highpass filter
% where the lowpass filter has a higher cut off frquency than the high pass filter.
%
% Usage GAUSSIANBPF(I,DO,D1)
% Example
% ima = imread('grass.jpg');
% ima = rgb2gray(ima);
% filtered_image = gaussianbpf(ima,30,120);
% Gaussian Bandpass Filter

nx=imsiz(1);
ny=imsiz(2);
% Initialize filter.
filter1 = ones(2*nx-1,2*ny-1);
filter2 = ones(2*nx-1,2*ny-1);
filter3 = ones(2*nx-1,2*ny-1);
for i = 1:2*nx-1
    for j =1:2*ny-1
        dist = ((i-(nx+1))^2 + (j-(ny+1))^2)^.5;
        % Use Gaussian filter.
        filter1(i,j) = exp(-dist^2/(2*d1^2));
        filter2(i,j) = exp(-dist^2/(2*d0^2));
        filter3(i,j) = 1.0 - filter2(i,j);
        filter3(i,j) = filter1(i,j).*filter3(i,j);
    end
end

filter3=ifftshift(filter3);
return