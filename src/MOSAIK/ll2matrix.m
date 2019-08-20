%====================================================================== 
%
% LL2MATRIX: Converts trajectories in the linked-list representation
%            to matrix form
%
% SYNTAX:  take = ll2matrix(peaks)
%
% INPUTS:  peaks   trajectories in linked list form as:
%
%         peaks{t}(:,1)    x (col)-positions at time t
%         peaks{t}(:,2)    y (row)-positions at time t
%         peaks{t}(:,3)    zero order intensity moments
%         peaks{t}(:,4)    second order intensity moments
%         peaks{t}(:,5)    area size
%         peaks{t}(:,6)    linked list index to the same
%                          particle at time t+1
%         trajLen          minimum trajectory length to be further
%                          processed
%
%
%
% OUTPUTS:  trajectories more than trajlen length :
%         take{i}[:,1] number frames
%         take{i}[:,2] x position
%         take{i}[:,3] y position
%         take{i}[:,4] moment 1
%         take{i}[:,5] moment 2
%         take{i}[:,6] moment segmented area size

% The function returns a cell list of matrices where matrices{i}
% is the i-th trajectory in the form of an N times 2 matrix
% where N is the length of the trajectory and each row contains
% a [x,y] vector of positions.
%
% Ivo Sbalzarini, 26.3.2003
% Institute of Computational Science, Swiss Federal
% Institute of Technology (ETH) Zurich. 
% E-mail: sbalzarini@inf.ethz.ch

%====================================================================== 

function take = ll2matrix(peaks,trajLen)

% convert linked list trajectories into a list of (x,y) matrices
matrices = [];
% matrices_m0=[];
% matrices_m2=[];
for ii=1:length(peaks)         % loop over all frames
    npart = length(peaks{ii}(:,1));
    for ipart=1:npart         % loop over all particles
	iframe = ii;
	next = peaks{iframe}(ipart,6);
	if (next > 0)         % if particle starts a trajectory,
	                       % follow it
	    matrix = [iframe,peaks{iframe}(ipart,1), peaks{iframe}(ipart,2),peaks{iframe}(ipart,3),peaks{iframe}(ipart,4),peaks{iframe}(ipart,5)];
	    peaks{iframe}(ipart,6) = -1;   % mark used
	    while (next > 0)  % convert to matrix form
		iframe = iframe + 1;
		matrix = [matrix; iframe, peaks{iframe}(next,1), peaks{iframe}(next,2),peaks{iframe}(next,3),peaks{iframe}(next,4),peaks{iframe}(next,5)];
		nextold = next;
		next = peaks{iframe}(next,6);    % mark used
		peaks{iframe}(nextold,6) = -1;
	    end
	    matrix = {matrix};
	    matrices = [matrices, matrix];
	end
    end;
end;

% determine which trajectories to include in the analysis.
% Only take those which have no outliers in the step length histogram (they
% usually correspond to wrong tracking assignments) and move by more than 1
% pixel per frame (others are considered stationary and excluded from
% motion analysis).


% choose trajectory more than 5 points in length
take = [];
for itraj=1:length(matrices),
    traj = matrices{itraj};
    tlen = size(traj,1);
    slen = sqrt((traj(2:tlen,2)-traj(1:tlen-1,2)).^2+(traj(2:tlen,3)-traj(1:tlen-1,3)).^2);
    % the tracker accuracy is about 0.2 pixel
    %if and(mean(slen) >= 0.3, std(slen) < 1),
    if tlen>trajLen
	take = [take, {matrices{itraj}}];
    end;
%     figure(1)
%     hist(slen,10);
%     pause 
end;
%disp(sprintf('%d trajectories excluded from analysis. Remaining: %d',length(matrices)-length(take),length(take)));    
return

