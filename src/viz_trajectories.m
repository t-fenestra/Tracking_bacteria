%====================================================================== 
%
% viz_trajectories: visualize trajectory of individual particles 
% SYNTAX:  viz_trajectories(NumberFiles,Image_stack,trajectories)
%
% INPUTS:   NumberFiles     files number in stack
%           Image_stack     stack with images
%           trajectories    cell list with trajectories
%           trajectories{t}(:,1)    frame number of traj. t
%           trajectories{t}(:,2)    x postitions of traj. t
%           trajectories{t}(:,3)    y positions of traj. t
%
%
% updated 21.12.2018
%
%====================================================================== 


function viz_trajectories(NumberFiles,Image_stack,trajectories)
    videofig(NumberFiles, @(frames) redraw_trajectories(frames, Image_stack,trajectories));
    % Display initial frame
    redraw_trajectories(1,Image_stack,trajectories)
end 