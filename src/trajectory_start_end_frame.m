%====================================================================== 
%
% trajectory_start_end_frame: find start and end frame of trajectory 
% SYNTAX:  trajectory_start_end_frame(trajectory)
%
% INPUTS:   trajectory     array
%           trajectory[1]  frame number
%           trajectory[2]  x coordinate
%           trajectory[3]  y coordinate
% updated 21.12.2018
%
%====================================================================== 

function [StartFrame,EndFrame]=trajectory_start_end_frame(trajectory)
        StartFrame=min(trajectory(:,1));
        EndFrame=max(trajectory(:,1));
end