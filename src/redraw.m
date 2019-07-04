%====================================================================== 
%
% redraw: visualize images per frame 
% 
% SYNTAX:  viz_trajectories(frame,Image_stack)
%
% INPUTS:   frame           frame number
%           Image_stack          stack with images
% updated 21.12.2018
%
%====================================================================== 

function redraw(frame,Image_stack)
        imshow(Image_stack(:,:,frame),[]);
end