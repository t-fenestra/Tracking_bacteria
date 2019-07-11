%====================================================================== 
%
% viz_image_stack: visualize trajectory of individual particles 
% SYNTAX:  viz_image_stack(NumberFiles,Image_stack)
%
% INPUTS:   NumberFiles     files number in stack
%           Image_stack     stack with images
% updated 21.12.2018
%
%====================================================================== 
function viz_image_stack(NumberFiles,Image_stack)
    videofig(NumberFiles, @(frames) redraw(frames, Image_stack));
    % Display initial frame
    redraw(NumberFiles,Image_stack);
end