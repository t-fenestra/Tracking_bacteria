%% Create movie from images
% code from https://de.mathworks.com/matlabcentral/answers/153925-how-to-make-a-video-from-images
function [] = create_movie(ImageStack,trajectories,OutputFile)
        % ImageStack multidimensional array of Images(:,:,frame)
        % OutputFile  video file frame
        % grey_map grey color map same for all images
        
        % fps video frame per second frame
        fps=1;
        % create the video writer with 1 fps
        
        writerObj = VideoWriter(strcat(OutputFile,'.avi'),'Motion JPEG AVI');
        writerObj.FrameRate = fps;
       

        % open the video writer
        open(writerObj);
 
        % write the frames to the video

        
        for frame=1:size(ImageStack,3)
            frame
            redraw_trajectories_movie(frame, ImageStack,trajectories)
            Videoframe = getframe(gcf);
           writeVideo(writerObj, Videoframe);  
        end
   close(writerObj);     
end
