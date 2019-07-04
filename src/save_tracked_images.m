function [] = save_tracked_images(ImageStack,trajectories,folder)
        
        for frame=1:size(ImageStack,3)
            frame
            redraw_trajectories_movie(frame, ImageStack,trajectories)
            
            F= getframe(gcf);
            [X, Map] = frame2im(F);
            imwrite(X,sprintf('%s/Result_frame-%i.tif',folder,frame))
        end
end
