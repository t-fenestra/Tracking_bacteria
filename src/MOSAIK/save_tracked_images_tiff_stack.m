function save_tracked_images_tiff_stack(ImageStack,trajectories,outputFileName)
        
        for frame=1:size(ImageStack,3)
            %frame
            redraw_trajectories_movie(frame, ImageStack,trajectories)
            
            F= getframe(gcf);
            [X, Map] = frame2im(F);
            imwrite(X, outputFileName,'WriteMode', 'append', 'Compression','none')
        end
end
