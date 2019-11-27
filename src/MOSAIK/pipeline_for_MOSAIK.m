function pipeline_for_MOSAIK(Prefix_file_writing,filename,init,final,w,trajLen,LinkingDistance)

%Step 1: Images preparation 
disp('set up filter')
images_tifthresh=read_image_stack(filename,init,final);
%viz_image_stack(NumberFiles,orig_images)
%viz_image_stack(NumberFiles,imagesFTT)
            
% Step2: Peaks segmentation (define peaks on image)
%  Peacks linking into trajectories across frames
[peaks,SegmentedImageStack]=tracker(images_tifthresh,w,LinkingDistance);
%ALI=approximate_ALI(images_tifthresh(:,:,1));
            
%%% check if peaks empty
if(~isempty(peaks))
    
    % discard trajectory less than trajLen
    trajectories =ll2matrix(peaks,trajLen);
    if (size(trajectories,1)>0)
        % write to file
        file_name=strcat(Prefix_file_writing,'Trajectories.txt');
        write2file_trajectory(file_name,trajectories);
        %plot_traj_in_one_plot(imagesFTT,trajectories)
        
        % save output to the data folder
        outputFileName = strcat(Prefix_file_writing,'_TRAJ.tif');
        
        save_tracked_images_tiff_stack(images_tifthresh,trajectories,outputFileName);
        
        % save segemented to the data folder
        outputFileName = strcat(Prefix_file_writing,'_SEG_TRAJ.tif');
        save_tracked_images_tiff_stack(SegmentedImageStack,trajectories,outputFileName);
        
        % save ALI coordinats
        %outputFileName = strcat(Prefix_file_writing,'_ALI.txt');
        %dlmwrite(outputFileName,ALI);
    end
    
                
    %Step3: Analyze trajectories
    %calculate diffusion coefficient and MSS
    %alysis_matrix = moments(trajectories,dx,dt);
    %file_name=strcat(Prefix_file_writing,'Analysis.txt');
    %write2file_analysis(file_name,alysis_matrix);
    
else
    disp('programme cannot find particles across the imagestack')
end
end