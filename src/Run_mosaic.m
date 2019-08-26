%======================================================================% 
%% Step 0: Setup
clear all
close all;


% Mosaik setup
addpath('MOSAIK')

% set up experiment folder and file 
% for MAc and Ubuntu
experiment_folder='/Volumes/mpistaff/Diaz_Pichugina_Pseudomona/Data/1-TIMELAPSES_2019_1-1/SM_1_03072019_FR';
%experiment_folder='/Volumes/mpistaff/Diaz_Pichugina_Pseudomona/Data/1-TIMELAPSES_2019_1-1/SM_1_03072019_FR';

%for Windows
experiment_folder='X:\Diaz_Pichugina_Pseudomona\Data\1-TIMELAPSES_2019_1-1\SM_1_03072019_FR';
filePattern = fullfile(experiment_folder, '*.tifFastRSegmentation.tif');
fileList = dir(filePattern);

Nfiles=size(fileList,1);

init=1;
final=100;
cd ../output/

%1:Nfiles


for i=1:Nfiles
    file_name_seg=fileList(i).name;
    file_name=split(file_name_seg,'_t1');
    file_name_native=strcat(file_name{1},'NATIVE.tif');
    disp(file_name_native)
    
    Prefix_file_writing=strsplit(file_name_native,'.');
    Prefix_file_writing=Prefix_file_writing{1}
    file_name_native=fullfile(experiment_folder,file_name_native);
    file_name_seg=fullfile(experiment_folder,file_name_seg);
    %----------------------------%
    % writing Log File
    diary (strcat(Prefix_file_writing,'_MosaikLog.txt'))
    %----------------------------%
    
    disp('Images properties')
    %pixel size mkm
    dx=0.160
    %time frame s
    dt=0.020
    
  
    %----------------------------%
    disp('Mosaik parameters')
    w =10          % size of circular mask to calculate moments
    trajLen=3     % minimum trajectory length in frames

    AreaLevel_top=1000 %select particals less than pixel in area
    AreaLevel_bottom=10 %select particals more than pixel in area
    LinkingDistance=25 %Linking distance in pixel 

    %----------------------------%
    
    %======================================================================%
    %% Step 1: Images preparation FTT filtering
    disp('set up filter')
    BoxFilter=10
    GausFilter_lambda=3
    %LowFreqBand=10;
    %HighFreqBand=500;
    NumberFiles=final-init+1;
    viz=0;
    %[imagesFTT,orig_images]=imagespreparation_FTT(file_name,init,final,viz,LowFreqBand,HighFreqBand);
    [images_restored,images_seg]=imagespreparation(file_name_native,file_name_seg,init,final,BoxFilter,GausFilter_lambda);
    %viz_image_stack(NumberFiles,orig_images)
    %viz_image_stack(NumberFiles,imagesFTT)
    
    %======================================================================%
    %% Step2: Peaks segmentation (define peaks on image)
    %  Peacks linking into trajectories across frames
    [peaks,SegmentedImageStack]=tracker(images_restored,w,AreaLevel_top,AreaLevel_bottom,LinkingDistance);

    
    % discard trajectory less than trajLen
    trajectories =ll2matrix(peaks,trajLen);
    if (size(trajectories,1)>0)
        % write to file
        file_name=strcat(Prefix_file_writing,'Trajectories.txt');
        write2file_trajectory(file_name,trajectories);
        %plot_traj_in_one_plot(imagesFTT,trajectories)
        
        % save output to the data folder
        outputFileName = strcat(Prefix_file_writing,'_TRAJ.tif');
        save_tracked_images_tiff_stack(images_restored,trajectories,outputFileName);
        
        % save segemented to the data folder
        outputFileName = strcat(Prefix_file_writing,'_SEG_TRAJ.tif');
        save_tracked_images_tiff_stack(SegmentedImageStack,trajectories,outputFileName);

    end
    
    
    %% Step3: Analyze trajectories
    %calculate diffusion coefficient and MSS
    alysis_matrix = moments(trajectories,dx,dt);
    file_name=strcat(Prefix_file_writing,'Analysis.txt');
    write2file_analysis(file_name,alysis_matrix);
    t1=toc
    diary off;
    clear orig_images
    clear imagesFTT
    clear SegmentedImageStack
    close all
end
