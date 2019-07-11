%======================================================================% 
%% Step 0: Setup
clear all
close all;


% Mosaik setup
addpath('MOSAIK')

% set up experiment folder and file 
% for MAc and Ubuntu
% experiment_folder='/Volumes/mpistaff/Diaz_Pichugina_Pseudomona/Data/1-TIMELAPSES_2019_1-1/SM_1_03072019_FR';

%for Windows
experiment_folder='X:\Diaz_Pichugina_Pseudomona\Data\1-TIMELAPSES_2019_1-1\SM_1_03072019_FR';
fast_regime=15;
file_name=sprintf('SM_1_03072019_FR_%dNATIVE.tif',fast_regime);

Prefix_file_writing=strsplit(file_name,'.');
Prefix_file_writing=Prefix_file_writing{1}
file_name=fullfile(experiment_folder,file_name);
init=1;
final=5;
%----------------------------%

% writing Log File
cd ../output
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
AreaLevel_top=200 %select particals less than 150 pixel in area
AreaLevel_bottom=w %select particals more than 10 pixel in area
%----------------------------%


%======================================================================% 
%% Step 1: Images preparation FTT filtering
disp('set up FTT filter')
LowFreqBand=10
HighFreqBand=500
NumberFiles=final-init+1;
viz=0;
[imagesFTT,orig_images]=imagespreparation_FTT(file_name,init,final,viz,LowFreqBand,HighFreqBand);
%viz_image_stack(NumberFiles,orig_images)
%viz_image_stack(NumberFiles,imagesFTT)

%======================================================================% 
%% Step2: Peaks segmentation (define peaks on image)
%  Peacks linking into trajectories across frames
[peaks,SegmentedImageStack]=tracker(imagesFTT,w,AreaLevel_top,AreaLevel_bottom);

% discard trajectory less than trajLen 
trajectories =ll2matrix(peaks,trajLen);

% write to file
file_name=strcat(Prefix_file_writing,'Trajectories.txt');
write2file_trajectory(file_name,trajectories);
%plot_traj_in_one_plot(imagesFTT,trajectories)

% save output to the data folder
outputFileName = strcat(Prefix_file_writing,'_TRAJ.tif');
save_tracked_images_tiff_stack(imagesFTT,trajectories,outputFileName);


%% Step3: Analyze trajectories
%calculate diffusion coefficient and MSS
alysis_matrix = moments(trajectories,dx,dt);  
file_name=strcat(Prefix_file_writing,'Analysis.txt');
write2file_analysis(file_name,alysis_matrix);
diary off;

