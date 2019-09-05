%======================================================================% 
%% Step 0: Setup
clear all
close all;


% Mosaik setup
addpath('MOSAIK')

%% Remote folder with images
% for Mac
%MainFolder='/Volumes/mpistaff/Diaz_Pichugina_Pseudomona/Data/1-TIMELAPSES_2019_1-1';
% for Windows
MainFolder='X:\Diaz_Pichugina_Pseudomona\Data\1-TIMELAPSES_2019_1-1\'

DirContent=dir(MainFolder);
% discard hidden folder and Dicarded trajectories
DirList=setdiff({DirContent.name},{'.','..','.DS_Store','Discarded'});

%% Folder to save results
%Mac
%ResultFolder='/Users/pichugina/Work/Data_Analysis/MOSAIK_traking/Traking_bacteria_git/output/'
%Windows
ResultFolder='C:\Work\output';

%mkdir(ResultFolder)

%% Parameters of the stack
init=1;
final=100;
% scan across experiments
for k=6:length(DirList)
    % Read file list in directory
    fprefix=DirList{k}
    ImageFolder=fullfile(MainFolder,fprefix);
    filePattern = fullfile(ImageFolder, '*NATIVE.tif');
    TifFiles = dir(filePattern);
    files={TifFiles.name};
    
    if (~isempty(files))     % check if the file list empty to avoid mistake
        % sort files according to the data frame
        files=natsortfiles(files);
        Nfiles=length(files);
        %Create folder to save results
        mkdir(ResultFolder,fprefix)
        cd(fullfile(ResultFolder,fprefix))
        
        
        % calculate till 15 hours
        for Sregime=1:15
            file_name=files(Sregime);
            file_name=file_name{1};
            Prefix_file_writing=strsplit(file_name,'.');
            Prefix_file_writing=Prefix_file_writing{1}
            disp(file_name)
            file_name=fullfile(ImageFolder,file_name);
            
            
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
            LinkingDistance=15 %Linking distance in pixel
            %----------------------------%
            
            %======================================================================%
            %% Step 1: Images preparation FTT filtering
            disp('set up filter')
            BoxFilter=20
            GausFilter_lambda=3
            NumberFiles=final-init+1;
            viz=0;
            images_restored=imagespreparation(file_name,init,final,BoxFilter,GausFilter_lambda);
            %viz_image_stack(NumberFiles,orig_images)
            %viz_image_stack(NumberFiles,imagesFTT)
            
            %% Step2: Peaks segmentation (define peaks on image)
            %  Peacks linking into trajectories across frames
            [peaks,SegmentedImageStack]=tracker(images_restored,w,AreaLevel_top,AreaLevel_bottom,LinkingDistance);
            
            %%% ??heck if peaks empty
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
            else
                disp('programme cannot find particles across the imagestack')
            end
            diary off;
            clear images_restored
            clear peaks
            clear SegmentedImageStack
            close all

        end;
        cd(ResultFolder)
    end;
    
  
end;