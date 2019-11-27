%======================================================================% 
%% Step 0: Setup
clear all
close all;


% Mosaik setup
addpath('MOSAIK')

%% Remote folder with images
% for Mac
MainFolder='/Volumes/mpistaff/Diaz/1-datatests-1/';
% for Windows
%MainFolder='X:\Diaz_Pichugina_Pseudomona\Data\1-TIMELAPSES_2019_1-1\'

DirContent=dir(MainFolder);
%discard not folder files
NotFolder=dir(strcat(MainFolder,'dataset*'))
% discard hidden folder and Dicarded trajectories

DirList=setdiff({DirContent.name},{'.','..','.DS_Store','Discarded',NotFolder.name});

%% Folder to save results
%Mac
ResultFolder='/Users/pichugina/Work/Data_Analysis/MOSAIK_traking/Traking_bacteria_git/output/'
%Windows
%ResultFolder='C:\Work\output';

%mkdir(ResultFolder)

%% Parameters of the stack
init=1;
final=100;
% scan across experiments
for k=14 %6:length(DirList)
    % Read file list in directory
    fprefix=DirList{k}
    ImageFolder=fullfile(MainFolder,fprefix);
    
    %for calculation moments
    filePattern = fullfile(ImageFolder, '*NATIVE.tifCALC.tifthresh.tif');
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
        for Sregime=10
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
            thresh=55 % Threshold for Andres
            %----------------------------%
            
            %======================================================================%
            %% Step 1: Images preparation 
            disp('set up filter')
            images_tifthresh=read_image_stack(file_name,init,final);
            ALI=approximate_ALI(images_tifthresh)
            
            outputFileName=strcat(Prefix_file_writing,'_ALI.tif');

            for frame=1:size(images_tifthresh,3)
                frame
                figure(1); imshow(images_tifthresh(:,:,frame),[]);
                set(gca, 'DataAspectRatioMode', 'auto')
                hold on;
                plot([1:2048]',ALI,'red','LineWidth',3)
                hold off;
                
                F= getframe(gcf);
                [X, Map] = frame2im(F);
                imwrite(X, outputFileName,'WriteMode', 'append', 'Compression','none')
            end      
            
            %viz_image_stack(NumberFiles,orig_images)
            %viz_image_stack(NumberFiles,imagesFTT)
            
           
            close all

        end;
        cd(ResultFolder)
    end;
    
  
end;