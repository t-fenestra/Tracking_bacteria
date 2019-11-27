%======================================================================% 
%% Step 0: Setup
clear all
close all;


% Mosaik setup
addpath('MOSAIK')

%% Remote folder with images
% for Mac
MainFolder='/Volumes/mpistaff/Diaz/1-staticmotile-1/';
% for Windows
%MainFolder='X:\Diaz_Pichugina_Pseudomona\Data\1-TIMELAPSES_2019_1-1\'

DirContent=dir(MainFolder);
%discard not folder files
NotFolder=dir(strcat(MainFolder,'dataset*'))
% discard hidden folder and Dicarded trajectories

DirList=setdiff({DirContent.name},{'.','..','.DS_Store','Discarded',NotFolder.name});
%DirList={'SM_1_04252019_FR'};
%% Folder to save results
%Mac
ResultFolder='/Users/pichugina/Work/Data_Analysis/MOSAIK_traking/Traking_bacteria_git/output/stats_motile/'
%Windows
%ResultFolder='C:\Work\output';

%mkdir(ResultFolder)

%% Parameters of the stack
init=1;
final=10;
% scan across experiments
for k=4:length(DirList)
    % Read file list in directory
    fprefix=DirList{k}
    ImageFolder=fullfile(MainFolder,fprefix);
    
    %extract all projection images
    filePattern = fullfile(ImageFolder, '*NATIVE.tifCALC.tifCALC.tifthresh.tif_projection.tif');
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
            %AreaLevel_top=1000 %select particals less than pixel in area
            %AreaLevel_bottom=10 %select particals more than pixel in area
            LinkingDistance=15 %Linking distance in pixel
           
            %----------------------------%
            % motile part
            Prefix_file_writing_motile=strcat(Prefix_file_writing,'_MOTILE')
            filename=strcat(MainFolder,fprefix,'/',Prefix_file_writing,'.tifCALC.tifCALC.tifthresh.tif_MOTILE.tif');
            pipeline_for_MOSAIK(Prefix_file_writing_motile,filename,init,final,w,trajLen,LinkingDistance)
            
            Prefix_file_writing_static=strcat(Prefix_file_writing,'_STATIC')
            filename=strcat(MainFolder,fprefix,'/',Prefix_file_writing,'.tifCALC.tifCALC.tifthresh.tif_STATIC.tif');
            pipeline_for_MOSAIK(Prefix_file_writing_static,filename,init,final,w,trajLen,LinkingDistance)
            
            diary off;
            clear images_tifthresh
            clear peaks
            clear SegmentedImageStack
            close all

        end;
        cd(ResultFolder)
    end;
    
  
end;