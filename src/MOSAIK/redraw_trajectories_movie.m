%====================================================================== 
%
% redraw_trajectories: visualize trajectory of individual particles on
% individual frame
% 
% SYNTAX:  viz_trajectories(frame,images,trajectories)
%
% INPUTS:   frame           frame number
%           images          stack with images
%           trajectories    cell list with trajectories
%           trajectories{t}(:,1)    frame number of traj. t
%           trajectories{t}(:,2)    x postitions of traj. t
%           trajectories{t}(:,3)    y positions of traj. t
%
%
% updated 21.12.2018
%
%====================================================================== 


function redraw_trajectories_movie(frame,images,trajectories,ALI)
       [StartFrame,EndFrame]=cellfun(@trajectory_start_end_frame,trajectories);
       
       % trajectories finished before frame
       IdxFinishedBefore=find(((EndFrame-frame)<0) & ((StartFrame-frame)<0));
       
       % trajectories is still running 
       IdxCurrent=find(((EndFrame-frame).*(StartFrame-frame)<=0));
       
       
       figure(1); imshow(images(:,:,frame),[]);
       set(gca, 'DataAspectRatioMode', 'auto')
       hold on;
       
       % visualize ALI
       %plot([1:2048]',ones(2048,1)*median(ALI),'cyan','LineWidth',2);
       
       %trajectory running in red
       for t=1:length(IdxCurrent)
           traj=trajectories{IdxCurrent(t)};
            traj=traj(traj(:,1)<=frame,2:3);
            plot(traj(1,1),traj(1,2),'x','MarkerSize',3,'Color','r');
            plot(traj(:,1),traj(:,2),'Color','r');
       end
        
       %trajectory finished  in green
       for t=1:length(IdxFinishedBefore)
           traj=trajectories{IdxFinishedBefore(t)};
           traj=traj(traj(:,1)<=frame,2:3);
           plot(traj(1,1),traj(1,2),'x','MarkerSize',3,'Color','g');
           plot(traj(:,1),traj(:,2),'Color','g');
       end
       
       set(gca, 'unit', 'normalize')
       set(gca, 'position', [0 0 1 1]);
       
       
      hold off;
end