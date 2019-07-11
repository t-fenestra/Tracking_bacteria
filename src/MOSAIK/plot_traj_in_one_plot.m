
function plot_traj_in_one_plot(images,trajectories)
       
     Nfiles=size(images,3);
     figure(1); imshow(images(:,:,Nfiles),[]);
     set(gca, 'DataAspectRatioMode', 'auto')
     hold on;
       
     %trajectory running in red
     for t=1:length(trajectories)
         traj=trajectories{t};
         traj=traj(:,2:3);
         plot(traj(:,1),traj(:,2),'r-','LineWidth',1);
     end
     hold off;

end