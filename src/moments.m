%====================================================================== 
%
% MOMENTS: Calculates the scaling coefficients of the moments of the
%          particle displacement over time.
%
% SYNTAX:  gamma = moments(trajectories,dx,dt)
%
% INPUTS:   trajectories     vector of cells defining the 
%                           trajectories of all particles in the
%                           following synatx:
%
%           trajectories{t}(:,1)    frame number of traj. t
%           trajectories{t}(:,2)    x postitions of traj. t
%           trajectories{t}(:,3)    y positions of traj. t
%
%           dx               length of a pixel in physical units of um
%           dt               time between frames in sec
% OUPUTS:
% The function returns a matrix TrajectoryData 
%           TrajectoryData[,1] trajectory ID
%           TrajectoryData[,2] trajectory frames length 
%           TrajectoryData[,3] trajectory Diffusion coefficient mkm^2/s
%           TrajectoryData[,4] trajectory RMS Diffusion
%           TrajectoryData[,5] trajectory MSS slope
%           TrajectoryData[,6] trajectory RMS MSS
%           
%
%
% based on matlab version of Mosaik by Ivo Sbalzarini, 30.7.2003
% updated 21.12.2018
%
%====================================================================== 

function  TrajectoryData = moments(trajs,dx,dt)
% vector of real values telling which
% moments of the displacement are to be
% computed
moments=1:3;

% optimal lagging time estimated from brownian trajectories analysis
n_lag=4;

% delt vector of integer delta t (in frames) values for 
% which the displacements are to be computed.

% determine maximum length in frames among trajectories 
% MaxFrame = max(cellfun('size',trajs,1)); 
% MaxDelt=floor(MaxFrame/3);
% delt=1:MaxDelt;



TrajectoryData=zeros(length(trajs),6);
gamma = zeros(length(moments),1);
Dcoef=zeros(length(moments),1);
MSS = zeros(length(moments),1);
RMS_D=zeros(length(moments),1);
RMS_MSS=zeros(length(moments),1);


for idx=1:length(trajs)
    traj = trajs{idx};
    tlen = size(traj,1);
    delt=1:n_lag;
    MSD = -1*ones(length(moments),length(delt));
    
    %----------------------------------------%
    % calculate MSD for each moment
    %----------------------------------------%
    for imoment=1:length(moments)
        for it=1:length(delt)
            td = delt(it);
            %if (3*td <= tlen)
                %%% calculate MSD for each moment
                dv = (traj([1+td:1:tlen],3)-traj([1:1:tlen-td],3)).^2 + ...
                (traj([1+td:1:tlen],4)-traj([1:1:tlen-td],4)).^2;
                dv = dx.*dx.*dv;   % convert to physical units
                dv = dv.^(imoment/2);
                MSD(imoment,it) = sum(dv)/length(dv);

            %end
        end
        
        %----------------------------------------%
        % aproximate Y=B1*X+B2 to find Diffusion and gamma
        dt_index= find(MSD(imoment,:)>-0.5);
        LogMSD_momement=log(MSD(imoment,dt_index));
        Logt=log(dt*delt(dt_index)); % convert to physical units
        X = [Logt', ones(size(Logt'))];
        Y=LogMSD_momement';
        B=pinv(X)*Y;
        gamma(imoment) = B(1);
        Dcoef(imoment) = 0.25*exp(B(2));
        RMS_D(imoment) =sqrt(sum((Y-(B(1)*Logt'+B(2))).^2))/length(Logt);
    end
    
    %----------------------------------------%
    % calculate MSS
    % aproximate Y=B1*X+B2
    X = [moments', ones(size(moments'))];
    Y=gamma;
    b = pinv(X)*Y;
    MSS(idx) = b(1);
    RMS_MSS(idx) = sqrt(sum((gamma'-(b(1)*moments+b(2))).^2)/length(moments));
    
    %------------------------------------------%
    TrajectoryData(idx,1)=idx;
    TrajectoryData(idx,2)=tlen;
    TrajectoryData(idx,3)=Dcoef(2);
    TrajectoryData(idx,4)=RMS_D(2);
    TrajectoryData(idx,5)=MSS(idx);
    TrajectoryData(idx,6)=RMS_MSS(idx);
    
    %------------------------------------------%
    % plot trajectory, Diffusion constant, MSS
    
   
    % trajectory
    % MSD per trajectory 
%     dv = (traj(2:tlen,2)-traj(1:tlen-1,2)).^2 + ...
%                 (traj(2:tlen,3)-traj(1:tlen-1,3)).^2;
%     dv = dx.*dx.*dv;   % convert to physical units mkm 
%     traj_MSD = sum(dv)/length(dv);
    
%     if(Dcoef(2)>100)
%         figure(idx)
%         subplot(1,3,1)
%         hold on;
%         X=traj(:,2)-mean(traj(:,2))
%         Y=traj(:,3)-mean(traj(:,3))
%         line(X,Y)
%         plot(X,Y,'o','MarkerFaceColor','b','MarkerSize',5)
%         title(sprintf('trajectory %d:',idx))
% 
% 
%         subplot(1,3,2)
%         dt_index=find(MSD(2,:)>-0.5)
%         moment=2;
%         t=dt*delt(dt_index) % convert to physical units
%         m=MSD(2,dt_index)
%         A = [log(t'), ones(size(t'))];
%         a = pinv(A)*log(m');
%         mm = min(log(t));
%         uu = max(log(t));
%         hold on;
%         plot(log(t),log(m),'o','MarkerSize',7,'MarkerFaceColor','r','MarkerEdgeColor','r')
%         xlabel('log(\Deltat)    log([s])');
%         ylabel(sprintf('log(MSD)    log([\\mum]^%d)',moment));
%         plot([mm; uu], [a(1)*mm+a(2); a(1)*uu+a(2)], 'b-','LineWidth',2)
%         title(sprintf('moment %d:D=%.2f, RMS=%f',moment,Dcoef(2),RMS_D(2)))
%         hold off;    
% 
%         subplot(1,3,3)
%         hold on
%         plot(moments,gamma,'o','MarkerSize',7,'MarkerFaceColor','r','MarkerEdgeColor','r')
%         maxmoment = moments(length(moments));
%         mm = min(moments);
%         uu = max(moments);
%         plot([mm; uu], [b(1)*mm+b(2); b(1)*uu+b(2)],'b-','LineWidth',2)
%         plot([0,maxmoment],[0,0.5*maxmoment],'g-','LineWidth',3)
%         xlabel('p')
%         ylabel('\gamma_p')
%         title(sprintf('MSS=%.2f',MSS(idx)))
%         hold off
%     
%         %set(gcf,'position',[x0,y0,width,height])
%         set(gcf,'position',[100,100,1500,500])
%         saveas(gcf,sprintf('Trajectory_%d.png',idx))
%     end  
% %     %pause;
    
 
end
end








