%====================================================================== 
% Write tajectories in matrix form to the file
%
% SYNTAX:  write2file_trajectory(filename,trajectories)
%
% INPUTS:  filename   trajectories in linked list form as:
%
%         trajectories{t}(:,1)    trajectory ID
%         trajectories{t}(:,2)    frame number
%         trajectories{t}(:,3)    x (col)-positions at time t
%         trajectories{t}(:,4)    y (row)-positions at time t
%         trajectories{t}(:,5)    first order intensity moments
%         trajectories{t}(:,6)    second order intensity moments
%
% Tatyana Pichugina, 12.12.2018 
% E-mail: pichugina@evolbio.mpg.de
%
%====================================================================== 
function write2file_trajectory(filename,trajectories)

% open a file for writing
fid = fopen(filename, 'w');
% first line
fprintf(fid, 'ID \t Timeframe \t X \t Y\t M1\t M2\t Area\t\n');
for nline=1:length(trajectories)
    ID=nline;
    tlen=size(trajectories{nline},1);
    for iframe=1:tlen
        fprintf(fid,'%d \t %d \t %f\t %f\t %f \t %f \t %f \t\n',ID,trajectories{nline}(iframe,:));
    end
end
fclose(fid);

return