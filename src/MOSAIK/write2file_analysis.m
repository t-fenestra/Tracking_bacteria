%====================================================================== 
% Write tajectories in matrix form to the file
%
% SYNTAX:  write2file_analysis(filename,trajectories)
%
% INPUTS:  filename   trajectories in linked list form as:
%
%         alysis_matrix(:,1)    trajectory ID
%         alysis_matrix(:,2)    frame length
%         alysis_matrix(:,3)    D mkm2/s
%         alysis_matrix(:,4)     RMS_D
%         alysis_matrix(:,5)     MSS
%         talysis_matrix(:,6)    MSS_D
%
% Tatyana Pichugina, 12.12.2018 
% E-mail: pichugina@evolbio.mpg.de
%
%====================================================================== 
function write2file_analysis(filename,alysis_matrix)
fid = fopen(filename, 'w');

% first line
fprintf(fid, 'ID \t Length \t D mkm^2/s \t RMS_D \t MSS \t MSS_D\t\n');

for nline=1:size(alysis_matrix,1)
    fprintf(fid,'%d \t %d \t %f \t %f\t %f \t %f \t \n',alysis_matrix(nline,:));
end
fclose(fid);

return;