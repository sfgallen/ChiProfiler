function writeTXTtables(chanDir, fileTag, kpMat, regMat)
% writeTXTtables takes four inputs:
%       1) chanDir: string of the directory where channel data is stored
%       2) fileTag = a short string to idenitfy your text files
%       3) kpMat = the knickpoint data generated from chi_profiler.m that
%       is run during chi_profiler_STREAMobj.m
%       4) regMat = the chi regression data generated from chi_profiler.m
%       that is run during chi_profiler_STREAMobj.m
%
% Author: Sean F. Gallen
% Data modified: 07/09/2015

[nr,nc] = size(kpMat);
if nr ~= 0
    % organize knickpoint data for ArcGIS
    % make a table that can be added to ArcGIS
    % get number of rows and columns
    [nrowsKP,ncolsKP] = size(kpMat);
    [nrowsR,ncolsR] = size(regMat);
    
    kpXYstrm = nan(nrowsKP,5);
    kpXYstrm(:,1) = kpMat(:,10);
    kpXYstrm(:,2) = kpMat(:,11);
    kpXYstrm(:,3) = kpMat(:,1);
    kpXYstrm(:,4) = kpMat(:,2);
    kpXYstrm(:,5) = kpMat(:,3);
    
    % Writing textfiles of the data
    file1 = [chanDir,'\',fileTag, '_kp_data_comp.txt'];
    file2 = [chanDir,'\',fileTag, '_kp_data_forArcGIS.txt'];
    % write headers to each column for the three text files to be produced
    
    fileID1 = fopen(file1,'w');
    fprintf(fileID1, '%10s\t%6s\t%8s\t%3s\t%4s\t%2s\t%3s\t%3s\t%5s\t%7s\t%7s\t%6s\t%6s\n',...
        'stream_num', 'kp_num', 'kp_class', 'chi', 'elev', 'DA', 'dfm',...
        'dfd', 'sm_el', 'x_coord', 'y_coord', 'xMatID', 'yMatID');
    
    fileID2 = fopen(file2,'w');
    fprintf(fileID2, '%7s\t%7s\t%10s\t%6s\t%6s\n',...
        'x_coord', 'y_coord','stream_num', 'kp_num', 'kp_class');
    
    % write the two knickpoint tables as tab delimited
    txt =['Writing text files "', fileTag, '_kp_data_comp.txt" and "', fileTag, '_kp_data_forArcGIS.txt".....'];
    h = waitbar(0,txt);
    
    for i = 1:nrowsKP
        fprintf(fileID1, '%g\t',kpMat(i,1:ncolsKP-1));
        fprintf(fileID1, '%g\n',kpMat(i,ncolsKP));
        
        fprintf(fileID2, '%g\t',kpXYstrm(i,1:4));
        fprintf(fileID2, '%g\n',kpXYstrm(i,5));
        
        f = i/nrowsKP;
        waitbar(f)
    end
    close(h)
    fclose(fileID1);
    fclose(fileID2);
end

[nr,nc] = size(regMat);
if nr ~= 0
    % Writing textfiles of the data
    file3 = [chanDir,'\',fileTag, '_chiFits_data_comp.txt'];
    
    % write headers to each column for the three text files to be produced
    fileID3 = fopen(file3,'w');
    fprintf(fileID3, '%10s\t%5s\t%7s\t%6s\t%9s\t%6s\t%2s\t%10s\t%10s\t%9s\t%9s\t%9s\t%6s\t%10s\n',...
        'stream_num', 'reach', 'chi_ksn', '2sigma', 'chi_slope', '2sigma',...
        'R2', 'chiReg_min', 'chiReg_max', 'elReg_min', 'elReg_max',...
        'proj_elev', '2sigma', 'mouth_elev');

    % write the on regression table as tab delimited
    txt =['Writing text file "', fileTag, '_chiFits_data_comp.txt".....'];
    h = waitbar(0,txt);
    
    for i = 1:nrowsR
        fprintf(fileID3, '%g\t',regMat(i,1:ncolsR-1));
        fprintf(fileID3, '%g\n',regMat(i,ncolsR));
        
        f = i/nrowsKP;
        waitbar(f)
    end
    close(h)
end
fclose(fileID3);