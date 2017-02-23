function [chiFits, kp_data, regV] = profile_chi(data, strmID, profile_fig, map_fig, Ao, mn, step, chanDir, minchi, maxchi)
% Function name: chi_profiler.m
% Author: Sean F. Gallen
% Date modified: 06/30/2015
% Purpose: Make regressions through chi-elevation data to determine ksn
% with the integral approach (e.g. Perron and Royden, 2013) using the
% function chiRegress.m and identify knickpoints
%
% Inputs: 
%       1) data: Stream channel data table produce with 
%       chi_profiler_STREAMobj.m
%       2) strmID: stream ID number
%       3) profile_fig: figure with profiles created by
%       chi_profiler_STREAMobj.m
%       4) map_fig: figure with profiles created by
%       chi_profiler_STREAMobj.m
%       5) Ao: reference drainage area for chi analysis, 
%       6) mn: the m/n ratio used for chi analysis.
%       7) step: size of binning increment for binned ksn plot
%       8) chandir: the directory where the channel data is stored
%       9) minchi: minimum chi value to scale plot axes
%       10) maxchi: maximum chi value to scale plot axes
%
% Outputs:
%       1) chiFits: table of the data related to chi regressions
%       2) kp_data: table of knickpoint data
%       3) regV: vector of ksn values for regressed channel segments
%
% Author: Sean F. Gallen
% date modified: 12/31/2015

dist = data(:,4);      % distance from mouth
zr = data(:,2);      % elevation
A = data(:,3);      % drainage area
%% other data in chandata
dfd = data(:,1);        % distance from divide
sel = data(:,5);        % smoothed elevation
xmat = data(:,6);       % x (row) index in matrix
ymat = data(:,7);       % y (column) index in matrix
x_coord = data(:,9);    % UTM latitude in meters
y_coord = data(:,10);   % UTM longitude in meters;
chi = data(:,11);

z = sel;
% regression vector
regV = nan(size(z));


% get bin data to get chi-ksn pattern
incs = floor(length(chi)/step);
spt = 1; ept = step;
mp = nan(length(incs),1);
binKsn = nan(length(incs),1);
for i = 1:incs
    mp(i) = nanmean(chi(spt:ept));
    chi_segMat = [ones(size(chi(spt:ept))) chi(spt:ept)];
    [b,bint,r,rint,stats] = regress(z(spt:ept),chi_segMat,0.05);
    binKsn(i) = b(2)*Ao^mn;
    spt = spt + step;
    ept = ept + step;
end
figure(profile_fig);
s3 = subplot(3,1,3); hold on
ksnBin_opt = 'n';
if length(binKsn(~isnan(binKsn))) >= 2
    ksnBin_opt = 'y';
    sK = plot(mp,binKsn,'bo');
    xlabel('\chi (m)'); ylabel(['k_{sn} (m^{',num2str(2*mn),'})']);
    axis([minchi maxchi+(maxchi-minchi)*.1 0 nanmax(binKsn)])
end

%set(s3, 'YScale','log');

txt = sprintf(['Do you want to make regressions through river channel',...
    ' segments on the chi-elevation plot?\n']);
disp(txt);
regress_opt = input('type "y" for yes or "n" for no:  ','s');

while ~strcmp(regress_opt,'y') && ~strcmp(regress_opt,'n'),
    %case where you didn't enter a, b, c, or d:
    disp(' ');
    regress_opt = input('YES (y) or NO (n) only!!!:  ','s');
end

reg_plots = 'n';
kp_plots = 'n';

if regress_opt == 'y'
    chiFits = [];
    rtextLrg = [];
    rp1Lrg = [];
    rp2Lrg = [];
    rp3Lrg = [];
    rpMLrg = [];
    seg = 0;
    bo = 1;
    
    while bo == 1;
        % Call chiRegess.m to make regressions through specified reach
        [chiKsn, ksnUC, chiSlope, UnCert, R2,...
            regBounds, rtext, reg_plot1, reg_plot2,...
            reg_plot3, reg_plotmap] =...
            chiRegress(chi, z, dist, x_coord, y_coord, Ao, mn, seg, profile_fig, map_fig);

        disp(' ');
        disp('Do you want to remember this fit?')
        fit_opt = input('type "y" for yes or "n" for no:  ','s');
        
        if fit_opt == 'y'
            % add data to chiFits table
            newdata = [strmID, seg+1, chiKsn, ksnUC, chiSlope, UnCert, R2,...
                regBounds, x_coord(end), y_coord(end)];
            chiFits = [chiFits; newdata];
            rtextLrg = [rtextLrg; rtext];
            rp1Lrg = [rp1Lrg; reg_plot1];
            rp2Lrg = [rp2Lrg; reg_plot2];
            rp3Lrg = [rp3Lrg; reg_plot3];
            rpMLrg = [rpMLrg; reg_plotmap];
            reg_plots = 'y';
        else
            delete([rtext; reg_plot1; reg_plot2; reg_plot3; reg_plotmap]);
            seg = seg-1;
        end
        
        disp(' ');
        disp('Do you want to fit another channel segment?')
        fit_opt2 = input('type "y" for yes or "n" for no:  ','s');
        
        while ~strcmp(fit_opt2,'y') && ~strcmp(fit_opt2,'n'),
            %case where you didn't enter a, b, c, or d:
            disp(' ');
            fit_opt2 = input('YES (y) or NO (n) only!!!:  ','s');
        end
        
        if fit_opt2 == 'y'
            bo = 1;
            seg = seg+1;
        elseif fit_opt2 == 'n'
            bo = 0;
        end
    end
    
    % put chi-ksn regression value in to regression vector.
    [nr,nc] = size(chiFits);
    if nr ~= 0
        for q = 1:length(nr);
            regV(z >= chiFits(nr,10) & z <= chiFits(nr,11))...
                = chiFits(nr,3);
        end
    else
        noDatVect = -9999.*ones(1,12);
        newdata = [strmID, noDatVect];
        chiFits = [newdata];
    end
elseif regress_opt == 'n'
    noDatVect = -9999.*ones(1,12);
    newdata = [strmID, noDatVect];
    chiFits = [newdata];
    %regV = [];
end

disp(' ');
disp('Do you want to MARK KNICKPOINTS on the chi-elevation (middle) plot?')
kp_opt = input('type "y" for yes or "n" for no:  ','s');

while ~strcmp(kp_opt,'y') && ~strcmp(kp_opt,'n'),
    %case where you didn't enter a, b, c, or d:
    disp(' ');
    kp_opt = input('Yes(y) or No(n)!!!!:  ','s');
end

if kp_opt == 'y'
    
    kp_data = [];
    kp = 0;
    bo = 1;
    
    while bo ==1;
        disp(' ');
        disp('SELECT a POINT on the chi-elevation (middle) plot.') ;
        
        figure(profile_fig)
        subplot(3,1,2);
        [chi_kp,elev_kp] = ginput(1);
        kp_ind = find(chi <= chi_kp,1,'first');
        while isempty(kp_ind);
            disp('Wrong plot or point not within bounds.');
            disp('SELECT a POINT on the chi-elevation (middle) plot.') ;
            [chi_kp,elev_kp] = ginput(1);
            kp_ind = find(chi <= chi_kp,1,'first');
        end
        % might need to add some error handling here
        kp_p2(kp+1) = plot(chi(kp_ind),z(kp_ind),'kv','MarkerFaceColor',[150/255 150/255 255/255]);
        
        subplot(3,1,1);
        kp_p1(kp+1) = plot(dist(kp_ind)/1000,z(kp_ind),'kv','MarkerFaceColor',[150/255 150/255 255/255]);
        
        figure(map_fig);
        kp_p3(kp+1) = plot(x_coord(kp_ind),y_coord(kp_ind),'ko','MarkerFaceColor',[0.8 0.8 0.8]);
        
        figure(profile_fig);
        
        disp(' ');
        disp('Do you want to REMEMBER this POINT?')
        fit_opt = input('type "y" for yes or "n" for no:  ','s');
        
        while ~strcmp(fit_opt,'y') && ~strcmp(fit_opt,'n'),
            %case where you didn't enter y or n:
            disp(' ');
            fit_opt = input('Yes(y) or No(n)!!!!:  ','s');
        end
        
        if fit_opt == 'y'
            disp(' ');
            disp('CLASSIFY the KNICKPOINT type with a number, for example,');
            kp_class = input('1 = major knickpoint, 2 = minor knickpoint, ect...');
            
            % add data to table
            newdata = [strmID, kp+1, kp_class, chi(kp_ind), zr(kp_ind), A(kp_ind),...
                dist(kp_ind), dfd(kp_ind), sel(kp_ind), x_coord(kp_ind),...
                y_coord(kp_ind), xmat(kp_ind), ymat(kp_ind),...
                x_coord(end), y_coord(end)];
            kp_data = [kp_data; newdata];
            kp = kp+1;
            kp_plots = 'y';
        elseif fit_opt == 'n'
            delete([kp_p1(kp+1), kp_p2(kp+1), kp_p3(kp+1)]);
        end
        
        disp(' ');
        disp('Do you want to select ANOTHER POINT?')
        fit_opt2 = input('type "y" for yes or "n" for no:  ','s');
        
        while ~strcmp(fit_opt2,'y') && ~strcmp(fit_opt2,'n'),
            %case where you didn't enter a, b, c, or d:
            disp(' ');
            fit_opt2 = input('Yes(y) or No(n)!!!!:  ','s');
        end
        
        if fit_opt2 == 'y'
            bo = 1;
        elseif fit_opt2 == 'n'
            bo = 0;
        end
    end
    %allkpDat = [allkpDat; kp_data];
elseif kp_opt == 'n'
    kp_data = [];
end

disp(' ');
disp('Do you want to SAVE this PLOT?');
save_opt = input('type "y" for yes or "n" for no:  ','s');

while ~strcmp(save_opt,'y') && ~strcmp(save_opt,'n'),
    %case where you didn't enter a, b, c, or d:
    disp(' ');
    save_opt = input('Yes(y) or No(n)!!!!:  ','s');
end

if save_opt == 'y'
    figure(profile_fig)
    plotname = [chanDir,'\', num2str(strmID), '_chi-plot'];
    eval ([' print -depsc ',plotname,'.eps'])
elseif save_opt == 'n'
end

if ksnBin_opt == 'y'
    delete(sK);
end

if reg_plots == 'y'
    delete([rtextLrg; rp1Lrg; rp2Lrg; rp3Lrg; rpMLrg]);
end

if kp_plots == 'y'
    delete([kp_p1,kp_p2,s3,kp_p3])
elseif kp_plots == 'n'
    delete(s3)
end

end
