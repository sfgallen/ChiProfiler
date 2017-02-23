function [chiKsn, ksnUC, chiSlope, UnCert, R2,...
    regBounds,rtext, reg_plot1, reg_plot2, reg_plot3, reg_plotmap] =...
    chiRegress(chi, elev, dist, x, y, Ao, mn, seg, profile_fig, map_fig)

% Function name: chiRegress.m
% Author: Sean F. Gallen
% Date Modified: 06/30/2015
% Purpose: Do a least-squares regression though chi-elevation data to
% determine the channel steepness index, ksn = (U/(K))^(1/n). This function
% also give the user the option to project a modeled river profile based on
% the user defined reference m/n (concavity index) and the derived ksn
% value.
%
% Inputs: chi = chi data (m) derived with the integral approach (Perron and
% Royden, 2012), elev = elevation in meters, dist = distance from mouth in
% meter, Ao = user defined reference drainage area, mn = user defined m/n
% (reference concavity), seg = segment number.
%
% Outputs: chiKsn = ksn derived from regression of chi data, ksnUC = +/- 
% 2 sigma uncertainty, chiSlope = slope of chi regression (M-chi = 
% (U/KAo^m)^(1/n)), UnCert = 2 sigma uncertainty in chi-slope, 
% R2 = the r-squared of the regression,regBounds = the x and y bounds of 
% the regression (this is a vector with 4 values x-min, x-max, y-min, ymax)
% proElm = elevation of the projected profile at the river mouth in meters
% proEl2sig = 2 sigma uncertainty in the elevation of the projected profile
% at the river mouth, mouthEl = elevation of the river mouth in meters,
% rtext = text that is place on plot with some of the regression stats
% reg_plot1 = plot 1 figure handle with regression, reg_plot2= plot 2 
% figure handle with regression, pro1 = plot 1 figure handle with projected
% profile, pro2 = plot 1 figure handle with projected profile + 1 sigma,
% pro3 = plot 1 figure handle with projected profile - 1 sigma.


disp(' ');
disp('Click on minimum THEN maximum bounds for chi from chi-elevation plot (fig1, plot 2; bottom plot)')
disp('Include at least 3 data points on chi-elevation plot')
figure(profile_fig)
subplot(3,1,2)
[chiP,elevP] = ginput(2);
while chiP(1,1) >= chiP(2,1),       %loop to make sure min, max values entered in correct order
    disp(' ');
    disp('Follow the directions!')
    disp('Click on MINIMUM CHI (left) then MAXIMUM CHI (right) bounds from figure 1, plot 2 (right plot)')
    disp('Include at least 3 data points on chi-elevation PLOT')
    [chiP,elevP] = ginput(2);
end
min_chi = chiP(1,1);
max_chi = chiP(2,1);

% Select range of data to fit.
ind = find(chi>=min_chi & chi<=max_chi);
chi_seg = chi(ind);
el_seg = elev(ind);
dist_seg = dist(ind);
regBounds = [min(chi_seg), max(chi_seg), min(el_seg), max(el_seg)];
chi_segMat = [ones(size(el_seg)) chi_seg];
x_seg = x(ind);
y_seg = y(ind);

lastwarn('dude');
% Calculate the least-squares linear fit and 95% confidence intervals.
[b,bint,r,rint,stats] = regress(el_seg,chi_segMat,0.05);
ymod = b(2).*chi_seg + b(1);

while ~strcmp(lastwarn, 'dude')
    disp(' ');
    disp('Regression error! Pick a longer stream segment');
    disp(' ');
    disp('Click on minimum THEN maximum bounds for chi from chi-elevation plot (fig1, plot 2; bottom plot)')
    disp('Include at least 3 data points on chi-elevation plot')
    figure (1)
    subplot(3,1,2)
    [chiP,elevP] = ginput(2);
    while chiP(1,1) >= chiP(2,1),       %loop to make sure min, max values entered in correct order
        disp(' ');
        disp('Follow the directions!')
        disp('Click on MINIMUM CHI (left) then MAXIMUM CHI (right) bounds from figure 1, plot 2 (right plot)')
        disp('Include at least 3 data points -- crosses on chi-elevation PLOT')
        [chiP,elevP] = ginput(2);
    end
    min_chi = chiP(1,1);
    max_chi = chiP(2,1);
    
    % Select range of drainage area to fit.
    ind = find(chi>=min_chi & chi<=max_chi);
    chi_seg = chi(ind);
    el_seg = elev(ind);
    dist_seg = dist(ind);
    regBounds = [min(chi_seg), max(chi_seg), min(el_seg), max(el_seg)];
    chi_segMat = [ones(size(el_seg)) chi_seg];
    x_seg = x(ind);
    y_seg = y(ind);
    
    lastwarn('dude');
    % Calculate the least-squares linear fit and 95% confidence intervals.
    [b,bint,r,rint,stats] = regress(el_seg,chi_segMat,0.05);
    ymod = b(2).*chi_seg + b(1);
    %while strcmp(msgid,'stats:regress:RankDefDesignMat')
end

% find least squares residual
SSresid = sum((el_seg - ymod).^2);
SStotal = sum((el_seg - nanmean(el_seg)).^2);
% calculate R^2
R2 = 1 - SSresid/SStotal;

% get other stats
chiSlope = b(2);
chiKsn = chiSlope*Ao^mn;
UnCert= (bint(4)-bint(2))/2;
ksnUC = UnCert*Ao^mn;
yInt = b(1);
%yUnCert= (bint(3)-bint(1))/2;

% plot regression on both plots.
figure(profile_fig);
subplot(3,1,1)
reg_plot1 = plot(dist_seg./1000,ymod,'r--',dist_seg(1)/1000,ymod(1),'rx',...
    dist_seg(end)/1000,ymod(end),'rx','Linewidth',1.5);
subplot(3,1,2)
reg_plot2 = plot(chi_seg,ymod,'r--',chi_seg(1),ymod(1),'rx',...
    chi_seg(end),ymod(end),'rx','Linewidth',1.5);
subplot(3,1,3)
ck = chiKsn.*ones(size(chi_seg));
reg_plot3 = plot(chi_seg,ck,'r--',chi_seg(1),chiKsn,'rx',...
    chi_seg(end),chiKsn,'rx','Linewidth',1.5);

figure(map_fig);
reg_plotmap = plot(x_seg,y_seg, 'r--', x_seg(1), y_seg(1),'rx',...
    x_seg(end), y_seg(end), 'rx','Linewidth',1.5);
figure(profile_fig);
%place theta, ksn on this plot, above the fit section:
minchipt=min(chi_seg);
maxchipt=max(chi_seg);
xloc=(minchipt+maxchipt)/2;
disp(' ');
disp('Click on plot TO SAVE DATA on the screen near your regression.')
rtext=gtext({['reach ',num2str(seg+1),' fit:'],['chi k_{sn}=',num2str(chiKsn,4),' \pm ',num2str(ksnUC,4)],['R^2: ',num2str(R2,3)],},'fontsize',8);

end
