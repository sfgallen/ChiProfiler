function [smoothEl] = smoothChannelZ(elev,win,cs)
% smoothChannelZ.m smooths channel data with filtfilt, or if too few data
% points with fastsmooth.m
%
% inputs: 
%       elev = elevation in meters
%       win = window size in meters,
%       cs = cellsize of grid in meters
%
% ouputs:
%       smoothedEl: smoothed elevation vector
%
% Author: Sean F. Gallen (modified from Scott R. Miller's SmoothChanZ.m)
% Date Modifed: 12/31/2015

win2 = 2.*round(((win/cs)+1)/2)-1; % round to the nearest odd integer
win2 = min(win2,length(elev));

if win2*3 >= length(elev)
    smoothEl = fastsmooth(elev,win2,3,1);
else
    smoothEl = filtfilt(ones(1,win2)/win2,1,elev);
end

smoothEl(1) = elev(1);
smoothEl(end) = elev(end);

if length(elev) >= 3;
    if smoothEl(end -1) < smoothEl(end)
        smoothEl(end-1) = (smoothEl(end) + smoothEl(end-2))/2;
    end
end

end