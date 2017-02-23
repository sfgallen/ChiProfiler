function [ksnSTREAMs] = binnedKsn(S,chi,Sz,win,cs,Ao,mn)
% binnedKsn.m calcuates ksn from chi elevation data using a moving window
%
% Inputs:
%       1) S: a TopoToolBox STREAMobj
%       2) chi: topologically ordered vector of chi values
%       3) Sz: topologically ordered vector of elevation values
%       4) win: window size in map units
%       5) cs: cellsize in map units
%       6) Ao: reference drainage area used for chi integration
%       7) mn: reference m/n (theta) value
%
% Outputs:
%       1) topologically ordered vector of ksn values
%
% Author: Sean F. Gallen
% Date Modified: 12/31/2015

step = (win/2);
indStep = ceil(step/cs);
% declare necessary variables from stream object
ordList = S.orderednanlist;         % ordered list of streams split by nans
strmBreaks = find(isnan(ordList));  % get position of nans
Six = S.ix;                         % doners
Sixc = S.ixc;                       % recievers
Sd = S.distance;                    % distance from mouth

% cast dumby vector to catch data
ksnSTREAMs = nan(size(Sd));

% from the donor position identifiy receiver nodes in stream object order
receiverID = nan(size(Sd));
donorID = nan(size(Sd));
for i = numel(Six):-1:1;
    receiverID(Six(i)) = Sixc(i);
    donorID(Six(i)) = Sixc(i);
end


id1 = 0;
h = waitbar(0,'Calculating ksn from chi for ksn map...');
for i = 1:length(strmBreaks);
    strmInds = ordList(id1+1:strmBreaks(i)-1);
    tChi = chi(strmInds);
    tZ = Sz(strmInds);
    tD = Sd(strmInds);
    trecID = receiverID(strmInds);
    tKsn = nan(length(strmInds),1);
    % if the stream is a trunk channel run smoothing window
    if isnan(trecID(end))
        for j = 1:length(tChi)
            % find within the window size.
            chiWin = tChi(tD > tD(j)-win & tD <= tD(j)+win);
            zWin = tZ(tD > tD(j)-win & tD <= tD(j)+win);
            
            % make sure there are enough points for the regression
            if length(chiWin) > 2
                chi_segMat = [ones(size(zWin)) chiWin];
                [b,bint,r,rint,stats] = regress(zWin,chi_segMat,0.05);
                tKsn(j) = b(2)*Ao^mn;
            else
                tKsn(j) = nan;
            end
        end
    % if the stream is a tributary allow smoothing window to continue down
    % stream by 1/2 of moving window width
    else
        addChi = nan(1,indStep);
        addZ = nan(1,indStep);
        addD = nan(1,indStep);
        nInd = trecID(end);
        n = 1;
        while n < length(addChi);
            addChi(n) = chi(nInd);
            addZ(n) = Sz(nInd);
            addD(n) = Sd(nInd);
            nInd = receiverID(nInd);
            if isnan(nInd) == 1 || isempty(addChi(isnan(addChi))) == 1
                n = length(addChi) + 1;
            else
                n = n+1;
            end
        end
        addChi = addChi(~isnan(addChi));
        addZ = addZ(~isnan(addZ));
        addD = addD(~isnan(addD));
        tChi2 = nan(length(tChi)+length(addChi),1);
        tZ2 = nan(length(tZ)+length(addZ),1);
        tD2 = nan(length(tD)+length(addD),1);
        
        tChi2(1:length(tChi)) = tChi;
        tChi2(length(tChi)+1:end) = addChi;
        tZ2(1:length(tZ)) = tZ;
        tZ2(length(tZ)+1:end) = addZ;
        tD2(1:length(tD)) = tD;
        tD2(length(tD)+1:end) = addD;
        
        for j = 1:length(tChi)
            % find within the window size.
            chiWin = tChi2(tD2 > tD2(j)-win & tD2 <= tD2(j)+win);
            zWin = tZ2(tD2 > tD2(j)-win & tD2 <= tD2(j)+win);
            
            % make sure there are enough points for the regression
            if length(chiWin) > 2
                chi_segMat = [ones(size(zWin)) chiWin];
                [b,bint,r,rint,stats] = regress(zWin,chi_segMat,0.05);
                tKsn(j) = b(2)*Ao^mn;
            else
                tKsn(j) = nan;
            end
        end
    end
    tKsn(tKsn < 0) = 0;
    ksnSTREAMs(strmInds) = tKsn;
    id1 = strmBreaks(i);
    f = i/length(strmBreaks);
    waitbar(f,h);
end
close(h)
end
        
        
        
        
        
