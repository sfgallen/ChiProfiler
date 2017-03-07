function chi_profiler(DEM, S, fileTag, varargin)
% chi_profiler.m allows for interative river profile analysis via
% the integral or chi method (e.g. Perron and Royden, 2013). Users can
% regress through channel segments to get ksn and select knickpoints along
% river profiles
%
% Inputs:
%       1) DEM: TopoToolBox DEM GRIDobj (required)
%       2) S: TopoToolBox STREAMobj (required)
%       3) fileTag: 'project' name as a string used for output files and 
%          folder names. (required)
%       4) crita: Critical drainage area for channel head initiation in map
%          units. (optional) {default --> 1e6}
%       5) mn: reference m/n (theta) value (optional) {default --> 0.45}
%       6) Ao: reference drainage area for chi integration (optional)
%          {default --> 1}. Note: it is recommended to always use an 'Ao'
%          of 1 as this results in chi versus elevation plots with a slope
%          that is equal to the normalized steepness index.
%       7) smoWin: size of window used to smooth elevation data (set this
%          to the cell size of the DEM if you don't want the data to be
%          smoothed) (optional) {default --> 250}
%       8) flowOption: string of 'fill' or 'carve' for flow routing.
%          (optional) {default --> empty assumes that DEM is already 
%          filled or carved}
%       
% Outputs:
%       No variable outputs to Matlab, but chi_profiler.m will produce a
%       series of files that can be imported into Matlab and ArcGIS.
%
%       - Data tables of each of the stream channels analyzed
%       - Shapefiles of chi, ksn, regressed ksn segments and knickpoints
%       - Tables of regressed ksn and knickpoint statistics
%       Optional:
%       - river profiler figures
%       - chi and ksn map of entire river network in DEM as a shapefile
%
%       See 'User_Guide.docx' for more information on the how to run
%       chi_profiler.m and the output file organization
%
% Example:
%  chi_profiler(DEM,S,'my_proj','crita',1e6,'mn',0.45,'Ao',1,'smoWin',250);
%
%
% Author: Sean F. Gallen
% Date Modified: 02/20/2017
% email: sean.gallen[at]erdw.ethz.ch

% Parse Inputs
p = inputParser;         
p.FunctionName = 'chi_profiler';

% required inputs
addRequired(p,'DEM', @(x) isa(x,'GRIDobj'));
addRequired(p,'S',@(x) isa(x,'STREAMobj'));
addRequired(p,'fileTag', @(x) ischar(x))

% optional inputs
addOptional(p,'crita', 1e6, @(x) isscalar(x));
addOptional(p,'mn', 0.45, @(x) isscalar(x));
addOptional(p,'Ao', 1, @(x) isscalar(x));
addOptional(p,'smoWin', 250, @(x) isscalar(x));
addOptional(p,'flowOption', []);

parse(p,DEM, S, fileTag, varargin{:});
DEM   = p.Results.DEM;
S     = p.Results.S;
fileTag     = p.Results.fileTag;
crita    = p.Results.crita;
mn = p.Results.mn;
Ao = p.Results.Ao;
smoWin = p.Results.smoWin;

% make folder to get all the data
folder = [fileTag, '_stream_data/'];
if ~exist(folder, 'dir')
  mkdir(cd,[fileTag, '_stream_data']);
end

%% create varables with topotoolbox functions
% set nan values if it hasn't already been done
DEM.Z(DEM.Z <= -9999) = NaN;
% declare cellsize
cs = DEM.cellsize;

% flow routing options
if isempty(p.Results.flowOption)
    FD = FLOWobj(DEM);
elseif strcmp(p.Results.flowOption, 'fill');
    DEM = fillsinks(DEM);
    FD = FLOWobj(DEM);
elseif strcmp(p.Results.flowOption, 'carve');
    FD = FLOWobj(DEM,'preprocess','carve');
    DEM = imposemin(FD,DEM);
else
    error('fillOption is not "fill" or "carve"');
end

% Calculate flow accumulation (A)
A   = flowacc(FD).*(cs^2);

% Calculate distance from the channel head
DFD = flowdistance(FD,'downstream'); % this is actually distance from channel head

% Save STEAMobj
chanDir = [cd,'/',fileTag, '_stream_data'];
%mkdir(cd,[fileTag, '_stream_data']);

fileName = [fileTag, '_pickedstreams.mat'];
save([chanDir,'/',fileName],'S');

% Declare STREAMobj variables for faster processing through forloop
ordList = S.orderednanlist;
strmBreaks = find(isnan(ordList));

ksnReg = DEM;
ksnReg.Z = nan(size(DEM.Z));

strmNumGrid = DEM;
strmNumGrid.Z = nan(size(DEM.Z));

GridID = S.IXgrid;

disp(' ');
disp(['you will be analyzing ' num2str(length(strmBreaks)) ' stream channels'])

%% calcuate chi
% declare variables for chi integration
Schi = zeros(size(S.distance));
Six = S.ix;                         % donors
Sixc = S.ixc;                       % recievers
Sd = S.distance;                    % distance from mouth
Sa = (Ao./(A.Z(GridID))).^mn;       % chi transformation variable

h = waitbar(0,'calculating \chi for user picked streams...');
% calculating chi for the entire river network
for lp = numel(Six):-1:1;
    Schi(Six(lp)) = Schi(Sixc(lp)) + (Sa(Sixc(lp))+(Sa(Six(lp))-Sa(Sixc(lp)))/2) *(abs(Sd(Sixc(lp))-Sd(Six(lp))));
    f = (numel(Six)+1 - lp)/numel(Six);
    waitbar(f,h);
end
close(h);

%% declare all other stream network variables
Sz = double(DEM.Z(GridID));                 % elevation
SmoZ = Sz;               % dumby vector to get smoothed data
Sx = S.x;                           % x_coordinate
Sy = S.y;                           % y_coordinate
Sdfd = double(DFD.Z(GridID));               % distance from 'divide' (acually channel head)
Sda = double(A.Z(GridID));                  % drainage area


%% plot stream network data on a map and as profiles

% plot a map with DEM and streams
map_fig = figure(2);
imageschs(DEM); hold on
plot(S,'-','LineWidth', 2,'Color', [0 0 0]);
set(map_fig,'units','centimeters', 'Position', [25 9 15 11])

% get axes limits for profile data
mindfm = nanmin(Sd)./1000; maxdfm = nanmax(Sd)./1000;
minel = nanmin(Sz); maxel = nanmax(Sz);
minchi = nanmin(Schi); maxchi = nanmax(Schi);

% plot all of the river river profile data as thin gray lines
h = waitbar(0,'Smoothing data and plotting all streams...');
profile_fig = figure(1);
id1 = 0;
for i = 1:length(strmBreaks);
    strmInds = ordList(id1+1:strmBreaks(i)-1);
    SmoZ(strmInds) = smoothChannelZ(SmoZ(strmInds),smoWin,cs);
    subplot(3,1,1);
    plot(Sd(strmInds)./1000,SmoZ(strmInds),'k-','lineWidth',0.5,'color',[0.5, 0.5, 0.5]); hold on
    xlabel('distance (km)'); ylabel('elevation (m)');
    subplot(3,1,2);
    plot(Schi(strmInds),SmoZ(strmInds),'k-','lineWidth',0.5,'color',[0.5, 0.5, 0.5]); hold on
    xlabel('\chi (m)'); ylabel('elevation (m)');
    id1 = strmBreaks(i);
    f = i/length(strmBreaks);
    waitbar(f,h);
end
close(h)
subplot(3,1,1);
axis([mindfm maxdfm+(maxdfm-mindfm)*.1 minel maxel+(maxel-minel)*0.1])
subplot(3,1,2);
axis([minchi maxchi+(maxchi-minchi)*.1 minel maxel+(maxel-minel)*0.1])
set(profile_fig,'units','centimeters', 'Position', [1 2 15 18])

%% run through all the river profiles individually and have the user
%% decide if they want to get ksn for river channel segments and pick kps

% define step size for ksn bins, this is based on true distance
step = round(smoWin./cs);

% create empty matrices to catch chi regression and kp data
chiFits = [];
kp_data = [];

id1 = 0;
strmNum = 1;
for i = 1:length(strmBreaks);
    strmInds = ordList(id1+1:strmBreaks(i)-1);
    dataMat = nan(length(strmInds),11);
    % save stream channel data in a matrix
    dataMat(:,1) = Sdfd(strmInds);
    dataMat(:,2) = Sz(strmInds);
    dataMat(:,3) = Sda(strmInds);
    dataMat(:,4) = Sd(strmInds);
    dataMat(:,5) = smoothChannelZ(dataMat(:,2),smoWin,cs);
    [mrows, mcols] = ind2sub(DEM.size,GridID(strmInds));    % row and column of data in matrix
    dataMat(:,6) = mrows;
    dataMat(:,7) = mcols;
    dataMat(:,8) = GridID(strmInds);
    dataMat(:,9) = Sx(strmInds);
    dataMat(:,10) = Sy(strmInds);
    dataMat(:,11) = Schi(strmInds);
    
    % save data for this stream
    dataFileName = [num2str(strmNum),'_', fileTag, '_chandata.mat'];
    save([chanDir,'/',dataFileName],'dataMat');
    
    % highlight river profile or interest
    figure(profile_fig)
    subplot(3,1,1);
    s1 = plot((dataMat(:,4))./1000,dataMat(:,5),'-','lineWidth',2,'color',[0/255 205/255 205/255]);
    xlabel('distance (km)'); ylabel('elevation (m)');
    title(['stream ', num2str(strmNum), ' of ', num2str(length(strmBreaks))]);
    subplot(3,1,2);
    s2 =  plot(dataMat(:,11),dataMat(:,5),'-','lineWidth',2,'color',[0/255 205/255 205/255]);
    xlabel('\chi (m)'); ylabel('elevation (m)');
    
    figure(map_fig)
    s3 = plot(dataMat(:,9),dataMat(:,10),'-','lineWidth',2,'color',[250/255 250/255 250/255]);
    title(['stream ', num2str(strmNum), ' of ', num2str(length(strmBreaks))]);
    figure(profile_fig)
    
    % run the chi profiler function
    [newCF, newKP, regV] = profile_chi(dataMat, strmNum, profile_fig, map_fig, Ao, mn, step, chanDir, minchi, maxchi);
    
    chiFits = [chiFits; newCF];
    kp_data = [kp_data; newKP];
    delete([s1,s2,s3]);
    
    ksnReg.Z(GridID(strmInds)) = regV;
    strmNumGrid.Z(GridID(strmInds)) = strmNum;
    
    chiFileName = [fileTag, '_chiFits.mat'];
    kpFileName = [fileTag, '_kpData.mat'];
    save([chanDir,'/',chiFileName],'chiFits');
    save([chanDir,'/',kpFileName],'kp_data');
    
    strmNum = strmNum + 1;
    id1 = strmBreaks(i);
end
    
if ~isempty(kp_data);
    % make knickpoint shapefile
    MP = struct('Geometry',{'Point'},...
        'X',num2cell(kp_data(:,10)),...
        'Y',num2cell(kp_data(:,11)),...
        'strm_num',num2cell(kp_data(:,1)),...
        'kp_num',num2cell(kp_data(:,2)),...
        'kp_type',num2cell(kp_data(:,3)),...
        'chi',num2cell(kp_data(:,4)),...
        'elev',num2cell(kp_data(:,5)),...
        'smo_el',num2cell(kp_data(:,9)),...
        'd_area',num2cell(kp_data(:,6)),...
        'dfm',num2cell(kp_data(:,7)),...
        'dfd',num2cell(kp_data(:,8)),...
        'GridX',num2cell(kp_data(:,12)),...
        'Gridy',num2cell(kp_data(:,13)),...
        'outletX',num2cell(kp_data(:,14)),...
        'outletY',num2cell(kp_data(:,15)));
    
    shapewrite(MP,[chanDir,'/',fileTag, '_kp_Data.shp']);
    
    % make knickpoint data excel table.
    X =kp_data(:,10);
    Y = kp_data(:,11);
    strm_num = kp_data(:,1);
    kp_num = kp_data(:,2);
    kp_type = kp_data(:,3);
    chi = kp_data(:,4);
    elev = kp_data(:,5);
    smo_el = kp_data(:,9);
    d_area = kp_data(:,6);
    dfm = kp_data(:,7);
    dfd = kp_data(:,8);
    GridX = kp_data(:,12);
    GridY = kp_data(:,13);
    outletX = kp_data(:,14);
    outletY = kp_data(:,15);
    
    T = table(X, Y, strm_num, kp_num, kp_type, chi, elev, smo_el,...
        d_area, dfm, dfd, GridX, GridY, outletX, outletY);
    
    filename = [chanDir,'/',fileTag, '_kp_Data.xlsx'];
    writetable(T,filename)   
end

ksnRegS = ksnReg.Z(GridID);
if ~isempty(ksnRegS(~isnan(ksnRegS)));
    
    % make shapefile of ksn regressions
    MP = struct('Geometry',{'Point'},...
        'X',num2cell(Sx(~isnan(ksnRegS))),...
        'Y',num2cell(Sy(~isnan(ksnRegS))),...
        'ksn_reg',num2cell(ksnRegS(~isnan(ksnRegS))));
    shapewrite(MP,[chanDir,'/',fileTag, '_ksn_regressions.shp']);
    
    chiFits = chiFits((chiFits(:,2) ~= -9999),:);
    save([chanDir,'/',chiFileName],'chiFits');
    
    % write excel table with regression data
    stream_ID = chiFits(:,1);
    segment_num = chiFits(:,2);
    ksn = chiFits(:,3);
    ksn_95uc = chiFits(:,4);
    r_squared = chiFits(:,7);
    min_chi = chiFits(:,8);
    max_chi = chiFits(:,9);
    min_elev = chiFits(:,10);
    max_elev = chiFits(:,11);
    outletX = chiFits(:,12);
    outletY = chiFits(:,13);
    
    T = table(stream_ID, segment_num, ksn, ksn_95uc, r_squared, min_chi,...
        max_chi, min_elev, max_elev, outletX, outletY);
    
    filename = [chanDir,'/',fileTag, '_ksn_regressions.xlsx'];
    writetable(T,filename)  
    
end

% if ~isempty(kp_data) && ~isempty(ksnRegS(~isnan(ksnRegS)));
%     % save data with fileTag
%     % make tab delimited table for knickpoint and chi data
%     writeTXTtables(chanDir, fileTag, kp_data, chiFits);
% end

ChiGrid = DEM;
ChiGrid.Z = nan(size(DEM.Z));
ChiGrid.Z(GridID) = Schi;

% making ksn map for streams analyzed
ksnStreams = binnedKsn(S,Schi,SmoZ,smoWin,cs,Ao,mn);
ksnG = DEM;
ksnG.Z = nan(size(DEM.Z));
ksnG.Z(GridID) = ksnStreams;

% Make a shapefile with the analyzed stream segments numbered (0 mean
% it was not analyzed
MS = STREAMobj2mapstruct(S,'seglength',smoWin,'attributes',...
    {'stream_num' strmNumGrid @min...
    'chi' ChiGrid @mean...
    'ksn' ksnG @mean...
    'DrArea' A @mean});
shapewrite(MS,[chanDir,'/',fileTag, '_strmData.shp']);

% make an ascii of the channel segment chi-ksn values. (I can figure
% out how to make a shape file from this using topotoolbox. Any
% suggestions would be great.)
GRIDobj2ascii(ksnReg,[chanDir,'/',fileTag, '_chiksn_segs.txt']);
delete([profile_fig,map_fig]);

%% finally ask the user if they would like to make a chi map and ksn map
%% written as a shapefile for the entire drainage network
txt = sprintf('\nWould you like to make a chi map and a ksn map for your full drainage network?');
disp(txt);
txt = sprintf('Note: this may take a little time.\n');
disp(txt);
cont_opt = input('\ntype "y" for yes or "n" for no:  ','s');

while ~strcmp(cont_opt,'y') && ~strcmp(cont_opt,'n'),
    %case where you didn't enter a, b, c, or d:
    disp(' ');
    disp('Yes(y) or No(n) only!!!')
    cont_opt = input('type "y" for yes or "n" for no:  ','s');
end

if cont_opt == 'y'
    
    S1 = STREAMobj(FD,'minarea',crita/(DEM.cellsize^2));
    
    % Declare STREAMobj variables for faster processing through forloop
    ordList = S1.orderednanlist;
    strmBreaks = find(isnan(ordList));
    
    GridID = S1.IXgrid;
    
    Sz = double(DEM.Z(GridID));                 % elevation
    SmoZ = Sz;               % dumby vector to get smoothed data
    
    % get variables ready for chi integration
    chis = zeros(size(S1.distance));
    Six = S1.ix;
    Sixc = S1.ixc;
    Sx = S1.distance;
    Sa = (Ao./(A.Z(S1.IXgrid))).^mn;
    
    h = waitbar(0,'calculating \chi for full stream network...');
    % calculating chi for the entire river network
    for lp = numel(Six):-1:1;
        chis(Six(lp)) = chis(Sixc(lp)) + (Sa(Sixc(lp))+(Sa(Six(lp))-Sa(Sixc(lp)))/2) *(abs(Sx(Sixc(lp))-Sx(Six(lp))));
        f = (numel(Six)+1 - lp)/numel(Six);
        waitbar(f,h);
    end
    close(h);
    
    ChiGrid = DEM;
    ChiGrid.Z = nan(size(DEM.Z));
    ChiGrid.Z(S1.IXgrid) = chis;
    
    % plot all of the river river profile data as thin gray lines
    h = waitbar(0,'Smoothing elevation data for full stream network...');
    id1 = 0;
    for i = 1:length(strmBreaks);
        strmInds = ordList(id1+1:strmBreaks(i)-1);
        SmoZ(strmInds) = smoothChannelZ(Sz(strmInds),smoWin,cs);
        id1 = strmBreaks(i);
        f = i/length(strmBreaks);
        waitbar(f,h);
    end
    close(h)
    
    % making ksn map for streams analyzed
    ksnStreams = binnedKsn(S1,chis,SmoZ,smoWin,cs,Ao,mn);
    ksnG = DEM;
    ksnG.Z = nan(size(DEM.Z));
    ksnG.Z(GridID) = ksnStreams;
    
    MS = STREAMobj2mapstruct(S1,'seglength',smoWin,'attributes',...
        {'chi' ChiGrid @mean...
        'ksn' ksnG @mean});
    fileName = [chanDir,'/',fileTag, '_chi_ksn_map.shp'];
    shapewrite(MS,fileName);
elseif cont_opt == 'n'
end
end


