function [DEM] = runFlowPathApp(DEM,crita,varargin)
% runFlowPathApp.m takes a string of a GEOtiff (demtxt) and a critical
% drainage area for channel head initiation (crita) and runs the
% TopoToolBox flowpathapp so that users can hand select river profilers
% from a DEM.
%
% The script will automatically export a DEM GRIDobj.
%
% The user should export the streams selected to the matlab workspace as a
% STREAMobj and then type in the name of the stream object into the matlab
% command window as prompted. One can then use chi_profiler_STREAMobj.m for
% river profile analysis.
%
%   Inputs:
%       1) DEM: TopoToolBox DEM GRIDobj (required)
%       2) crita: Critical drainage area for channel head initiation in map
%          units. (reguired
%       3) flowOption: string of 'fill' or 'carve' for flow routing.
%          (optional) {default --> empty assumes DEM is filled or carved}
%
% Example:
%       DEM = runFlowPathApp(DEM,crita,'flowOption','carve');
%
% Author: Sean F. Gallen
% Date Modified: 02/20/2017


% Parse Inputs
p = inputParser;         
p.FunctionName = 'runFlowPathApp';
addRequired(p,'DEM', @(x) isa(x,'GRIDobj') | ischar(x));
addRequired(p,'crita', @(x) isscalar(x))

% optional inputs
addOptional(p,'flowOption', 'fill');

parse(p,DEM,crita,varargin{:});
DEM = p.Results.DEM;
crita = p.Results.crita;

if ischar(DEM)
    if ~isempty(strfind(DEM,'.tif') || strfind(DEM,'.txt'))
        % make GRIDobj and remove nan values from DEM
        DEM = GRIDobj(DEM);
        DEM.Z(DEM.Z<=-9999) = nan;
    else
        error('input file name does not have .tif or .txt file extension');
    end
end

% flow routing
if strcmp(p.Results.flowOption, 'fill');
    DEM = fillsinks(DEM);
    FD = FLOWobj(DEM);
elseif strcmp(p.Results.flowOption, 'carve');
    % Calculate flow object (FD) and distance to divide (dfd)
    FD = FLOWobj(DEM,'preprocess','carve');
    DEM = imposemin(FD,DEM);
else
    error('fillOption is not "fill" or "carve"');
end

S1 = STREAMobj(FD,'minarea', crita/(DEM.cellsize^2));

txt = sprintf(['\nPick streams from your DEM using the flowpathapp.\n'...
    'Export your streams to the workspace as STREAMobj.\n']);
disp(txt);

flowpathapp(FD,DEM,S1);
end