% All you need to do is modify the variable in the line bracked by percent
% (%) signs and hit run. Directions will pop up in the command window as
% you go. For information on input and output files simply type "help
% kp_picker_db" in the matlab command window.
%
% Author: Sean F. Gallen (sean.gallen[at]erdw.ethz.ch)
% Date modified: 02/20/2017

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
demtxt = 'dem.tif';     % name of your DEM tiff file with extension
fileTag = 'proj1';      % tag used to identify specific files
crita = 1e6;            % threshold drainage area for channel head initiation
mn = 0.45;              % m/n or concavity index
Ao = 1;                 % reference drainage area for chi-analysis
smoWin = 250;           % size of window (in map units) used to smooth elevation data
flowOption = 'fill';    % Option for flow routing. Either 'carve' or 'fill'
                        % See help FLOWobj for more info on these options

% file path to your topotoolbox folder
addpath(genpath('C:\Users\sgallen\Documents\topo_toolbox\topotoolbox-master'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEM = runFlowPathApp(demtxt,crita,'flowOption',flowOption);
figure(1)
cont_opt = input('When you are finished type in the name of your STREAMobj variable here: ', 's');
disp(' ');

varTest = 0;
while varTest == 0;
    varTest = exist(cont_opt,'var');
    if varTest ~= 1
        txt = sprintf(['\nYour STREAMobj variable name and the name that you\n'...
            'just input do not match.\n']);
        disp(txt);
        cont_opt = input('Retype your STREAMobj variable name: ', 's');
    end
end
% close flowpathapp figures
close('Main');
close('Profiles');

% declare user picked stream object.
S = eval(cont_opt);
fileTag = [fileTag,'_', cont_opt];
chi_profiler(DEMc, S, fileTag, 'crita', crita, 'mn', mn, 'Ao', Ao,...
    'smoWin', smoWin);