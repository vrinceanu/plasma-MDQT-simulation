function [directory,f,ax] = analyzeMDQT(directory,useErrorBars)
%% Program Notes
    % This program's features have all been divided into subsections as found below. The user specifies via
    % 'directory' which data folders are to be analyzed. Then for each data folder the program will load the information
    % from it. In the last section that information will be plotted. Please see the notes in each individual section for 
    % more information.

    % Inputs
        % directory: (1xn cell) each element contains strings corresponding to the path of each selected data folder
            % Note: The user may select a single folder that contains more than one MDQT data folder or can select more
            % than one MDQT data folder when prompted by the 'uigetfile_n_dir' function.

    % Outputs
        % directory: (nx1 cell) elements contain strings corresponding to the full path of MDQT data folders to be analyzed
        % f: (1x1 fig obj) figure object for the figure plotted below in 'plotLoadedData'
        % ax: (1x1 axes obj) axes object contained by f

%% Process Directory Input
% This step ensures that the user selected valid folders for analysis and formats them appropriately. See inside the
% function for more details. 
[directory] = processDirectoryInput(directory);

%% Request User Input
% This is where the user is prompted to select program options, which are used in many other subfunctions below. The user
% will be asked what they want to plot, what units they want, what quantities to use in the legend, etc.
[flags] = askUserForInput(directory);

%% Hard-Coded Quantities
% All hard-coded quantities are contained within structure 'hc'
hc.eps = 8.854e-12; % Electric permitivity constant in SI units
hc.mi = 1.455e-25; % Sr+ mass in kg
hc.e = 1.602e-19; % Electric charge in C
hc.gammaSP = 1.41e8; % Natural linewidth of the 2S_(1/2)->2P_(3/2) transition in s^-1
hc.kb = 1.381e-23; % Boltzmann constant in SI units

%% Load Simulation Parameters from Foldere Names
% This step loads the initial simulation parameter information contained within
% each simulation folder name.
simParams = cell(size(directory));
for i = 1:length(directory)
    [simParams{i}] = loadSimParamsFromFolderName(directory{i},hc);
end

%% Obtain job info
% Inside each sim folder (e.g. 'Ge...' file) there are simulation results
% for each 'job' to be averaged over. This step determines which job
% numbers are present within each directory
jobs = cell(size(directory));
for i = 1:length(directory)
    [jobs{i}] = obtainJobInfo(directory{i});
end

%% Load information for plotting
% This step loads user-selected information from 'directory' and places it in the following matrices, which will be plotted
% in the 'plot data' section.
xData = cell(size(directory));
yData = cell(size(directory));
yDataErr = cell(size(directory));
plotInfo = cell(size(directory));
for i = 1:length(directory)
    [xData{i},yData{i},yDataErr{i},plotInfo{i}] = loadPlotInfo(directory{i},jobs{i},simParams{i},hc,flags);
end

%% Plot data loaded in previous section
[f,ax] = plotLoadedData(xData,yData,yDataErr,plotInfo,simParams,flags,useErrorBars);

end

