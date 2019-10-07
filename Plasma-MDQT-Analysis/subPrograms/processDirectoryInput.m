function [directory] = processDirectoryInput(directory)
%% Program Notes
    % In 'mainSimAnalysis', the user is prompted to select one or more
    % directories. The user should either select a single folder that contains
    % one or more simulation folders or the user should select one or more
    % simulation folders. By simulation folder, I refer to a folder whose name
    % begins with 'Ge...', which was produced by our MDQT c++ program.

    % The input 'directory' is a 1xm cell whose elements are full paths to folder(s) selected by the user, where m is the
    % number of folders selected. Note that the output 'directory' will be a nx1 cell whose elements are full paths to MDQT
    % data folders, where n just needs to be greater than zero. All the if statements and nonsense below just ensures that
    % 'directory' is outputted with this criteria in mind.

%% Format 'directory' properly when one folder supplied
% Check to make sure user has either selected one or more simulation
% folders or a single folder that contains one or more simulation folders
numOfFolders = length(directory);
if numOfFolders == 1
    % Check whether only supplied folder is sim folder, if not look for sim
    % folders within it.
    if ~(contains(directory{1},'Ge') && contains(directory{1},'Density') && contains(directory{1},'Sig'))
        % Import files from supplied directory
        files = struct2table(dir(directory{1}));
        files(1:2,:) = [];
        % Search for sim folders within contents of supplied directory
        ind = zeros(size(files.name));
        for i = 1:length(files.name)
            if (contains(files.name{i},'Ge') && contains(files.name{i},'Density') && contains(files.name{i},'Sig'))
                ind(i) = 1;
            end
        end
        % Keep sim folders found within supplied directory
        geFiles = files.name(logical(ind));
        % Check to make sure simulation folders were found
        if size(geFiles,1) == 0
            error('No simulation folders were found within user''s chosen directory.')
        end
        % Add found sim folders to 'directory' with full path
        for i = 1:length(geFiles)
            geFiles{i} = [directory{1} '\' geFiles{i}];
        end
        directory = geFiles;
    end
end

%% Format 'directory' properly when more than one folder supplied
% When 'directory' contains more than one entry, ensure that all specified
% folders are valid simulation folders.
if numOfFolders > 1
    % For each supplied folder, check whether it's a simulation folder
    ind = zeros(size(directory));
    for i = 1:length(directory)
        if (contains(directory{i},'Ge') && contains(directory{i},'Density') && contains(directory{i},'Sig'))
            ind(i) = 1;
        end
    end
    % Throw an error if any non-simulation folders were supplied
    if ~logical(min(ind))
        error('When specifying multiple directories, all folders must be valid simulation folders.')
    end
    % Transpose 'directory' so that it is a column vector
    directory = directory';
end

end

