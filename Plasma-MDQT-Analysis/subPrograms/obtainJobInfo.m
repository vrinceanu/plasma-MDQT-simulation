function [jobs] = obtainJobInfo(directory)
%% FUNCTION NOTES
% This function figures out how many 'job' files are contained within 'directory'. 
% Many different "jobs" (or instances) of the simulation program may be
% submitted to the cluster at the same time that all have the same exact
% simulation parameters. Inside 'directory', there will be a series of
% 'job' folders that contain the simulation data for each instance of the
% program that run. It's important to know how many jobs are contained
% within 'directory' prior to beginning analysis.

%   Outputs
%       - jobs: (nx1 double, n is total number of jobs) each element contains a unique job number


%% OBTAIN JOB INFORMATION
% Import job folders contained within simulation data folder
files = struct2table(dir(directory));
files(1:2,:) = [];

% Iterate through each job file and extract the job number into the 'jobs' vector
jobs = zeros(size(files.name));
for ii = 1:length(files.name)
    jobs(ii) = str2double(extractAfter(files.name{ii},'job'));
end

% Sort the jobs vector in ascending order
jobs = sort(jobs);

end

