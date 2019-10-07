function [timeStep,timeForPlotInUs] = timeStepClosestToSpecifiedTime(specifiedTimeInUs,availableTimesInUs)
%% Notes
% This function outputs the timeStep (1x1 double in dimensionless time units of the 
% simulation program) closest to specifiedTimeInUs (1x1 double in us) by searching 
% through the list of available times in availableTimesInUs (should be a
% vector double in us). Note that the rows of availableTimeInUs specify the
% time step for a particular time (e.g. if 1 us is located in the 100th row of availableTimeInUs,
% that means it's the 100th time step).

%% Function
% Calculate difference between available and specified time
tDif = abs(availableTimesInUs-specifiedTimeInUs);
% Sort the differences in ascending order
[~,index] = min(tDif);
% The first value of 'index' corresponds to the timeStep because the values
% of 'tDif' were sorted in ascending order. Subtract 1 because the file
% names start at 000000 instead of 000001
timeStep = index(1)-1;
timeForPlotInUs = availableTimesInUs(timeStep);

if specifiedTimeInUs > max(availableTimesInUs)
    warningText = ['Specified time outside available range: [' num2str(min(availableTimesInUs)) ' us, ' num2str(max(availableTimesInUs)) ' us]'];
    warning(warningText)
end

end

