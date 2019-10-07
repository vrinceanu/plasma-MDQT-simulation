function [xData,yData,yDataErr,plotInfo] = plotStatePopForSpecificV(xData,yData,yDataErr,simParams,plotInfo,flags)
%% Filter data based on user's preference
% For each data set, find index of xData that is closest to user-specified velocity
statePop = zeros(length(xData),size(yData{1},2));
for i = 1:length(xData)
    % Find index of velocity within xData that is closest to user specified velocity
    vDiff = abs(xData{i}-flags.plotForSpecificV);
    [~,minInd] = min(vDiff);
    % Extract the value of yData at that index
    statePop(i,:) = yData{i}(minInd,:);
    statePopErr(i,:) = yDataErr{i}(minInd,:);
end
yData = [];
yDataErr = [];
yData{1} = statePop;
yDataErr{1} = statePopErr;
% Assign values to 'xData' based on user's selection of initial parameter.
xData = [];
xData{1} = zeros(size(plotInfo));
baseLegendOffWhat = flags.baseLegendOffWhat;
for i = 1:length(plotInfo)
    for j = 1:length(xData)
        %   Load simulation parameters
        [currSimParam] = simParams{i};
        if baseLegendOffWhat == 1 % Use density in the legend
            xData{1}(i) = currSimParam.n*1e-14;
            plotInfo{i}.xlabel = 'Density (1x10^1^4 m^-^3)';
        elseif baseLegendOffWhat == 2 % Use electron temperature
            xData{1}(i) = currSimParam.Te;
            plotInfo{i}.xlabel = 'T_e (K)';
        elseif baseLegendOffWhat == 3% Use initial plasma size
            xData{1}(i) = currSimParam.sigma*1e3;
            plotInfo{i}.xlabel = '\sigma (mm)';
        elseif baseLegendOffWhat == 4% Use SP rabi frequency
            xData{1}(i) = currSimParam.omegaSP*1e-6;
            plotInfo{i}.xlabel = '\Omega_s_p (Hz)';
        end
    end
    plotInfo{i}.title = [plotInfo{i}.title 'and v = ' num2str(flags.plotForSpecificV) ' m/s'];
    plotInfo{i}.markerSize = 20;
end

end

