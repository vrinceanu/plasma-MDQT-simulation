function [f,ax] = plotLoadedData(xData,yData,yDataErr,plotInfo,simParams,flags,useErrorBars)
%% Create string for legend
% Ask user which simParam to use for legend if more than one data set
% supplied
if length(plotInfo) > 1
    %   Ask user what quantity to use in the legend
    baseLegendOffWhat = flags.baseLegendOffWhat;
    newLegend = {};
    for i = 1:length(xData)
        %   Load simulation parameters
        [currSimParam] = simParams{i};
        if baseLegendOffWhat == 1 % Use density in the legend
            legendNum(i) = currSimParam.n*1e-14;
            legendStr = num2str(legendNum(i));
            newLegend{length(newLegend)+1} = ['n = ' legendStr 'x10^1^4 m^-^3'];
        elseif baseLegendOffWhat == 2 % Use electron temperature
            legendNum(i) = currSimParam.Te;
            legendStr = num2str(legendNum(i));
            newLegend{length(newLegend)+1} = ['T_e = ' legendStr ' K'];
        elseif baseLegendOffWhat == 3% Use initial plasma size
            legendNum(i) = currSimParam.sigma*1e3;
            legendStr = num2str(legendNum(i));
            newLegend{length(newLegend)+1} = ['\sigma(0) = ' legendStr ' mm'];
        elseif baseLegendOffWhat == 4% Use SP rabi frequency
            legendNum(i) = currSimParam.omegaSP*1e-6;
            legendStr = num2str(legendNum(i));
            newLegend{length(newLegend)+1} = ['\Omega_S_P = ' legendStr ' MHz'];
        end
    end
    % Sort everything by what's in the legend
    [~,sortInd] = sort(legendNum,'descend');
    xData = xData(sortInd);
    yData = yData(sortInd);
    yDataErr = yDataErr(sortInd);
    plotInfo = plotInfo(sortInd);
    simParams = simParams(sortInd);
    newLegend = newLegend(sortInd);
    % Generate string for legend
    legendStr = cell(1,length(plotInfo)*length(plotInfo{1}.legend));
    ind = 0;
    for i = 1:length(plotInfo)
        for j = 1:length(plotInfo{i}.legend)
            ind = ind+1;
            legendStr{ind} = [plotInfo{i}.legend{j} ' for ' newLegend{i}];
        end
    end
else
    newLegend = [];
end

%% Option for plotting state population vs initial parameter
%
if ~isempty(flags.plotForSpecificV)
    [xData,yData,yDataErr,plotInfo] = plotStatePopForSpecificV(xData,yData,yDataErr,simParams,plotInfo,flags);
end


%% Do plotting
colorInd = [0 0.4470 0.7410;0.8500 0.3250 0.0980;0.9290 0.6940 0.1250;0.4940 0.1840 0.5560;0.4660 0.6740 0.1880;0.3010 0.7450 0.9330;0.6350 0.0780 0.1840;1 0 0;0 1 0;0 0 1;0 1 1;1 0 1;1 1 0;0 0 0];
% Create figure - edit some properties
f = figure;
f.Units = 'inches';
f.Position = [2.5 2.5 7 6];
hold on
% Plot data within figure
ind = 0;
for i = 1:length(xData)
    if useErrorBars == 0
        plot(xData{i},yData{i},'Linewidth',2)
    elseif useErrorBars == 1
        for j = 1:size(yData{i},2)
            ind = ind+1;
            lineProps.col{j} = colorInd(ind,:);
        end
        if isempty(flags.plotForSpecificV)  % use shaded error bars when plotting multiple data sets
            mseb(xData{i}',yData{i}',yDataErr{i}',lineProps)
        else    % otherwise, use standard error bars
            errorbar(xData{i}',yData{i}',yDataErr{i}')
        end
    end
end
% Edit axes properties of figure
ax = gca;
ax.XLabel.String = plotInfo{1}.xlabel;
ax.YLabel.String = plotInfo{1}.ylabel;
ax.Title.String = plotInfo{1}.title;
ax.FontSize = 12;

if ~isempty(newLegend) && length(xData)>1
    legend(legendStr)
else
    legend(plotInfo{1}.legend)
end

if isfield(plotInfo{1},'markerSize')
    ax.Children.Marker = '.';
    ax.Children.MarkerSize = plotInfo{1}.markerSize;
    ax.Children.LineStyle = 'none';
end

end
