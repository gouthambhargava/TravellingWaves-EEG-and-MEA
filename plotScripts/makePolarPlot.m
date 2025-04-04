function makePolarPlot(phaseValues,binWidth,plotAxis,color)

% specify phase values and width of bin for plotting polar histogram

if nargin<4
    color = [0.8500 0.3250 0.0980];
end

phaseValues = reshape(phaseValues,[1,numel(phaseValues)]);
phaseValues(isnan(phaseValues)) = [];

%bin the data and get counts
nBins = 360/binWidth;
allBins = linspace(0,360,nBins);
wrappedAngles = wrapTo360(rad2deg(phaseValues));
meanAngle = circ_mean(deg2rad(wrappedAngles)');
meanAngle = deg2rad(mod(rad2deg(meanAngle),360));
binnedVals = discretize(wrappedAngles,allBins);
binnedVals(isnan(binnedVals)) = [];
uniqueVals = unique(binnedVals);
counts = zeros(1,numel(uniqueVals));
for i = 1:numel(uniqueVals)
    counts(uniqueVals(i)) = numel(find(binnedVals==uniqueVals(i)));
end
counts = counts/max(counts); %normalize counts 

%plot the counts
%get bin co-ordinates
binCord = exp(1i*deg2rad([allBins,allBins(2)]));

for i = 1:numel(counts)
    xcords = counts(i)*[real(binCord(i)),real(binCord(i+1))];
    ycords = counts(i)*[imag(binCord(i)),imag(binCord(i+1))];
    patch([0 xcords 0],[0 ycords 0],color,'parent',plotAxis)
    hold(plotAxis,'on')
end

arrowColor = color+(1-color)*0.2;
quiver(0,0,cos(meanAngle),sin(meanAngle),'AutoScaleFactor',0.9,'LineWidth',3,'Color',arrowColor,'parent',plotAxis)
hold(plotAxis,'on')

line([-1 1],[0 0],'Color',[0.811764705882353,0.811764705882353,0.811764705882353],'parent',plotAxis)
line([0 0],[-1 1],'Color',[0.811764705882353,0.811764705882353,0.811764705882353],'parent',plotAxis)
circle(0,0,1,[0,0,0],plotAxis);
circle(0,0,0.5,[0.811764705882353,0.811764705882353,0.811764705882353],plotAxis);
axis(plotAxis,[-1,1,-1,1])
% text(1.01,0,['0',char(176)])
% text(0,1.05,['90',char(176)])
% text(-1.2,0,['180',char(176)])
% text(0,-1.05,['-270',char(176)])
axis(plotAxis,'square')
axis(plotAxis,'off')


    function circle(x,y,r,colorVal,plotAxis)
        hold on
        th = 0:pi/50:2*pi;
        xunit = r * cos(th) + x;
        yunit = r * sin(th) + y;
        plot(plotAxis,xunit, yunit,'Color',colorVal);
        hold off
    end
end