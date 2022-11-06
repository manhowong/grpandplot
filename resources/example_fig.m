load carbig.mat
cars = table(org,when,cyl4,Horsepower,Weight);
filter = [cars.when~='M'];
carsNoMid = cars(filter(:,1),:);

close all

fig = figure("Position",[0 0 700 550]);
layout = tiledlayout(2,2,Padding="tight",TileSpacing="tight");

colors = [1 0 1; 0 0 1; 0 .5 .5];
t1And2 = grpandplot(carsNoMid,"Weight",yTitle='Weight (lbs)',xFactor="cyl4",cFactor="org",tFactor="when",...
                showXLine=true,showVln=true,showBox=false,showNum=true,numYPos=700,pntSize=5, vlnEdgeC='k',...
                cmap=colors,pntFillC='k',w=0.3,gap=0.5,parent=layout);

whenOrder = {'Early','Mid','Late'};
t3 = grpandplot(cars,"Weight",yTitle='Weight (lbs)',xFactor="when",cFactor="org",xOrder=whenOrder,...
                showXLine=true,showVln=true,showBox=false,showNum=true,numYPos=500,pntSize=2, ...
                pntFillC='k',vlnAlpha=0,xAxisMode=1,showLegend=false,parent=layout);
t3.XTickLabelRotation = 90;

t4 = grpandplot(cars,"Horsepower",yTitle='Horsepower',xFactor="when",cFactor="cyl4",xOrder=whenOrder,...
                pntSize=15,pntEdgeC='w',boxAlpha=0,boxEdgeC='k',pntOnTop=false,w=0.3,gap=0.8,parent=layout);
leg = legend(t4,Location="southoutside");
title(leg, "cyl4");