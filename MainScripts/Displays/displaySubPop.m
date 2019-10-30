


function displaySubPop(PARAMS,sidesList,allVertices,fullLegs,dispPARAMS,dataCurv)
%% Display subpopulations


figure
hold on

for graphPt = 1:length(sidesList)
    h{graphPt} = plotSeg(sidesList{graphPt},allVertices,PARAMS);
    h{graphPt}.EdgeColor = dispPARAMS{graphPt}.EdgeColor;
    h{graphPt}.LineWidth = dispPARAMS{graphPt}.LineWidth;
end

% Overlay with mesh in blue only if the dataCurv arg is provided
if nargin == 6
    h{length(sidesList)+1} = plotMesh(dataCurv,PARAMS);
    h{length(sidesList)+1}.FaceColor = [0;0.45;0.74];
end

axis equal
box

legend(fullLegs);

end