
function displaySubplotCellContour(cellListError, cellListContour, dim)
% Display the cell contour of a cell list into a subplot

% Display the cells with a fit error
if numel(cellListError) > 20
    fprinf('Too many cells were not fitted properly. Only fitting the first 20 cells\n');
    nSub = 20;
else
    nSub = numel(cellListError);
end
[subPlotSize,~] = numSubplots(nSub);

figure
hold on
for bioCell = 1:nSub
   subplot(subPlotSize(1),subPlotSize(2),bioCell)
   plot(cellListContour{bioCell}(1,:),cellListContour{bioCell}(2,:))
   legend(sprintf('cell %d',cellListError(bioCell)))
end

supAxes=[.08 .08 .84 .84];
[~,~]=suplabel(sprintf('Not fitted cells (in %d dimension)\n',dim) ,'t', supAxes); 

end

