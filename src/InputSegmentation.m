

close all
clear
clc

mask_path = 'GT.tif';

% ImageJ mask.
I = imread( mask_path );

tic 
[ objects, junction_graph ] = mask_to_objects( I );
fprintf('Analyzed a %d x %d mask in %.1f seconds.\n', size(I,1), size(I,2), toc )

%% Plot everything.

figure
imshow( ~I , [ 0 2 ], 'Border', 'tight' )
hold on

plot( junction_graph, ...
    'XData', junction_graph.Nodes.Centroid(:,1), ...
    'YData', junction_graph.Nodes.Centroid(:,2), ...
    'LineWidth', 2, ...
    'EdgeColor', 'b', ...
    'EdgeAlpha', 1, ...
    'Marker', 'o', ...
    'MarkerSize', 4, ...
    'NodeColor', 'r' )
