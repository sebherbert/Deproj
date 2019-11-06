

function [sides, vertices] = segNetwork( boundaries, dataSeg )

% This function recreates the cell junctions network for lighter display
% from a list of closed polygons based on their junctions (triple or more).
% any contour vertex that is simly part of a cell limit (not a junction)
% will be lost in the process.It will not work if the objects are not
% contiguous. In that case, use a polygon decimation function, but here we
% want to keep the exact position of the polygon contour.
% 
% % INPUT:
%   - boundaries: N-by-1 structure array each containing at least: a cellCt field
%   containing the ordered contour list of 2D vertices for a single polygon
%   in spatial unit, a cellCtPix field containing the same contour in pixel
%   units
% 
% OUTPUT:
%   - sides: All the polygon junctions between 2 vertices
%   - verticesPix: All the triple (or more) junctions between the polygons.
%   - vertices: All the triple (or more) junctions between the polygons.
%
% ALGORITHM:
% (1) overlay of all the contours in a single matrix. Any position >2 is a
% network vertex
% (2) Connect all the relevant pairs of vertices 
%
% AUTHOR: Sebastien HERBERT (herbert.sebastien@gmail.com)
% 

imSize = size(dataSeg);

[ verticesPix, vertices ] = findVertices( boundaries, imSize );


end


function [ verticesPix, vertices ] = findVertices( boundaries, imSize )
% find all the vertices of the segmentation network.

% Set the whole image to 0 occurences
allContours = zeros( imSize, 'uint8' );

% Find the occurence of each pixel in the contours list
for bioCell = 1:numel(boundaries)
    for pix = 1:length(boundaries{bioCell}.cellCtPix)
    allContours( boundaries{bioCell}.cellCtPix( pix , 1), ...
        boundaries{bioCell}.cellCtPix( pix , 2) ) = ...
        allContours( boundaries{bioCell}.cellCtPix( pix , 1), ...
        boundaries{bioCell}.cellCtPix( pix , 2) ) + 1 ;
    end
end

figure; imshow( allContours, [] )

end


function sides = findSides( verticesPix, boundaries )
%% for each polygon contour find the included vertices and recreate the sides in a 
% CW manner by following the contour




end