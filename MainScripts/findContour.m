

function [boundaries, cellIdx] = findContour( labelIm, PARAMS )
%% finds the contour around a set of pixels originating from a previous segmentation 

% find the correct order of the points defining the polygon contour of each
% individual cell in a CW order

% INPUT
% - dataSeg: a mask image of the segmentation ridges
% - PARAMS: a structure containing image specs
% OUTPUT
% - cellContout2D: structure with 3 fields
% - cellContour2D.Line: paired vertices matrix to infer the pixel
% connections -> deleted
% - cellContour2D.Vertices: XY sorted vertices -> deleted
% - cellContour2D.cellCt: XY sorted vertices in "CW" order

fprintf('Recreating the cellular contours\n');

% Find all the objects in the image
cc = bwconncomp(labelIm, 4);

% Prepare the morphological structuring element to recover the ridge element loss
se =strel('disk', 1);

% initialize output structure
boundaries = cell( cc.NumObjects, 1 );

% Extend each boundary to the center of the ridge
parfor bioCell = 1 : cc.NumObjects

    % Exclude cells touching the image border.
    % TODO elsewhere?
    
    emptyImage = zeros( cc.ImageSize, 'logical' );
    emptyImage( cc.PixelIdxList{bioCell} ) = 1;
    emptyImage = imdilate( emptyImage, se );
    
    b = bwboundaries( emptyImage, 8, 'noholes');
    
    % Only save objects that are not connected to the border
    if (min(b{1}(:,1)) == 1) ||...                        % Top
            (min(b{1}(:,2))== 1) ||...                    % Left
            (max(b{1}(:,1)) == cc.ImageSize(1)) ||...     % Bottom
            (max(b{1}(:,2)) == cc.ImageSize(2))           % Right
       continue 
    end
        
    boundaries{ bioCell }.cellCtPix = b{ 1 };
        
end

% Clean the cell structure of its empty positions
for bioCell = cc.NumObjects : -1 : 1 
    if isempty(boundaries{ bioCell })
        boundaries( bioCell ) = [];
    end
end

% figure % => To display all the segmented cells
% imshow( labelIm, [] )
% hold on
% for bioCell = 1 : numel(boundaries)
%     line( boundaries{ bioCell }.cellCt(:,2), boundaries.cellCt{ bioCell }(:,1) , 'Marker', '.' )  
% end

% Rescale to physical size
for bioCell = 1 : numel(boundaries)
    boundaries{ bioCell }.cellCt = boundaries{ bioCell }.cellCtPix * PARAMS.imSettings.latPixSize;
end

cellIdx = 1:numel(boundaries);

end












