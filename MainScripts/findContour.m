

function boundaries = findContour( labelIm, PARAMS )
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



% %% Create a temporary binary image of each individual cell before
% % applying a boundary finding method == OLD VERSION
% for bioCell = 1:length(contour2Dpos) % for all cells
%     clear miniContour miniImg
%     contour2Dpos{bioCell} = round(contour2Dpos{bioCell}/PARAMS.imSettings.latPixSize); % rescale from µm to pix for image creation
%     miniContour=[contour2Dpos{bioCell}(:,1) contour2Dpos{bioCell}(:,2)]-(min(contour2Dpos{bioCell})-1);
%     miniImg = zeros(max(miniContour(:,1)),max(miniContour(:,2))); % columns and lines are switched
%     for i = 1:size(miniContour,1)
%         miniImg(miniContour(i,1),miniContour(i,2))=1;
%     end
%     
%     % append 1 layer in all directions to avoid border effect
%     % could be merged with the previous paragraph drawing the contour but simpler like this
%     miniImg = [zeros(1,size(miniImg,2)) ; miniImg ; zeros(1,size(miniImg,2))];
%     miniImg = [zeros(size(miniImg,1),1) , miniImg , zeros(size(miniImg,1),1)];
%     %     figure;imshow(testIm)
%     
%     % Fill contour
%     miniImg = imfill(miniImg,'holes');
%     
%     switch bmethod
%         
%         case 'msquare'
%             % Find contour
%             [cellContour2D{bioCell}.Lines,cellContour2D{bioCell}.Vertices]=isocontour(miniImg,0.5);
%             % Replace in full image % -1 for start table at 1 and -1 for
%             % appened lines
%             cellContour2D{bioCell}.Vertices = cellContour2D{bioCell}.Vertices+...
%                 min(contour2Dpos{bioCell})-2; % replace in full image % -1 for start table at 1 and -1 for appened lines
%             
%         case 'msquareImproved'
%             % Enlarge and erode
%             miniImg = imresize(miniImg, 2);
%             selem = ones(3,3); %// square, 8-Negihbours
%             miniImg = imerode(miniImg, selem);
%             % Find contour
%             [cellContour2D{bioCell}.Lines,cellContour2D{bioCell}.Vertices]=isocontour(miniImg,0.5);
%             % Replace in full image % -1 for start table at 1 and -1 for
%             % appened lines and /2 for the enlargment             
%             cellContour2D.Vertices = (cellContour2D{bioCell}.Vertices)/2-2+...
%                 min(contour2Dpos{bioCell});
%             
%         case 'bwboundary'
%             [boundaries,~] = bwboundaries(miniImg,'noholes'); % simpler call is boundaries = bwboundaries(miniImg,'noholes');
%             cellContour2D{bioCell}.cellCt = boundaries{1}+...
%                 min(contour2Dpos{bioCell})-2; % replace in full image % -1 for start table at 1 and -1 for appened lines
%             
%         case 'mask2poly'
%             polyCell=mask2poly(miniImg,'Inner','MINDIST');
%             % switch columns 1 and 2
%             polyCell(:,3) = polyCell(:,1);
%             polyCell(:,1) = polyCell(:,2);
%             polyCell(:,2) = polyCell(:,3);
%             polyCell(:,3) = [];
%             cellContour2D{bioCell}.cellCt = polyCell+...
%                 min(contour2Dpos{bioCell})-2; % replace in full image % -1 for start table at 1 and -1 for appened lines
%         
%         otherwise
%             fprintf('Error: Unknown contour method: %s\n', bmethod);
%     end
%     
%     switch bmethod % this step to reorganize the msquare output
%         case {'msquare','msquareImproved'}
%             Ori = 1;
%             tempVert = zeros(length(cellContour2D{bioCell}.Lines),2);
%             for lines = 1:length(cellContour2D{bioCell}.Lines)
%                 tempVert(lines,:) = cellContour2D{bioCell}.Vertices(cellContour2D{bioCell}.Lines(Ori,1),:);
%                 Dest = cellContour2D{bioCell}.Lines(Ori,2);
%                 if length(find(cellContour2D{bioCell}.Lines(:,1)==Dest))~=1
%                     fprintf('Error: a contour line is open. Check cell %d at Origin = %d\n',...
%                         bioCell,Dest);
%                     break
%                 end
%                 Ori = find(cellContour2D{bioCell}.Lines(:,1)==Dest);
%             end
%             cellContour2D{bioCell}.cellCt = tempVert;
%             clear tempVert
%     end

%%


% Exclude cells touching the image border.
labelIm = excludeBorderObjects( labelIm );

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
    
    boundaries{ bioCell } = b{ 1 } * PARAMS.imSettings.latPixSize;
        
end

% figure % => To display all the segmented cells
% imshow( labelIm, [] )
% hold on
% for bioCell = 1 : cc.NumObjects
%     line( boundaries{ bioCell }(:,2), boundaries{ bioCell }(:,1) , 'Marker', '.' )  
% end

end



function labelIm_NoBorder = excludeBorderObjects( labelIm )
% Exclude objects that are connected to the border
% Keep in mind that objects 1 pixel away from the border will not be
% rejected!
% Only works for label image

% Copy original image
labelIm_NoBorder = labelIm;

% Find all the objects of interest
objectsList2reject = unique( labelIm( : , 1) );
objectsList2reject = unique( [ objectsList2reject ; labelIm( : , end) ] );
objectsList2reject = unique( [ objectsList2reject ; labelIm( 1 , :)' ] );
objectsList2reject = unique( [ objectsList2reject ; labelIm( end , :)' ] );

% Set all objects to reject to 0
for object = 1:numel(objectsList2reject)
    labelIm_NoBorder( labelIm_NoBorder == objectsList2reject(object) ) = 0;
end


end













