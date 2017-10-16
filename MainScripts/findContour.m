
%% finds the contour around a set of pixels originating from a previous segmentation 

function cellContour2D = findContour(contour2Dpos,bmethod,PARAMS)
% find the correct order of the points defining the polygon contour of the
% cell instead of sorted by XY position
% Also create the connectivity pairs
% INPUT
% - contour2Dpos: a cell array with N cells each representing the biological
% cells each of dim = nx2 for the n pixels composing the segmented cell in
% 2D.
% - bmethod: 1 of 4 possibilities of boundary contour method
%        - msquare
%        - msquareImproved: method improved with dilatation and erosion
%        - bwboundary: default Matlab method
%        - mask2poly: fex boundary finding method, very close to bwboundary
%        although with more options
% OUTPUT
% - cellContout2D: structure with 3 fields
% - cellContour2D.Line: paired vertices matrix to infer the pixel
% connections -> deleted
% - cellContour2D.Vertices: XY sorted vertices -> deleted
% - cellContour2D.cellCt: XY sorted vertices in "CW" order

fprintf('Recreating the cellular contours\n');



%% Create a temporary binary image of each individual cell before applying a boundary finding method
for bioCell = 1:length(contour2Dpos) % for all cells
    clear miniContour miniImg
    contour2Dpos{bioCell} = round(contour2Dpos{bioCell}/PARAMS.imSettings.latPixSize); % rescale from µm to pix for image creation
    miniContour=[contour2Dpos{bioCell}(:,1) contour2Dpos{bioCell}(:,2)]-(min(contour2Dpos{bioCell})-1);
    miniImg = zeros(max(miniContour(:,1)),max(miniContour(:,2))); % columns and lines are switched
    for i = 1:size(miniContour,1)
        miniImg(miniContour(i,1),miniContour(i,2))=1;
    end
    
    % append 1 layer in all directions to avoid border effect
    % could be merged with the previous paragraph drawing the contour but simpler like this
    miniImg = [zeros(1,size(miniImg,2)) ; miniImg ; zeros(1,size(miniImg,2))];
    miniImg = [zeros(size(miniImg,1),1) , miniImg , zeros(size(miniImg,1),1)];
    %     figure;imshow(testIm)
    
    % Fill contour
    miniImg = imfill(miniImg,'holes');
    
    switch bmethod
        
        case 'msquare'
            % Find contour
            [cellContour2D{bioCell}.Lines,cellContour2D{bioCell}.Vertices]=isocontour(miniImg,0.5);
            % Replace in full image % -1 for start table at 1 and -1 for
            % appened lines
            cellContour2D{bioCell}.Vertices = cellContour2D{bioCell}.Vertices+...
                min(contour2Dpos{bioCell})-2; % replace in full image % -1 for start table at 1 and -1 for appened lines
            
        case 'msquareImproved'
            % Enlarge and erode
            miniImg = imresize(miniImg, 2);
            selem = ones(3,3); %// square, 8-Negihbours
            miniImg = imerode(miniImg, selem);
            % Find contour
            [cellContour2D{bioCell}.Lines,cellContour2D{bioCell}.Vertices]=isocontour(miniImg,0.5);
            % Replace in full image % -1 for start table at 1 and -1 for
            % appened lines and /2 for the enlargment             
            cellContour2D.Vertices = (cellContour2D{bioCell}.Vertices)/2-2+...
                min(contour2Dpos{bioCell});
            
        case 'bwboundary'
            [boundaries,~] = bwboundaries(miniImg,'noholes');
            cellContour2D{bioCell}.cellCt = boundaries{1}+...
                min(contour2Dpos{bioCell})-2; % replace in full image % -1 for start table at 1 and -1 for appened lines
            
        case 'mask2poly'
            polyCell=mask2poly(miniImg,'Inner','MINDIST');
            % switch columns 1 and 2
            polyCell(:,3) = polyCell(:,1);
            polyCell(:,1) = polyCell(:,2);
            polyCell(:,2) = polyCell(:,3);
            polyCell(:,3) = [];
            cellContour2D{bioCell}.cellCt = polyCell+...
                min(contour2Dpos{bioCell})-2; % replace in full image % -1 for start table at 1 and -1 for appened lines
        
        otherwise
            fprintf('Error: Unknown contour method: %s\n', bmethod);
    end
    
    switch bmethod % this step to reorganize the msquare output
        case {'msquare','msquareImproved'}
            Ori = 1;
            tempVert = zeros(length(cellContour2D{bioCell}.Lines),2);
            for lines = 1:length(cellContour2D{bioCell}.Lines)
                tempVert(lines,:) = cellContour2D{bioCell}.Vertices(cellContour2D{bioCell}.Lines(Ori,1),:);
                Dest = cellContour2D{bioCell}.Lines(Ori,2);
                if length(find(cellContour2D{bioCell}.Lines(:,1)==Dest))~=1
                    fprintf('Error: a contour line is open. Check cell %d at Origin = %d\n',...
                        bioCell,Dest);
                    break
                end
                Ori = find(cellContour2D{bioCell}.Lines(:,1)==Dest);
            end
            cellContour2D{bioCell}.cellCt = tempVert;
            clear tempVert
    end
    cellContour2D{bioCell}.cellCt = cellContour2D{bioCell}.cellCt*PARAMS.imSettings.latPixSize; % scale back to µm
end
cellContour2D = cellContour2D';

end