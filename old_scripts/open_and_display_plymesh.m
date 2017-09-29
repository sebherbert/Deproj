

% open and parse the ply file produced by morphographX
[FVOri.vertices,FVOri.faces] = read_ply('processedMesh_bin.ply');
FV = FVOri;


% Display the reduced mesh
% Reduce the resolution of the faces (1% of the surface) => ~smoothing
FV = reducepatch(FVOri,0.01);
FV = FVOri;
figure; 
patchedMesh = patch(FV, 'FaceColor','none','LineWidth',1);
axis equal;


% Display mesh with z info encoded in colormap without transparency
figure; 
patchedMesh = patch(FV, 'FaceVertexCData',FV.vertices(:,3),'FaceColor','interp','LineStyle','none');
axis equal;


% Display mesh with z info encoded in colormap and transparency
figure; 
patchedMesh = patch(FV, 'FaceVertexCData',FV.vertices(:,3),'FaceColor','interp','LineStyle','none','FaceVertexAlphaData',0.3,...
    'FaceAlpha','flat');
axis equal;

% % Display mesh with normal info encoded in colormap
% figure; 
% patchedMesh = patch(FV, 'FaceVertexCData',N,'FaceColor','interp','LineStyle','none','FaceNormal',N);
% axis equal;

% Display vertex normals only
figure;
[FV.vertNorm,FV.faceNorm] = compute_normal(FV.vertices,FV.faces);
FV.vertNorm = FV.vertNorm';
FV.faceNorm = FV.faceNorm';
for face = 1:size(FV.faces,1)
    FV.faceCentroid(face,:) = mean(FV.vertices(FV.faces(face,:),:));
end
quiver3(FV.vertices(:,1),FV.vertices(:,2),FV.vertices(:,3),FV.vertNorm(:,1),FV.vertNorm(:,2),FV.vertNorm(:,3));
axis equal;
hold on
quiver3(FV.faceCentroid(:,1),FV.faceCentroid(:,2),FV.faceCentroid(:,3),FV.faceNorm(:,1),FV.faceNorm(:,2),FV.faceNorm(:,3));













%% % % 
% % % %% test normal from web
% % % theta = gallery('uniformdata',[100,1],0)*2*pi;
% % % phi = gallery('uniformdata',[100,1],1)*pi;
% % % x = cos(theta).*sin(phi);
% % % y = sin(theta).*sin(phi);
% % % z = cos(phi);
% % % DT = delaunayTriangulation(x,y,z);
% % % [T,Xb] = freeBoundary(DT);
% % % TR = triangulation(T,Xb);
% % % 
% % % figure
% % % trisurf(T,Xb(:,1),Xb(:,2),Xb(:,3), ...
% % %      'FaceColor', 'cyan', 'faceAlpha', 0.8);
% % % axis equal;
% % % hold on;
% % % 
% % % % Calculate the incenters and face normals. 
% % % P = incenter(TR);
% % % fn = faceNormal(TR);  
% % % 
% % % % Display the normal vectors on the surface. 
% % % quiver3(P(:,1),P(:,2),P(:,3), ...
% % %      fn(:,1),fn(:,2),fn(:,3),0.5, 'color','r');
% % % hold off;



