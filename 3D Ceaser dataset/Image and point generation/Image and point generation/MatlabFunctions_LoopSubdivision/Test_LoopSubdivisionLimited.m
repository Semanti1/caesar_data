% Test: Mesh subdivision using the Loop scheme.
%
% Author: Zexi Liu

%% Example: Face
clc;
clear;

num = 0;
model = 'nemo';

ppt = load( [ 'person_pt_HD_', num2str(num), '_plus' ] );
Ft = cell2mat(struct2cell(load( [model, '_feature_pts_125'] )));

pose = 1;

Idx = GetKeyIdx(ppt.uv,1,1);
u = ppt.uv(:,2*pose-1);
v = ppt.uv(:,2*pose);
X = u(Idx);
Y = v(Idx);
Tri = delaunay(X,Y);

figure;
subplot(1,2,1);
plotMesh( [X, Y], Tri, 'k' );
[ Res, ~, ~ ] = MeshResolution( [X, Y, zeros(numel(X),1) ], Tri );
subplot(1,2,2);
vertices = [X, Y, zeros(numel(X),1) ];
faces = Tri;
for i = 1:3
    [ vertices, faces ] = LoopSubdivisionLimited( vertices, faces, 1 );
end

plotMesh( vertices(:,1:2), faces, 'k' );
 

%% Example: Box
vertices = [10 10 10; -10 10 10; 10 -10 10; -10 -10 10; 10 10 -10; -10 10 -10; 10 -10 -10; -10 -10 -10];
faces = [1 2 3; 4 3 2; 1 3 5; 7 5 3; 1 5 2; 6 2 5; 8 6 7; 5 7 6; 8 7 4; 3 4 7; 8 4 6; 2 6 4];

figure(1);
subplot(1,4,1);
plotMesh(vertices, faces);
for i=2:4
 subplot(1,4,i);
 [vertices, faces] = loopSubdivision(vertices, faces); 
 plotMesh(vertices, faces);
end


%% Example: Tetrahedron
vertices = [10 10 10; -100 10 -10; -100 -10 10; 10 -10 -10]';
faces = [1 2 3; 1 3 4; 1 4 2; 4 3 2]';

figure(2);
subplot(1,4,1);
plotMesh(vertices, faces);
for i=2:4
 subplot(1,4,i);
 [vertices, faces] = loopSubdivision(vertices, faces); 
 plotMesh(vertices, faces);
end


%% Example: Cylinder
vertices = [0 -5 0; 0 5 0; 10 -5 0; 9.65 -5 2.58; 8.66 -5 5; 7.07 -5 7.07; 5 -5 8.66; 2.58 -5 9.65; 0 -5 10; -2.58 -5 9.65; -5 -5 8.66; -7.07 -5 7.07; -8.66 -5 5; -9.65 -5 2.58; -10 -5 0; -9.65 -5 -2.58; -8.66 -5 -5; -7.07 -5 -7.07; -5 -5 -8.66; -2.58 -5 -9.65; -0 -5 -10; 2.58 -5 -9.65; 5 -5 -8.66; 7.07 -5 -7.07; 8.66 -5 -5; 9.65 -5 -2.58; 10 5 0 ; 9.65 5 2.58; 8.66 5 5; 7.07 5 7.07; 5 5 8.66; 2.58 5 9.65; 0 5 10; -2.58 5 9.65; -5 5 8.66; -7.07 5 7.07; -8.66 5 5; -9.65 5 2.58; -10 5 0; -9.65 5 -2.58; -8.66 5 -5; -7.07 5 -7.07; -5 5 -8.66; -2.58 5 -9.65; -0 5 -10; 2.58 5 -9.65; 5 5 -8.66; 7.07 5 -7.07; 8.66 5 -5; 9.65 5 -2.58]';
faces = [1 3 4; 1 4 5; 1 5 6; 1 6 7; 1 7 8; 1 8 9; 1 9 10; 1 10 11; 1 11 12; 1 12 13; 1 13 14; 1 14 15; 1 15 16; 1 16 17; 1 17 18; 1 18 19; 1 19 20; 1 20 21; 1 21 22; 1 22 23; 1 23 24; 1 24 25; 1 25 26; 1 26 3; 2 28 27; 2 29 28; 2 30 29; 2 31 30; 2 32 31; 2 33 32; 2 34 33; 2 35 34; 2 36 35; 2 37 36; 2 38 37; 2 39 38; 2 40 39; 2 41 40; 2 42 41; 2 43 42; 2 44 43; 2 45 44; 2 46 45; 2 47 46; 2 48 47; 2 49 48; 2 50 49; 2 27 50; 3 27 28; 3 28 4; 4 28 29; 4 29 5; 5 29 30; 5 30 6; 6 30 31; 6 31 7; 7 31 32; 7 32 8; 8 32 33; 8 33 9; 9 33 34; 9 34 10; 10 34 35; 10 35 11; 11 35 36; 11 36 12; 12 36 37; 12 37 13; 13 37 38; 13 38 14; 14 38 39; 14 39 15; 15 39 40; 15 40 16; 16 40 41; 16 41 17; 17 41 42; 17 42 18; 18 42 43; 18 43 19; 19 43 44; 19 44 20; 20 44 45; 20 45 21; 21 45 46; 21 46 22; 22 46 47; 22 47 23; 23 47 48; 23 48 24; 24 48 49; 24 49 25; 25 49 50; 25 50 26; 26 50 27; 26 27 3]';

figure(3);
subplot(1,4,1);
plotMesh(vertices, faces);
for i=2:4
 subplot(1,4,i);
 [vertices, faces] = loopSubdivision(vertices, faces); 
 plotMesh(vertices, faces);
end


%% Example: Grid
vertices = [-1 -1 0; 1 -1 0; -1 1 0; 1 1 0;];
faces = [ 2 1 4; 1 3 4 ];

figure(4);
subplot(1,4,1);
trimesh( faces, vertices(:,1), vertices(:,2), 'Color', 'k' );
axis tight;
axis square; 
axis off;
view(3);
% plotMesh(vertices, faces,'k');
for i=2:4
 subplot(1,4,i);
 [vertices, faces] = LoopSubdivision(vertices, faces); 
 plotMesh(vertices, faces,'k');
end
