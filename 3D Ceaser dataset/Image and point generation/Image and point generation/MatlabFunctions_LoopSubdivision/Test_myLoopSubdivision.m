% Test: Mesh subdivision using the Loop scheme.
%
% Author: Jesus Mena


%% Example: Face
clc;
clear;

num = 0;
model = 'nemo';
pose = 1;

ppt = load( [ 'person_pt_HD_', num2str(num), '_plus' ] );
Ft = cell2mat(struct2cell(load( [model, '_feature_pts_125'] )));
gmt = cell2mat(struct2cell(load( [model,'_125'], 'Tri' )));
gmv = cell2mat(struct2cell(load( [model,'_125'], 'Ver' )));

% x = gmv(:,1);
% y = gmv(:,2);
% z = gmv(:,3);
% xlin = linspace(min(x),max(x),1280);
% ylin = linspace(min(y),max(y),1280);
% 
% [X,Y] = meshgrid(xlin,ylin);
% tic
% Z = griddata(x,y,z,X,Y,'cubic');toc
% Z = (Z + min(min(Z)))/max(max(Z));
% figure(1)
% % mesh(X,Y,Z) %interpolated
% % axis tight; hold on
% % h2=plot3(X,Y,Z,'.','MarkerSize',15) %nonuniform
% % %set(h2,'MarkerSize',[3.0]);
% imshow(flipud(Z))

% return;

Idx = GetKeyIdx(ppt.uv,1,1);
u = ppt.uv(:,2*pose-1);
v = ppt.uv(:,2*pose);
X = u(Idx);
Y = v(Idx);
Tri = delaunay(X,Y);
Ft = Ft(Idx',:);

figure;
ri1 = randi(size(Tri,1),1,6);
subplot(2,2,1);
plotMesh([X, -Y], Tri, 'k');
for i = 1 : length(ri1)
    cc(i) = ChooseColors(i);
    triplot(Tri(ri1(i),:), X, -Y, cc(i),'lineWidth', 2);
    for j = 1 : 3
        text( X(Tri(ri1(i),j)), -Y(Tri(ri1(i),j)), num2str(j) );
    end
end
subplot(2,2,2);
plotMesh(Ft(:,1:2), Tri, 'k');
for i = 1 : length(ri1)
    triplot(Tri(ri1(i),:), Ft(:,1), Ft(:,2), cc(i),'lineWidth', 2);
    for j = 1 : 3
        text( Ft(Tri(ri1(i),j),1), Ft(Tri(ri1(i),j),2), num2str(j) );
    end
end

ver1 = [X, Y, zeros(numel(X),1) ];
ver2 = [ Ft(:,1:2), zeros(numel(X),1) ];
% ver2 = Ft(:,1:3);
faces1 = Tri;
faces2 = Tri;
tic;
for i = 1:5
    [ ver1, faces1, sync ] = myLoopSubdivision( ver1, faces1 );
    [ ver2, faces2 ] = myLoopSubdivision( ver2, faces2, sync );
end
toc;

subplot(2,2,3);
ri2 = randi(size(faces1,1),1,6);
plotMesh( [ver1(:,1),-ver1(:,2)], faces1, 'k');
for i = 1 : length(ri2)
    cc(i) = ChooseColors(i);
    triplot(faces1(ri2(i),:), ver1(:,1),-ver1(:,2), cc(i),'lineWidth', 2);
    for j = 1 : 3
        text( ver1(faces1(ri2(i),j),1), -ver1(faces1(ri2(i),j),2), num2str(j) );
    end
end
subplot(2,2,4);
plotMesh( [ver2(:,1),ver2(:,2)], faces2, 'k');
for i = 1 : length(ri2)
    triplot( faces2(ri2(i),:), ver2(:,1),ver2(:,2), cc(i),'lineWidth', 2);
    for j = 1 : 3
        text( ver2(faces2(ri2(i),j),1), ver2(faces2(ri2(i),j),2), num2str(j) );
    end
end

% return;

[ ~, Z ] = VisibleVertices( ver2, gmv, gmt, [0,0,1], 0 );

I = imread( 'P0_1.JPG' );

figure;
imshow( I, 'InitialMagnification', 'fit' );hold on;
plot3(ver1(:,1),ver1(:,2),Z*500,'.b')

L2i(:,1) = round( ver1(:,1) );
L2i(:,2) = round( ver1(:,2) );

Ind = sub2ind( [ 1280, 1280 ], L2i(:,2), L2i(:,1) );

R = I(:,:,1);
G = I(:,:,2);
B = I(:,:,3);
VCData_R( :, 1 ) = double( R(Ind) ) / 255;
VCData_G( :, 1 ) = double( G(Ind) ) / 255;
VCData_B( :, 1 ) = double( B(Ind) ) / 255;

campos_t = [ 0  0 20 ]';
figure;
ph = patch( 'Faces', faces1, 'Vertices', [ver1(:,1),-Z*400,-ver1(:,2)] ); axis image; hold on;
set( ph, 'EdgeColor', 'none' );
campos( [campos_t(1) campos_t(2) campos_t(3)] );    % camera pose
xlabel('X');ylabel('Y');zlabel('Z');
view([0 0]);
set(ph, 'FaceVertexCData', [VCData_R, VCData_G, VCData_B] );
set(ph, 'FaceColor', 'flat');

 return;
 

%% Example: Box
vertices = [10 10 10; -10 10 10; 10 -10 10; -10 -10 10; 10 10 -10; -10 10 -10; 10 -10 -10; -10 -10 -10];
faces = [1 2 3; 4 3 2; 1 3 5; 7 5 3; 1 5 2; 6 2 5; 8 6 7; 5 7 6; 8 7 4; 3 4 7; 8 4 6; 2 6 4];

figure;
subplot(1,4,1);
plotMesh(vertices, faces,'k');
for i = 2:4
 subplot(1,4,i);
 [vertices, faces] = LoopSubdivision(vertices, faces); 
 plotMesh(vertices, faces,'k');
end


%% Example: Tetrahedron
vertices = [10 10 10; -100 10 -10; -100 -10 10; 10 -10 -10];
faces = [1 2 3; 1 3 4; 1 4 2; 4 3 2];

figure;
subplot(1,4,1);
plotMesh(vertices, faces,'k');
for i=2:4
 subplot(1,4,i);
 [vertices, faces] = LoopSubdivision(vertices, faces); 
 plotMesh(vertices, faces,'k');
end


%% Example: Cylinder
vertices = [0 -5 0; 0 5 0; 10 -5 0; 9.65 -5 2.58; 8.66 -5 5; 7.07 -5 7.07; 5 -5 8.66; 2.58 -5 9.65; 0 -5 10; -2.58 -5 9.65; -5 -5 8.66; -7.07 -5 7.07; -8.66 -5 5; -9.65 -5 2.58; -10 -5 0; -9.65 -5 -2.58; -8.66 -5 -5; -7.07 -5 -7.07; -5 -5 -8.66; -2.58 -5 -9.65; -0 -5 -10; 2.58 -5 -9.65; 5 -5 -8.66; 7.07 -5 -7.07; 8.66 -5 -5; 9.65 -5 -2.58; 10 5 0 ; 9.65 5 2.58; 8.66 5 5; 7.07 5 7.07; 5 5 8.66; 2.58 5 9.65; 0 5 10; -2.58 5 9.65; -5 5 8.66; -7.07 5 7.07; -8.66 5 5; -9.65 5 2.58; -10 5 0; -9.65 5 -2.58; -8.66 5 -5; -7.07 5 -7.07; -5 5 -8.66; -2.58 5 -9.65; -0 5 -10; 2.58 5 -9.65; 5 5 -8.66; 7.07 5 -7.07; 8.66 5 -5; 9.65 5 -2.58];
faces = [1 3 4; 1 4 5; 1 5 6; 1 6 7; 1 7 8; 1 8 9; 1 9 10; 1 10 11; 1 11 12; 1 12 13; 1 13 14; 1 14 15; 1 15 16; 1 16 17; 1 17 18; 1 18 19; 1 19 20; 1 20 21; 1 21 22; 1 22 23; 1 23 24; 1 24 25; 1 25 26; 1 26 3; 2 28 27; 2 29 28; 2 30 29; 2 31 30; 2 32 31; 2 33 32; 2 34 33; 2 35 34; 2 36 35; 2 37 36; 2 38 37; 2 39 38; 2 40 39; 2 41 40; 2 42 41; 2 43 42; 2 44 43; 2 45 44; 2 46 45; 2 47 46; 2 48 47; 2 49 48; 2 50 49; 2 27 50; 3 27 28; 3 28 4; 4 28 29; 4 29 5; 5 29 30; 5 30 6; 6 30 31; 6 31 7; 7 31 32; 7 32 8; 8 32 33; 8 33 9; 9 33 34; 9 34 10; 10 34 35; 10 35 11; 11 35 36; 11 36 12; 12 36 37; 12 37 13; 13 37 38; 13 38 14; 14 38 39; 14 39 15; 15 39 40; 15 40 16; 16 40 41; 16 41 17; 17 41 42; 17 42 18; 18 42 43; 18 43 19; 19 43 44; 19 44 20; 20 44 45; 20 45 21; 21 45 46; 21 46 22; 22 46 47; 22 47 23; 23 47 48; 23 48 24; 24 48 49; 24 49 25; 25 49 50; 25 50 26; 26 50 27; 26 27 3];

figure;
subplot(1,4,1);
plotMesh(vertices, faces,'k');
for i=2:4
 subplot(1,4,i);
 [vertices, faces] = myLoopSubdivision(vertices, faces); 
 plotMesh(vertices, faces,'k');
end


%% Example: Grid
vertices = [-4 -4 0; -2 -4 0; 0 -4 0; 2 -4 0; 4 -4 0; -4 -2 0; -2 -2 0; 0 -2 0; 2 -2 0; 4 -2 0; -4 0 0; -2 0 0; 0 0 0; 2 0 0; 4 0 0; -4 2 0; -2 2 0; 0 2 0; 2 2 0; 4 2 0; -4 4 0; -2 4 0; 0 4 0; 2 4 0; 4 4 0];
faces = [7 2 1; 1 6 7; 8 3 2; 2 7 8; 9 4 3; 3 8 9; 10 5 4; 4 9 10; 12 7 6; 6 11 12; 13 8 7; 7 12 13; 14 9 8; 8 13 14; 15 10 9; 9 14 15; 17 12 11; 11 16 17; 18 13 12; 12 17 18; 19 14 13; 13 18 19; 20 15 14; 14 19 20; 22 17 16; 16 21 22; 23 18 17; 17 22 23; 24 19 18; 18 23 24; 25 20 19; 19 24 25];

figure(4);
subplot(1,4,1);
plotMesh(vertices, faces,'k');
for i=2:4
 subplot(1,4,i);
 [vertices, faces] = myLoopSubdivision(vertices, faces); 
 plotMesh(vertices, faces,'k');
end

