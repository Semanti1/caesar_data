% --------------------------- Descriptions --------------------------------
% mesh subdivision based on Loop's algorithm
% 
% File name:    myLoopSubdivision.m
% Date created: 03/21/2013
% Last revise:  03/22/2013
% -------------------------------------------------------------------------

function [ new_Ver, new_Tri, sync ] = myLoopSubdivision( Ver, Tri, varargin )


%% init
NoV_old = size(Ver,1);
NoT = size(Tri,1);
TR = TriRep( Tri, Ver );
B_old = freeBoundary( TR );
if size( B_old ~= 0, 1 )   % open surface
    tmp = [ B_old(:,2),B_old(:,1) ];
    B_old = [ B_old; tmp ];
    B_old = sortrows(B_old,[1,2]);
%     [ ~, Bfo_old, ~ ] = unique( B_old(:,1), 'first' );    % first occurrence
%     [ ~, Blo_old, ~ ] = unique( B_old(:,1), 'last' );     % last occurrence
end

SI = reshape( repmat((1:NoT),3,1), NoT*3, 1 );  % triangles indices in the form of: 1 1 1 2 2 2 3 3 3 4 4 4 ...
M = [ .5, .5, 0; 0, .5, .5; .5, 0, .5 ];        % barycentric coordinates of midpoints
Bar = repmat( M, NoT, 1 );                        % in accordance with SI: m1 m2 m3 m1 m2 m3 m1 m2 m3 ...


%% 1st pass
%
%      1                   1
%     / \                 /  \
%    /   \               1'   3'
%   /     \             /      \
%  2 ----- 3    ->     2-- 2' --3
%   \     /             \      /
%    \   /               \    /
%     \ /                 \  /
%      4                   4
% 
Car = baryToCart( TR, SI, Bar );    % cartesian coordinates of midpoints
new_Ver = [ Ver; Car ];           % add to new vertices, note that there are many duplicated points in C


%% 2nd pass 
% new_Tri = zeros( NoT*4, 3 );

%      1                   1
%     / \                 /  \
%    /   \               1'-- 3'
%   /     \             / \  / \
%  2 ----- 3    ->     2-- 2' --3
%   \     /             \      /
%    \   /               \    /
%     \ /                 \  /
%      4                   4
% 
% organize new triangles with above pattern
% see the for loop below:
% 1' is "n+1", 1 is Tri(i,1)
% 2' is "n+2", 2 is Tri(i,2)
% 3' is "n+3", 3 is Tri(i,3)

% for loop
% m = 0;
% n = NoV_old;
% for i = 1:NoT
%     new_Tri(m+1,1:3) = [ Tri(i,1), n+1, n+3 ];
%     new_Tri(m+2,1:3) = [ n+1, Tri(i,2), n+2 ];
%     new_Tri(m+3,1:3) = [ n+3, n+2, Tri(i,3) ];
%     new_Tri(m+4,1:3) = [ n+1, n+2, n+3 ];
%     m = m + 4;
%     n = n + 3;
% end

% vectorized for loop
n = ( NoV_old : 3 : (NoV_old+(NoT-1)*3) )';
nT1 = [ Tri(:,1), n+1, n+3 ];
nT2 = [ n+1, Tri(:,2), n+2 ];
nT3 = [ n+3, n+2, Tri(:,3) ];
nT4 = [ n+1, n+2, n+3 ];
new_Tri = [ nT1; nT2; nT3; nT4 ];

if size(varargin,1) == 0
    [ ~, ia, ic ] =  unique( new_Ver, 'rows' );    % patch slim: ignore duplicated vertices
    new_Ver = new_Ver( ia, : );
    new_Tri = ic(new_Tri);
    sync.ia = ia;
    sync.ic = ic;
else
    ia = varargin{1}.ia;
    ic = varargin{1}.ic;
    new_Ver = new_Ver( ia, : );
    new_Tri = ic(new_Tri);
end

new_TR = TriRep( new_Tri, new_Ver );


%% refine (method 1)
NoV = size(new_Ver,1);
VB = zeros(NoV,1);      % if VB is all zero, this is a closed surface
B = freeBoundary(new_TR);
if size( B ~= 0, 1 )	% open surface
    VB(B(:,1),1) = 1;
end
E = edges(new_TR);          % edge list (unidirection), this is to make unidirection edge list
tmp = [ E(:,2),E(:,1) ];    % become bidirection edge list for example, we'll have both 1->2 and 2->1 
E = [ E; tmp ];             % (same edge, different direction). We need both directions in order to find adjacent vertices
E = sortrows(E,[1,2]);      % sort this list so that adjacent vertices are easy to index
[ ~, Efo, ~ ] = unique( E(:,1), 'first' );    % first occurrence
[ ~, Elo, ~ ] = unique( E(:,1), 'last' );     % last occurrence
V = zeros(NoV,3);           % new vertices (smoothed vertices)
ALPHA = zeros(max(diff(Efo)),1);   % pre-calculate alpha(n) to increase speed
for i = 1:length(ALPHA)
    ALPHA(i) = alpha(i);
end
% plotMesh(new_Ver, new_Tri,'k');
% k = 1;
for i = 1 : NoV
    if VB(i,1) == 1
        if ia(i) <= NoV_old
            nei = B_old( B_old(:,1) == ia(i), 2 );
            V(i,:) = 6/8*Ver(ia(i),:) + 1/8*sum(Ver(nei,:),1);        % boundary vertices
%             k = k + 1;
%             h1 = plot(Ver(ia(i),1),Ver(ia(i),2),'sr');
%             h2 = plot(Ver(nei,1),Ver(nei,2),'or');
%             delete(h1);
%             delete(h2);
        else
            V(i,:) = new_Ver(i,:);
        end
    else
        nei = E( Efo(i):Elo(i), 2 );
        n = length(nei);
%         V(i,:) = ( alpha(n)*new_Ver(i,:) + sum(new_Ver(nei,:)) ) / ( alpha(n) + n );	% inner vertices
        V(i,:) = ( ALPHA(n)*new_Ver(i,:) + sum(new_Ver(nei,:),1) ) / ( ALPHA(n) + n );	% inner vertices
    end
end

new_Ver = V;


%% refine (method 2)
% NoV = size(new_Ver,1);
% VB = zeros(NoV,1);      % if VB is all zero, this is a closed surface
% B = freeBoundary(new_TR);
% if size( B ~= 0, 1 )   % open surface
%     VB(B(:,1),1) = 1;
%     tmp = [ B(:,2),B(:,1) ];
%     B = [ B; tmp ];
%     B = sortrows(B,[1,2]);
%     [ ~, Bfo, ~ ] = unique( B(:,1), 'first' );    % first occurrence
%     [ ~, Blo, ~ ] = unique( B(:,1), 'last' );     % last occurrence
% end
% E = edges(new_TR);          % edge list (unidirection), this is to make unidirection edge list
% tmp = [ E(:,2),E(:,1) ];    % become bidirection edge list for example, we'll have both 1->2 and 2->1 
% E = [ E; tmp ];             % (same edge, different direction). We need both directions in order to find adjacent vertices
% E = sortrows(E,[1,2]);      % sort this list so that adjacent vertices are easy to index
% [ ~, Efo, ~ ] = unique( E(:,1), 'first' );    % first occurrence
% [ ~, Elo, ~ ] = unique( E(:,1), 'last' );     % last occurrence
% V = zeros(NoV,3);           % new vertices (smoothed vertices)
% ALPHA = zeros(max(diff(Efo)),1);   % pre-calculate alpha(n) to increase speed
% for i = 1:length(ALPHA)
%     ALPHA(i) = alpha(i);
% end
% % plotMesh(new_Ver, new_Tri,'k');
% k = 1;
% for i = 1 : NoV
%     if VB(i,1) == 1
%         nei = B( Bfo(k):Blo(k), 2 );
%         V(i,:) = 6/8*new_Ver(i,:) + 1/8*sum(new_Ver(nei,:),1);        % boundary vertices
%         k = k + 1;
%     else
%         nei = E( Efo(i):Elo(i), 2 );
%         n = length(nei);
% %         V(i,:) = ( alpha(n)*new_Ver(i,:) + sum(new_Ver(nei,:)) ) / ( alpha(n) + n );	% inner vertices
%         V(i,:) = ( ALPHA(n)*new_Ver(i,:) + sum(new_Ver(nei,:),1) ) / ( ALPHA(n) + n );	% inner vertices
%     end
% end
% 
% new_Ver = V;


function [ a ] = alpha(n)
a = n*(1-beta(n))/beta(n);

function [ b ] = beta(n)
b = 5/4 - ((3+2*cos(2*pi/n))^2)/32;
% b = 5/8 - (3/8 + 1/4*cos(2*pi/n))^2;

