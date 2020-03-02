% --------------------------- Descriptions --------------------------------
% mesh subdivision based on Loop's algorithm
% 
% File name:    myLoopSubdivision2.m
% Date created: 03/21/2013
% Last revise:  04/10/2013
% -------------------------------------------------------------------------

function [ new_Ver, new_Tri ] = myLoopSubdivision2( Ver, Tri, varargin )


Inner_Smooth = 0;
Outer_Smooth = 0;
switch max( size( varargin ) )
    case 0
        Inner_Smooth = 1;
        Outer_Smooth = 1;
    case 1
        Inner_Smooth = 1;
    otherwise
        % do nothing
end


%% init
NoV_old = size(Ver,1);
% NoT_old = size(Tri,1);
TR = TriRep( Tri, Ver );
E = edges(TR);
B_old = freeBoundary( TR );
if size( B_old ~= 0, 1 )   % open surface
    tmp = [ B_old(:,2),B_old(:,1) ];
    B_old = [ B_old; tmp ];
    B_old = sortrows(B_old,[1,2]);
    [ ~, Bfo_old, ~ ] = unique( B_old(:,1), 'first' );    % first occurrence
    [ ~, Blo_old, ~ ] = unique( B_old(:,1), 'last' );     % last occurrence
end
NoE = size(E,1);


%% 1st pass
Head = Ver( E(:,1), 1:3 );
Tail = Ver( E(:,2), 1:3 );
M = ( Head + Tail ) / 2;

EMAT = sparse( [E(:,1),E(:,2)], [E(:,2),E(:,1)], [(1:NoE)';(1:NoE)'], NoV_old, NoV_old, 2*NoE );
% EMAT = full(EMAT);
new_Ver = [ Ver; M ];
% return;


%% 2nd pass 
Ind1 = sub2ind( [NoV_old, NoV_old], Tri(:,1), Tri(:,2) );
Ind2 = sub2ind( [NoV_old, NoV_old], Tri(:,2), Tri(:,3) );
Ind3 = sub2ind( [NoV_old, NoV_old], Tri(:,3), Tri(:,1) );

newP1 = EMAT(Ind1);
newP2 = EMAT(Ind2);
newP3 = EMAT(Ind3);

nT1 = [ Tri(:,1), newP1+NoV_old, newP3+NoV_old ];
nT2 = [ newP1+NoV_old, Tri(:,2), newP2+NoV_old ];
nT3 = [ newP3+NoV_old, newP2+NoV_old, Tri(:,3) ];
nT4 = [ newP1+NoV_old, newP2+NoV_old, newP3+NoV_old ];
new_Tri = [ nT1; nT2; nT3; nT4 ];

% [ new_Ver, new_Tri ] = PatchSlim( new_Ver, new_Tri );
new_TR = TriRep( new_Tri, new_Ver );


%% refine
if Inner_Smooth
    % NoV = size(new_Ver,1);
    % VB = zeros(NoV,1);      % if VB is all zero, this is a closed surface
    B = freeBoundary(new_TR);
    % if size( B ~= 0, 1 )	% open surface
    %     VB(B(:,1),1) = 1;
    % end
    E = edges(new_TR);          % edge list (unidirection), this is to make unidirection edge list
    tmp = [ E(:,2),E(:,1) ];    % become bidirection edge list for example, we'll have both 1->2 and 2->1 
    E = [ E; tmp ];             % (same edge, different direction). We need both directions in order to find adjacent vertices
    E = sortrows(E,[1,2]);      % sort this list so that adjacent vertices are easy to index
    [ ~, Efo, ~ ] = unique( E(:,1), 'first' );    % first occurrence
    [ ~, Elo, ~ ] = unique( E(:,1), 'last' );     % last occurrence
    % V = zeros(NoV,3);           % new vertices (smoothed vertices)
    ALPHA = zeros(max(diff(Efo)),1);   % pre-calculate alpha(n) to increase speed
    for i = 1:length(ALPHA)
        ALPHA(i) = alpha(i);
    end

    % plotMesh(new_Ver, new_Tri,'k');
    % k = 1;
    % for i = 1 : NoV
    %     if VB(i,1) == 1
    %         if i <= NoV_old
    % %             nei = B_old( B_old(:,1) == i, 2 );
    %             nei = B_old( Bfo_old(k):Blo_old(k), 2 );
    %             V(i,:) = 6/8*Ver(i,:) + 1/8*sum(Ver(nei,:),1);        % boundary vertices
    %             k = k + 1;
    % %             h1 = plot(Ver(ia(i),1),Ver(ia(i),2),'sr');
    % %             h2 = plot(Ver(nei,1),Ver(nei,2),'or');
    % %             delete(h1);
    % %             delete(h2);
    %         else
    %             V(i,:) = new_Ver(i,:);
    %         end
    %     else
    % %         n = size( Tv{i}', 1) ;
    % %         nei = myUnique( sort( reshape( new_Tri( Tv{i}',:), n*3, 1 ) ) );
    % %         nei = nei(nei~=i);
    %         nei = E( Efo(i):Elo(i), 2 );
    %         n = length(nei);
    % %         V(i,:) = ( alpha(n)*new_Ver(i,:) + sum(new_Ver(nei,:)) ) / ( alpha(n) + n );	% inner vertices
    %         V(i,:) = ( ALPHA(n)*new_Ver(i,:) + sum(new_Ver(nei,:),1) ) / ( ALPHA(n) + n );	% inner vertices
    %     end
    % end

    [ SUM_Ver ] = regMask( new_Ver(E(:,2),:), Efo, Elo );
    n = Elo - Efo + 1;
    V = ( [ALPHA(n),ALPHA(n),ALPHA(n)].*new_Ver + SUM_Ver ) ./ ( [ALPHA(n) + n,ALPHA(n) + n,ALPHA(n) + n] );    % inner vertices

    if size( B ~= 0, 1 )	% open surface
        [ SUM_VerB ] = regMask( Ver(B_old(:,2),:), Bfo_old, Blo_old );
        V(B,:) = new_Ver(B,:);
        if Outer_Smooth
            V(unique(B_old(:,1)),:) = 6/8*Ver( unique(B_old(:,1)), : ) + 1/8*SUM_VerB;    % boundary vertices
        end
    end

    new_Ver = V;
end


function [ a ] = alpha(n)
a = n*(1-beta(n))/beta(n);

function [ b ] = beta(n)
b = 5/4 - ((3+2*cos(2*pi/n))^2)/32;


