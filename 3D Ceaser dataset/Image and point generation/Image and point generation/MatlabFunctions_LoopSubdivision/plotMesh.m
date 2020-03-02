function plotMesh(vertices, faces, color )
    hold on;
    
    mv = size( vertices, 1 );
    nv = size( vertices, 2 );
    mf = size( faces, 1 );
    nf = size( faces, 2 );
    if mf < nf
        faces = faces';
    end
    if mv < nv
        vertices = vertices';
        nv = size( vertices, 2 );
    end
    if nv == 3
%         trimesh( faces, vertices(:,1), vertices(:,2), vertices(:,3), 'EdgeColor', color );
        patch( 'Faces', faces(:,1:3), 'Vertices', vertices, ...
               'FaceColor', 'none', 'FaceAlpha', 1 ,...
               'EdgeColor', color ); hold on
%         axis tight;
        axis equal; 
%         axis off;
        view(3);
    else
%         trimesh( faces, vertices(:,1), vertices(:,2), 'Color', color );
        patch( 'Faces', faces(:,1:3), 'Vertices', [ vertices, zeros(mv,1) ], ...
               'FaceColor', 'none', 'FaceAlpha', 1,...
               'EdgeColor', color ); hold on
        axis equal; 
    end

end