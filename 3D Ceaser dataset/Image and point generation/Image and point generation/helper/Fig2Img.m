function [matName, LM_FR_PM_major, I]=Fig2Img(filePath,fileName, idx_boundary, idx_main, ia, FR_IDX, pos)
    t = load('triConnection');
    P_Connection = t.t;
    
    [PM] = importGM(filePath+fileName);
    
    for i = 1:3
    [ PM , P_Connection ] = myLoopSubdivision2( PM, P_Connection, 1);
    end
    PM_dense = PM;
    
    if pos == 1
        PM_FR_boundary = [PM_dense(idx_boundary,1), PM_dense(idx_boundary,3)];
        pos = "f_";
    else
        PM_FR_boundary = [PM_dense(idx_boundary,2), PM_dense(idx_boundary,3)];
        pos = "s_";
    end
    
    
%     [idx, ~] = NNSearch2DFEX(PM_FR_boundary, [PM(FR_IDX,1), PM(FR_IDX,3)]);
    figure,
    fill(PM_FR_boundary(:,1), PM_FR_boundary(:,2), 'k','edgeColor', 'none')
    axis equal
    ylim([0, 2000])
    axis off
    F = getframe(gca);
    [X, Map] = frame2im(F);
    close
    [rows, cols, ~] = size(X);
    height = max(PM_FR_boundary(:,2)) - min(PM_FR_boundary(:,2));
    height = 2000;
    PM_FR_boundary(:,2) = PM_FR_boundary(:,2);
    PM_FR_boundary = PM_FR_boundary /height*rows;
    PM_FR_boundary(:,1) = PM_FR_boundary(:,1) + (cols + 1)/2;
    PM_FR_boundary = round(PM_FR_boundary);
    PM_FR_boundary(:,2) = rows - PM_FR_boundary(:,2);
    
    PM_FR_boundary = PM_FR_boundary(ia, :);
    LM_FR_PM_major = PM_FR_boundary(idx_main, 1:2);
    
%     LM_FR_PM_major(:,2) = rows - LM_FR_PM_major(:,2);
%     LM_FR_PM_major(:,1) = cols - LM_FR_PM_major(:,1);

    p.n = length(PM_FR_boundary);
    p.x = PM_FR_boundary(:,2)';
    p.y = PM_FR_boundary(:,1)';
    
    t = ones(1,length(PM_FR_boundary))*2;
    t(idx_main) = 0;
    p.t = t;
%     p.t = zeros(1,length(LM_FR_PM_major));
    p.I = rgb2gray(im2double(X));
    if exist('ceasar_mat') 
        save_file = 'ceasar_mat/train_'+pos+fileName;
    else
        mkdir ceasar_mat
        save_file = 'ceasar_mat/train_'+pos+fileName;
    end
    save(save_file, 'p')

    I=im2double(X);  
    matName=save_file;