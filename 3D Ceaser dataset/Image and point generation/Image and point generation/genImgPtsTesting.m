%genImg / genPts

P_dense = importdata('meanShape_dense.mat');
Data = importdata('Data.mat');
FR_CTRL_PT = Data.FR_CTRL_PT;
FR_LIST = Data.FR_LIST;
FR_IDX = Data.FR_IDX;

idx_boundary = boundary(P_dense(:,1), P_dense(:,3), 1);
P_FR_boundary = [P_dense(idx_boundary,1), P_dense(idx_boundary,3)];
[~, i_max] = max(P_FR_boundary(:,2));
idx_boundary = [idx_boundary(i_max:end);idx_boundary(1:i_max-1)];

P_FR_boundary = [P_dense(idx_boundary,1), P_dense(idx_boundary,3)];
% P_FR_boundary = [P_FR_boundary(i_max:end,:);P_FR_boundary(1:i_max-1,:)];
[idx_main, ~] = NNSearch2DFEX(P_FR_boundary, [FR_CTRL_PT(:,1), FR_CTRL_PT(:,3)]);

figure,
fill(P_FR_boundary(:,1), P_FR_boundary(:,2), 'k','edgeColor', 'none')
axis equal
axis off
F = getframe(gca);
[X, Map] = frame2im(F);

%% Frontal image/pts generation
[rows, cols, ~] = size(X);
height = max(P_FR_boundary(:,2)) - min(P_FR_boundary(:,2));
P_FR_boundary(:,2) = P_FR_boundary(:,2) - min(P_FR_boundary(:,2));
P_FR_boundary = P_FR_boundary /height*rows;
P_FR_boundary(:,1) = P_FR_boundary(:,1) + (cols - 1)/2;
P_FR_boundary = round(P_FR_boundary);
P_FR_boundary(:,2) = rows - P_FR_boundary(:,2);
LM_FR_major = P_FR_boundary(idx_main, 1:2);

[~, ia, ~] = unique(round(P_FR_boundary),'rows','stable');
P_FR_boundary = P_FR_boundary(ia,:);
[idx_main, ~] = NNSearch2DFEX(P_FR_boundary, LM_FR_major);

fprintf("\n Generating frontal images and point sets.")
for i = 101:150
    
    filePath = "C:\Users\Semanti Basu\Documents\OneDrive_2020-02-19\3D Ceaser dataset\newModel\newModel\";
    fileName = "csr"+num2str(i+4000)+"b.mat";
    fprintf("\n ... Computing model " + fileName);
    pos = 1;
    [matName, PM_CTRL_PT_resized, I] = Fig2Img(filePath,fileName, idx_boundary, idx_main, ia, FR_IDX, pos);

    if pos == 1
        pos_i = "f_";
    else
        pos_i = "s_";
    end

    if ~exist('ceasar_mat_test') 
        mkdir ceasar_mat_test
    end
    save_file = "ceasar_mat_test/test_"+pos_i+fileName;

    points = PM_CTRL_PT_resized;

    save(save_file, 'points');
    imwrite(I, save_file+".jpg");
end

%% Profile image/pts generation

fprintf("\n Generating profile images and point sets.")
SD_CTRL_PT = Data.SD_CTRL_PT;
SD_LIST = Data.SD_LIST;
SD_IDX = Data.SD_IDX;

idx_boundary = boundary(P_dense(:,2), P_dense(:,3), 1);
P_SD_boundary = [P_dense(idx_boundary,2), P_dense(idx_boundary,3)];
[~, i_max] = max(P_SD_boundary(:,2));
idx_boundary = [idx_boundary(i_max:end);idx_boundary(1:i_max-1)];

P_SD_boundary = [P_dense(idx_boundary,2), P_dense(idx_boundary,3)];

[idx_main, ~] = NNSearch2DFEX(P_SD_boundary, [SD_CTRL_PT(:,2), SD_CTRL_PT(:,3)]);
figure,
fill(P_SD_boundary(:,1), P_SD_boundary(:,2), 'k','edgeColor', 'none')
axis equal
axis off
F = getframe(gca);
[X, Map] = frame2im(F);
figure,
imshow(X)

[rows, cols, ~] = size(X);
height = max(P_SD_boundary(:,2)) - min(P_SD_boundary(:,2));
P_SD_boundary(:,2) = P_SD_boundary(:,2) - min(P_SD_boundary(:,2));
P_SD_boundary = P_SD_boundary /height*rows;
P_SD_boundary(:,1) = P_SD_boundary(:,1) + (cols + 1)/2;
P_SD_boundary = round(P_SD_boundary);
P_SD_boundary(:,2) = rows - P_SD_boundary(:,2);
LM_SD_major = P_SD_boundary(idx_main, 1:2);

[~, ia, ~] = unique(round(P_SD_boundary),'rows','stable');
P_SD_boundary = P_SD_boundary(ia,:);
[idx_main, ~] = NNSearch2DFEX(P_SD_boundary, LM_SD_major);

for i = 101:150
    
    filePath = "C:\Users\Semanti Basu\Documents\OneDrive_2020-02-19\3D Ceaser dataset\newModel\newModel\";
    fileName = "csr"+num2str(i+4000)+"b.mat";
    fprintf("\n ... Computing model " + fileName);
    pos = 2;
    [matName, PM_CTRL_PT_resized, I] = Fig2Img(filePath,fileName, idx_boundary, idx_main, ia, FR_IDX, pos);

    if pos == 1
        pos_i = "f_";
    else
        pos_i = "s_";
    end

    if ~exist('ceasar_mat_test') 
        mkdir ceasar_mat_test
    end
    save_file = "ceasar_mat_test/test_"+pos_i+fileName;

    points = PM_CTRL_PT_resized;

    save(save_file, 'points');
    imwrite(I, save_file+".jpg");
end
