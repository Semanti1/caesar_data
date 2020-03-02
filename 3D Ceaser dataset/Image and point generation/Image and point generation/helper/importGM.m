function [GM] = importGM(filePath)
% filePath = "C:\Users\Chenxi\Documents\MATLAB\caesar\meanShape.mat";
genericModel = load(filePath);
genericModel = genericModel.points;

y_mean = mean(genericModel(:,2));
x_mean = mean(genericModel(:,1));
genericModel(:,2) = genericModel(:,2) - y_mean;
genericModel(:,1) = genericModel(:,1) - x_mean;

GM = genericModel;

% figure,
% plot3(GM(:,1), GM(:,2), GM(:,3),'k.', 'MarkerSize', 3)
% xlabel('(Width, x-axis, mm)')
% ylabel('(Depth, y-axis, mm)')
% zlabel('(Height, z-axis, mm)')
% axis equal
% title('Generic model')
% axis off

