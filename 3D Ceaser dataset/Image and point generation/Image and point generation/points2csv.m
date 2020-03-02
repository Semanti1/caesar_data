% c={'f1',reshape(points',1,112);'f2',[1,2,3,4]};
% T=cell2table(c);
% writetable(T,'testsimmy2.csv','WriteVariableNames',false);
clear;clc;
filePath = "C:\Users\Semanti Basu\Documents\OneDrive_2020-02-19\3D Ceaser dataset\Image and point generation\Image and point generation\ceasar_mat\";
% Files=dir(fullfile(filePath,'*.mat'));
filePattern = fullfile(filePath, '*.mat');
matFiles = dir(filePattern);
matData={};
for k = 1:100
% for k=1:3
  baseFileName = matFiles(k).name;
  fullFileName = fullfile(filePath, baseFileName);
  fprintf(1, 'Now reading %s\n', fullFileName);
  temp = load(fullFileName);
%   temp.points
  matData(end+1,:)={baseFileName+".jpg",reshape((temp.points)',1,112)};
%   matData{baseFileName}=reshape((temp.points)',1,112)
end
T=cell2table(matData);
writetable(T,'frontalpoints.csv','WriteVariableNames',false);
%%
matData={};
for k = 101:200
% for k=1:3
  baseFileName = matFiles(k).name;
  fullFileName = fullfile(filePath, baseFileName);
  fprintf(1, 'Now reading %s\n', fullFileName);
  temp = load(fullFileName);
%   temp.points
  matData(end+1,:)={baseFileName+".jpg",reshape((temp.points)',1,52)};
%   matData{baseFileName}=reshape((temp.points)',1,112)
end
T2=cell2table(matData);
writetable(T2,'sideviewpoints.csv','WriteVariableNames',false);