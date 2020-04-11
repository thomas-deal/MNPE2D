clc
clear
close all
path(pathdef)
%% Path Settings
TEXPath = pwd;
cd(fullfile('..','..'))
MNPEPath = pwd;
DOCPath = fullfile(MNPEPath,'doc');
IMGPath = fullfile(DOCPath,'images');
%% Update Images
cd(IMGPath)
ExampleImages
close all
%% Update Documentation
texFile = 'MNPEintro';
cd(TEXPath)
latexCMD = ['pdflatex -shell-escape ' texFile];
system(latexCMD);
system(['bibtex ' texFile]);
system(latexCMD);
system(latexCMD);
if ~exist(DOCPath,'dir')
    mkdir(DOCPath)
end
copyfile([texFile '.pdf'],DOCPath)
open(fullfile(DOCPath,[texFile '.pdf']))
%% Clean Up
CleanupList = {'.aux' '.bbl' '.blg' '.lof' '.lot' '.out' '.pdf' '.toc'};
for i=1:length(CleanupList)
    try
        delete(['*' CleanupList{i}])
    catch me
    end
end