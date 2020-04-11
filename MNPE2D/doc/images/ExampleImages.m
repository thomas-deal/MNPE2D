close all
%% Settings
RunDiaryFile = 'RunOutput.txt';
peout1DiaryFile = 'peout1Output.txt';
PromptFile = 'prompt.txt';
pressbin = 'press.bin';
apvrbin = 'apvr.bin';
apvzbin = 'apvz.bin';
%% Examples
% Flat
Example(1).Name = 'Flat';
Example(1).pefile = {pressbin};
Example(1).options = {'opt1=2;'};
Example(1).pos = [10 10 800 400];
Example(1).cax = [20 80];
Example(1).imgfile = {'FlatP.eps'};
% ASA Wedge
Example(2).Name = 'ASAWedge';
Example(2).pefile = {pressbin};
Example(2).options = {'opt1=2;'};
Example(2).pos = [10 10 800 400];
Example(2).cax = [20 80];
Example(2).imgfile = {'ASAWedgeP.eps'};
% Monopole
Example(3).Name = 'Monopole';
Example(3).pefile = {pressbin, apvzbin};
Example(3).options = {'opt1=2;', 'opt1=2;'};
Example(3).pos = [10 10 800 350];
Example(3).cax = [20 80];
Example(3).imgfile = {'MonopoleP.eps','MonopoleVZ.eps'};
% Dipole
Example(4).Name = 'Dipole';
Example(4).pefile = {pressbin, apvzbin};
Example(4).options = {'opt1=2;', 'opt1=2;'};
Example(4).pos = [10 10 800 350];
Example(4).cax = [20 80];
Example(4).imgfile = {'DipoleP.eps','DipoleVZ.eps'};
% Rough Surface
Example(5).Name = 'RoughSurface';
Example(5).pefile = {pressbin};
Example(5).options = {'opt1=2;'};
Example(5).pos = [10 10 800 400];
Example(5).cax = [20 80];
Example(5).imgfile = {'RoughSurfaceP.eps'};
% Travel Time
Example(6).Name = 'TravelTime';
Example(6).pefile = {pressbin};
Example(6).options = {'opt1=6;optt=1;optt1=1;rngout=1;'};
Example(6).pos = [10 10 800 400];
Example(6).cax = [60 140];
Example(6).imgfile = {'TravelTimeP.eps'};
%
Nexamples = length(Example);
%% Setup Paths
IMGPath = pwd;
cd(fullfile('..','..'))
MNPEPath = pwd;
BINPath = fullfile(MNPEPath,'bin');
addpath(BINPath)
EXPath = fullfile(MNPEPath,'examples');
RunDiary = fullfile(IMGPath,RunDiaryFile);
peout1Diary = fullfile(IMGPath,peout1DiaryFile);
%% Log Outputs for MNPE Example
cd(fullfile(EXPath,'Monopole'))
% Log executable outputs
if exist(RunDiary,'file')
    delete(RunDiary)
end
diary(RunDiary)
% Write temporary file to answer MMPE command line prompt for accuracy/efficiency
fid=fopen(PromptFile,'w');fprintf(fid,'1\n');fclose(fid);
system([fullfile(BINPath,'MNPE2D') ' < ' PromptFile]);
delete(PromptFile);
diary off
% Log header info
doPlots = true; %#ok<*NASGU>
scriptedinput = true;
pefile = pressbin;
opt1 = 2;
savedat = 'n';
if exist(peout1Diary,'file')
    delete(peout1Diary)
end
diary(peout1Diary)
peout1
diary off
%% Run Examples and Save Plots
doPlots = true;
scriptedinput = true;
savedat = 'n';
for iex = 1:Nexamples
    cd(fullfile(EXPath,Example(iex).Name))
    % Write temporary file to answer MNPE command line prompt for accuracy/efficiency
    fid=fopen(PromptFile,'w');fprintf(fid,'1\n');fclose(fid);
    system([fullfile(BINPath,'MNPE2D') ' < ' PromptFile]);
    delete(PromptFile);
    for jex = 1:length(Example(iex).pefile)
        pefile = Example(iex).pefile{jex};
        eval(Example(iex).options{jex});
        peout1
        peout2
        caxis(Example(iex).cax)
        set(gcf,'position',Example(iex).pos, ...
            'paperpositionmode','auto', ...
            'renderer','painters');
        try
            print(gcf,'-depsc2',fullfile(IMGPath,Example(iex).imgfile{jex}))
        catch me
            disp(['Unable to save ' fullfile(IMGPath,Example(iex).imgfile{jex}) ', check write status.'])
        end
    end
end
%% Finish
cd(IMGPath)