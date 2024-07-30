%% This code creates Fig. 4, which describes the metrics computed from the 
% psychometric curves
project_path = "C:\Users\hannako\Downloads\ExploringPerceptionandAction\ExploringPerceptionandAction";
addpath(genpath(fullfile(project_path, 'Code')));

%% Blue curves
%%
PsychometricData_F=85-[140:-10:30]';
PsychometricData_F(:,2)=[0; 0; 0; 0; 1; 2; 6; 6; 7; 8; 8; 8];
PsychometricData_F(:,3)=8;
PsychometricData_SS=85-[140:-10:30]';
PsychometricData_SS(:,2)=[0; 0; 0; 0; 0; 1; 0; 0; 3; 5; 4; 6];
PsychometricData_SS(:,3)=8;

[PSE1,PSE2]=PsychometricFitting_ForMetricDescription(PsychometricData_F,PsychometricData_SS, 'Real');
set(gca, 'fontname', 'Times New Roman')
set(gca, 'fontsize', 10)

fig=gcf;
fig.Position=[0, 0, 250, 210];

%% Red curves
%%
PsychometricData_F=85-[140:-10:30]';
PsychometricData_F(:,2)=[0; 0; 2; 5; 5; 5; 6; 6; 7; 8; 8; 8];
PsychometricData_F(:,3)=8;
PsychometricData_SS=85-[140:-10:30]';
PsychometricData_SS(:,2)=[0; 0; 0; 2; 2; 1; 3; 5; 7; 5; 6; 5];
PsychometricData_SS(:,3)=8;


[PSE1,PSE2]=PsychometricFitting_ForMetricDescription(PsychometricData_F,PsychometricData_SS, 'Predicted');
set(gca, 'fontname', 'Times New Roman')
set(gca, 'fontsize', 10)

fig=gcf;
fig.Position=[0, 0, 250, 210];
