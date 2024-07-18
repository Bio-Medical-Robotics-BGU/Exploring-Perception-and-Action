%% This code analyses the network predictions - creates the real and predicted curves (for 
% example, those shown in Fig. 6), and computes the regression of the
% predicted augmentation against the real one, and compares the real and
% predicted curves in each condition. (Fig. 6, 8, 9, and 10).
% This also contains the code for the percentage of erros (model or participant
% per comparison stiffness level (Fig. 7(a-f))).
clear all; clc; 
cd 'D:\OneDrive\PerceptionActionReview'
addpath('D:\OneDrive\PerceptionActionReview')
addpath(genpath('D:\OneDrive\PerceptionActionReview'))

%replace following directory with the location of the saved network
%predictions
saved_path = 'D:\OneDrive\PerceptionActionReview\saved_predictions';

%replace following directory with the location of the saved preprocessed
%data
data_path = 'D:\OneDrive\PerceptionActionReview\Preprocessed';

% replace following directory with the location of the codes
project_path = 'D:\OneDrive\PerceptionActionReview\Code\Analyses';

%load the splits
folds = load('TestIndsSplit_AllParticipants.mat');

%getting the field names of the struct
fns = fieldnames(folds);

%getting the model type from user
% model_name = input('What network would you like (e.g., which signals or parts of the network are included?) \n' ,'s');
model_name = 'LogisticRegression';
runs = input('Which predictions would you like to use? \n 1. Model \n 2. Network \n');

if runs == 2
    Run = input('Which number run would you like? \n', 's');
end

%for PSE and JND errors, and regression
RealPSEs = [];
PredPSEs = [];
RealJNDs = [];
PredJNDS = [];

%for accuracy
Accuracies = zeros(10, 1);

%for percentage of errors per comparison stiffness level
AllLabels = [];
AllPredictions = [];
AllTdgains = [];
AllKComps = [];

%for finding the two large error participants
count = 0;
%% Creating psychometric curves, and getting PSE and JND values
% The psychometric curves in Fig. 6(a-d) and Fig. 8(a) were created here.

for f = 1:10 %run over the 10 folds
    %getting the participants in this fold
    participants = folds.(fns{f});
    
    
    FoldLabels = [];
    FoldPreds = [];
    %running over the participants of the fold
    for p = 1:length(participants)

        count = count + 1;

        cd(data_path)
        Labels = load(['Labels_SN', num2str(participants(p)), '.mat']);
        Labels = Labels.AllPlabels(1:192);
        
        KComps = load(['Kcomps_SN', num2str(participants(p)), '.mat']);
        KComps = KComps.AllKcomps(1:192);
        
        TdGains = load(['TdGains_SN', num2str(participants(p)), '.mat']);
        TdGains = TdGains.AllTdGains(1:192);

        %getting predictions
        cd(saved_path)
        if runs == 2
            preds = load(['Preds_SN', num2str(participants(p)), '_', model_name, '_Run', Run, '.mat']);
        else
            preds = load(['Preds_SN', num2str(participants(p)), '_', model_name, '.mat']);
        end
        preds = preds.Preds;

        cd(project_path)
        [pses, pred_pses, jnds, pred_jnds] = PsychometricMats(KComps, TdGains, Labels, preds);

        RealPSEs = [RealPSEs; pses];
        PredPSEs = [PredPSEs; pred_pses];
        RealJNDs = [RealJNDs; jnds];
        PredJNDS = [PredJNDS; pred_jnds];

        AllLabels = [AllLabels; Labels];
        AllPredictions = [AllPredictions; preds];
        AllTdgains = [AllTdgains; TdGains];
        AllKComps = [AllKComps; KComps];

        FoldLabels = [FoldLabels; Labels];
        FoldPreds = [FoldPreds; preds];

    end %end of participants loop
    Accuracies(f) = (length(find(FoldLabels == FoldPreds))) / length(FoldLabels);
end %end of folds loop

%% Metrics
%% Accuracy
MeanRunAcc = mean(Accuracies);
disp(MeanRunAcc)
%% PSE Error
PSE_ForceErrors = RealPSEs(:, 1) - PredPSEs(:, 1);
PSE_StretchErrors = RealPSEs(:, 2) - PredPSEs(:, 2);

mes = ['PSE average errors are: \n Force condition ', num2str(mean(PSE_ForceErrors))...
    '\n Stretch condition ', num2str(mean(PSE_StretchErrors))];

fprintf(mes)

% Plotting Fig. 6(g), 8(c) and 10(b)
hold on
b = bar([0.5, 1, 2, 2.5], [mean(RealPSEs(:, 1)), mean(PredPSEs(:, 1)), mean(RealPSEs(:, 2)), mean(PredPSEs(:, 2))]);
b.FaceColor = 'flat';
b.CData(1,:) = [113 198 255]./255;
b.CData(2,:) = [255 159 159]./255;

b.CData(3,:) = [0 121 204]./255;
b.CData(4,:) = [192 0 0]./255;

errorbar(0.5, mean(RealPSEs(:, 1)), std(RealPSEs(:, 1))/2, 'color', 'k', 'linewidth', 0.5)
errorbar(1, mean(PredPSEs(:, 1)), std(PredPSEs(:, 1))/2, 'color', 'k', 'linewidth', 0.5)
errorbar(2, mean(RealPSEs(:, 2)), std(RealPSEs(:, 2))/2, 'color', 'k', 'linewidth', 0.5)
errorbar(2.5, mean(PredPSEs(:, 2)), std(PredPSEs(:, 2))/2, 'color', 'k', 'linewidth', 0.5)

set(gca, 'fontname', 'Times New Roman', 'fontsize', 14)

set(gca, 'Xtick', [0.5, 1, 2, 2.5], 'XtickLabels',...
    {'Real Force', 'Predicted Force', 'Real Stretch', 'Predicted Stretch'}, 'fontname', 'Times New Roman', 'fontsize', 14)
ylabel('Average PSE [N/m]', 'fontname', 'Times New Roman', 'fontsize', 16)

fig = gcf;
fig.Position = [200 200 440 415];

box on

%% JND Error
clc;
JND_ForceErrors = RealJNDs(:, 1) - PredJNDS(:, 1);
JND_StretchErrors = RealJNDs(:, 2) - PredJNDS(:, 2);

mes2 = ['JND average errors are: \n Force condition ', num2str(mean(JND_ForceErrors))...
    '\n Stretch condition ', num2str(mean(JND_StretchErrors))];

fprintf(mes2)

disp(mean(PredJNDS(:, 1)))
disp(mean(PredJNDS(:, 2)))

%% Regresssion
%% Regresssion
% Creates the plots for Fig. 6(f), 8(b), 9(a-b), 9(d), 10(a)
AllReals = RealPSEs(:, 2) - RealPSEs(:, 1);
AllPreds = PredPSEs(:, 2) - PredPSEs(:, 1);


% Predicted PSE vs PSE 
[sorted_reals1, sorted_inds1] = sort(AllReals);
sorted_preds1 = AllPreds(sorted_inds1);

%regressing
[sim1, ~, ~, ~, stats1] = regress(AllPreds, [ones(length(AllReals),1) AllReals]);
y_fit1 = sim1(1) + sim1(2)*sorted_reals1;


figure
hold on
plot(sorted_reals1, sorted_preds1, '*', 'linewidth', 0.5, 'markersize', 6, 'color', [0.5 0.5 0.5])
p2 = plot([0; sorted_reals1], [sim1(1); y_fit1], 'linewidth', 1.5, 'color', 'k');

plot([-25, 85], [-25, 85], 'linewidth', 0.5, 'color', [0.7 0.7 0.7], 'linestyle', '--', 'HandleVisibility', 'Off');

set(gca, 'fontname', 'Times New Roman', 'fontsize', 12)
xlabel('Real \DeltaPSE', 'fontname', 'Times New Roman', 'fontsize', 14)
ylabel('Predicted \DeltaPSE', 'fontname', 'Times New Roman', 'fontsize', 14)

l = legend(p2,['Predicted \DeltaPSE = ', num2str(round(sim1(1)*100)/100), ' + ', num2str(round(sim1(2)*100)/100), '\cdotReal \DeltaPSE'],...
    'Location', 'northwest');
l.FontSize = 12;

set(gca, 'fontname', 'Times New Roman', 'fontsize', 12)

% ylim([-25 85])

xlim(get(gca, 'Ylim'))
axis square

plot([0, 0], get(gca, 'Ylim'), 'k--', 'HandleVisibility','off')
plot(get(gca, 'Xlim'), [0, 0], 'k--', 'HandleVisibility','off')
box on

disp(stats1)
%% To plot the regressions of no position, velocity and acceleration on the same graph (Fig. 9(c)):
%Run code above 3 times and save the predicted values using the following
%lines:
%for saving in each run:
%No position:
AllRealsPos = RealPSEs(:, 2) - RealPSEs(:, 1);
AllPredsPos = PredPSEs(:, 2) - PredPSEs(:, 1);

%No velocity:
AllRealsVel = RealPSEs(:, 2) - RealPSEs(:, 1);
AllPredsVel = PredPSEs(:, 2) - PredPSEs(:, 1);

%No acceleration:
AllRealsAcc = RealPSEs(:, 2) - RealPSEs(:, 1);
AllPredsAcc = PredPSEs(:, 2) - PredPSEs(:, 1);

%create plot
[sim11, bint, r, rint, stats] = regress(AllPredsPos, [ones(length(AllRealsPos),1) AllRealsPos]);
[sorted_reals_pos1, sorted_inds_pos] = sort(AllRealsPos);
y_fit_pos1 = sim11(1) + sim11(2)*sorted_reals_pos1;

%no vel
[sorted_reals_vel2, sorted_inds_vel] = sort(AllRealsVel);
[sim22, bint, r, rint, stats] = regress(AllPredsVel, [ones(length(AllRealsVel),1) AllRealsVel]);
y_fit_vel2 = sim22(1) + sim22(2)*sorted_reals_vel2;

% No acc
[sorted_reals_acc3, sorted_inds_acc] = sort(AllRealsAcc);
[sim33, bint, r, rint, stats] = regress(AllPredsAcc, [ones(length(AllRealsAcc),1) AllRealsAcc]);
y_fit_acc3 = sim33(1) + sim33(2)*sorted_reals_acc3;

figure
p1 = plot([0; sorted_reals_pos1], [sim11(1); y_fit_pos1], 'linewidth', 1.5, 'color', [171 183 201]./255);
hold on

p2 = plot([0; sorted_reals_vel2], [sim22(1); y_fit_vel2], 'linewidth', 1.5, 'color', [94 117 148]./255);

p3 = plot([0; sorted_reals_acc3], [sim33(1); y_fit_acc3], 'linewidth', 1.5, 'color', [34 42 53]./255);

p4 = plot([-25 85], [-25 85], 'linewidth', 0.5, 'color', [0.7 0.7 0.7], 'linestyle', '--', 'HandleVisibility', 'Off');

set(gca, 'fontname', 'Times New Roman', 'fontsize', 12)
xlabel('Real Perceptual Augmentation [N/m]', 'fontname', 'Times New Roman', 'fontsize', 14)
ylabel('Predicted Perceptual Augmentation [N/m]', 'fontname', 'Times New Roman', 'fontsize', 14)
l = legend(['{\bfNo Position:\rm} Predicted \DeltaPSE = ', num2str(round(sim11(1)*100)/100), ' + ', num2str(round(sim11(2)*100)/100), '\cdotReal \DeltaPSE'], ...
    ['{\bfNo Velocity:\rm} Predicted \DeltaPSE = ', num2str(round(sim22(1)*100)/100), ' + ', num2str(round(sim22(2)*100)/100), '\cdotReal \DeltaPSE'],...
    ['{\bfNo Acceleration:\rm} Predicted \DeltaPSE = ', num2str(round(sim33(1)*100)/100), ' + ', num2str(round(sim33(2)*100)/100), '\cdotReal \DeltaPSE']);
l.FontSize = 8.4;
ylim([-25 85])
xlim([-25 85])

axis square

plot([0, 0], get(gca, 'Ylim'), 'k--', 'HandleVisibility','off')
plot(get(gca, 'Xlim'), [0, 0], 'k--', 'HandleVisibility','off')

%% Percentage of errors made by model per comparison stiffness level (Fig. 7(a-c))
UniqueK = unique(AllKComps);

%Force
TotalCounter1 = zeros(size(UniqueK));
ErrorCounter1 = zeros(size(UniqueK));
CorrectCounter1 = zeros(size(UniqueK));

AllPreds1 = AllPredictions(find(AllTdgains == 0));
AllLabels1 = AllLabels(find(AllTdgains == 0));
Kcomps1 = AllKComps(find(AllTdgains == 0));

for i = 1:length(Kcomps1)
    ind = find(UniqueK == Kcomps1(i));
    TotalCounter1(ind) = TotalCounter1(ind) + 1;
    
    if AllPreds1(i) == AllLabels1(i) %correct
        CorrectCounter1(ind) = CorrectCounter1(ind) + 1;
    else %mistake
        ErrorCounter1(ind) = ErrorCounter1(ind) + 1;
    end

end

% Skin Stretch
UniqueK = unique(AllKComps);
TotalCounter2 = zeros(size(UniqueK));
ErrorCounter2 = zeros(size(UniqueK));
CorrectCounter2 = zeros(size(UniqueK));

AllPreds2 = AllPredictions(find(AllTdgains == 80));
AllLabels2 = AllLabels(find(AllTdgains == 80));
Kcomps2 = AllKComps(find(AllTdgains == 80));

for i = 1:length(Kcomps2)
    ind = find(UniqueK == Kcomps2(i));
    TotalCounter2(ind) = TotalCounter2(ind) + 1;
    
    if AllPreds2(i) == AllLabels2(i) %correct
        CorrectCounter2(ind) = CorrectCounter2(ind) + 1;
    else %mistake
        ErrorCounter2(ind) = ErrorCounter2(ind) + 1;
    end

end

% All trials together
ErrorCounter = ErrorCounter1 + ErrorCounter2;
TotalCounter = TotalCounter1 + TotalCounter2;
CorrectCounter = CorrectCounter1 + CorrectCounter2;

% plotting
%All trials (Fig. 7(d))
figure
bar(UniqueK, 100*(ErrorCounter./TotalCounter), 'facecolor', [0.7 0.7 0.7])
xlabel('Comparison Stiffness Level [N/m]', 'fontsize', 14)
ylabel('Percentage of Errors [%]', 'fontsize', 14)
set(gca, 'fontname', 'Times New Roman', 'fontsize', 14)
ylim([0, 80])

%Force trials (Fig. 7(e))
figure
bar(UniqueK, 100*(ErrorCounter1./TotalCounter1), 'facecolor', [255 159 159]./255)
xlabel('Comparison Stiffness Level [N/m]', 'fontsize', 14)
ylabel('Percentage of Errors [%]', 'fontsize', 14)
set(gca, 'fontname', 'Times New Roman', 'fontsize', 14)
ylim([0, 80])

%Skin stretch trials (Fig. 7(f))
figure
bar(UniqueK, 100*(ErrorCounter2./TotalCounter2), 'facecolor', [192 0 0]./255)
xlabel('Comparison Stiffness Level [N/m]', 'fontsize', 14)
ylabel('Percentage of Errors [%]', 'fontsize', 14)
set(gca, 'fontname', 'Times New Roman', 'fontsize', 14)
ylim([0, 80])


