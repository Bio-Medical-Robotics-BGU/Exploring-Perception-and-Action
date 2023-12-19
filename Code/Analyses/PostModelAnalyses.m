%% This code analyses the network predictions - creates the real and predicted curves (for 
% example, those shown in Fig. 6), and computes the regression of the
% predicted augmentation against the real one, and compares the real and
% predicted curves in each condition. 

%replace following directory with the location of the saved network
%predictions
saved_path = 'C:\Users\hanna\OneDrive\lab\PhD\python\PerceptionActionInSurface\Final\saved_predictions';

%replace following directory with the location of the saved preprocessed
%data
data_path = 'C:\Users\hanna\OneDrive\MATLAB\lab\PhD\Perception_and_GF_prediction\DL_perception\OnlyInSurface\ParticipantMatsAndVecs\final';

% replace following directory with the location of the codes
project_path = 'C:\Users\hanna\OneDrive\MATLAB\lab\PhD\Perception_and_GF_prediction\DL_perception\OnlyInSurface\final';

%load the splits
folds = load('TestIndsSplit.mat');

%getting the field names of the struct
fns = fieldnames(folds);

%getting the model type from user
model_name = input('What network would you like (e.g., which signals or parts of the network are included?) \n' ,'s');

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
LargeErrors = []; %to save indices of large error participants
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

        if (participants(p) == 10 || participants(p) == 17)
            LargeErrors = [LargeErrors; count];
        end
        cd(data_path)
        Labels = load(['Labels_SN', num2str(participants(p)), '.mat']);
        Labels = Labels.AllPlabels(1:192);
        
        KComps = load(['Kcomps_SN', num2str(participants(p)), '.mat']);
        KComps = KComps.AllKcomps(1:192);
        
        TdGains = load(['TdGains_SN', num2str(participants(p)), '.mat']);
        TdGains = TdGains.AllTdGains(1:192);

        %getting predictions
        cd(saved_path)
        preds = load(['Preds_SN', num2str(participants(p)), model_name, '.mat']);
        preds = preds.Preds;

        cd(project_path)
        [pses, pred_pses, jnds, pred_jnds] = PsychometricMats(KComps, TdGains, Labels, preds);

        RealPSEs = [RealPSEs; pses];
        PredPSEs = [PredPSEs; pses];
        RealJNDs = [RealJNDs; pses];
        PredJNDS = [PredJNDS; pses];

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

%% PSE Error
PSE_ForceErrors = RealPSEs(:, 1) - PredPSEs(:, 1);
PSE_StretchErrors = RealPSEs(:, 2) - PredPSEs(:, 2);

mes = ['PSE average errors and standard deviations are: \n Force condition ', num2str(mean(PSE_ForceErrors))...
    ', ' num2str(std(PSE_ForceErrors)), '\n Stretch condition ', num2str(mean(PSE_StretchErrors)),...
    ', ' num2str(std(PSE_StretchErrors))];

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
JND_ForceErrors = RealJNDs(:, 1) - PredJNDS(:, 1);
JND_StretchErrors = RealJNDs(:, 2) - PredJNDS(:, 2);

mes2 = ['JND average errors and standard deviations are: \n Force condition ', num2str(mean(JND_ForceErrors))...
    ', ' num2str(std(JND_ForceErrors)), '\n Stretch condition ', num2str(mean(JND_StretchErrors)),...
    ', ' num2str(std(JND_StretchErrors))];

fprintf(mes2)

%% Regresssion
% Creates the plots for Fig. 6(f), 8(b), 9(a-b), 9(d), 10(a)
AllReals = RealPSEs(:, 2) - RealPSEs(:, 1);
AllPreds = PredPSEs(:, 2) - PredPSEs(:, 1);

AllReals_part = AllReals;
AllPreds_part = AllPreds;

AllReals_part(LargeErrors) = []; %eliminating the big error participants
AllPreds_part(LargeErrors) = [];

% Predicted PSE vs PSE - No HV
[sorted_reals1, sorted_inds1] = sort(AllReals);
sorted_preds1 = AllPreds(sorted_inds1);

%regressing
[sim1, ~, ~, ~, stats1] = regress(AllPreds, [ones(length(AllReals),1) AllReals]);
y_fit1 = sim1(1) + sim1(2)*sorted_reals1;

[sorted_reals2, sorted_inds2] = sort(AllReals_part);
sorted_preds2 = AllPreds_part(sorted_inds2);

%regressing|
[sim2, ~, ~, ~, stats2] = regress(AllPreds_part, [ones(length(AllReals_part),1) AllReals_part]);
y_fit2 = sim2(1) + sim2(2)*sorted_reals2;

figure
hold on
plot(sorted_reals1, sorted_preds1, '*', 'linewidth', 0.5, 'markersize', 6, 'color', [0.5 0.5 0.5])
plot(AllReals(13), AllPreds(13), '*', 'linewidth', 0.5, 'markersize', 6, 'color', [109 202 16]./255)
plot(AllReals(14), AllPreds(14), '*', 'linewidth', 0.5, 'markersize', 6, 'color', [109 202 16]./255)
p1 = plot([0; sorted_reals2], [sim2(1); y_fit2], 'linewidth', 1.5, 'color', 'k');
p2 = plot([0; sorted_reals1], [sim1(1); y_fit1], 'linewidth', 1.5, 'color', [109 202 16]./255, 'linestyle', '--');

plot([-25, 85], [-25, 85], 'linewidth', 0.5, 'color', [0.7 0.7 0.7], 'linestyle', '--', 'HandleVisibility', 'Off');

set(gca, 'fontname', 'Times New Roman', 'fontsize', 12)
xlabel('Real \DeltaPSE', 'fontname', 'Times New Roman', 'fontsize', 14)
ylabel('Predicted \DeltaPSE', 'fontname', 'Times New Roman', 'fontsize', 14)

l = legend([p1, p2], ['Predicted \DeltaPSE = ', num2str(round(sim2(1)*100)/100), ' + ', num2str(round(sim2(2)*100)/100), '\cdotReal \DeltaPSE'], ...
    ['Predicted \DeltaPSE = ', num2str(round(sim1(1)*100)/100), ' + ', num2str(round(sim1(2)*100)/100), '\cdotReal \DeltaPSE'],...
    'Location', 'northwest');
l.FontSize = 12;

set(gca, 'fontname', 'Times New Roman', 'fontsize', 12)

ylim([-25 85])

xlim(get(gca, 'Ylim'))
axis square

plot([0, 0], get(gca, 'Ylim'), 'k--', 'HandleVisibility','off')
plot(get(gca, 'Xlim'), [0, 0], 'k--', 'HandleVisibility','off')
box on

%% Error sizes (Fig. 6(e))
diffs = abs(AllReals - AllPreds);
hist(diffs);
h = findobj(gca,'Type','patch');
h.FaceColor = [0.5 0.5 0.5];
xlabel('Error [N/m]')
ylabel('Number of Participants')
set(gca, 'fontname', 'Times New Roman', 'fontsize', 14)

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

AllRealsPos(LargeErrors) = [];
AllPredsPos(LargeErrors) = [];

[sorted_reals_pos, sorted_inds_pos] = sort(AllRealsPos);
sorted_preds_pos = AllPredsPos(sorted_inds_pos);

%regressing
[sim1, bint, r, rint, stats] = regress(AllPredsPos, [ones(length(AllRealsPos),1) AllRealsPos]);
y_fit_pos = sim1(1) + sim1(2)*sorted_reals_pos;

%no vel
[sorted_reals_vel2, sorted_inds_vel] = sort(AllRealsVel);
[sim22, bint, r, rint, stats] = regress(AllPredsVel, [ones(length(AllRealsVel),1) AllRealsVel]);
y_fit_vel2 = sim22(1) + sim22(2)*sorted_reals_vel2;

AllRealsVel(LargeErrors) = [];
AllPredsVel(LargeErrors) = [];

[sorted_reals_vel, sorted_inds_vel] = sort(AllRealsVel);
sorted_preds_vel = AllPredsVel(sorted_inds_vel);

%regressing
[sim2, bint, r, rint, stats] = regress(AllPredsVel, [ones(length(AllRealsVel),1) AllRealsVel]);
y_fit_vel = sim2(1) + sim2(2)*sorted_reals_vel;

% No acc
[sorted_reals_acc3, sorted_inds_acc] = sort(AllRealsAcc);
[sim33, bint, r, rint, stats] = regress(AllPredsAcc, [ones(length(AllRealsAcc),1) AllRealsAcc]);
y_fit_acc3 = sim33(1) + sim33(2)*sorted_reals_acc3;

AllRealsAcc(LargeErrors) = [];
AllPredsAcc(LargeErrors) = [];

[sorted_reals_acc, sorted_inds_acc] = sort(AllRealsAcc);
sorted_preds_acc = AllPredsPos(sorted_inds_acc);

%regressing
[sim3, bint, r, rint, stats] = regress(AllPredsAcc, [ones(length(AllRealsAcc),1) AllRealsAcc]);
y_fit_acc = sim3(1) + sim3(2)*sorted_reals_acc;


figure
p1 = plot([0; sorted_reals_pos], [sim1(1); y_fit_pos], 'linewidth', 1.5, 'color', [171 183 201]./255);
hold on
p11 = plot([0; sorted_reals_pos1], [sim11(1); y_fit_pos1], 'linewidth', 1.5, 'color', [171 183 201]./255, 'linestyle', '--');

p2 = plot([0; sorted_reals_vel], [sim2(1); y_fit_vel], 'linewidth', 1.5, 'color', [94 117 148]./255);
p22 = plot([0; sorted_reals_vel2], [sim22(1); y_fit_vel2], 'linewidth', 1.5, 'color', [94 117 148]./255, 'linestyle', '--');

p3 = plot([0; sorted_reals_acc], [sim3(1); y_fit_acc], 'linewidth', 1.5, 'color', [34 42 53]./255);
p33 = plot([0; sorted_reals_acc3], [sim33(1); y_fit_acc3], 'linewidth', 1.5, 'color', [34 42 53]./255, 'linestyle', '--');

p4 = plot([-25 85], [-25 85], 'linewidth', 0.5, 'color', [0.7 0.7 0.7], 'linestyle', '--', 'HandleVisibility', 'Off');

set(gca, 'fontname', 'Times New Roman', 'fontsize', 12)
xlabel('Real Perceptual Augmentation [N/m]', 'fontname', 'Times New Roman', 'fontsize', 14)
ylabel('Predicted Perceptual Augmentation [N/m]', 'fontname', 'Times New Roman', 'fontsize', 14)
l = legend(['{\bfNo Position:\rm} Predicted \DeltaPSE = ', num2str(round(sim1(1)*100)/100), ' + ', num2str(round(sim1(2)*100)/100), '\cdotReal \DeltaPSE'], ...
    ['Predicted \DeltaPSE = ', num2str(round(sim11(1)*100)/100), ' + ', num2str(round(sim11(2)*100)/100), '\cdotReal \DeltaPSE'], ...
    ['{\bfNo Velocity:\rm} Predicted \DeltaPSE = ', num2str(round(sim2(1)*100)/100), ' + ', num2str(round(sim2(2)*100)/100), '\cdotReal \DeltaPSE'],...
    ['Predicted \DeltaPSE = ', num2str(round(sim22(1)*100)/100), ' + ', num2str(round(sim22(2)*100)/100), '\cdotReal \DeltaPSE'],...
    ['{\bfNo Acceleration:\rm} Predicted \DeltaPSE = ', num2str(round(sim3(1)*100)/100), ' + ', num2str(round(sim3(2)*100)/100), '\cdotReal \DeltaPSE'],...
    ['Predicted \DeltaPSE = ', num2str(round(sim33(1)*100)/100), ' + ', num2str(round(sim33(2)*100)/100), '\cdotReal \DeltaPSE']);
l.FontSize = 8.8;
ylim([-25 85])
xlim([-25 85])

axis square

plot([0, 0], get(gca, 'Ylim'), 'k--', 'HandleVisibility','off')
plot(get(gca, 'Xlim'), [0, 0], 'k--', 'HandleVisibility','off')

%% Percentage of errors made by model per comparison stiffness level (Fig. 7(d-f))
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
UniqueK = unique(Kcomps);
TotalCounter2 = zeros(size(UniqueK));
ErrorCounter2 = zeros(size(UniqueK));
CorrectCounter2 = zeros(size(UniqueK));

AllPreds2 = AllPredictions(find(AllTdgains == 80));
AllLabels2 = AllLabels(find(AllTdgains == 80));
Kcomps2 = AllKComps(find(AllTdgains == 80));

for i = 1:length(Kcomps2)
    ind = find(UniqueK == Kcomps(i));
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

%% Percentage of errors made by participants per comparison stiffness level (Fig. 7(a-c))
AllRlabels = zeros(size(AllKComps)); 
%find all the places where the stiffer object is standard - so the true
%answer is 1
stands = find(AllKComps < 85);
AllRlabels(stands) = 1; %otherwise, the correct answer is comparison

%Force
TotalCounter11 = zeros(size(UniqueK));
ErrorCounter11 = zeros(size(UniqueK));
CorrectCounter11 = zeros(size(UniqueK));

AllLabels11 = AllRlabels(find(AllTdgains == 0));

for i = 1:length(Kcomps1)
    ind = find(UniqueK == Kcomps1(i));
    TotalCounter11(ind) = TotalCounter11(ind) + 1;
    
    if AllLabels11(i) == AllLabels1(i) %correct
        CorrectCounter11(ind) = CorrectCounter11(ind) + 1;
    else %mistake
        ErrorCounter11(ind) = ErrorCounter11(ind) + 1;
    end

end

% Skin Stretch
TotalCounter22 = zeros(size(UniqueK));
ErrorCounter22 = zeros(size(UniqueK));
CorrectCounter22 = zeros(size(UniqueK));

AllLabels22 = AllRlabels(find(AllTdgains == 80));

for i = 1:length(Kcomps2)
    ind = find(UniqueK == Kcomps(i));
    TotalCounter22(ind) = TotalCounter22(ind) + 1;
    
    if AllLabels22(i) == AllLabels2(i) %correct
        CorrectCounter22(ind) = CorrectCounter22(ind) + 1;
    else %mistake
        ErrorCounter22(ind) = ErrorCounter22(ind) + 1;
    end

end

% All trials together
ErrorCounter12 = ErrorCounter11 + ErrorCounter22;
TotalCounter12 = TotalCounter11 + TotalCounter22;
CorrectCounter12 = CorrectCounter11 + CorrectCounter22;

%% plotting
%All trials (Fig. 7(a))
figure
bar(UniqueK, 100*(ErrorCounter12./TotalCounter12), 'facecolor', [1 1 1])
xlabel('Comparison Stiffness Level [N/m]', 'fontsize', 14)
ylabel('Percentage of Errors [%]', 'fontsize', 14)
set(gca, 'fontname', 'Times New Roman', 'fontsize', 14)
ylim([0, 80])

%Force trials (Fig. 7(b))
figure
bar(UniqueK, 100*(ErrorCounter11./TotalCounter11), 'facecolor', [113 198 255]./255)
xlabel('Comparison Stiffness Level [N/m]', 'fontsize', 14)
ylabel('Percentage of Errors [%]', 'fontsize', 14)
set(gca, 'fontname', 'Times New Roman', 'fontsize', 14)
ylim([0, 80])

%Skin stretch trials (Fig. 7(c))
figure
bar(UniqueK, 100*(ErrorCounter22./TotalCounter22), 'facecolor', [0 121 204]./255)
xlabel('Comparison Stiffness Level [N/m]', 'fontsize', 14)
ylabel('Percentage of Errors [%]', 'fontsize', 14)
set(gca, 'fontname', 'Times New Roman', 'fontsize', 14)
ylim([0, 80])
