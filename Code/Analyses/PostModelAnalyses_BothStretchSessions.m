%% This code analyses the network predictions - creates the real and predicted curves and computes the 
% regression of the predicted augmentation against the real one, and compares the real and
% predicted curves in each condition for both stretch sessions (Fig. 10).

%replace following directory with the project path
clear all; clc; 
project_path = "C:\Users\hannako\Downloads\ExploringPerceptionandAction\ExploringPerceptionandAction";

addpath(genpath(project_path))

%load the splits
folds = load('TestIndsSplit_AllParticipants.mat');

Run = input('Which number run would you like? \n', 's');

%getting the field names of the struct
fns = fieldnames(folds);

%for PSE and JND errors, and regression
RealPSEs = [];
PredPSEs = [];
RealJNDs = [];
PredJNDS = [];

%for accuracy
Accuracies = zeros(10, 1);
count = 0;
%% Creating psychometric curves, and getting PSE and JND values
% The psychometric curves in Fig. 11(a-c) were created here.

for f = 1:10 %run over the 10 folds
    %getting the participants in this fold
    participants = folds.(fns{f});
    
    FoldLabels = [];
    FoldPreds = [];
    %running over the participants of the fold
    for p = 1:length(participants)

        count = count + 1; 

       
        Labels = load(['Labels_BothStretch_SN', num2str(participants(p)), '.mat']);
        Labels = Labels.AllPlabels;
        
        KComps = load(['Kcomps_BothStretch_SN', num2str(participants(p)), '.mat']);
        KComps = KComps.AllKcomps;
        
        TdGains = load(['TdGains_BothStretch_SN', num2str(participants(p)), '.mat']);
        TdGains = TdGains.AllTdGains;

        %getting predictions
        preds = load(['Preds_SN', num2str(participants(p)), '_BothStretch_Run', Run, '.mat']);
        preds = preds.Preds;

        [pses, pred_pses, jnds, pred_jnds] = PsychometricMatsBothStretch(KComps, TdGains, Labels, preds);


        RealPSEs = [RealPSEs; pses];
        PredPSEs = [PredPSEs; pred_pses];
        RealJNDs = [RealJNDs; jnds];
        PredJNDS = [PredJNDS; pred_jnds];

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
PSE_ForceErrors2 = RealPSEs(:, 3) - PredPSEs(:, 3);
PSE_NegStretchErrors = RealPSEs(:, 4) - PredPSEs(:, 4);

mes = ['PSE average errors are: \n Force condition in positive stretch session ',...
    num2str(mean(PSE_ForceErrors)), '\n Positive stretch condition ',...
    num2str(mean(PSE_StretchErrors)),...
    '\n Force condition in negative stretch session ', num2str(mean(PSE_ForceErrors2)),...
    '\n Negative stretch condition ', num2str(mean(PSE_NegStretchErrors))];

fprintf(mes)

% Plotting Fig. 11(d)
hold on
b = bar([0.5, 1, 2, 2.5, 3.5, 4, 5, 5.5],...
    [mean(RealPSEs(:, 1)), mean(PredPSEs(:, 1)), mean(RealPSEs(:, 2)), mean(PredPSEs(:, 2)), mean(RealPSEs(:, 3)), mean(PredPSEs(:, 3)), mean(RealPSEs(:, 4)), mean(PredPSEs(:, 4))]);
b.FaceColor = 'flat';
CIr1 = bootci(2000, @mean, RealPSEs(:, 1));
CIp1 = bootci(2000, @mean, PredPSEs(:, 1));
b.CData(1,:) = [113 198 255]./255;
errorbar(0.5, mean(RealPSEs(:, 1)), mean(RealPSEs(:, 1)) - CIr1(1), CIr1(2) - mean(RealPSEs(:, 1)), 'color', 'k', 'linewidth', 0.5)
b.CData(2,:) = [255 159 159]./255;
errorbar(1, mean(PredPSEs(:, 1)), mean(PredPSEs(:, 1)) - CIp1(1), CIp1(2) - mean(PredPSEs(:, 1)), 'color', 'k', 'linewidth', 0.5)

CIr2 = bootci(2000, @mean, RealPSEs(:, 2));
CIp2 = bootci(2000, @mean, PredPSEs(:, 2));
b.CData(3,:) = [0 121 204]./255;
errorbar(2, mean(RealPSEs(:, 2)), mean(RealPSEs(:, 2)) - CIr2(1), CIr2(2)-  mean(RealPSEs(:, 2)), 'color', 'k', 'linewidth', 0.5)
b.CData(4,:) = [192 0 0]./255;
errorbar(2.5, mean(PredPSEs(:, 2)), mean(PredPSEs(:, 2))- CIp2(1), CIp2(2) - mean(PredPSEs(:, 2)), 'color', 'k', 'linewidth', 0.5)

CIr3 = bootci(2000, @mean, RealPSEs(:, 3));
CIp3 = bootci(2000, @mean, PredPSEs(:, 3));
b.CData(5,:) = [208 144 248]./255;
errorbar(3.5, mean(RealPSEs(:, 3)), mean(RealPSEs(:, 3)) - CIr3(1), CIr3(2)-  mean(RealPSEs(:, 3)), 'color', 'k', 'linewidth', 0.5)
b.CData(6,:) = [248 203 173]./255;
errorbar(4, mean(PredPSEs(:, 3)), mean(PredPSEs(:, 3))- CIp3(1), CIp3(2) - mean(PredPSEs(:, 3)), 'color', 'k', 'linewidth', 0.5)

CIr4 = bootci(2000, @mean, RealPSEs(:, 4));
CIp4 = bootci(2000, @mean, PredPSEs(:, 4));
b.CData(7,:) = [144 14 224]./255;
errorbar(5, mean(RealPSEs(:, 4)), mean(RealPSEs(:, 4)) - CIr4(1), CIr4(2)-  mean(RealPSEs(:, 4)), 'color', 'k', 'linewidth', 0.5)
b.CData(8,:) = [197 90 17]./255;
errorbar(5.5, mean(PredPSEs(:, 4)), mean(PredPSEs(:, 4))- CIp4(1), CIp4(2) - mean(PredPSEs(:, 4)), 'color', 'k', 'linewidth', 0.5)

set(gca, 'fontname', 'Times New Roman', 'fontsize', 14)

set(gca, 'Xtick', [0.5, 1, 2, 2.5, 3.5, 4, 5, 5.5], 'XtickLabels',...
    {'Real FP', 'Predicted FP', 'Real SP', 'Predicted SP', 'Real FN', 'Predicted FN', 'Real SN', 'Predicted SN'}, 'fontname', 'Times New Roman', 'fontsize', 14)
ylabel('Average PSE [N/m]', 'fontname', 'Times New Roman', 'fontsize', 16)

fig = gcf;
fig.Position = [200 200 800 400];

box on

%% JND Error
clc;

JND_ForceErrors = RealJNDs(:, 1) - PredJNDS(:, 1);
JND_StretchErrors = RealJNDs(:, 2) - PredJNDS(:, 2);
JND_ForceErrors2 = RealJNDs(:, 3) - PredJNDS(:, 3);
JND_NegStretchErrors = RealJNDs(:, 4) - PredJNDS(:, 4);

mes2 = ['JND average errors are: \n Force condition in positive stretch session ',...
    num2str(mean(JND_ForceErrors)), '\n Positive stretch condition ',...
    num2str(mean(JND_StretchErrors)),...
    '\n Force condition in negative stretch session ', num2str(mean(JND_ForceErrors2)),...
    '\n Negative stretch condition ', num2str(mean(JND_NegStretchErrors))];

fprintf(mes2)

%% Regresssion
% Creates the plot for Fig. 11(e)

AllRealsPos = RealPSEs(:, 2) - RealPSEs(:, 1);
AllPredsPos = PredPSEs(:, 2) - PredPSEs(:, 1);

AllRealsNeg = RealPSEs(:, 4) - RealPSEs(:, 3);
AllPredsNeg = PredPSEs(:, 4) - PredPSEs(:, 3);

AllReals = [AllRealsPos; AllRealsNeg];
AllPreds = [AllPredsPos; AllPredsNeg];

% Predicted PSE vs PSE 
[sorted_reals1, sorted_inds1] = sort(AllReals);
sorted_preds1 = AllPreds(sorted_inds1);

%regressing
[sim1, ~, ~, ~, stats1] = regress(AllPreds, [ones(length(AllReals),1) AllReals]);
y_fit1 = sim1(1) + sim1(2)*sorted_reals1;


Marks = ['o', '+', '.', 'x', 's', 'd', '<', '>', 'v', '^', 'p', 'h'];
    
figure
hold on
for i = 1:12
    plot(AllReals(i), AllPreds(i), 'marker', Marks(i), 'linewidth', 0.5, 'markersize', 6, 'color', [0 121 204]./255)
    plot(AllReals(length(AllRealsNeg) + i), AllPreds(length(AllRealsNeg) + i), 'marker', Marks(i), 'linewidth', 0.5, 'markersize', 6, 'color', [144 14 224]./255)
end

Marks2 = ['o', 's', 'd', '<', '>', 'v', '^', 'p', 'h'];
for i = 13:21
    plot(AllReals(i), AllPreds(i), 'marker', Marks2(i-12), 'linewidth', 0.5, 'markersize', 6, 'color', [0 121 204]./255, 'markerfacecolor', [0 121 204]./255)
    plot(AllReals(length(AllRealsNeg) + i), AllPreds(length(AllRealsNeg) + i), 'marker', Marks2(i-12), 'linewidth', 0.5, 'markersize', 6, 'color', [144 14 224]./255, 'markerfacecolor', [144 14 224]./255)
end

font = 'ZapfDingbats';
m = char([0x2721, 0x2740, 0x2722, 0x2737, 0x2713, 0x2744, 0x2764, 0x273A, 0x274B, 0x2768, 0x273F, 0x2765, 0x27A4, 0x270F, 0x2708, 0x2756, 0x2746]);
for i = 22:38
    text(AllReals(i), AllPreds(i), m(i-21),'fontname',font,'fontsize',10,'color', [0 121 204]./255)
    text(AllReals(length(AllRealsNeg) + i), AllPreds(length(AllRealsNeg) + i), m(i-21),'fontname',font,'fontsize',10,'color', [144 14 224]./255)
end

p2 = plot([0; sorted_reals1], [sim1(1); y_fit1], 'linewidth', 1.5, 'color', 'k');

plot([-40, 120], [-40, 120], 'linewidth', 0.5, 'color', [0.7 0.7 0.7], 'linestyle', '--', 'HandleVisibility', 'Off');

set(gca, 'fontname', 'Times New Roman', 'fontsize', 12)
xlabel('Real \DeltaPSE', 'fontname', 'Times New Roman', 'fontsize', 14)
ylabel('Predicted \DeltaPSE', 'fontname', 'Times New Roman', 'fontsize', 14)

l = legend(p2, ['Predicted \DeltaPSE = ', num2str(round(sim1(1)*100)/100), ' + ', num2str(round(sim1(2)*100)/100), '\cdotReal \DeltaPSE'],...
    'Location', 'northwest');
l.FontSize = 12;

set(gca, 'fontname', 'Times New Roman', 'fontsize', 12)

ylim([-20 100])

xlim(get(gca, 'Ylim'))
axis square

plot([0, 0], get(gca, 'Ylim'), 'k--', 'HandleVisibility','off')
plot(get(gca, 'Xlim'), [0, 0], 'k--', 'HandleVisibility','off')
box on

disp(stats1)

%% plotting Delta PSEs
Real_Pos_Minus_Neg = AllRealsPos - AllRealsNeg; %difference between perceptual augmentation caused 
%by positive stretch and that caused by negative stretch.
Pred_Pos_Minus_Neg = AllPredsPos - AllPredsNeg; %the predicted difference

sign_real = sign(Real_Pos_Minus_Neg);
sign_pred = sign(Pred_Pos_Minus_Neg);
allerrs = find(sign_real ~= sign_pred);
allcors = setdiff([1:38], allerrs);


