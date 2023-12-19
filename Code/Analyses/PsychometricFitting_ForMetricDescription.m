function [PSE1,PSE3]=PsychometricFitting_ForMetricDescription(dat1 ,dat3, title_str)

%blue - real
if strcmp(title_str, 'Real')
    color1 = [113 198 255]./255;
    color2 = [0 121 204]./255;
else
    %red - predicted
    color1 = [255 159 159]./255;
    color2 = [192 0 0]./255;
end


figure('position',[100 100 500 300]);
hp1=plotpd(dat1,'marker','o','color',[255 255 255]./255,'MarkerFaceColor',[255 255 255]./255,'MarkerSize',8);
hold on;
hp3=plotpd(dat3,'marker','d','color',[255 255 255]./255,'MarkerFaceColor',[255 255 255]./255,'MarkerSize',8);

shape = 'logistic';                                            % [0.25 0.5 0.75]
prefs = batch('shape', shape, 'n_intervals', 1, 'runs', 999, 'cuts',[0.25, 0.5, 0.75],'conf',[0.025 0.975]);
prefs3 = batch('shape', shape, 'n_intervals', 1, 'runs', 999, 'cuts',[0.25, 0.5, 0.75],'conf',[0.025 0.975]);
outputPrefs1 = batch('write_pa', 'pa1', 'write_th', 'th1');
outputPrefs3 = batch('write_pa', 'pa3', 'write_th', 'th3');
psignifit(dat1, [prefs outputPrefs1]);
psignifit(dat3, [prefs outputPrefs3]);


h1=plotpf(shape, pa1.est,'color',color1,'LineStyle','-','linewidth',2);
h3=plotpf(shape, pa3.est,'color',color2,'LineStyle','-','linewidth',2);

ax = gca;
ax.FontSize = 18;
ax.YTick=[0, 0.25, 0.5, 0.75, 1];
fig = gcf;
fig.Position = [0 0 250 400];

h = legend([h1 h3],'Force','Stretch','Location','northwest');
title(title_str)

legend('Boxoff');
h.FontSize = 10;
h.FontName = 'Times New Roman';
PSE1=th1.est;
PSE3=th3.est;
