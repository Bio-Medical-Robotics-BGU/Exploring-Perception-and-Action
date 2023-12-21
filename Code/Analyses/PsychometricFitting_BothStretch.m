function [PSE1,PSE2,PSE3,PSE4]=PsychometricFitting_BothStretch(dat1,dat2,dat3,dat4, title_str)


if strcmp(title_str, 'Real')
    color1 = [113 198 255]./255;
    color2 = [0 121 204]./255;
    color3 = [208 144 248]./255;
    color4 = [144 14 224]./255;
else
    %red - predicted
    color1 = [255 159 159]./255;
    color2 = [192 0 0]./255;
    color3 = [248 203 173]./255;
    color4 = [197 90 17]./255;

end



figure('position',[100 100 500 300]);
hp1=plotpd(dat1,'marker','o','color',color1,'MarkerFaceColor',color1,'MarkerSize',10);
hold on;
hp2=plotpd(dat2,'marker','p','color',color2,'MarkerFaceColor',color2,'MarkerSize',10);
hp3=plotpd(dat3,'marker','d','color',color3,'MarkerFaceColor',color3,'MarkerSize',10);
hp4=plotpd(dat4,'marker','>','color',color4, 'MarkerFaceColor', color4,'MarkerSize',10);


shape = 'logistic';                                           
prefs = batch('shape', shape, 'n_intervals', 1, 'runs', 999, 'cuts',[0.25, 0.5, 0.75],'conf',[0.025 0.975]);
prefs3 = batch('shape', shape, 'n_intervals', 1, 'runs', 999, 'cuts',[0.25, 0.5, 0.75],'conf',[0.025 0.975]);
outputPrefs1 = batch('write_pa', 'pa1', 'write_th', 'th1');
outputPrefs2 = batch('write_pa', 'pa2', 'write_th', 'th2');
outputPrefs3 = batch('write_pa', 'pa3', 'write_th', 'th3');
outputPrefs4 = batch('write_pa', 'pa4', 'write_th', 'th4');
psignifit(dat1, [prefs outputPrefs1]);
psignifit(dat2, [prefs outputPrefs2]);
psignifit(dat3, [prefs outputPrefs3]);
psignifit(dat4, [prefs outputPrefs4]);

h1=plotpf(shape, pa1.est,'color',color1,'LineStyle','-','linewidth',2);
h2=plotpf(shape, pa2.est,'color',color2,'LineStyle','-','linewidth',2);
h3=plotpf(shape, pa3.est,'color',color3,'LineStyle','-','linewidth',2);
h4=plotpf(shape, pa4.est,'color',color4,'LineStyle','-','linewidth',2);


%%
drawHeights1 = psi(shape, pa1.est, th1.est);
line(th1.lims(:, 2),ones(size(th1.lims(:, 2),1),1)*drawHeights1(2),'color',color1,'LineStyle','-');
drawHeights2 = psi(shape, pa2.est, th2.est);
line(th2.lims(:, 2),ones(size(th2.lims(:, 2),1),1)*drawHeights2(2),'color',color2,'LineStyle','-');
drawHeights3 = psi(shape, pa3.est, th3.est);
line(th3.lims(:, 2),ones(size(th3.lims(:, 2),1),1)*drawHeights3(2),'color',color3,'LineStyle','-');
drawHeights4 = psi(shape, pa4.est, th4.est);
line(th4.lims(:, 2),ones(size(th4.lims(:, 2),1),1)*drawHeights4(2),'color',color4,'LineStyle','-');

    ax = gca;
    ax.FontSize = 18;
    ax.YTick=[0, 0.25, 0.5, 0.75, 1];
    fig = gcf;
    fig.Position = [0 0 250 400];

h = legend([h1 h2 h3 h4],'FP','SP','FN','SN','Location','northwest');

title(title_str)
legend('Boxoff');
h.FontSize = 10;
h.FontName = 'Times New Roman';
PSE1=th1.est;
PSE2=th2.est;
PSE3=th3.est;
PSE4=th4.est;
