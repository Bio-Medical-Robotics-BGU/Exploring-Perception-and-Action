function [pses, pred_pses, jnds, pred_jnds] = PsychometricMats(comps, gains, responses, predictions)

PsychometricData_F=85-[140:-10:30]';
PsychometricData_F(:,2:3)=0;
PsychometricData_SS=85-[140:-10:30]';
PsychometricData_SS(:,2:3)=0;

PsychometricData_F_pred=85-[140:-10:30]';
PsychometricData_F_pred(:,2:3)=0;
PsychometricData_SS_pred=85-[140:-10:30]';
PsychometricData_SS_pred(:,2:3)=0;


for i=1:length(comps) 
    ind=comps(i)/10-2;
    if (gains(i) == 0)
        %real
        if (responses(i) == 1)
            PsychometricData_F(ind,3)=PsychometricData_F(ind,3)+1;
        else
            PsychometricData_F(ind,2)=PsychometricData_F(ind,2)+1;
            PsychometricData_F(ind,3)=PsychometricData_F(ind,3)+1;
        end
        
        %predicted
        if (predictions(i) == 1)
            PsychometricData_F_pred(ind,3)=PsychometricData_F_pred(ind,3)+1;
        else
            PsychometricData_F_pred(ind,2)=PsychometricData_F_pred(ind,2)+1;
            PsychometricData_F_pred(ind,3)=PsychometricData_F_pred(ind,3)+1;
        end
    end
    
    if (gains(i) == 80)
        %real
        if (responses(i) == 1)
            PsychometricData_SS(ind,3)=PsychometricData_SS(ind,3)+1;
        else
            PsychometricData_SS(ind,2)=PsychometricData_SS(ind,2)+1;
            PsychometricData_SS(ind,3)=PsychometricData_SS(ind,3)+1;
        end
        
        %predicted
        if (predictions(i) == 1)
            PsychometricData_SS_pred(ind,3)=PsychometricData_SS_pred(ind,3)+1;
        else
            PsychometricData_SS_pred(ind,2)=PsychometricData_SS_pred(ind,2)+1;
            PsychometricData_SS_pred(ind,3)=PsychometricData_SS_pred(ind,3)+1;
        end
    end
    
end

%%
[PSE1,PSE2]=PsychometricFitting(PsychometricData_F,PsychometricData_SS, 'Real');

set(gca, 'fontname', 'Times New Roman')
set(gca, 'fontsize', 10)

fig=gcf;
fig.Position=[0, 0, 250, 210];


[PSE1_p,PSE2_p]=PsychometricFitting(PsychometricData_F_pred,PsychometricData_SS_pred, 'Predicted');

set(gca, 'fontname', 'Times New Roman')
set(gca, 'fontsize', 10)

fig=gcf;
fig.Position=[100, 100, 250, 210];

pses = [PSE1(2), PSE2(2)];
pred_pses = [PSE1_p(2), PSE2_p(2)];
jnds = [(PSE1(3) - PSE1(1))/2, (PSE2(3) - PSE2(1))/2];
pred_jnds = [(PSE1_p(3) - PSE1_p(1))/2, (PSE2_p(3) - PSE2_p(1))/2];

end

