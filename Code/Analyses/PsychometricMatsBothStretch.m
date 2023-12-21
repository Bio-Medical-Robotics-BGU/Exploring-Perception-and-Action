function [pses, pred_pses, jnds, pred_jnds] = PsychometricMatsBothStretch(comps, gains, responses, predictions)

%positive stretch session
PsychometricData_F=85-[140:-10:30]';
PsychometricData_F(:,2:3)=0;
PsychometricData_SS=85-[140:-10:30]';
PsychometricData_SS(:,2:3)=0;

PsychometricData_F_pred=85-[140:-10:30]';
PsychometricData_F_pred(:,2:3)=0;
PsychometricData_SS_pred=85-[140:-10:30]';
PsychometricData_SS_pred(:,2:3)=0;

%negative stretch session
PsychometricData_F2=85-[140:-10:30]';
PsychometricData_F2(:,2:3)=0;
PsychometricData_SS_Neg=85-[140:-10:30]';
PsychometricData_SS_Neg(:,2:3)=0;

PsychometricData_F2_pred=85-[140:-10:30]';
PsychometricData_F2_pred(:,2:3)=0;
PsychometricData_SS_Neg_pred=85-[140:-10:30]';
PsychometricData_SS_Neg_pred(:,2:3)=0;

%positive stretch session
for i=1:192 
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

%negative stretch session
for i=193:length(comps) 
    ind=comps(i)/10-2;
    if (gains(i) == 0)
        %real
        if (responses(i) == 1)
            PsychometricData_F2(ind,3)=PsychometricData_F2(ind,3)+1;
        else
            PsychometricData_F2(ind,2)=PsychometricData_F2(ind,2)+1;
            PsychometricData_F2(ind,3)=PsychometricData_F2(ind,3)+1;
        end
        
        %predicted
        if (predictions(i) == 1)
            PsychometricData_F2_pred(ind,3)=PsychometricData_F2_pred(ind,3)+1;
        else
            PsychometricData_F2_pred(ind,2)=PsychometricData_F2_pred(ind,2)+1;
            PsychometricData_F2_pred(ind,3)=PsychometricData_F2_pred(ind,3)+1;
        end
    end
    
    if (gains(i) == -80)
        %real
        if (responses(i) == 1)
            PsychometricData_SS_Neg(ind,3)=PsychometricData_SS_Neg(ind,3)+1;
        else
            PsychometricData_SS_Neg(ind,2)=PsychometricData_SS_Neg(ind,2)+1;
            PsychometricData_SS_Neg(ind,3)=PsychometricData_SS_Neg(ind,3)+1;
        end
        
        %predicted
        if (predictions(i) == 1)
            PsychometricData_SS_Neg_pred(ind,3)=PsychometricData_SS_Neg_pred(ind,3)+1;
        else
            PsychometricData_SS_Neg_pred(ind,2)=PsychometricData_SS_Neg_pred(ind,2)+1;
            PsychometricData_SS_Neg_pred(ind,3)=PsychometricData_SS_Neg_pred(ind,3)+1;
        end
    end
end
%%
[PSE1, PSE2, PSE3, PSE4]=PsychometricFitting_BothStretch(PsychometricData_F,PsychometricData_SS, PsychometricData_F2, PsychometricData_SS_Neg,'Real');

set(gca, 'fontname', 'Times New Roman')
set(gca, 'fontsize', 10)

fig=gcf;
fig.Position=[0, 0, 250, 210];

[PSE1_p, PSE2_p, PSE3_p, PSE4_p]=PsychometricFitting_BothStretch(PsychometricData_F_pred,PsychometricData_SS_pred, PsychometricData_F2_pred, PsychometricData_SS_Neg_pred, 'Predicted');

set(gca, 'fontname', 'Times New Roman')
set(gca, 'fontsize', 10)

fig=gcf;
fig.Position=[100, 100, 250, 210];

pses = [PSE1(2), PSE2(2), PSE3(2), PSE4(2)];
pred_pses = [PSE1_p(2), PSE2_p(2), PSE3_p(2), PSE4_p(2)];

jnds = [(PSE1(3) - PSE1(1))/2, (PSE2(3) - PSE2(1))/2, (PSE3(3) - PSE3(1))/2, (PSE4(3) - PSE4(1))/2];
pred_jnds = [(PSE1_p(3) - PSE1_p(1))/2, (PSE2_p(3) - PSE2_p(1))/2, (PSE3_p(3) - PSE3_p(1))/2, (PSE4_p(3) - PSE4_p(1))/2];


end

