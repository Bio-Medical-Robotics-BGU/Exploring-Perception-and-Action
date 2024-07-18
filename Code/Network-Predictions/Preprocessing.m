%% This code does the preprocessing of the participants' data
% Data extraction
% Filtering
% Interpolation
% Creating the matrices and vectors containing each participants
% preprocessed data

project_path = "C:\Users\hannako\Downloads\ExploringPerceptionAction\ExploringPerceptionandAction";
%%
forvector = [1:21, 23:39]; %which participants to include (two outlier participants were removed)
ForVector = [25:216, 241:432]; %the relevant trials (excluding training trials)

ds = input ('Which dataset would you like to create? \n 1. Without Negative Stretch \n 2. With Negative Stretch\n');
for k = 1:length(forvector) %loop that runs over all the participants
    h = forvector(k);
    
    err = 0;
    
    AllCompSigs = [];
    AllRefSigs = [];
    AllPlabels = [];
    AllKcomps = [];
    AllTdGains = [];
    
    
    
    filename = ['SN', num2str(h), '.mat'];
    
    %data location
    cd(fullfile(project_path, "Data\Participant-Data"));
    
    data = load(filename);
    data = data.M;
    for i = 1:length(ForVector) %loop that runs over all the trials
        j = ForVector(i);
        trial_data = data{j}; %get trial data
        try %this is for the case of data not saved properly during the
            %stiffness discrimination experiment. This happened for one
            %trial for one participant (Participant 8, j = 253)
            
            dt_gain = trial_data.Gain;
            
            %The next line defines that only the data from force or positive stretch trials
            % are included.
            if ds == 1
                if strcmp(dt_gain, '-80')
                    continue
                end
            end
            k_standard = trial_data.RefStiffnessVal;
            k_comp = trial_data.CompStiffnessVal;
            
            %getting participant answer
            Panswer_str = trial_data.Answer;
            if strcmp(Panswer_str, 'Comp')
                Panswer = 0;
            elseif strcmp(Panswer_str, 'Ref')
                Panswer = 1;
            end
            
            %getting signals:
            %time
            t_ref = trial_data.DataRef(:, 1);
            t_comp = trial_data.DataComp(:, 1);
            allt = [t_ref; t_comp];
            [t, inds] = sort(allt);
            
            
            %position
            y_ref = trial_data.DataRef(:, 3);
            y_comp = trial_data.DataComp(:, 3);
            ally = [y_ref; y_comp];
            y = ally(inds);
            
            %velocity
            yvel_ref = trial_data.DataRef(:, 6);
            yvel_comp = trial_data.DataComp(:, 6);
            allyvel = [yvel_ref; yvel_comp];
            yvel = allyvel(inds);
            
            %grip force
            gf_ref = -trial_data.DataRef(:, 13);
            gf_comp = -trial_data.DataComp(:, 13);
            allgf = [gf_ref; gf_comp];
            gf = allgf(inds);
            
            %force
            lf_ref = trial_data.DataRef(:, 9);
            lf_comp = trial_data.DataComp(:, 9);
            alllf = [lf_ref; lf_comp];
            lf = alllf(inds);
            
            %finding indices of intercation with the comparison and
            %standard objects
            [~, comp_inds] = intersect(t, t_comp, 'stable');
            [~, standard_inds] = intersect(t, t_ref, 'stable');
            
            %creating 0/1 vector
            surface_label = zeros(length(t), 1);
            surface_label(standard_inds) = 1;
            
            
            %filtering
            %getting average sampling frequency
            diffed = diff(t);
            if length(find(diffed == 0))~= 0
                badinds = find(diffed == 0);
                goodinds = setdiff(length(diffed), badinds);
                diffed = diffed(goodinds);
            end
            fs = mean(1./diffed);
            if length(find(diffed == 0))~= 0
                diffed(badinds) = 1/fs;
            end
            
            %computing acceleration
            acc = diff(yvel)./diffed;
            acc = [acc; acc(end)];
            
            
            %filtering signals
            d1 = designfilt('lowpassiir', 'FilterOrder', 2, 'HalfPowerFrequency', 12, 'DesignMethod', 'butter', 'SampleRate', fs);
            gfd = filtfilt(d1, gf);
            d2 = designfilt('lowpassiir', 'FilterOrder', 2, 'HalfPowerFrequency', 10, 'DesignMethod', 'butter', 'SampleRate', fs);
            yd = filtfilt(d2, y);
            yveld = filtfilt(d2, yvel);
            accd = filtfilt(d2, acc);
            
            
            %getting relevant indices - those of the interactions with
            %the objects (in which force was applied)
            rel_inds = find(lf > 0);
            
            y_rel = yd(rel_inds);
            yvel_rel = yveld(rel_inds);
            gf_rel = gfd(rel_inds);
            acc_rel = accd(rel_inds);
            
            surface_label_rel = surface_label(rel_inds);
            
            standard_inds_rel = find(surface_label_rel == 1);
            comp_inds_rel = find(surface_label_rel == 0);
            
            %seperating into the interactions with standard and
            %comparison
            
            y_ref = y_rel(standard_inds_rel);
            yvel_ref = yvel_rel(standard_inds_rel);
            acc_ref = acc_rel(standard_inds_rel);
            gf_ref = gf_rel(standard_inds_rel);
            
            y_comp = y_rel(comp_inds_rel);
            yvel_comp = yvel_rel(comp_inds_rel);
            acc_comp = acc_rel(comp_inds_rel);
            gf_comp = gf_rel(comp_inds_rel);
            
            sigs_comp = [y_comp, yvel_comp, acc_comp, gf_comp];
            sigs_ref = [y_ref, yvel_ref, acc_ref, gf_ref];
            
            %Interpolation and normalization
            interp_lenght = 150;
            
            taken_comp = linspace(1, size(sigs_comp, 1), interp_lenght);
            taken_ref = linspace(1, size(sigs_ref, 1), interp_lenght);
            interped_comp = zeros(interp_lenght, size(sigs_comp, 2));
            interped_ref = zeros(interp_lenght, size(sigs_ref, 2));
            
            for kk = 1:size(sigs_ref, 2)
                temp_comp = interp1(1:size(sigs_comp, 1), sigs_comp(:, kk), taken_comp, 'linear');
                temp_ref = interp1(1:size(sigs_ref, 1), sigs_ref(:, kk), taken_ref, 'linear');
                
                combined = [temp_comp, temp_ref];
                
                meani = mean(combined);
                stdi = std(combined);
                
                interped_comp(:, kk) = (temp_comp - meani)/stdi;
                interped_ref(:, kk) = (temp_ref - meani)/stdi;
            end
            
            %final setup into matrices and vectors
            AllCompSigs = cat(1, AllCompSigs, reshape(interped_comp, 1, size(interped_comp, 1), size(interped_comp, 2)));
            AllRefSigs = cat(1, AllRefSigs, reshape(interped_ref, 1, size(interped_ref, 1), size(interped_ref, 2)));
            AllPlabels = [AllPlabels; Panswer];
            AllKcomps = [AllKcomps; k_comp];
            AllTdGains = [AllTdGains; str2num(trial_data.Gain)];
        catch
            err = err + 1;
        end
        disp(i)
        
    end%end of running of k's trials
    
    %location for saving data
    mkdir(fullfile(project_path, 'Preprocessed'));
    cd(fullfile(project_path, 'Preprocessed'));
    if ds == 1
        save(['CompSignals_SN', num2str(h)], 'AllCompSigs')
        save(['StandardSignals_SN', num2str(h)], 'AllRefSigs')
        save(['Labels_SN', num2str(h)], 'AllPlabels')
        save(['Kcomps_SN', num2str(h)], 'AllKcomps')
        save(['TdGains_SN', num2str(h)], 'AllTdGains')
        
    else
        %to save for all trials - including both positive and negative stretch
        save(['CompSignals_BothStretch_SN', num2str(h)], 'AllCompSigs')
        save(['StandardSignals_BothStretch_SN', num2str(h)], 'AllRefSigs')
        save(['Labels_BothStretch_SN', num2str(h)], 'AllPlabels')
        save(['Kcomps_BothStretch_SN', num2str(h)], 'AllKcomps')
        save(['TdGains_BothStretch_SN', num2str(h)], 'AllTdGains')
    end
    disp(k)
    
end%end of loop on all participants
