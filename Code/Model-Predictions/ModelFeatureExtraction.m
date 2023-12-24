%% This code computes the features for the models
% Data extraction 
% Filtering
% Computes the features for the models
%%
forvector = [1:21, 23:39];
ForVector = [25:216, 241:432];
NumFeatures = 4;

for k = 1:length(forvector)
    h = forvector(k);
    
    AllFeatures = zeros(288, 2, NumFeatures);
    
    count = 1;
    err = 0;
    
    filename = ['SN', num2str(h), '.mat'];
    
    %in following line - replace directory to data location
    cd C:\Users\hanna\OneDrive\MATLAB\lab\First_Degree\NegativeStretch;

    data = load(filename);
    data = data.M;

    for i = 1:length(ForVector) %loop that runs over all the trials
        j = ForVector(i);
        trial_data = data{j}; %get trial data
        try %this is for the case of data not saved properly during the
            %stiffness discrimination experiment. This happened for one
            %trial for one participant (Participant 8, i = 205)

            dt_gain = trial_data.Gain;
            %The next line defines that only the data from force or positive stretch trials
            % are included.

            if (strcmp(dt_gain, '0') || strcmp(dt_gain, '80')) 

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
                
                %filtering signals
                d1 = designfilt('lowpassiir', 'FilterOrder', 2, 'HalfPowerFrequency', 12, 'DesignMethod', 'butter', 'SampleRate', fs);
                gfd = filtfilt(d1, gf);
                lfd = filtfilt(d1, lf);
                d2 = designfilt('lowpassiir', 'FilterOrder', 2, 'HalfPowerFrequency', 10, 'DesignMethod', 'butter', 'SampleRate', fs);
                yd = filtfilt(d2, y);
                yveld = filtfilt(d2, yvel);
                
                
                %splitting filtered signals back into standard and comparison
                %grip force
                filtered_gf_comp = gfd(comp_inds);
                filtered_gf_ref = gfd(standard_inds);
                
                %load force
                filtered_lf_comp = lfd(comp_inds);
                filtered_lf_ref = lfd(standard_inds);
                
                %position
                filtered_y_comp = yd(comp_inds);
                filtered_y_ref = yd(standard_inds);
                
                %velocity
                filtered_yvel_comp = yveld(comp_inds);
                filtered_yvel_ref = yveld(standard_inds);
                
                %relevant indices
                rel_refs = find(lf_ref > 0);
                rel_comps = find(lf_comp > 0);
                
                %getting features
                
                %1. Min position (max penetration)
                Min_pos_ref = -min(filtered_y_ref(rel_refs));
                Min_pos_comp = -min(filtered_y_comp(rel_comps));
                
                %2. Max velocity
                Max_vel_ref = max(filtered_yvel_ref(rel_refs));
                Max_vel_comp = max(filtered_yvel_comp(rel_comps));
                
                %3. Mean velocity
                Mean_vel_ref = mean(filtered_yvel_ref(rel_refs));
                Mean_vel_comp = mean(filtered_yvel_comp(rel_comps));
                
                %4. Max grip force
                Max_gf_ref = max(filtered_gf_ref(rel_refs));
                Max_gf_comp = max(filtered_gf_comp(rel_comps));
                
                Ref_features = [Min_pos_ref, Max_vel_ref, Mean_vel_ref, Max_gf_ref];
                
                Comp_features = [Min_pos_comp, Max_vel_comp, Mean_vel_comp, Max_gf_comp];
                
                AllFeatures(count, 1, :) = Ref_features;
                AllFeatures(count, 2, :) = Comp_features;
                
                count = count + 1;
                
            end%end of no noise condition
        catch
            err = err + 1;
        end
        
        disp(i)
    end%end of running of k's trials
    
    AllFeatures = AllFeatures(1:end - err, :, :);
    
    %replace following directory with desired location for saving data
    cd C:\Users\hanna\OneDrive\MATLAB\lab\PhD\Perception_and_GF_prediction\DL_perception\OnlyInSurface\ParticipantMatsAndVecs\final;
    save(['Features_', filename(1:end-4)], 'AllFeatures');
   
    disp(k)
end%end of loop on all participants