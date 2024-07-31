%% This code creates Fig. 3 in the manuscript - the figure that demonstrates the 
% preprocessing of the signals
%% Choose participant for plotting and get data
ForVector = [25:216, 241:432]; %all the trials

k = 13; %take a signal from participant 13
filename = ['SN', num2str(k), '.mat'];

%in following line - replace directory to data location
cd D:\OneDrive\PerceptionActionReview\Data\Data;

data = load(filename);
data = data.M;

i = 100; %choose a trial to plot

j = ForVector(i);
trial_data = data{j};
%% Time data
%getting signals
t_ref = trial_data.DataRef(:, 1);
t_comp = trial_data.DataComp(:, 1);
allt = [t_ref; t_comp];
[t, inds] = sort(allt);

%% Interactions with each object
%finding comparison and standard indices
[~, comp_inds] = intersect(t, t_comp, 'stable');
[~, standard_inds] = intersect(t, t_ref, 'stable');

%creating 0/1 vector
surface_label = zeros(length(t), 1);
surface_label(standard_inds) = 1;

%% Action Signals
%position
y_ref = 1000*trial_data.DataRef(:, 3); %times by 1000 to convert from [m] to [mm]
y_comp = 1000*trial_data.DataComp(:, 3);
ally = [y_ref; y_comp];
y = ally(inds);

%velocity
v_ref = trial_data.DataRef(:, 6); 
v_comp = trial_data.DataComp(:, 6);
allv = [v_ref; v_comp];
v = allv(inds);

%grip force
gf_ref = -trial_data.DataRef(:, 13); 
gf_comp = -trial_data.DataComp(:, 13);
allgf = [gf_ref; gf_comp];
gf = allgf(inds);

%% Filtered Versions
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

d1 = designfilt('lowpassiir', 'FilterOrder', 2, 'HalfPowerFrequency', 12, 'DesignMethod', 'butter', 'SampleRate', fs);
gfd = filtfilt(d1, gf);
d2 = designfilt('lowpassiir', 'FilterOrder', 2, 'HalfPowerFrequency', 10, 'DesignMethod', 'butter', 'SampleRate', fs);
yd = filtfilt(d2, y);
vd = filtfilt(d2, v);

%% Raw
PreprocessingHelper(y, t, surface_label, y, 1, 'Position [mm]', 1)
PreprocessingHelper(v, t, surface_label, y, 0, 'Velocity [mm/s]', 1)
PreprocessingHelper(gf, t, surface_label, y,0, 'Grip Force [N]', 1)

%% Filtered
PreprocessingHelper(yd, t, surface_label, y, 1, 'Position [mm]', 0)
PreprocessingHelper(vd, t, surface_label, y, 0, 'Velocity [mm/s]', 0)
PreprocessingHelper(gfd, t, surface_label, y, 0, 'Grip Force [N]', 0)

