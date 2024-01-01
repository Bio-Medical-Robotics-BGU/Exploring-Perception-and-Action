%% This code creates Fig. 3 in the manuscript - the figure that demonstrates the 
% preprocessing of the signals
%%
ForVector = [25:216, 241:432]; %all the trials

k = 8; %take a signal from participant 8
filename = ['SN', num2str(k), '.mat'];

%in following line - replace directory to data location
cd D:\OneDrive\PerceptionAction\Data\Data;

data = load(filename);
data = data.M;

i = 100; %choose a trial to plot

j = ForVector(i);
trial_data = data{j};

%getting signals
t_ref = trial_data.DataRef(:, 1);
t_comp = trial_data.DataComp(:, 1);
allt = [t_ref; t_comp];
[t, inds] = sort(allt);

%finding comparison and standard indices
[~, comp_inds] = intersect(t, t_comp, 'stable');
[~, standard_inds] = intersect(t, t_ref, 'stable');

%creating 0/1 vector
surface_label = zeros(length(t), 1);
surface_label(standard_inds) = 1;

%position
y_ref = 1000*trial_data.DataRef(:, 3); %times by 1000 to convert from [m] to [mm]
y_comp = 1000*trial_data.DataComp(:, 3);
ally = [y_ref; y_comp];
y = ally(inds);

%interpolating signals for smoother plotting due to the breaking up of the
%signals into seperate parts
N = length(t);
taken = linspace(1, N, N*100);
t2 = interp1(1:N, t, taken, 'linear')';
y2 = interp1(1:N, y, taken, 'linear')';
surface_label2 = interp1(1:N, surface_label, taken, 'nearest')';


t_ref = t2(surface_label2 == 1);
t_comp = t2(surface_label2 == 0);

y_ref = y2(surface_label2 == 1);
y_comp = y2(surface_label2 == 0);

%in and out of surface indices
in_inds_comp = find(y_comp < 0);
out_inds_comp = find(y_comp >= 0);

in_inds_ref = find(y_ref < 0);
out_inds_ref = find(y_ref >= 0);

%in and out of surface signals
t_comp_in = t_comp(in_inds_comp);
y_comp_in = y_comp(in_inds_comp);

t_comp_out = t_comp(out_inds_comp);
y_comp_out = y_comp(out_inds_comp);

t_ref_in = t_ref(in_inds_ref);
y_ref_in = y_ref(in_inds_ref);

t_ref_out = t_ref(out_inds_ref);
y_ref_out = y_ref(out_inds_ref);

%inserting NaNs - for plotting in seperate parts
%comp in 
jumps_comp_in = find(diff(in_inds_comp) > 1);
for i = 1:length(jumps_comp_in) - 1
    t_comp_in = [t_comp_in(1:jumps_comp_in(i) + i - 1); NaN; t_comp_in((jumps_comp_in(i) + i): end)];
    y_comp_in = [y_comp_in(1:jumps_comp_in(i) + i - 1); NaN; y_comp_in((jumps_comp_in(i) + i): end)];
end
t_comp_in = [t_comp_in(1:jumps_comp_in(end) + length(jumps_comp_in) - 1); NaN; t_comp_in((jumps_comp_in(end) + length(jumps_comp_in)): end)];
y_comp_in = [y_comp_in(1:jumps_comp_in(end) + length(jumps_comp_in) - 1); NaN; y_comp_in((jumps_comp_in(end) + length(jumps_comp_in)): end)];


%comp out
jumps_comp_out1 = find(diff(out_inds_comp) > 1);
jumps_comps = find(diff(t_comp_out) > 1);
jumps_comp_out = [jumps_comp_out1; jumps_comps];
jumps_comp_out = sort(jumps_comp_out);
for i = 1:length(jumps_comp_out) - 1
    t_comp_out = [t_comp_out(1:jumps_comp_out(i) + i - 1); NaN; t_comp_out((jumps_comp_out(i) + i): end)];
    y_comp_out = [y_comp_out(1:jumps_comp_out(i) + i - 1); NaN; y_comp_out((jumps_comp_out(i) + i): end)];
end
t_comp_out = [t_comp_out(1:jumps_comp_out(end) + length(jumps_comp_out) - 1); NaN; t_comp_out((jumps_comp_out(end) + length(jumps_comp_out)): end)];
y_comp_out = [y_comp_out(1:jumps_comp_out(end) + length(jumps_comp_out) - 1); NaN; y_comp_out((jumps_comp_out(end) + length(jumps_comp_out)): end)];

%ref in 
jumps_ref_in = find(diff(in_inds_ref) > 1);
for i = 1:length(jumps_ref_in) - 1
    t_ref_in = [t_ref_in(1:jumps_ref_in(i) + i - 1); NaN; t_ref_in((jumps_ref_in(i) + i): end)];
    y_ref_in = [y_ref_in(1:jumps_ref_in(i) + i - 1); NaN; y_ref_in((jumps_ref_in(i) + i): end)];
end
t_ref_in = [t_ref_in(1:jumps_ref_in(end) + length(jumps_ref_in) - 1); NaN; t_ref_in((jumps_ref_in(end) + length(jumps_ref_in)): end)];
y_ref_in = [y_ref_in(1:jumps_ref_in(end) + length(jumps_ref_in) - 1); NaN; y_ref_in((jumps_ref_in(end) + length(jumps_ref_in)): end)];

%ref out
jumps_ref_out1 = find(diff(out_inds_ref) > 1);
jumps_refs = find(diff(t_ref_out) > 1);
jumps_ref_out = [jumps_ref_out1; jumps_refs];
jumps_ref_out = sort(jumps_ref_out);
for i = 1:length(jumps_ref_out) - 1
    t_ref_out = [t_ref_out(1:jumps_ref_out(i) + i - 1); NaN; t_ref_out((jumps_ref_out(i) + i): end)];
    y_ref_out = [y_ref_out(1:jumps_ref_out(i) + i - 1); NaN; y_ref_out((jumps_ref_out(i) + i): end)];
end
t_ref_out = [t_ref_out(1:jumps_ref_out(end) + length(jumps_ref_out) - 1); NaN; t_ref_out((jumps_ref_out(end) + length(jumps_ref_out)): end)];
y_ref_out = [y_ref_out(1:jumps_ref_out(end) + length(jumps_ref_out) - 1); NaN; y_ref_out((jumps_ref_out(end) + length(jumps_ref_out)): end)];


fig = figure;
hold on
plot(t_ref_in, y_ref_in, 'color', [209, 1, 36]./255, 'linewidth', 2);
plot(t_comp_in, y_comp_in, 'k', 'linewidth', 2);

plot(t_ref_out, y_ref_out, 'color', [255, 183, 195]./255, 'linewidth', 2, 'linestyle', '--');
plot(t_comp_out, y_comp_out, 'color', [0.6, 0.6, 0.6], 'linewidth', 2, 'linestyle', '--');

ys = get(gca, 'ylim');
xs = get(gca, 'xlim');

patch([xs(1), xs(2), xs(2), xs(1)], [ys(1), ys(1), 0, 0], 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.2)

set(gca, 'fontname', 'Times New Roman', 'fontsize', 14)
xlabel('Time [s]', 'fontname', 'Times New Roman', 'fontsize', 18)
ylabel('Position [mm]', 'fontname', 'Times New Roman', 'fontsize', 18)
xlim([min(t), 6.85])
l = legend('Standard', 'Comparison');
set(l,'FontSize', 12, 'fontname', 'Times New Roman');
fig.Position = [100 100 700 550];
box on
