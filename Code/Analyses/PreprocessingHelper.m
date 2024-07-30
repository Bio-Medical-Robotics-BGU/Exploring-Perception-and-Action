function PreprocessingHelper(sig, t, surface_label, y, p, yl)
%This function plots the data for the raw and preprocessed signals
%figure.
%Inputs:
%   sig - the action signal
%   t - the time vector
%   surface_label - a vector indicating when the participant interacted
%   with each of the two objects (1 for standard, 0 for comparison)
%   y - position signal
%   p - 0 or 1. Defines if we want to add the gray shaded area
%   describing where the object is. In the paper, we did this only for the
%   position signal.
%   yl - ylabel (string)

%% interpolating signals for smoother plotting due to the breaking up of the
%signals into seperate parts
N = length(t);
taken = linspace(1, N, N*100);
t2 = interp1(1:N, t, taken, 'linear')';
sig2 = interp1(1:N, sig, taken, 'linear')';

%% Seperating into the interactions with the two objects
surface_label2 = interp1(1:N, surface_label, taken, 'nearest')';

t_ref = t2(surface_label2 == 1);
t_comp = t2(surface_label2 == 0);

sig_ref = sig2(surface_label2 == 1);
sig_comp = sig2(surface_label2 == 0);

%% Finding when participant interacted with the object
%in and out of surface indices
y2 = interp1(1:N, y, taken, 'linear')';
y_ref = y2(surface_label2 == 1);
y_comp = y2(surface_label2 == 0);

in_inds_comp = find(y_comp < 0);
out_inds_comp = find(y_comp >= 0);

in_inds_ref = find(y_ref < 0);
out_inds_ref = find(y_ref >= 0);

%in and out of surface signals
t_comp_in = t_comp(in_inds_comp);
sig_comp_in = sig_comp(in_inds_comp);

t_comp_out = t_comp(out_inds_comp);
sig_comp_out = sig_comp(out_inds_comp);

t_ref_in = t_ref(in_inds_ref);
sig_ref_in = sig_ref(in_inds_ref);

t_ref_out = t_ref(out_inds_ref);
sig_ref_out = sig_ref(out_inds_ref);

%% inserting NaNs - for plotting in seperate parts
%comp in
jumps_comp_in = find(diff(in_inds_comp) > 1);
for i = 1:length(jumps_comp_in) - 1
    t_comp_in = [t_comp_in(1:jumps_comp_in(i) + i - 1); NaN; t_comp_in((jumps_comp_in(i) + i): end)];
    sig_comp_in = [sig_comp_in(1:jumps_comp_in(i) + i - 1); NaN; sig_comp_in((jumps_comp_in(i) + i): end)];
end
t_comp_in = [t_comp_in(1:jumps_comp_in(end) + length(jumps_comp_in) - 1); NaN; t_comp_in((jumps_comp_in(end) + length(jumps_comp_in)): end)];
sig_comp_in = [sig_comp_in(1:jumps_comp_in(end) + length(jumps_comp_in) - 1); NaN; sig_comp_in((jumps_comp_in(end) + length(jumps_comp_in)): end)];


%comp out
jumps_comp_out1 = find(diff(out_inds_comp) > 1);
jumps_comps = find(diff(t_comp_out) > 1);
jumps_comp_out = [jumps_comp_out1; jumps_comps];
jumps_comp_out = sort(jumps_comp_out);
for i = 1:length(jumps_comp_out) - 1
    t_comp_out = [t_comp_out(1:jumps_comp_out(i) + i - 1); NaN; t_comp_out((jumps_comp_out(i) + i): end)];
    sig_comp_out = [sig_comp_out(1:jumps_comp_out(i) + i - 1); NaN; sig_comp_out((jumps_comp_out(i) + i): end)];
end
t_comp_out = [t_comp_out(1:jumps_comp_out(end) + length(jumps_comp_out) - 1); NaN; t_comp_out((jumps_comp_out(end) + length(jumps_comp_out)): end)];
sig_comp_out = [sig_comp_out(1:jumps_comp_out(end) + length(jumps_comp_out) - 1); NaN; sig_comp_out((jumps_comp_out(end) + length(jumps_comp_out)): end)];

%ref in
jumps_ref_in = find(diff(in_inds_ref) > 1);
for i = 1:length(jumps_ref_in) - 1
    t_ref_in = [t_ref_in(1:jumps_ref_in(i) + i - 1); NaN; t_ref_in((jumps_ref_in(i) + i): end)];
    sig_ref_in = [sig_ref_in(1:jumps_ref_in(i) + i - 1); NaN; sig_ref_in((jumps_ref_in(i) + i): end)];
end
t_ref_in = [t_ref_in(1:jumps_ref_in(end) + length(jumps_ref_in) - 1); NaN; t_ref_in((jumps_ref_in(end) + length(jumps_ref_in)): end)];
sig_ref_in = [sig_ref_in(1:jumps_ref_in(end) + length(jumps_ref_in) - 1); NaN; sig_ref_in((jumps_ref_in(end) + length(jumps_ref_in)): end)];

%ref out
jumps_ref_out1 = find(diff(out_inds_ref) > 1);
jumps_refs = find(diff(t_ref_out) > 1);
jumps_ref_out = [jumps_ref_out1; jumps_refs];
jumps_ref_out = sort(jumps_ref_out);
for i = 1:length(jumps_ref_out) - 1
    t_ref_out = [t_ref_out(1:jumps_ref_out(i) + i - 1); NaN; t_ref_out((jumps_ref_out(i) + i): end)];
    sig_ref_out = [sig_ref_out(1:jumps_ref_out(i) + i - 1); NaN; sig_ref_out((jumps_ref_out(i) + i): end)];
end
t_ref_out = [t_ref_out(1:jumps_ref_out(end) + length(jumps_ref_out) - 1); NaN; t_ref_out((jumps_ref_out(end) + length(jumps_ref_out)): end)];
sig_ref_out = [sig_ref_out(1:jumps_ref_out(end) + length(jumps_ref_out) - 1); NaN; sig_ref_out((jumps_ref_out(end) + length(jumps_ref_out)): end)];

%% Plotting
fig = figure;
hold on
plot(t_ref_in, sig_ref_in, 'color', [209, 1, 36]./255, 'linewidth', 2);
plot(t_comp_in, sig_comp_in, 'k', 'linewidth', 2);

plot(t_ref_out, sig_ref_out, 'color', [255, 183, 195]./255, 'linewidth', 2, 'linestyle', '--');
plot(t_comp_out, sig_comp_out, 'color', [0.6, 0.6, 0.6], 'linewidth', 2, 'linestyle', '--');

fig.Position = [100 100 700 200];

ys = get(gca, 'ylim');
xs = get(gca, 'xlim');

if p == 1
    patch([xs(1), xs(2), xs(2), xs(1)], [ys(1), ys(1), 0, 0], 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.2)
    l = legend('Standard', 'Comparison');
    set(l,'FontSize', 12, 'fontname', 'Times New Roman');
end

set(gca, 'fontname', 'Times New Roman', 'fontsize', 12)
xlabel('Time [s]', 'fontname', 'Times New Roman', 'fontsize', 14)
ylabel(yl, 'fontname', 'Times New Roman', 'fontsize', 14)
xlim([min(t), 6.85])

box on

end

