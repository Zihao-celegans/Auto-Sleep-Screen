clear; clc; close all;

load('data2.mat');

%% Keep only strains in strainInfo
resultsTable2 = resultsTable2(ismember(resultsTable2.Strain, strainInfo.Strain_Name), :);

%% Mean per strain
[G, groupNames] = findgroups(resultsTable2.Strain);

mean_post = splitapply(@mean, resultsTable2.Qf_4hr_postUV, G);
mean_pre  = splitapply(@mean, resultsTable2.Qf_30min_preUV, G);

%% Pull N2 values
n2_idx = groupNames == "N2 UV";

if ~any(n2_idx)
    error('N2 UV not found in groupNames.');
end

n2_pre_values  = resultsTable2.Qf_30min_preUV(resultsTable2.Strain == "N2 UV");
n2_post_values = resultsTable2.Qf_4hr_postUV(resultsTable2.Strain == "N2 UV");

n2_mean_pre  = mean(n2_pre_values);
n2_mean_post = mean(n2_post_values);

%% Histogram of the mean quiescence

cOrder = colororder;
Color_postUV = cOrder(4,:);

figure();
hold on;
histogram(mean_pre, 15, ...
    'FaceColor', 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'k', 'EdgeAlpha', 0.3)
histogram(mean_post, 15, ...
    'FaceColor', Color_postUV, 'FaceAlpha', 0.3, 'EdgeColor', Color_postUV, 'EdgeAlpha', 0.3)

plot([n2_mean_pre n2_mean_pre], ylim, '--k', 'LineWidth', 2)

plot([n2_mean_post n2_mean_post], ylim, '--', 'Color', Color_postUV, 'LineWidth', 2)

hold off;


xlabel('Mean quiescence fraction')
ylabel('Strain count')

ax = gca; ax.FontSize = 25;
pbaspect([1 1 1]);

xticks(0 : 0.2 : 1);

xlim([0 1])

f = gcf;
f.Units = 'pixels';
f.Position = [100 100 800 800];
box off;

%exportgraphics(gca, "Fig3A.svg");