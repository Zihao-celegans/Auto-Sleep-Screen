clear; clc; close all;

load('data.mat');

%% Filtered table only for strain-level histogram
% Convert to string first (safe for manipulation)
strainInfo.StrainName = string(strainInfo.StrainName);

% Add " UV" to each strain
strainInfo.StrainName = strainInfo.StrainName + " UV";

% Remove any accidental spaces
strainInfo.StrainName = strtrim(strainInfo.StrainName);

% Convert to categorical
strainInfo.StrainName = categorical(strainInfo.StrainName);


strainInfo.Genotype = categorical(cellstr(strainInfo.Genotype));
resultsTable2 = resultsTable(ismember(resultsTable.Strain, strainInfo.StrainName), :);

%% Mean per strain (for histogram)
[G, groupNames] = findgroups(resultsTable2.Strain);
mean_pre  = splitapply(@mean, resultsTable2.Qf_pre, G);
mean_post = splitapply(@mean, resultsTable2.Qf_post, G);

%% Reconstruct N2 UV wormotels from ORIGINAL unfiltered table
dates = unique(resultsTable.Date);

n2_mean_pre_by_wormotel  = [];
n2_mean_post_by_wormotel = [];
n2_block_sizes           = [];
n2_block_date            = strings(0,1);

for d = 1:numel(dates)

    thisDate = dates(d);
    Tdate = resultsTable(resultsTable.Date == thisDate, :);

    % Find N2 UV rows within this date in the ORIGINAL table
    n2_rows = find(Tdate.Strain == "N2 UV");

    if isempty(n2_rows)
        continue;
    end

    % Split into contiguous blocks
    breaks = [1; find(diff(n2_rows) > 1) + 1];
    numBlocks = numel(breaks);

    for b = 1:numBlocks

        startIdx = breaks(b);

        if b < numBlocks
            endIdx = breaks(b+1) - 1;
        else
            endIdx = numel(n2_rows);
        end

        blockRows = n2_rows(startIdx:endIdx);
        blockSize = numel(blockRows);

        n2_block_sizes(end+1,1) = blockSize;
        n2_block_date(end+1,1)  = string(thisDate);

        n2_mean_pre_by_wormotel(end+1,1) = ...
            mean(Tdate.Qf_pre(blockRows), 'omitnan');

        n2_mean_post_by_wormotel(end+1,1) = ...
            mean(Tdate.Qf_post(blockRows), 'omitnan');
    end
end

%% Inspect block sizes
% disp('N2 UV wormotel block sizes:')
% disp(n2_block_sizes)
% 
% disp(table(n2_block_date, n2_block_sizes, ...
%     n2_mean_pre_by_wormotel, n2_mean_post_by_wormotel))

%% Plot
cOrder = colororder;
Color_postUV = cOrder(4,:);

figure;
hold on;
edges= [0:0.03:1];

%% Histogram of N2 wormotel means (foreground)

histogram(n2_mean_pre_by_wormotel,...
    edges, 'Normalization', 'probability', ...
    'FaceColor', [0.8 0.2 0], 'FaceAlpha', 0.4, 'EdgeColor', [0.8 0.2 0], 'EdgeAlpha', 1, ...
    'LineStyle','none')

histogram(n2_mean_post_by_wormotel, ...
    edges, 'Normalization', 'probability', ...
    'FaceColor', 'g', 'FaceAlpha', 0.4, 'EdgeColor', 'g', 'EdgeAlpha', 1, ...
    'LineStyle','none')


%% Histogram across strains (background)
histogram(mean_pre, edges, 'Normalization', 'probability', ...
    'FaceColor', 'k', 'FaceAlpha', 0, 'EdgeColor', 'k', 'EdgeAlpha', 1,  ...
    'LineStyle','--', 'LineWidth', 3);

histogram(mean_post, edges, 'Normalization', 'probability',...
    'FaceColor', Color_postUV, 'FaceAlpha', 0, 'EdgeColor', Color_postUV, 'EdgeAlpha', 1, ...
    'LineStyle','-', 'LineWidth', 3)
%% Add mean lines
yl = ylim;

% plot([mean(n2_mean_pre_by_wormotel,'omitnan') mean(n2_mean_pre_by_wormotel,'omitnan')], ...
%     yl, '--k', 'LineWidth', 2.5);
% 
% plot([mean(n2_mean_post_by_wormotel,'omitnan') mean(n2_mean_post_by_wormotel,'omitnan')], ...
%     yl, '--', 'Color', Color_postUV, 'LineWidth', 2.5);

hold off;

xlabel('Mean quiescence fraction');
ylabel('Relative frequency');

ax = gca;
ax.FontSize = 25;
ax.LineWidth = 2;
pbaspect([1 1 1]);

xticks(0:0.2:1);
xlim([0 1]);

f = gcf;
f.Units = 'pixels';
f.Position = [100 100 800 800];

box off;

exportgraphics(gca, "Fig3A.svg");

%% Test of variance and mean

% test of variance
disp('test of variance for N2 and MMP')
var_test_data_postUV = [mean_post; n2_mean_post_by_wormotel];
groups_data_postUV = [ones(size(mean_post)); 2 * ones(size(n2_mean_post_by_wormotel))];
vartestn(var_test_data_postUV, groups_data_postUV, 'TestType', 'BrownForsythe');

var_test_data_preUV = [mean_pre; n2_mean_pre_by_wormotel];
groups_data_preUV = [ones(size(mean_pre)); 2 * ones(size(n2_mean_pre_by_wormotel))];
vartestn(var_test_data_preUV, groups_data_preUV, 'TestType', 'BrownForsythe');

% test of mean
disp('test of mean for N2 and MMP')
p = ranksum(mean_post, n2_mean_post_by_wormotel, 'tail', 'right')
p = ranksum(mean_pre, n2_mean_pre_by_wormotel, 'tail', 'right')