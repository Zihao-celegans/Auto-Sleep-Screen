%% Script for ploting Figure 1 Sleep Screen Paper

clear all;
close all;
clc;

load('data.mat');
load('robot_MMP_Qf.mat')
load('sleep_data_cleaned.mat')
stressed_quiescence_mean_kerry = stressed_quiescence_mean'/60;
unstressed_quiescence_mean_kerry = unstressed_quiescence_mean'/60;
load('sleep_data_N2_Kerry')
load('robot_N2_Qf_byWM.mat')

%% remove duplicate strains in Kerry's dataset

[unique_MMP_Names_kerry, ~, idx] = unique(lower(strain_names'), 'stable');
unique_stressed_quiescence_mean_kerry = accumarray(idx(:), stressed_quiescence_mean_kerry(:), [], @mean);
unique_unstressed_quiescence_mean_kerry = accumarray(idx(:), unstressed_quiescence_mean_kerry(:), [], @mean);

%% ratio for pbaspect ratio
ratio = 1;

%% === Fig1B ===
close all;

gpLab = strings(1,141);                 % 1×48 of ""
gpLab([24 72 122]) = ["N2","\it aptf-1(tm3287)","\it ceh-17(np1)"];  % set only those ticks
% gpLab = strings(1,45);                 % 1×48 of ""
% gpLab([7 23 39]) = ["N2","\it aptf-1","\it ceh-17"];  % set only those ticks

figure
qPRE = datPRE';
qPOST = datPOST';
qALL=[qPRE(1:10:end-1,:); qPOST(1:10:end,:)];
imagesc(qALL') %make heatmap showing fracQ data
title('Activity')
xticks(1:36:size(qALL,1))         % ticks every hour
xticklabels(0:4)
xlabel('Time (hour)')
% ylim([0.5,143])
% yticks([1, size(qALL,2)])
yticks([])
% yticklabels()
ax = gca;
ax.FontSize = 30;

% ytickangle(90)
colorbar
colormap hot

% === draw 3 white lines every 16 rows ===
hold on
ax.Layer = 'top';                  % keep lines on top of image
nRows = size(qALL,2);              % = 48 (rows after transpose)
blk   = 48;
edges = blk:blk:nRows-1;             % [16 32 48]
for k = edges
    y = k + 0.5;                   % boundary between row k and k+1
    if y <= nRows + 0.5
        plot(xlim, [y y], 'w-', 'LineWidth', 5);
    end
end

hold off

pbaspect([ratio 1 1]);

f = gcf;
f.Units = 'pixels';
f.Position = [100 100 800 800/ratio];

exportgraphics(gca, "Fig1B.svg");

%% === Fig1C ===
close all;
figure
[qPRE,t3PRE,p3PRE]=fracQ(datPRE',1,10,60);
[qPOST,t3POST,p3POST]=fracQ(datPOST',1,10,60);
% qALL=[qPRE(600:10:end,:); qPOST(600:10:end,:)];
qALL=[qPRE(1:10:end,:); qPOST(1:10:end,:)];

imagesc(qALL')                     % heatmap
title('Quiescence fraction')
xticks(1:360:size(qALL,1))         % ticks every hour
xticklabels(0:4)
xlabel('Time (hour)')
% ylim([0.5,143]);
% yticks([1, size(qALL,2)])
yticks([])
% yticks(1:143)
ax = gca; ax.FontSize = 30;
% yticklabels(gpLab); 
% ytickangle(90)
colorbar; colormap jet

% === draw 3 white lines every 16 rows ===
hold on
ax.Layer = 'top';                  % keep lines on top of image
nRows = size(qALL,2);              % = 48 (rows after transpose)
blk   = 48;
edges = blk:blk:nRows-1;             % [16 32 48]
for k = edges
    y = k + 0.5;                   % boundary between row k and k+1
    if y <= nRows + 0.5
        plot(xlim, [y y], 'w-', 'LineWidth', 5);
    end
end
hold off

pbaspect([ratio 1 1]);

f = gcf;
f.Units = 'pixels';
f.Position = [100 100 800 800/ratio];

exportgraphics(gca, "Fig1C.svg");

%% === Fig1D ===
close all;

% Scatterplot of total quiescence fraction

%ENTER wells to be censored
% cens=[1, 16, 23];
cens = [14, 17, 32, 34, 48, 50, 53, 61, 71, 108, 109, 131, 144]; % N2, N2, N2, aptf-1, aptf-1, aptf-1, ceh-17, ceh-17, ceh-17

%ENTER THE HOW MANY STRAINS
N_STRAINS = 3;

%condition names will be entered based on strains provided above
%(can be copy/pasted directly to prism
gpNames=cell(N_STRAINS,1);
% gpNames{1}=["D3-N2"];
% gpNames{2}=["D3-HBR232"];
% gpNames{3}=["D3-IB16"];
gpNames{1}=["N2"];
gpNames{2}=["\it aptf-1(tm3287)"];
gpNames{3}=["\it ceh-17(np1)"];

%ENTER wells for corresponding conditions
%list which animals go in which group

% a vector store the begin/end well number for each group
% the length of the vector must be an even number, since there are two
% number for each group, begin well number and end well number.
% the pattern: [begin gp 1, end gp 1; begin gp 2, end gp 2; ...]
% well_num_begin_end = [1,16;17,32;33,48];
well_num_begin_end = [1,48;49,96;97,144];

gp = cell(1, N_STRAINS);
gp{1}=[well_num_begin_end(1,1) : well_num_begin_end(1,2)]; %STR A
gp{2}=[well_num_begin_end(2,1) : well_num_begin_end(2,2)]; %STR B
gp{3}=[well_num_begin_end(3,1) : well_num_begin_end(3,2)]; %STR C
gpN=length(gp); %calculate number of groups
gpLens=zeros(1,gpN); %N per group (filled below)
nOrig=zeros(1,gpN); %number of animals in original groups
nCens=zeros(1,gpN); %number of animals after removing censored animals
for i = 1:length(gp)
    nOrig(i)=length(gp{i});
    
    gp{i}=setdiff(gp{i},cens); %remove censored wells from groups
    gpLens(i)=length(gp{i}); %get number of animals in each group
    
    nCens(i)=nOrig(i)-gpLens(i);
end
clear i
gpMax=max(gpLens); %find largest group
TimepImg = 10;
EPM= 60 / TimepImg; %epochs per min (1 image/10 seconds)
groupSz=48; %size of groups for input into prism (standard size makes it easier to copy onto preexisting files)

%out3 qTotal - sum of quiescence over all 4 hours
%use prism grouped input for 2-way anova design
%so that genotype and treatment are matched appropriately across groups

%same set up as out1, except sum over 4 hours

ATp=0.01; %activity threshold: fraction of 95th percentile of activity
QT=0.25; %qui threshold: fraction of moving window needed to be quiescent
qWin=60; %qWin is moving window size (in number of frames)
binTime = 60; % binTime time duration (min) of each bin window
binSz = binTime * 4 * EPM; %binSize here refers to number of frames in each bin

[qTime,~,~,qTimePRE,~,durPre] = fracQbin230817(datPRE',datPOST',gp,ATp,QT,qWin,binSz,EPM);

nCol=length(gp); %number of groups of columns in output

out3_qTot=NaN(2,groupSz*nCol);
for iCols = 1:nCol
    colSt=(iCols-1) * groupSz + 1;
    
    
    quiPRE=qTimePRE{iCols} * (EPM/binSz);
    quiPOST=qTime{iCols} * (EPM/binSz);
    colEn=colSt+gpLens(iCols)-1;
    
    out3_qTot(1,colSt:colEn)=quiPRE;
    out3_qTot(2,colSt:colEn)=quiPOST;
    
end

figure

% flatten the data for beeswarm plot

% category vector w/o merge
cate_vec = [];

% pre-UV
QfTot_vec_preUV = [];

% post-UV
QfTot_vec_postUV = [];


% all
cate_vec_all = [];
QfTot_vec_all = [];

% color for plotting different groups
cmap = [];
cOrder = colororder;
Color_preUV = [0 0 0];
Color_postUV = cOrder(4,:);

for i = 1 : size(gp, 2)
    
    % data frame
    gp_QfTot = out3_qTot(:, well_num_begin_end(i, 1): well_num_begin_end(i, 2));

    % remove NaN
    gp_QfTot_preUV = rmmissing(gp_QfTot(1,:));
    gp_QfTot_postUV = rmmissing(gp_QfTot(2,:));
    
    % category vector
    cate_vec = [cate_vec, i * ones(1, size(gp_QfTot_postUV, 2))];

    % pre-UV
    QfTot_vec_preUV = [QfTot_vec_preUV, gp_QfTot_preUV];

    
    % add pre-UV to merged vector
    cate_vec_all = [cate_vec_all, (2 * i - 1) * ones(1, size(gp_QfTot_preUV, 2))];
    QfTot_vec_all = [QfTot_vec_all, gp_QfTot_preUV];

    cmap = [cmap; Color_preUV]; % plot color for pre-UV
    

    % post-UV
    QfTot_vec_postUV = [QfTot_vec_postUV, gp_QfTot_postUV];
    
    % add post-UV to merged vector
    QfTot_vec_all = [QfTot_vec_all, gp_QfTot_postUV];
    cate_vec_all = [cate_vec_all, (2 * i) * ones(1, size(gp_QfTot_postUV, 2))];

    cmap = [cmap; Color_postUV]; % use a different plot color with pre-UV

end

dot_sz = 3;

xbee = beeswarm(cate_vec_all',QfTot_vec_all', ...
    'dot_size', dot_sz, ...
    'overlay_style','sd',...
    'corral_style','random',...
    'colormap',cmap, ...
    'MarkerFaceAlpha', 0.5, ...
    'MarkerEdgeColor', 'w');


N2_pre_UV = QfTot_vec_all(cate_vec_all == 1);
N2_post_UV = QfTot_vec_all(cate_vec_all == 2);
aptf_1_pre_UV = QfTot_vec_all(cate_vec_all == 3);
aptf_1_post_UV = QfTot_vec_all(cate_vec_all == 4);
ceh_17_pre_UV = QfTot_vec_all(cate_vec_all == 5);
ceh_17_post_UV = QfTot_vec_all(cate_vec_all == 6);

fprintf('N2 vs aptf-1\n')
ranksum(N2_pre_UV, aptf_1_pre_UV, 'tail', 'right')
ranksum(N2_post_UV, aptf_1_post_UV, 'tail', 'right')

fprintf('N2 vs ceh-17\n')
ranksum(N2_pre_UV, ceh_17_pre_UV, 'tail', 'both')
ranksum(N2_post_UV, ceh_17_post_UV, 'tail', 'right')

ylim([0, 1]);
xlim([0, 2 * size(gp, 2)+0.8]);

ylabel('Quiescence fraction')

xticks(1.5 : 2 : 2 * size(gp, 2));
xticklabels(gpNames)

y_leg_preUV = 0.90;
y_leg_postUV = 0.83;

% h1 = text(2 * size(gp, 2) + 0.7, y_leg_preUV,'Pre-UV','Color',Color_preUV);
% h1.FontSize = 20;
% h2 = text(2 * size(gp, 2) + 0.7, y_leg_postUV,'Post-UV','Color',Color_postUV);
% h2.FontSize = 20;

% draw the legend
hold on
% scatter(2 * size(gp, 2) + 0.4, y_leg_preUV, dot_sz*36, 'filled', ...
%     'MarkerFaceColor', Color_preUV, 'MarkerFaceAlpha', 0.2, 'MarkerEdgeColor', 'k')
% scatter(2 * size(gp, 2) + 0.4, y_leg_postUV, dot_sz*36, 'filled', ...
%     'MarkerFaceColor', Color_postUV, 'MarkerFaceAlpha', 0.2, 'MarkerEdgeColor', 'k')
% rectangle('Position',[2 * size(gp, 2) + 0.2 0.84 1.4 0.15])
hold off

ax = gca; ax.FontSize = 33;
ax.LineWidth = 2;
xtickangle(20)

pbaspect([ratio 1 1]);

f = gcf;
f.Units = 'pixels';
f.Position = [100 100 800 800/ratio];

exportgraphics(gca, "Fig1D.svg");

%% === Fig1E : Plot the histogram of the mean quiescence for MMPs ===

close all;

% calculate the pre and post-UV for N2
N_N2_kerry_preUV = size(unstressed_quiescence_mean_N2_kerry, 2) * 8;
N_N2_kerry_postUV = size(stressed_quiescence_mean_N2_kerry, 2) * 8;

N_N2_john_preUV = size(Qf_all_ctrl_preUV, 1);
N_N2_john_postUV = size(Qf_all_ctrl, 1);

% Calculate the mean for pool N2s
meanQf_N2_preUV_pool = (N_N2_kerry_preUV * mean(unstressed_quiescence_mean_N2_kerry) + ...
    N_N2_john_preUV * mean(Qf_all_ctrl_preUV))/(N_N2_kerry_preUV + N_N2_john_preUV);

meanQf_N2_postUV_pool = (N_N2_kerry_postUV * mean(stressed_quiescence_mean_N2_kerry) + ...
    N_N2_john_postUV * mean(Qf_all_ctrl))/(N_N2_kerry_postUV + N_N2_john_postUV);


% sperate N2
% N2_idx_kerry = find(strcmpi(strain_names_kerry, 'N2'));
% N2_unstressed_quiescence_mean_kerry = unstressed_quiescence_mean_kerry(N2_idx_kerry);
% N2_stressed_quiescence_mean_kerry = stressed_quiescence_mean_kerry(N2_idx_kerry);
% 
% unstressed_quiescence_mean_kerry_N2removed = unstressed_quiescence_mean_kerry([1:N2_idx_kerry-1, N2_idx_kerry+1:end]);
% stressed_quiescence_mean_kerry_N2removed = stressed_quiescence_mean_kerry([1:N2_idx_kerry-1, N2_idx_kerry+1:end]);

merge_mean_Qf_preUV = [unique_unstressed_quiescence_mean_kerry;ranked_mean_Qf_preUV];
merge_mean_Qf_postUV = [unique_stressed_quiescence_mean_kerry;ranked_mean_Qf];


% remove strain having both 0 pre and post UV
idx_to_remove = (merge_mean_Qf_preUV == 0) .* (merge_mean_Qf_postUV == 0);

merge_mean_Qf_preUV_filtered = merge_mean_Qf_preUV(~idx_to_remove);
merge_mean_Qf_postUV_filtered = merge_mean_Qf_postUV(~idx_to_remove);

% Histogram of the mean quiescence for MMPs
edges= [0:0.03:1];
figure();
hold on;

N2_postUV_Qf_all = [stressed_quiescence_mean_N2_kerry'; robot_N2_Qf_postUV_byWM'];
histogram(N2_postUV_Qf_all, ...
    edges, 'Normalization', 'probability', ...
    'FaceColor', 'g', 'FaceAlpha', 0.4, 'EdgeColor', 'g', 'EdgeAlpha', 1, ...
    'LineStyle','none')


N2_preUV_Qf_all = [unstressed_quiescence_mean_N2_kerry'; robot_N2_Qf_preUV_byWM'];
histogram(N2_preUV_Qf_all, ...
    edges, 'Normalization', 'probability', ...
    'FaceColor', [0.8 0.2 0], 'FaceAlpha', 0.4, 'EdgeColor', [0.8 0.2 0], 'EdgeAlpha', 1, ...
    'LineStyle','none')


histogram(merge_mean_Qf_preUV_filtered, edges, 'Normalization', 'probability', ...
    'FaceColor', 'k', 'FaceAlpha', 0, 'EdgeColor', 'k', 'EdgeAlpha', 1,  ...
    'LineStyle','--', 'LineWidth', 3)
histogram(merge_mean_Qf_postUV_filtered, edges, 'Normalization', 'probability',...
    'FaceColor', Color_postUV, 'FaceAlpha', 0, 'EdgeColor', Color_postUV, 'EdgeAlpha', 1, ...
    'LineStyle','-', 'LineWidth', 3)


% plot([meanQf_N2_preUV_pool meanQf_N2_preUV_pool], ylim, '--k', 'LineWidth', 2)
% 
% plot([meanQf_N2_postUV_pool meanQf_N2_postUV_pool], ylim, '--', 'Color', Color_postUV, 'LineWidth', 2)

hold off;



xlabel('Mean quiescence fraction')
ylabel('Relative frequency')

ax = gca; ax.FontSize = 25;
ax.LineWidth = 2;
pbaspect([1 1 1]);

xticks(0 : 0.2 : 1);

xlim([0 1])

f = gcf;
f.Units = 'pixels';
f.Position = [100 100 800 800];
box off;

exportgraphics(gca, "Fig2A.svg");

% test of variance
disp('test of variance for N2 and MMP')
var_test_data_postUV = [merge_mean_Qf_postUV_filtered; N2_postUV_Qf_all];
groups_data_postUV = [ones(size(merge_mean_Qf_postUV_filtered)); 2 * ones(size(N2_postUV_Qf_all))];
vartestn(var_test_data_postUV, groups_data_postUV, 'TestType', 'BrownForsythe');

var_test_data_preUV = [merge_mean_Qf_preUV_filtered; N2_preUV_Qf_all];
groups_data_preUV = [ones(size(merge_mean_Qf_preUV_filtered)); 2 * ones(size(N2_preUV_Qf_all))];
vartestn(var_test_data_preUV, groups_data_preUV, 'TestType', 'BrownForsythe');

% test of mean
disp('test of mean for N2 and MMP')
p = ranksum(merge_mean_Qf_postUV_filtered,N2_postUV_Qf_all, 'tail', 'right')
p = ranksum(merge_mean_Qf_preUV_filtered,N2_preUV_Qf_all, 'tail', 'right')

% disp(strcat(strcat(" Population Mean of mean quiescence = ", num2str(mean(ranked_mean_Qf)))));
% disp(strcat(strcat(" Population SD of mean quiescence = ", num2str(std(ranked_mean_Qf)))));
% disp(strcat(strcat(" Population Mean - 2 * SD = ", num2str(mean(ranked_mean_Qf) - 2 * std(ranked_mean_Qf)))));

% Scatter plot of the mean quiescence for MMPs
figure();
hold on;
scatter(merge_mean_Qf_preUV_filtered, merge_mean_Qf_postUV_filtered, ...
    'filled', 'MarkerFaceColor', 0.5*Color_postUV + 0.5*[0 0 0], 'MarkerFaceAlpha', 0.6);
scatter(N2_preUV_Qf_all, N2_postUV_Qf_all, 'filled', 'MarkerFaceColor', 'g', 'MarkerFaceAlpha', 0.2);

xlabel('Mean pre-UV quiescence fraction')
ylabel('Mean post-UV quiescence fraction')
ax = gca; ax.FontSize = 25;
ax.LineWidth = 2;
xlim([0 1])
ylim([0 1])

xl = xlim;
plot(xl, xl, '--k', 'LineWidth', 2)

f = gcf;
f.Units = 'pixels';
f.Position = [100 100 800 800];
pbaspect([1 1 1]);

hold off
exportgraphics(gca, "Fig2B.svg");


%% === Fig1F : Dot Plot — Quiescence ===

% === Aggregate per strain ===
[G, strainList] = findgroups(Strain_names_all);

% --- Quiescence ---
preUV_Qf_mean = splitapply(@mean, preUV_Qf_all, G);
postUV_Qf_mean = splitapply(@mean, postUV_Qf_all, G);
SE_pre_Qf_mean = splitapply(@(x) std(x)/sqrt(numel(x)), preUV_Qf_all, G);
SE_post_Qf_mean = splitapply(@(x) std(x)/sqrt(numel(x)), postUV_Qf_all, G);

% === Identify N2 and reorder so it plots last ===
isN2 = contains(string(strainList), 'N2', 'IgnoreCase', true);
order = [find(~isN2); find(isN2)];

strainList = strainList(order);

% Quiescence
preUV_Qf_mean  = preUV_Qf_mean(order);
postUV_Qf_mean = postUV_Qf_mean(order);
SE_pre_Qf_mean  = SE_pre_Qf_mean(order);
SE_post_Qf_mean = SE_post_Qf_mean(order);


% === Define gene-group colors ===
geneColors = containers.Map( ...
    {'N2', 'strd-1', 'cla-1', 'egl-8','aptf-1', 'ceh-17'}, ...
    {[0.5 0.5 0.5], ...   % N2  (grey)
     [0.70, 0.00, 0.70], ...   % strd-1 (purple)
     [0.20, 0.90, 0.20], ...   % cla-1 (green)
     [1.00, 0.00, 0.00], ...   % egl-8 (red)
     [1.00, 0.50, 0.00], ...   % aptf-1 (orange)
     [0.15, 0.50, 0.85]} ...   % ceh-17 (blue)
);

% === Define shapes per gene ===
geneShapes = containers.Map( ...
    {'N2', 'strd-1','cla-1', 'egl-8', 'aptf-1', 'ceh-17'}, ...
    {'p', 's', '^', 'd', 'v', 'o'} ...
);

% === Strain → Gene mapping ===
strainToGene = containers.Map( ...
    {'N2 UV', ...
     'RB1771 UV', 'RB1775 UV', ...
     'VC387 UV', 'RB778 UV', 'RB1777 UV', 'VC596 UV', 'VC631 UV', ...
     'CB6614 UV', 'RM2221 UV', 'IK777 UV', 'JT47 UV', 'MT1083 UV', ...
     'HBR232 UV', ...
     'IB16 UV'}, ...
    {'N2', ...
     'strd-1', 'strd-1', ...
     'cla-1', 'cla-1', 'cla-1', 'cla-1', 'cla-1', ...
     'egl-8', 'egl-8', 'egl-8', 'egl-8', 'egl-8', ...
     'aptf-1', ...
     'ceh-17'} ...
);
%     'cysl-1', 'cysl-1', 'cysl-1', 'cysl-1', ...
%     'AR10 UV', 'AR4 UV', 'AR11 UV', 'RB899 UV', ...

%     'Y37A1A.4/C27H2.2', 'Y37A1A.4/C27H2.2', ...
%     'RB1723 UV', 'tm10921 UV', ...

%, 'RB1012 UV'
%, 'egl-8'
% === Appearance ===
defaultColor   = [0.7 0.7 0.7];
baseMarkerSize = 35;
geneMarkerSize = 70;

% === Legend containers ===
legendEntries = [];
legendLabels  = {};

figure('Color','w'); hold on; box on;

% ============================================================
%   STRAIN ORDER: plot non-N2 first, N2 last
% ============================================================
strainList = categories(Strain_names_all);

isN2 = false(size(strainList));

for i = 1:numel(strainList)
    sName = strainList{i};

    if isKey(strainToGene, sName) && strcmp(strainToGene(sName), 'N2')
        isN2(i) = true;
    end
end

% N2 first → plotted first → appears underneath
strainList = [strainList(isN2); strainList(~isN2)];

% === legend containers ===
legendEntries = gobjects(0);
legendLabels  = {};

% ellipse smoothness
theta = linspace(0, 2*pi, 120);

% ============================================================
%   LOOP — one SEM ellipse per wormmotel trial directory
% ============================================================

N_each_trial = [];
for s = 1:numel(strainList)

    strainName = strainList{s};

    % all worms for this strain
    idxStrain = Strain_names_all == strainName;

    x_all   = preUV_Qf_all(idxStrain);
    y_all   = postUV_Qf_all(idxStrain);
    dir_all = trialDir_all(idxStrain);

    if isempty(x_all)
        continue;
    end

    % gene mapping
    if isKey(strainToGene, strainName)
        geneName  = strainToGene(strainName);
        geneColor = geneColors(geneName);
    else
        continue;
    end

    % unique wormmotel trial folders for this strain
    uniqueDirs = categories(removecats(dir_all));

    for k = 1:numel(uniqueDirs)

        thisDir = categorical(uniqueDirs(k));
        idxDir  = dir_all == thisDir;

        xDir = x_all(idxDir);
        yDir = y_all(idxDir);

        valid = ~isnan(xDir) & ~isnan(yDir);
        xDir = xDir(valid);
        yDir = yDir(valid);

        if isempty(xDir)
            continue;
        end

        % mean per trial
        xMean = mean(xDir);
        yMean = mean(yDir);
        
        % SEM per trial
        if numel(xDir) > 1
            N_each_trial(end+1) = numel(xDir);
            xSEM = std(xDir) / sqrt(numel(xDir));
        else
            xSEM = 0;
        end
        
        if numel(yDir) > 1
            ySEM = std(yDir) / sqrt(numel(yDir));
        else
            ySEM = 0;
        end
        
        % --- vertical errorbar ---
        h = errorbar(xMean, yMean, ySEM, ...
            'Color', geneColor, ...
            'LineWidth', 1.5, ...
            'CapSize', 0);
        
        % --- horizontal errorbar ---
        errorbar(xMean, yMean, xSEM, ...
            'horizontal', ...
            'Color', geneColor, ...
            'LineWidth', 1.5, ...
            'CapSize', 0);
        
        % small center point
        plot(xMean, yMean, '.', ...
            'Color', geneColor, ...
            'MarkerSize', 14);

        % add legend once per gene
        if ~any(strcmp(legendLabels, geneName))
            legendEntries(end+1) = h;
            legendLabels{end+1}  = geneName;
        end
    end
end

fprintf('N max = %d\n', max(N_each_trial));
fprintf('N min = %d\n', min(N_each_trial));
fprintf('N trial = %d\n', size(N_each_trial, 2));

% ============================================================
%   AXES / DIAGONAL / STYLE
% ============================================================
allVals = [preUV_Qf_all(:); postUV_Qf_all(:)];
allVals = allVals(~isnan(allVals));

lims = [min(allVals), max(allVals)];

% % add a little padding
% pad = 0.03 * range(lims);
% if pad == 0
%     pad = 0.01;
% end
% lims = [lims(1)-pad, lims(2)+pad];

plot(lims, lims, 'k--', 'LineWidth', 2);

% axis equal;
% xlim(lims);
xlim([0 0.4]);
ylim(lims);

xlabel('Pre-UV quiescence fraction');
ylabel('Post-UV quiescence fraction');
% title('Pre vs Post-UV Quiescence');

set(gca, ...
    'FontSize', 11, ...
    'LineWidth', 1, ...
    'TickDir', 'in', ...
    'Layer', 'top');

grid off;

% legend(legendEntries, legendLabels, ...
%     'Location', 'eastoutside', ...
%     'Box', 'off');

f = gcf;
f.Units = 'pixels';
f.Position = [100 100 800 800];
pbaspect([1 1 1]);


ax = gca; ax.FontSize = 25;
ax.LineWidth = 2;

hold off;

box off;
exportgraphics(gca, "Fig3B.svg");

%% Output a list of MMP strain name

% Convert both arrays to string arrays
names1 = string(ranked_MMP_name_list(:));      % categorical -> string column
names2 = string(strain_names(:));     % cell -> string column

% Drop " UV" from categorical-array strain names only
names1 = erase(names1, " UV");

% Merge
allNames = [names1; names2];

% Clean formatting
allNames = strtrim(allNames);        % remove leading/trailing spaces
allNames = upper(allNames);          % capitalize everything

% Optional: only use this if apostrophes are actually part of the cell text
allNames = erase(allNames, "'");

% Remove missing/empty entries
allNames = allNames(~ismissing(allNames));
allNames = allNames(allNames ~= "");

% Remove a specific strain name
allNames = allNames(allNames ~= "VC40165");

% Remove duplicates, preserving original order
allNames = unique(allNames, "stable");

% Create extra columns
Genotype = repmat("MMP", size(allNames));

Identifier = "https://cgc.umn.edu/strain/" + allNames;

Source = repmat("Caenorhabditis Genetics Center", size(allNames));

% Create table
T_mmp_info = table(allNames, Genotype, Identifier, Source, ...
    'VariableNames', {'StrainName', 'Genotype', 'Identifier', 'Source'});

% Write to Excel
f_tableS1 = "C:\Users\lizih\Dropbox\Publication\SleepScreen\Supplemental\Table_S1.xlsx";
% writetable(T_mmp_info, "MMP_strain_info.xlsx");

sheetName = "MMP_screen";

writetable(T_mmp_info, f_tableS1, ...
    'Sheet', sheetName, ...
    'WriteMode', 'overwritesheet');

%%
%fracQbin function (normalizes post-treatment data to pre-treatment data)
%identifies periods of quiescence based on two criteria:
%for a given epoch to be considered quiescent it must meet the following
%criteria:
%(1) ActVal must be below ACTIVITY THRESHOLD
%(2) fracQ must be above QUIESCENCE THRESHOLD
%in this script the activity threshold is established relative to an
%untreated control group. the activity threshold is 1% of the 95th
%percentile of ActVal for control animals (median of 95th percentiles from each animal)
%fracQ is the fraction of epochs in a 10 minute moving window
%around a given epoch (5 mins before/5mins after), where ActVal falls below
%the activity threshold
%after calculating fracQ, this fxn will sum up the time spent quiescent
%into bins, as well as average quiescence (normalized)

%inputs:
%pre and post (MxN matrices): Activity Matrices w/ #pixels changed in
%subtracted frames (from before/after UV respectively)
%%M (rows) are successive time frames, N (cols) are wells
%%ActVal'

%gp (cell array): each cell contains vector of well numbers pertaining to
%%each group

%ATp (scalar): activity threshold percentage. this is the percentage of the
%%'max' activity in untreated groups that will be used as activity threshold
%%note ATp is a fraction (if you want 1%, then ATp=0.001)

%QT (scalar): quiescence threshold. this is the percentage of time in a
%%moving window where an animal must be inactive to be considered quiescent
%%0.25 is a good staring point

%qWin (scalar): quiescence window. this is the size of the moving average
%%window used to calculate fracQ at each point
%%use number of epochs, eg for 10min window w/ images every 10secs, qWin=60 

%binSz (scalar): size of bins (in # of epochs) to add up qui and act into

%EPM (scalar): epochs per minute (eg 1 image per 10sec = 6 epochs per min)

%outputs:
%qTime (cell array): each cell contains matrix of time spent quiescent
%%values are minutes
%%rows are binned time points and cols are wells

%ActNorm (cell array): each cell contains activity values divided by the
%%normalization factor (95th percentile of control activity)
%%rows are time points, cols are wells

%ActAve (cell array): average ActNorm value for each bin

function [qTime,ActNorm,ActAve,qTimePRE,ActAvePRE,durPre] = fracQbin230817(pre,post,gp,ATp,QT,qWin,binSz,EPM)

nGp=length(gp); %get number of groups used for analysis

mRowPRE=size(pre,1); %get size of pre data
mRow=size(post,1); %get size of post data

%if qWin is an even number, add 1
if mod(qWin,2)>0
    qWin=qWin+1;
end

%since we are using a moving window average
%there will be a period of time at the beginning and end of the recording
%that is not covered by a sufficiently large time window to give the data
winLag=floor(qWin/2); %find the lagging window size (number of epochs before or after)

stRanPRE=winLag + 1; %find the first position covered by moving window
enRanPRE=mRowPRE - winLag; %find the last position covered by moving window
RanPRE=[stRanPRE:enRanPRE]; %find all positions covered by moving window
bigRanPRE=[stRanPRE-winLag:enRanPRE+winLag]; %find all positions used for moving mean
durPre=length(RanPRE) * (1/EPM);

stRan=winLag + 1; %find the first position covered by moving window
enRan=mRow - winLag; %find the last position covered by moving window
Ran=[stRan:enRan]; %find all positions covered by moving window
bigRan=[stRan-winLag:enRan+winLag]; %find all positions used for moving mean

%find normalization-factor (NF)
%median of 95th percentile of ActVal in preData for each group
NF=zeros(1,nGp);
for i = 1:nGp
    wellsI=gp{i};
    pct95i=prctile(pre(:,wellsI)',95);
    NF(i)=median(pct95i,'omitnan');
end
clear i gpI wellsI pct95i

%make a matrix for each group that is normalized to its baseline
ActNormPre=cell(1,nGp);
QuiNormPre=cell(1,nGp);
fracQpre=cell(1,nGp);
qOutPre=cell(1,nGp);
ActNorm=cell(1,nGp);
QuiNorm=cell(1,nGp);
fracQ=cell(1,nGp);
qOut=cell(1,nGp);
for j = 1:nGp
    NFj=NF(j); %get the normalization factor for this group
    wellsJ=gp{j}; %get the wells for this group
    
    %make binary for inactivity based on activity threshold for this group
    ATj = NFj * ATp; %activity threshold specific to this group

    %pre-treatment analysis
    ActJpre=pre(bigRanPRE,wellsJ); %get actval for these wells in the time range
    ActNormPre{j}=ActJpre(RanPRE,:) * (1/NFj); %divide all activity values by norm factor
    QuiNormPre{j}=ActNormPre{j} <= ATp; %find all times where normalized activity is below activity threshold
    QuiJpre=ActJpre; 
    QuiJpre=QuiJpre <= ATj; %QuiJ is essentially the same as QuiNorm (all values less than ATp percent of normalization factor)
    fracQjPRE=movmean(QuiJpre,qWin); %fracQ is a moving average of time spent quiescent
    fracQpre{j}=fracQjPRE(RanPRE,:);
    qnJpre=QuiNormPre{j};
    frJpre=fracQpre{j} >= QT;
    qOutPre{j}= qnJpre .* frJpre; %qui = times when: qnJpre=1 (activity below activity threshold); AND frJpre=1 (fracQ above quiescence threshold)
    
    %post-treatment analysis
    ActJ=post(bigRan,wellsJ); %get actval for these wells in the time range
    ActNorm{j}=ActJ(Ran,:) * (1/NFj); %divide all activity values by norm factor
    QuiNorm{j}=ActNorm{j} <= ATp; %find all times where normalized activity is below activity threshold
    QuiJ=ActJ;
    QuiJ=QuiJ <= ATj; %QuiJ is essentially the same as QuiNorm (all values less than ATp percent of normalization factor)
    fracQj=movmean(QuiJ,qWin); %fracQ is a moving average of time spent quiescent
    fracQ{j}=fracQj(Ran,:); %take only range covered by window size
    qnJ=QuiNorm{j}; %=1 if activity is below ATp 
    frJ=fracQ{j} >= QT; %=1 if fracQ is above QT
    qOut{j} = qnJ .* frJ; %multiply (pairwise) so that only QUI if both <ATp and >QT
end
clear j

%sum quiescent time and activity according to binSize
RanLen=length(Ran); %get number of epochs covered by range
Nbins = floor(RanLen/binSz); %find number of bins to cut this range into
binStarts=[1:binSz:Nbins*binSz+1]; %find start times for each bin
binEnds=binStarts + (binSz-1); %find end times for each bin
qTime=cell(1,nGp); %set output variable for total quiescent time
ActAve=cell(1,nGp); %set output variable for activity in each bin
%for pre-treatment, sum up the entire time frame (minus the lag window)
%then transform the time frame so that it equals bin size of post-treatment
RanLenPRE=length(RanPRE);
qTimePRE=cell(1,nGp);
ActAvePRE=cell(1,nGp);
for i = 1:nGp %loop for each group
    quiI=qOut{i}; %get qui data for this group
    actI=ActNorm{i}; %get act data for this group
    
    %pre
    qTimePRE{i}=sum(qOutPre{i},'omitnan') * (1/EPM); %convert qui time in minutes
    qTimePRE{i}=qTimePRE{i} * (binSz/RanLenPRE); %scale qui-time total to per hour
    ActAvePRE{i}=sum(ActNormPre{i},'omitnan') * (1/RanLenPRE); %average quiescent value over all epochs
    
    %post
    for j = 1:Nbins
        ranJ=binStarts(j):binEnds(j);
        qTime{i}(j,:)=sum(quiI(ranJ,:),'omitnan') * (1/EPM);
        ActAve{i}(j,:)=sum(actI(ranJ,:),'omitnan') * (1/binSz);
    end
    
end


end

%%
%turn matrix from XY data into prism input (eg timecourse from fracQbin)

%cellArray: cell array of MxN matrices, M (subjects), N (time points),
%could be output from fracQbin
%N (scalar): maximum size of group in output
%ord (vector): order to arrange output(correspond to indices of cellArray)
%X (vector): x-values corresponding to rows of matrices in cellArray

%prismOut (matrix): input for prism XY graph, each group is N-members long
%with NaN filling in blank spaces. first column is X-vector

function prismOut=prismXY0915(cellArray,N,ord,X)

nGps=length(ord); %how many groups are in new matrix
nRows=length(X); %how many rows are in new matrix
newMat=NaN(nRows,nGps*N); %make new dummy matrix
for i = 1:nGps
    gpI=ord(i); %get group number for this loop
    matI=cellArray{gpI}; %get matrix for this loop
    [mI,nI]=size(matI);
    colSt=(i-1)*N + 1; %get first column for input
    colEn=colSt + nI - 1; %get last column for input
    newMat(:,colSt:colEn)=matI(1:nRows,:);
end

if nRows==size(X,1) %if X is row vector, transpose it
    prismOut=[X,newMat];
else
    prismOut=[X',newMat];
end

end

%%
%fracQ.m 
% Calculates the fraction of quiescence from activity data in p. p is the transpose of the activity trace.  qt is the
% number of pixels required to move to call the frame not-quiescent.  1 is
% a good value.  spi is the number of seconds each image represents.  ps is
% the moving average size. ps must be odd.  If it is not, the program
% automatically adds 1 to ps.

% fracQ returns PAV, a fraction of quiescence vector equal in size to p,
% and t1 which is the time vector associated with PAV.
function [PAV t1 pnew] = fracQ(p,qt,spi,ps)

if mod(ps,2)==0
    ps=ps+1;
end

pq=p<qt;
[sizex sizey]=size(p);
t1=[1:1/spi:sizex]/3600*spi;
rois=1:sizey;
pstart=ceil(ps/2);
pf=floor(ps/2);

pav=zeros(ceil(sizex/ps),sizey);
for i=pstart:ps:sizex-ps
    for j=1:sizey
        %pav((i+pf)/ps,j)=nanmean(pq((i-pf):(i+pf),j));
        pav((i+pf)/ps,j)=mean(pq((i-pf):(i+pf),j), "omitnan");
    end
end
pax=[1:ceil(sizex/ps)]*ps/3600*spi;
pax2=[1:sizex]/3600*spi;
PAV=interp2(rois,pax',pav,rois,t1');
pnew=interp2(rois,pax2',p,rois,t1')/spi;

end

%%
function x = beeswarm(x,y,varargin)
%function xbee = beeswarm(x,y)
%
% Input arguments:
%   x               column vector of groups (only tested for integer)
%   y               column vector of data
%
% Optional input arguments:
%   sort_style      ('nosort' - default | 'up' | 'down' | 'fan' | 'rand' | 'square' | 'hex')
%   corral_style    ('none' default | 'gutter' | 'omit' | 'rand')
%   dot_size        relative. default=1
%   overlay_style   (false default | 'box' | 'sd' | 'ci')
%   use_current_axes (false default | true)
%   colormap        (lines default | 'jet' | 'parula' | 'r' | Nx3 matrix of RGB values]
%
% Output arguments:
%   xbee            optimized layout positions
%
% Known Issues:
%       x locations depend on figure aspect ratio. resizing the figure window and rerunning may give different results
%       setting corral to 'none' still has a gutter when the width is large
%
% Usage example:
% 	x = round(rand(150,1)*5);
%   y = randn(150,1);
%   beeswarm(x,y,3,'sort_style','up','overlay_style','ci')
%
% % Ian Stevenson, CC-BY 2019

p = inputParser;
addRequired(p,'x')
addRequired(p,'y')
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
addOptional(p,'sort_style','nosort')
addOptional(p,'corral_style','none')
addOptional(p,'dot_size',11/sqrt(length(x)),validScalarPosNum)
addOptional(p,'overlay_style',false)
addOptional(p,'use_current_axes',false)
addOptional(p,'colormap','lines')
addOptional(p,'MarkerFaceColor','')
addOptional(p,'MarkerFaceAlpha',.5)
addOptional(p,'MarkerEdgeColor','none')
parse(p,x,y,varargin{:});

% extra parameters
rwid = .05; % width of overlay box/dash

dcut=8; % spacing factor
nxloc=512; % resolution for optimization
chanwid = .9; % percent width of channel to use
yl = [min(y) max(y)]; % default y-limits
asp_rat = 1;
keep_hold = false;

% get aspect ratio for a figure window
if isfinite(p.Results.dot_size)
    if ~p.Results.use_current_axes
        % make new axes
        s=scatter(x,y);
        xl=[min(x)-.5 max(x)+.5];
    else
        xl=xlim();
    end
    yl=ylim();
    pasp_rat = get(gca,'PlotBoxAspectRatio');
    dasp_rat = get(gca,'DataAspectRatio');
    asp_rat = pasp_rat(1)/pasp_rat(2);
    
    % pix-scale
    pf = get(gcf,'Position');
    pa = get(gca,'Position');
    as = pf(3:4).*pa(3:4); % width and height of panel in pixels
    dcut = dcut*sqrt(p.Results.dot_size)/as(1)*(range(unique(x))+1);
    if ~ishold
        cla
    else
        keep_hold = true;
    end
end

% sort/round y for different plot styles
yorig=y;
switch lower(p.Results.sort_style)
    case 'up'
        [y,sid]=sort(y);
    case 'fan'
        [~,sid]=sort(abs(y-mean(y)));
        sid=[sid(1:2:end); sid(2:2:end)];
        y=y(sid);
    case 'down'
        [y,sid]=sort(y,'descend');
    case 'rand'
        sid=randperm(length(y));
        y=y(sid);
    case 'square'
        nxloc=.9/dcut;
%         [~,e,b]=histcounts(y,ceil((range(x)+1)*chanwid*nxloc/2/asp_rat));
        edges = linspace(min(yl),max(yl),ceil((range(x)+1)*chanwid*nxloc/asp_rat));
        [~,e,b]=histcounts(y,edges);
        y=e(b)'+mean(diff(e))/2;
        [y,sid]=sort(y);
    case 'hex'
        nxloc=.9/dcut;
%         [~,e,b]=histcounts(y,ceil((range(x)+1)*chanwid*nxloc/2/sqrt(1-.5.^2)/asp_rat));
        edges = linspace(min(yl),max(yl),ceil((range(x)+1)*chanwid*nxloc/sqrt(1-.5.^2)/asp_rat));
        [n,e,b]=histcounts(y,edges);
        oddmaj=0;
        if sum(mod(n(1:2:end),2)==1)>sum(mod(n(2:2:end),2)==1),
            oddmaj=1;
        end
        y=e(b)'+mean(diff(e))/2;
        [y,sid]=sort(y);
        b=b(sid);
    otherwise
        sid=1:length(y);
end
x=x(sid);
yorig=yorig(sid);
[ux,~,ic] = unique(x);
% rmult=(range(ux)+1)*2;
rmult=5;

% for each group...
for i=1:length(ux)
    fid = find(ic==i);   
    
    % set of possible x locations
    xi = linspace(-chanwid/2*rmult,chanwid/2*rmult,nxloc*rmult+(mod(nxloc*rmult,2)==0))'+ux(i);

    % rescale y to that things are square visually
    zy=(y(fid)-min(yl))/(max(yl)-min(yl))/asp_rat*(range(ux)+1)*chanwid;
    
    % precalculate y distances so that we only worry about nearby points
    D0=squareform(pdist(zy))<dcut*2;    
    
    if length(fid)>1
        % for each data point in the group sequentially...
        for j=1:length(fid)
            if strcmp(lower(p.Results.sort_style),'hex')
                xi = linspace(-chanwid/2*rmult,chanwid/2*rmult,nxloc*rmult+(mod(nxloc*rmult,2)==0))'+ux(i);
                if mod(b(fid(j)),2)==oddmaj
                    xi = linspace(-chanwid/2*rmult,chanwid/2*rmult,nxloc*rmult+(mod(nxloc*rmult,2)==0))'+ux(i)+mean(diff(xi))/2;
                end
            end
            zid = D0(j,1:j-1);
            e = (xi-ux(i)).^2; % cost function
            if ~strcmp(lower(p.Results.sort_style),'hex') && ~strcmp(lower(p.Results.sort_style),'square')
                if sum(zid)>0
                    D = pdist2([xi ones(length(xi),1)*zy(j)], [x(fid(zid)) zy(zid)]);
                    D(D<=dcut)=Inf;
                    D(D>dcut & isfinite(D))=0;
                    e = e + sum(D,2) + randn(1)*10e-6; % noise to tie-break
                end
            else
                if sum(zid)>0
                    D = pdist2([xi ones(length(xi),1)*zy(j)], [x(fid(zid)) zy(zid)]);
                    D(D==0)=Inf;
                    D(D>dcut & isfinite(D))=0;
                    e = e + sum(D,2) + randn(1)*10e-6; % noise to tie-break
                end
            end

            if strcmp(lower(p.Results.sort_style),'one')
                e(xi<ux(i))=Inf;
            end
            [~,mini] = min(e);
            if mini==1 && rand(1)>.5, mini=length(xi); end
            x(fid(j)) = xi(mini);
        end
    end
%     x(fid)=x(fid)-median(x(fid))+ux(i); % center x locations by median
end

if strcmp(lower(p.Results.sort_style),'randn')
    x=ux(ic)+randn(size(ic))/4;
end

% corral any points outside of the channel
out_of_range = abs(x-ux(ic))>chanwid/2;
switch lower(p.Results.corral_style)
    case 'gutter'
        id = (x-ux(ic))>chanwid/2;
        x(id)=chanwid/2+ux(ic(id));
        id = (x-ux(ic))<-chanwid/2;
        x(id)=-chanwid/2+ux(ic(id));
    case 'omit'
        x(out_of_range)=NaN;
    case 'random'
        x(out_of_range)=ux(ic(out_of_range))+rand(sum(out_of_range),1)*chanwid-chanwid/2;
end

% plot groups and add overlay
if isfinite(p.Results.dot_size)
    if isnumeric(p.Results.colormap)
        cmap=p.Results.colormap;
    else
        cmap = feval(p.Results.colormap,length(ux));
    end
    for i=1:length(ux)
        if isempty(p.Results.MarkerFaceColor')
            scatter(x(ic==i),y(ic==i),p.Results.dot_size*36,'filled','MarkerFaceAlpha', ...
                p.Results.MarkerFaceAlpha,'MarkerEdgeColor',p.Results.MarkerEdgeColor,'MarkerFaceColor',cmap(i,:))
        else
            scatter(x(ic==i),y(ic==i),p.Results.dot_size*36,'filled','MarkerFaceAlpha', ...
                p.Results.MarkerFaceAlpha,'MarkerEdgeColor',p.Results.MarkerEdgeColor,'MarkerFaceColor',p.Results.MarkerFaceColor)
        end
        hold on
        iqr = prctile(yorig(ic==i),[25 75]);
        switch lower(p.Results.overlay_style)
            case 'box'
                rectangle('Position',[ux(i)-rwid iqr(1) 2*rwid iqr(2)-iqr(1)],'EdgeColor','k','LineWidth',2)
                line([ux(i)-rwid ux(i)+rwid],[1 1]*median(yorig(ic==i)),'LineWidth',3,'Color',cmap(i,:))
            case 'sd'
                % erbar_color = cmap(i,:);
                erbar_color = 'b';
                line([1 1]*ux(i),mean(yorig(ic==i))+[-1 1]*std(yorig(ic==i)),'Color',erbar_color,'LineWidth',4)
                line([ux(i)-4*rwid ux(i)+4*rwid],[1 1]*mean(yorig(ic==i)),'LineWidth',4,'Color',erbar_color)
            case 'se'
                line([1 1]*ux(i),mean(yorig(ic==i))+[-1 1]*std(yorig(ic==i))/sqrt(sum(ic==i)),'Color',cmap(i,:),'LineWidth',4)
                line([ux(i)-2*rwid ux(i)+2*rwid],[1 1]*mean(yorig(ic==i)),'LineWidth',4,'Color',cmap(i,:))
            case 'ci'
                line([1 1]*ux(i),mean(yorig(ic==i))+[-1 1]*std(yorig(ic==i))/sqrt(sum(ic==i))*tinv(0.975,sum(ic==i)-1),'Color',cmap(i,:),'LineWidth',2)
                line([ux(i)-2*rwid ux(i)+2*rwid],[1 1]*mean(yorig(ic==i)),'LineWidth',3,'Color',cmap(i,:))
        end
        
    end
    hold off
    if keep_hold
        hold on
    end
    xlim(xl)
    ylim(yl)
end

% unsort so that output matches the original y data
x(sid)=x;

end
