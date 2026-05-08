clear; clc; close all;

load('data.mat');


%% === Fig2A ===

% === FILTER strains: keep only those in strainInfo ===
resultsTable2 = resultsTable2(ismember(resultsTable2.Strain, strainInfo.Strain_Name), :);

% Summaries by strain
[G, groupNames] = findgroups(resultsTable2.Strain);
strain_means = splitapply(@mean, resultsTable2.Qf_4hr_postUV, G);
strain_sems  = splitapply(@(x) std(x)/sqrt(numel(x)), resultsTable2.Qf_4hr_postUV, G);
strain_counts = splitapply(@numel, resultsTable2.Qf_4hr_postUV, G); % count worms

% Identify N2 (wild type)
wt_idx = find(groupNames == 'N2 UV');

wt_data = resultsTable2.Qf_4hr_postUV(resultsTable2.Strain == groupNames(wt_idx));

% Compute p-values vs WT
pvals = nan(size(groupNames));
for i = 1:numel(groupNames)
    if i == wt_idx, continue; end
    strain_data = resultsTable2.Qf_4hr_postUV(resultsTable2.Strain == groupNames(i));
    [~, p] = ttest2(wt_data, strain_data);
    pvals(i) = p;
end

% === FDR correction ===
p_fdr = mafdr(pvals, 'BHFDR', true);

% Sort all strains by descending mean (including N2)
[~, sortOrder] = sort(strain_means, 'descend');
groupNames = groupNames(sortOrder);
strain_means = strain_means(sortOrder);
strain_sems  = strain_sems(sortOrder);
strain_counts = strain_counts(sortOrder);
pvals = pvals(sortOrder);
p_fdr = p_fdr(sortOrder);

% Plot bar + SEM
fig = figure('Color','w'); hold on;
set(fig, 'WindowState', 'maximized');

b = bar(strain_means, 'FaceColor','flat','EdgeColor','none');

% Colors: N2 = green, significant = red, others = blue
b.CData = repmat([0 0 1], numel(groupNames), 1);
b.CData(groupNames == 'N2 UV', :) = [0 0.7 0];
sig_idx = find(p_fdr < 0.05 & ~isnan(p_fdr));
b.CData(sig_idx,:) = repmat([1 0 0], numel(sig_idx), 1);

% Error bars
errorbar(1:numel(groupNames), strain_means, strain_sems, ...
    'k', 'LineStyle','none', 'LineWidth',1.5);

% ===== Build one-line x tick labels: Strain ; Genotype =====
strainNames = cellstr(groupNames);

infoNames  = cellstr(strainInfo.Strain_Name);
infoGenos  = cellstr(strainInfo.Genotype);

genotypes = strings(numel(strainNames),1);
for ii = 1:numel(strainNames)
    idx = strcmp(infoNames, strainNames{ii});
    if any(idx)
        g = infoGenos(idx);
        genotypes(ii) = string(g{1});
    else
        genotypes(ii) = "";
    end
end

combinedLabels = strcat(strainNames, {') '}, genotypes);
combinedLabels = strcat(string(1:numel(combinedLabels))', {'- '}, combinedLabels);

set(gca, 'XTick', 1:numel(groupNames), ...
         'XTickLabel', combinedLabels, ...
         'TickLabelInterpreter','none', ...
         'FontSize', 7, 'LineWidth', 1.2);
xtickangle(-45);

ylabel('Quiescence fraction (4h post UV)', 'FontSize', 12);
title('UV-induced sleep: N2 vs other SKAT strains', 'FontSize', 13);

% Annotate only significant p-values
ylims = ylim;
for ii = sig_idx'
    yPos = strain_means(ii) + strain_sems(ii) + 0.1*range(ylims);
    text(ii, yPos, sprintf('q = %.2g', p_fdr(ii)), ...
        'HorizontalAlignment','center','VerticalAlignment','middle', ...
        'Rotation', -90, ...
        'FontSize', 9,'FontWeight','bold');
end

% Annotate number of worms (dots) per strain
for ii = 1:numel(groupNames)
    yPos = strain_means(ii) - strain_sems(ii) - 0.1*range(ylims);
    text(ii, yPos, sprintf('%d', strain_counts(ii)), ...
        'HorizontalAlignment','center','VerticalAlignment','bottom', ...
        'FontSize', 8, 'Color','k', ...
        'BackgroundColor','w','Margin',1,'EdgeColor','k','LineWidth',0.8);
end

ylim([0, ylims(2)*1.08]);
box off; ylim([0 1]);
set(gca,'TickLength',[0 0]);

% Small legend
hN2  = bar(nan, 'FaceColor', [0 0.7 0]);
hSig = bar(nan, 'FaceColor', [1 0 0]);
hNS  = bar(nan, 'FaceColor', [0 0 1]);
hDot = scatter(nan, nan, 40, 'square', 'MarkerFaceColor','w', 'MarkerEdgeColor','k');

legend([hN2 hSig hNS hDot], {'N2','Significant vs N2','Not significant','number of worms'}, ...
    'FontSize', 8, 'Location', 'northeastoutside', 'Box', 'off');

close all;

%% === Fig2B, Fig 2C, Fig2D, Fig2E, Fig2F, Fig2G ===

%% ============================================================
%   STRAIN → GENE MAP
%% ============================================================

strainToGene = containers.Map( ...
    {'N2 UV', ...
     'RB1771 UV', 'RB1775 UV', ...
     'VC387 UV','RB778 UV','RB1777 UV','VC596 UV','VC631 UV', ...
     'CB6614 UV','RM2221 UV','IK777 UV','JT47 UV','MT1083 UV'}, ...
    {'N2', ...
     'strd-1', 'strd-1', ...
     'cla-1','cla-1','cla-1','cla-1','cla-1', ...
     'egl-8','egl-8','egl-8','egl-8','egl-8'} ...
     );

geneList = {'strd-1', 'cla-1', 'egl-8'};

alpha_pre  = 0.05;
alpha_post = 0.05;

% =========================
% CGC-based display labels
% =========================
strainDisplay = containers.Map( ...
    {'N2 UV', ...
     'RB1771 UV', 'RB1775 UV', ...
     'VC387 UV', 'RB778 UV', 'RB1777 UV', 'VC596 UV', 'VC631 UV', ...
     'CB6614 UV', 'RM2221 UV', 'IK777 UV', 'JT47 UV', 'MT1083 UV'}, ...
    {'N2', ...
     'ok2275', 'ok2283', ...
     'ok618', 'ok560', 'ok2285', 'gk352', 'ok937', ...
     'e2917', 'md1971', 'nj77', 'sa47', 'n488'} ...
    );
% , 'RB1012 UV'
% , 'egl-8(ok934)'
%% ============================================================
%   PATHS
%% ============================================================

baseDirs = { ...
    'C:\Users\lizih\Dropbox\Shared_folders\Helya_John_Shared_folder\Writing\Raw_data\SleepScreenSelectedMOMO', ...
    'C:\Users\lizih\Dropbox\Shared_folders\Helya_John_Shared_folder\Writing\Raw_data\SleepScreenSelectedKAWKAB'};

outDir = 'C:\Users\jl200\Dropbox\Publication\SleepScreen\Figures\Helya_figures\Figure_2\ResultsImages';
% if ~exist(outDir,'dir'), mkdir(outDir); end

%% ============================================================
%   LOAD DATA
%% ============================================================

allStrains = strings(0,1);
allPost = [];
allPre = [];
allDates = strings(0,1);
alldatPRE = {};
alldatPOST = {};

for b = 1:numel(baseDirs)

    root = baseDirs{b};
    folders = dir(root);
    folders = folders([folders.isdir] & ~ismember({folders.name},{'.','..'}));

    for i = 1:numel(folders)

        sub = folders(i).name;
        matFile = fullfile(root,sub,'Analysis',[sub '_ForPlot.mat']);
        if ~exist(matFile,'file'), continue; end

        S = load(matFile,'Qf_output_struct');

        analysisDir = fullfile(root, sub, 'Analysis');

        FilePRE  = dir(fullfile(analysisDir, '*_NS.mat'));
        FilePOST = dir(fullfile(analysisDir, '*_UVC.mat'));

        if isempty(FilePRE) || isempty(FilePOST)
            continue;
        end

        Spre  = load(fullfile(analysisDir, FilePRE(1).name));
        datPRE = Spre.ActVal;

        Spost = load(fullfile(analysisDir, FilePOST(1).name));
        datPOST = Spost.ActVal;

        strains = string(S.Qf_output_struct.strain_category);
        strains = regexprep(strtrim(strains),'^D\d+-','');

        post = S.Qf_output_struct.Qf_4hr_postUV(:);
        pre  = S.Qf_output_struct.Qf_30min_preUV(:);

        n = min([numel(post), numel(pre), size(datPRE,1), size(datPOST,1), numel(strains)]);

        strains = strains(1:n);
        post    = post(1:n);
        pre     = pre(1:n);
        datPRE  = datPRE(1:n,:);
        datPOST = datPOST(1:n,:);

        dateToken = string(regexp(sub,'^\d{2}-\d{2}-\d+','match','once'));

        allStrains = [allStrains; strains];
        allPost = [allPost; post];
        allPre  = [allPre; pre];
        allDates = [allDates; repmat(dateToken,n,1)];

        for k = 1:n
            alldatPRE{end+1,1}  = datPRE(k,:);
            alldatPOST{end+1,1} = datPOST(k,:);
        end
    end
end

resultsTable = table( ...
    categorical(allStrains), ...
    allPost, allPre, categorical(allDates), ...
    'VariableNames',{'Strain','Qf_post','Qf_pre','Date'});

traceTable = table( ...
    categorical(allStrains), ...
    categorical(allDates), ...
    alldatPRE(:), ...
    alldatPOST(:), ...
    'VariableNames',{'Strain','Date','datPRE','datPOST'});

%% ============================================================
%   FILTER TO GENE STRAINS + SAME-DATE N2
%% ============================================================

allowedStrains = string(keys(strainToGene));
allowedCat = categorical(allowedStrains, categories(resultsTable.Strain));

resultsTable = resultsTable(ismember(resultsTable.Strain, allowedCat), :);

mutDates = unique(resultsTable.Date(resultsTable.Strain ~= "N2 UV"));

resultsTable = resultsTable( resultsTable.Strain ~= "N2 UV" | ...
                            (resultsTable.Strain=="N2 UV" & ismember(resultsTable.Date,mutDates)), :);

allowedCat_trace = categorical(allowedStrains, categories(traceTable.Strain));

traceTable = traceTable(ismember(traceTable.Strain, allowedCat_trace), :);

mutDates_trace = unique(traceTable.Date(traceTable.Strain ~= "N2 UV"));

traceTable = traceTable( traceTable.Strain ~= "N2 UV" | ...
                        (traceTable.Strain=="N2 UV" & ismember(traceTable.Date,mutDates_trace)), :);

%% ============================================================
%   PER-GENE GROUP ANALYSIS
%% ============================================================

for g = 1:numel(geneList)

    gene = geneList{g};

    %% Find strains for this gene
    strainCats = categories(resultsTable.Strain);
    geneStrains = strings(0,1);

    for i = 1:numel(strainCats)
        s = string(strainCats{i});
        if isKey(strainToGene,s) && strcmp(strainToGene(s),gene)
            geneStrains(end+1,1) = s;
        end
    end

    geneStrains = unique(geneStrains,'stable');

    strainsToKeep = [geneStrains; "N2 UV"];
    strains_cat = categorical(strainsToKeep, categories(resultsTable.Strain));

    subT = resultsTable(ismember(resultsTable.Strain,strains_cat), :);

    strains_cat_trace = categorical(strainsToKeep, categories(traceTable.Strain));
    subTrace = traceTable(ismember(traceTable.Strain, strains_cat_trace), :);

    %% Same-date control
    mutDates = unique(subT.Date(subT.Strain ~= "N2 UV"));
    subT = subT(subT.Strain ~= "N2 UV" | ismember(subT.Date, mutDates), :);

    mutDates_trace = unique(subTrace.Date(subTrace.Strain ~= "N2 UV"));
    subTrace = subTrace(subTrace.Strain ~= "N2 UV" | ismember(subTrace.Date, mutDates_trace), :);

    % ============================================================
    %   TRACE HEATMAP FOR THIS GENE GROUP
    % ============================================================

    % Force N2 rows to appear first
    isN2 = subTrace.Strain == "N2 UV";
    subTrace.isN2 = ~isN2;   % N2 = 0, others = 1
    subTrace = sortrows(subTrace, {'isN2','Strain','Date'});
    subTrace.isN2 = [];

    if isempty(subTrace)
        warning('subTrace is empty for gene %s', gene);
        continue;
    end

    preCells  = subTrace.datPRE;
    postCells = subTrace.datPOST;

    % force all traces to row vectors
    preCells  = cellfun(@(x) x(:)', preCells,  'UniformOutput', false);
    postCells = cellfun(@(x) x(:)', postCells, 'UniformOutput', false);

    lenPre  = cellfun(@numel, preCells);
    lenPost = cellfun(@numel, postCells);

    % keep all nonempty worms only
    keepIdx = (lenPre > 0) & (lenPost > 0);
    traceKeep = subTrace(keepIdx,:);

    preCells  = preCells(keepIdx);
    postCells = postCells(keepIdx);

    if isempty(traceKeep)
        warning('No nonempty traces for gene %s', gene);
        continue;
    end

    % IMPORTANT: preserve sorted row order for grouping
    strainOrder = unique(string(traceKeep.Strain), 'stable');
    traceKeep.Strain = categorical(string(traceKeep.Strain), strainOrder);

    [GtraceKeep, sTraceKeep] = findgroups(traceKeep.Strain);
    nPerStrain = splitapply(@numel, traceKeep.Strain, GtraceKeep);
    edges = cumsum(nPerStrain);
    if ~isempty(edges)
        edges(end) = [];
    end

    % pad PRE to max PRE length
    maxPreLen = max(cellfun(@numel, preCells));
    alldatPRE_gene = nan(numel(preCells), maxPreLen);

    for ii = 1:numel(preCells)
        x = preCells{ii};
        alldatPRE_gene(ii,1:numel(x)) = x;
    end

    % pad POST to max POST length
    maxPostLen = max(cellfun(@numel, postCells));
    alldatPOST_gene = nan(numel(postCells), maxPostLen);

    for ii = 1:numel(postCells)
        x = postCells{ii};
        alldatPOST_gene(ii,1:numel(x)) = x;
    end

    [qPRE,~,~]  = fracQ(alldatPRE_gene',1,10,60);
    [qPOST,~,~] = fracQ(alldatPOST_gene',1,10,60);

    % safe indexing
    startPre  = min(600, size(qPRE,1));
    startPost = min(600, size(qPOST,1));

    qALL = [qPRE(startPre:10:end,:); qPOST(startPost:10:end,:)];

    if isempty(qALL)
        warning('qALL is empty for gene %s', gene);
        continue;
    end

    figure('Color','w','WindowState','maximized');
    imagesc(qALL');

    ax = gca;
    ax.FontSize = 12;
    colorbar
    colormap jet

    starts  = [1; edges(:) + 1];
    stops   = [edges(:); sum(nPerStrain)];
    centers = (starts + stops) / 2;

    strainLabels = string(sTraceKeep);
    displayLabels = strings(numel(strainLabels),1);

    for ii = 1:numel(strainLabels)
        s = strainLabels(ii);
        if isKey(strainDisplay, s)
            displayLabels(ii) = strainDisplay(s);
        else
            displayLabels(ii) = s;
        end
    end

    y_tick_labels = cellstr(displayLabels);

    for k = 2:size(y_tick_labels, 1)
        y_tick_labels{k} = ['\it{' y_tick_labels{k} '}'];  % italic allele names
    end


    ax = gca;
    ax.TickLabelInterpreter = 'tex';
    set(ax, 'YTick', centers, ...
             'YTickLabel', y_tick_labels);

    % put ticks every hour
    xticks_steps = 0:360:size(qALL, 1);   % 0h, 1h, 2h, 3h, 4h
    set(gca, 'XTick', xticks_steps + 1);   % +1 because MATLAB indices start at 1
    set(gca, 'XTickLabel', string(xticks_steps/360));

    hold on
    for k2 = edges(:)'
        plot(xlim, (k2 + 0.5)*[1 1], 'Color', [1 1 1], 'LineWidth', 5);
    end
    hold off

    % title(['Quiescence fraction - ' gene], 'FontSize',16);
    xlabel('Time (hour)');

    f = gcf;
    f.Units = 'pixels';
    f.Position = [100 100 800 800];

    ax = gca; ax.FontSize = 25;
    pbaspect([1 1 1]);

    exportgraphics(gca, [ gene,'.svg']);

    % Per-strain stats
    % preserve display/order consistency in subT too
    subT.Strain = categorical(string(subT.Strain), unique(string(subT.Strain), 'stable'));

    [G,sNames] = findgroups(subT.Strain);
    [Gtrace,sTraceNames] = findgroups(subTrace.Strain);

    preTraceByStrain  = splitapply(@(x) {x}, subTrace.datPRE, Gtrace);
    postTraceByStrain = splitapply(@(x) {x}, subTrace.datPOST, Gtrace);

    meanPost = splitapply(@mean, subT.Qf_post, G);
    semPost  = splitapply(@(x) std(x)/sqrt(numel(x)), subT.Qf_post, G);
    meanPre  = splitapply(@mean, subT.Qf_pre, G);
    semPre   = splitapply(@(x) std(x)/sqrt(numel(x)), subT.Qf_pre, G);
    nWorms   = splitapply(@numel, subT.Qf_post, G);

    % Pre + Post P-values
    pPost = nan(size(sNames));
    pPre  = nan(size(sNames));

    for i = 1:numel(sNames)

        if sNames(i) == "N2 UV", continue; end

        idxS = subT.Strain == sNames(i);
        dates_i = unique(subT.Date(idxS));

        idxN2 = subT.Strain=="N2 UV" & ismember(subT.Date,dates_i);

        wt_post = subT.Qf_post(idxN2);
        wt_pre  = subT.Qf_pre(idxN2);
        dat_post = subT.Qf_post(idxS);
        dat_pre  = subT.Qf_pre(idxS);

        [~,pPost(i)] = ttest2(wt_post, dat_post);
        [~,pPre(i)]  = ttest2(wt_pre, dat_pre, 'Tail', 'left');
    end

    qPost = mafdr(pPost,'BHFDR',true);
    qPre  = mafdr(pPre,'BHFDR',true);

    % Put N2 first
    n2_idx = find(sNames=="N2 UV");
    if ~isempty(n2_idx)
        order = [n2_idx, setdiff(1:numel(sNames),n2_idx)];
        sNames   = sNames(order);
        meanPost = meanPost(order);
        semPost  = semPost(order);
        meanPre  = meanPre(order);
        semPre   = semPre(order);
        nWorms   = nWorms(order);
        qPost    = qPost(order);
        qPre     = qPre(order);
    end

    % ============================================================
    %   Scatterplot: pre/post colors fixed + qPre/qPost + n worms
    % ============================================================

    figure('Color','w','WindowState','maximized'); hold on;

    col_pre  = [0 0 0];
    col_post = [0.4940 0.1840 0.5560];

    x_pre  = 2*(1:numel(sNames)) - 1;
    x_post = x_pre + 1;

    yl = [0 1];

    cOrder = colororder;
    Color_preUV = [0 0 0];
    Color_postUV = cOrder(4,:);

    for i = 1:numel(sNames)

        idx = subT.Strain == sNames(i);
        dat_pre  = subT.Qf_pre(idx);
        dat_post = subT.Qf_post(idx);

        beeswarm( x_pre(i)*ones(numel(dat_pre),1), dat_pre, ...
            'corral_style','random', ...
            'overlay_style','sd', ...
            'colormap', Color_preUV);

        beeswarm( x_post(i)*ones(numel(dat_post),1), dat_post, ...
            'corral_style','random', ...
            'overlay_style','sd', ...
            'colormap', Color_postUV);

    end

    ylim([yl(1)-0.15 yl(2)+0.20]);
    xlim([0 max(x_post)+1]);

    ylabel('Quiescence Fraction');

    x_mid = (x_pre + x_post)/2;

    strainLabels = string(sNames);
    displayLabels = strings(numel(strainLabels),1);

    for ii = 1:numel(strainLabels)
        s = strainLabels(ii);
        if isKey(strainDisplay, s)
            displayLabels(ii) = strainDisplay(s);
        else
            displayLabels(ii) = s;
        end
    end

    set(gca,'XTick',x_mid,'XTickLabel',cellstr(displayLabels));
    xtickangle(45);

    for i = 1:numel(sNames)
        text(x_mid(i), yl(1)-0.08, sprintf('n=%d', nWorms(i)), ...
            'HorizontalAlignment','center', 'FontSize',12);
    end

    for i = 1:numel(sNames)

        if ~isnan(qPre(i))
            labelPre = sprintf('qpre=%.3f', qPre(i));
            if qPre(i) < 0.05 && sNames(i) ~= "N2 UV"
                labelPre = [labelPre ' *'];
            end
            text(x_pre(i), yl(2)+0.05, labelPre, ...
                'HorizontalAlignment','center','FontSize',12);
        end

        if ~isnan(qPost(i))
            labelPost = sprintf('qpost=%.3f', qPost(i));
            if qPost(i) < 0.05 && sNames(i) ~= "N2 UV"
                labelPost = [labelPost ' *'];
            end
            text(x_post(i), yl(2)+0.10, labelPost, ...
                'HorizontalAlignment','center','FontSize',12);
        end
    end

    hPre  = plot(nan,nan,'o','MarkerFaceColor',col_pre,  'MarkerEdgeColor','none');
    hPost = plot(nan,nan,'o','MarkerFaceColor',col_post,'MarkerEdgeColor','none');

    legend([hPre hPost], {'Pre-UV','Post-UV'}, ...
        'Location','northeast','Box','off');

    allSc = findobj(gca,'Type','Scatter');
    set(allSc, ...
        'MarkerEdgeColor','none', ...
        'MarkerFaceAlpha',0.7, ...
        'SizeData',30, ...
        'Marker','o');

    box off;
    title(['Gene Group: ' gene],'FontSize',16,'FontWeight','bold');

end

%%
% fracQ.m
function [PAV, t1, pnew] = fracQ(p,qt,spi,ps)

if mod(ps,2)==0
    ps=ps+1;
end

pq=p<qt;
[sizex, sizey]=size(p);
t1=[1:1/spi:sizex]/3600*spi;
rois=1:sizey;
pstart=ceil(ps/2);
pf=floor(ps/2);

pav=zeros(ceil(sizex/ps),sizey);
for i=pstart:ps:sizex-ps
    for j=1:sizey
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

rwid = .05;
dcut=8;
nxloc=512;
chanwid = .9;
yl = [min(y) max(y)];
asp_rat = 1;
keep_hold = false;

if isfinite(p.Results.dot_size)
    if ~p.Results.use_current_axes
        scatter(x,y);
        xl=[min(x)-.5 max(x)+.5];
    else
        xl=xlim();
    end
    yl=ylim();
    pasp_rat = get(gca,'PlotBoxAspectRatio');
    asp_rat = pasp_rat(1)/pasp_rat(2);

    pf = get(gcf,'Position');
    pa = get(gca,'Position');
    as = pf(3:4).*pa(3:4);
    dcut = dcut*sqrt(p.Results.dot_size)/as(1)*(range(unique(x))+1);
    if ~ishold
        cla
    else
        keep_hold = true;
    end
end

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
        edges = linspace(min(yl),max(yl),ceil((range(x)+1)*chanwid*nxloc/asp_rat));
        [~,e,b]=histcounts(y,edges);
        y=e(b)'+mean(diff(e))/2;
        [y,sid]=sort(y);
    case 'hex'
        nxloc=.9/dcut;
        edges = linspace(min(yl),max(yl),ceil((range(x)+1)*chanwid*nxloc/sqrt(1-.5.^2)/asp_rat));
        [n,e,b]=histcounts(y,edges);
        oddmaj=0;
        if sum(mod(n(1:2:end),2)==1)>sum(mod(n(2:2:end),2)==1)
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
rmult=5;

for i=1:length(ux)
    fid = find(ic==i);

    xi = linspace(-chanwid/2*rmult,chanwid/2*rmult,nxloc*rmult+(mod(nxloc*rmult,2)==0))'+ux(i);
    zy=(y(fid)-min(yl))/(max(yl)-min(yl))/asp_rat*(range(ux)+1)*chanwid;
    D0=squareform(pdist(zy))<dcut*2;

    if length(fid)>1
        for j=1:length(fid)
            if strcmp(lower(p.Results.sort_style),'hex')
                xi = linspace(-chanwid/2*rmult,chanwid/2*rmult,nxloc*rmult+(mod(nxloc*rmult,2)==0))'+ux(i);
                if mod(b(fid(j)),2)==oddmaj
                    xi = linspace(-chanwid/2*rmult,chanwid/2*rmult,nxloc*rmult+(mod(nxloc*rmult,2)==0))'+ux(i)+mean(diff(xi))/2;
                end
            end
            zid = D0(j,1:j-1);
            e = (xi-ux(i)).^2;
            if ~strcmp(lower(p.Results.sort_style),'hex') && ~strcmp(lower(p.Results.sort_style),'square')
                if sum(zid)>0
                    D = pdist2([xi ones(length(xi),1)*zy(j)], [x(fid(zid)) zy(zid)]);
                    D(D<=dcut)=Inf;
                    D(D>dcut & isfinite(D))=0;
                    e = e + sum(D,2) + randn(1)*10e-6;
                end
            else
                if sum(zid)>0
                    D = pdist2([xi ones(length(xi),1)*zy(j)], [x(fid(zid)) zy(zid)]);
                    D(D==0)=Inf;
                    D(D>dcut & isfinite(D))=0;
                    e = e + sum(D,2) + randn(1)*10e-6;
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
end

if strcmp(lower(p.Results.sort_style),'randn')
    x=ux(ic)+randn(size(ic))/4;
end

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

if isfinite(p.Results.dot_size)
    if isnumeric(p.Results.colormap)
        cmap=p.Results.colormap;
    else
        cmap = feval(p.Results.colormap,length(ux));
    end
    for i=1:length(ux)
        if isempty(p.Results.MarkerFaceColor')
            scatter(x(ic==i),y(ic==i),p.Results.dot_size*36,'filled', ...
                'MarkerFaceAlpha',p.Results.MarkerFaceAlpha, ...
                'MarkerEdgeColor',p.Results.MarkerEdgeColor, ...
                'MarkerFaceColor',cmap(i,:))
        else
            scatter(x(ic==i),y(ic==i),p.Results.dot_size*36,'filled', ...
                'MarkerFaceAlpha',p.Results.MarkerFaceAlpha, ...
                'MarkerEdgeColor',p.Results.MarkerEdgeColor, ...
                'MarkerFaceColor',p.Results.MarkerFaceColor)
        end
        hold on
        iqr = prctile(yorig(ic==i),[25 75]);
        switch lower(p.Results.overlay_style)
            case 'box'
                rectangle('Position',[ux(i)-rwid iqr(1) 2*rwid iqr(2)-iqr(1)],'EdgeColor','k','LineWidth',2)
                line([ux(i)-rwid ux(i)+rwid],[1 1]*median(yorig(ic==i)),'LineWidth',3,'Color',cmap(i,:))
            case 'sd'
                line([1 1]*ux(i),mean(yorig(ic==i))+[-1 1]*std(yorig(ic==i)),'Color',cmap(i,:),'LineWidth',2)
                line([ux(i)-2*rwid ux(i)+2*rwid],[1 1]*mean(yorig(ic==i)),'LineWidth',3,'Color',cmap(i,:))
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

x(sid)=x;

end