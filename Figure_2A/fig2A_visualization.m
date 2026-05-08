clear; clc;

load('data2.mat');

%% =========================
% Filter strains
% ==========================
resultsTable2 = resultsTable2(ismember(resultsTable2.Strain, strainInfo.Strain_Name), :);

%% =========================
% Group by strain
% ==========================
[G, groupNames] = findgroups(resultsTable2.Strain);

% Pre UV
mean_pre = splitapply(@mean, resultsTable2.Qf_30min_preUV, G);
sd_pre   = splitapply(@std,  resultsTable2.Qf_30min_preUV, G);

% Post UV
mean_post = splitapply(@mean, resultsTable2.Qf_4hr_postUV, G);
sd_post   = splitapply(@std,  resultsTable2.Qf_4hr_postUV, G);

% Counts
N = splitapply(@numel, resultsTable2.Qf_4hr_postUV, G);

%% =========================
% Add genotype info
% ==========================
strainNames = cellstr(groupNames);
infoNames = cellstr(strainInfo.Strain_Name);
infoGenos = cellstr(strainInfo.Genotype);

genotypes = strings(numel(strainNames),1);

for i = 1:numel(strainNames)
    idx = strcmp(infoNames, strainNames{i});
    if any(idx)
        genotypes(i) = string(infoGenos{find(idx,1)});
    else
        genotypes(i) = "";
    end
end

%% =========================
% Create table
% ==========================
summaryTable = table( ...
    strainNames, ...
    genotypes, ...
    N, ...
    mean_pre, sd_pre, ...
    mean_post, sd_post, ...
    'VariableNames', { ...
        'Strain', ...
        'Genotype', ...
        'N', ...
        'Mean_PreUV', 'SD_PreUV', ...
        'Mean_PostUV', 'SD_PostUV'});

%% =========================
% Display
% ==========================
%disp(summaryTable);

%% =========================
% Save outputs
% ==========================
save('Quiescence_Summary.mat', 'summaryTable');   % MATLAB file
