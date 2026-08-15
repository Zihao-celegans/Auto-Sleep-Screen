%% Combine Fig2 and Fig3 supplemental CSVs into one workbook (like Table_S1.xlsx).
% Sheet "MMP_screen"       <- figure2/Fig2_MMP_quiescence_data.csv
% Sheet "candidate_screen" <- figure3/Fig3_mutants_quiescence_data.csv
%
% Run this script from anywhere; paths below are resolved relative to
% this file's location (tableS1/), one level under the repo root.

clear; clc;

repoRoot = fileparts(fileparts(mfilename('fullpath')));

fig2Csv = fullfile(repoRoot, 'figure2', 'Fig2_MMP_quiescence_data.csv');
fig3Csv = fullfile(repoRoot, 'figure3', 'Fig3_mutants_quiescence_data.csv');
outFile = fullfile(fileparts(mfilename('fullpath')), 'Table_S1.xlsx');

if isfile(outFile)
    delete(outFile); % avoid stale sheets from a previous run
end

mmpScreen = readtable(fig2Csv);
candidateScreen = readtable(fig3Csv);

% readtable treats the "NA" text in SKAT_Rank_Percentile as a missing
% numeric value (NaN), which writetable then exports as a blank cell.
% Restore "NA" as literal text so it survives into the xlsx.
pct = candidateScreen.SKAT_Rank_Percentile;
isMissingPct = isnan(pct);
pctStr = strings(numel(pct), 1);
pctStr(~isMissingPct) = string(pct(~isMissingPct));
pctStr(isMissingPct) = "NA";
candidateScreen.SKAT_Rank_Percentile = pctStr;

writetable(mmpScreen, outFile, 'Sheet', 'MMP_screen');
writetable(candidateScreen, outFile, 'Sheet', 'candidate_screen');

fprintf('Wrote %s\n', outFile);
fprintf('  MMP_screen: %d rows\n', height(mmpScreen));
fprintf('  candidate_screen: %d rows\n', height(candidateScreen));
