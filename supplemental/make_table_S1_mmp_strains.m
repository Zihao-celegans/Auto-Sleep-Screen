%% Build Supplemental Table S1 (MMP strain list).
% Output: Table_S1.xlsx in this folder.
%
% Sources:
%   ../figure2/robot_MMP_Qf.mat   (ranked_MMP_name_list)
%   ../figure2/sleep_data_cleaned.mat (strain_names)

clear; clc;

S1 = load(fullfile('..', 'figure2', 'robot_MMP_Qf.mat'), 'ranked_MMP_name_list');
S2 = load(fullfile('..', 'figure2', 'sleep_data_cleaned.mat'), 'strain_names');

ranked_MMP_name_list = S1.ranked_MMP_name_list;
strain_names         = S2.strain_names;

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

% Write to Excel (relative path; output stays next to the script)
f_tableS1 = "Table_S1.xlsx";

sheetName = "MMP_screen";

writetable(T_mmp_info, f_tableS1, ...
    'Sheet', sheetName, ...
    'WriteMode', 'overwritesheet');
