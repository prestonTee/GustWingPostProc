function [vDat] = extractVartInfo(vFileLoc, vPrefix, vSuffix, delta_ts)

%% Varts Preprocessing
list = dir([vFileLoc, vPrefix, '*', vSuffix]);
for i = 1:length(list)
    vIndexList{i} = strrep(list(i).name, vPrefix, '');
    vIndexList{i} = str2num(strrep(vIndexList{i}, vSuffix, ''));
end
vIndexList = sort(cell2mat(vIndexList));

%% Read In vart Data
vDat = zeros(delta_ts, 7);
for i = 1:length(list)
    dat = load([vFileLoc, vPrefix, num2str(vIndexList(i)), vSuffix]);
    if vDat(1, 1) == 0
        vDat(1:size(dat, 1), :) = dat(:, :);
    else
        index = find(vDat(:, 1) == dat(1, 1)-1, 1)+1;
        vDat(index:index+size(dat, 1)-1, :) = dat(:, :);
    end
end

end
