function trl = remove_overlaps(trl)

% remove any mistaken trials
kickOut = [];
if size(trl, 1) > 1
    kickOut(1, 1) = false;
    for i = 2:size(trl, 1)
        if trl(i, 1) <= trl(i-1, 1) + 12*1024
            kickOut(i, 1) = true; %#ok<*AGROW>
        else
            kickOut(i, 1) = false;
        end
    end
    trl(logical(kickOut), :) = [];
end

end