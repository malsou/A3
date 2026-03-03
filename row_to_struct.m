function s = row_to_struct(row, varNames)
s = struct;
for j = 1:numel(varNames)
    s.(varNames(j)) = row(j);
end
end
