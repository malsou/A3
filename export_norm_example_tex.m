function export_norm_example_tex(varNames, varUnits, xRaw, lb, ub, xN, outTex)
fid = fopen(outTex,'w');
fprintf(fid,"%% Auto-generated normalization example table\n");
fprintf(fid,"\\begin{table}[H]\n\\centering\n");
fprintf(fid,"\\caption{Example pre- and post-normalization values (min--max).}\\label{tab:norm_example}\n");
fprintf(fid,"\\renewcommand{\\arraystretch}{1.15}\n");
fprintf(fid,"\\begin{tabular}{lccccc}\n\\toprule\n");
fprintf(fid,"\\textbf{Var} & \\textbf{Units} & \\textbf{Raw} & \\textbf{Min} & \\textbf{Max} & \\textbf{Norm}\\\\\n\\midrule\n");

for i = 1:numel(varNames)
    fprintf(fid,"%s & %s & %.4g & %.4g & %.4g & %.4g\\\\\n", ...
        varNames(i), varUnits(i), xRaw(i), lb(i), ub(i), xN(i));
end

fprintf(fid,"\\bottomrule\n\\end{tabular}\n\\end{table}\n");
fclose(fid);
end
