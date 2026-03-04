function export_norm_example_tex(varNames, varUnits, Xraw, lb, ub, Xn, outTex)
%EXPORT_NORM_EXAMPLE_TEX Export multi-row raw/normalized samples as LaTeX.
% Xraw and Xn can be 1 x k or nRows x k for the same samples.

if isvector(Xraw); Xraw = reshape(Xraw,1,[]); end
if isvector(Xn);   Xn   = reshape(Xn,1,[]);   end

nRows = size(Xraw,1);
k = numel(varNames);

if size(Xraw,2) ~= k || size(Xn,2) ~= k || size(Xn,1) ~= nRows
    error('Dimension mismatch in export_norm_example_tex inputs.');
end

fid = fopen(outTex,'w');
fprintf(fid,"%% Auto-generated normalization example tables\n");
fprintf(fid,"\\begin{table}[H]\n\\centering\n");
fprintf(fid,"\\caption{Pre/post min-max normalization examples for %d sample rows. Variables shown with units in headers.}\\label{tab:norm_example_multi}\n", nRows);
fprintf(fid,"\\renewcommand{\\arraystretch}{1.15}\n");

for r = 1:nRows
    fprintf(fid,"\\textbf{Sample %d}\\\\[2pt]\n", r);
    fprintf(fid,"\\begin{tabular}{lccccc}\\toprule\n");
    fprintf(fid,"\\textbf{Var} & \\textbf{Units} & \\textbf{Raw} & \\textbf{Min} & \\textbf{Max} & \\textbf{Norm}\\\\\\midrule\n");
    for i = 1:k
        fprintf(fid,"%s & %s & %.4g & %.4g & %.4g & %.4g\\\\\n", ...
            varNames(i), varUnits(i), Xraw(r,i), lb(i), ub(i), Xn(r,i));
    end
    fprintf(fid,"\\bottomrule\\end{tabular}\n");
    if r < nRows
        fprintf(fid,"\\vspace{6pt}\\\\\n");
    end
end

fprintf(fid,"\\end{table}\n");
fclose(fid);
end
