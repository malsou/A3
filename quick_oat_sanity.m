% Quick one-at-a-time sanity check around a baseline design.
clear; clc;

op = operating_points();
opts = engine_defaults();

base = struct;
base.BPR = 8.0;
base.FPR = 1.60;
base.LPR = 1.8;
base.HPR = 18.0;
base.Tt4 = 1600;
base.SMc = 0.15;
base.kSt = 1.0;
base.porosity = 0.010;
base.Cd = 0.80;
base.alpha = 45;
base.tTBC = 250;
base.nozzle = "mixed";
base.cooling = "film";

vars = ["BPR","FPR","LPR","HPR","Tt4","SMc","kSt","porosity","Cd","alpha","tTBC"];

y0 = eval_engine(base, op, opts);
fprintf('Baseline TSFC=%.6e, MT_worst=%.6f\n', y0.TSFC_cruise, y0.MT_worst);

for i = 1:numel(vars)
    vn = vars(i);
    yp = base;
    ym = base;
    yp.(vn) = 1.10 * base.(vn);
    ym.(vn) = 0.90 * base.(vn);

    ypOut = eval_engine(yp, op, opts);
    ymOut = eval_engine(ym, op, opts);

    dTSFC = ypOut.TSFC_cruise - ymOut.TSFC_cruise;
    dMT = ypOut.MT_worst - ymOut.MT_worst;
    fprintf('%-10s : DeltaTSFC=%+.6e   DeltaMT=%+.6f\n', vn, dTSFC, dMT);
end

function op = operating_points()
op(1).name = "SLS_takeoff"; op(1).alt_ft = 0;      op(1).M0 = 0.00; op(1).Freq = 120e3;
op(2).name = "climb";       op(2).alt_ft = 19685;  op(2).M0 = 0.40; op(2).Freq = 60e3;
op(3).name = "cruise";      op(3).alt_ft = 36089;  op(3).M0 = 0.78; op(3).Freq = 25e3;
end
