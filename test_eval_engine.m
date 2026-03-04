% Quick baseline smoke-test for eval_engine
clear; clc;

des = struct;
des.BPR = 8.0;
des.FPR = 1.60;
des.LPR = 1.8;
des.HPR = 18.0;
des.Tt4 = 1600;
des.nozzle = "mixed";
des.cooling = "film";
des.porosity = 0.010;
des.tTBC = 250;
des.kSt = 1.0;

op(1) = struct('alt_ft',0,'M0',0.00,'Freq',120e3);
op(2) = struct('alt_ft',19685,'M0',0.40,'Freq',60e3);
op(3) = struct('alt_ft',36089,'M0',0.78,'Freq',25e3);

opts = engine_defaults();
opts.SMc = 0.15;
opts.Cd = 0.80;
opts.alpha = 45;

y = eval_engine(des, op, opts);

fprintf('TSFC_cruise = %.6f lbm/hr/lbf\n', y.TSFC_cruise);
fprintf('phi_crit    = %.6f\n', y.phi_crit);
fprintf('MT_worst    = %.6f K\n', y.MT_worst);
fprintf('MF_worst    = %.6f\n', y.MF_worst);
