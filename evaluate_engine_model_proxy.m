function y = evaluate_engine_model_proxy(x)
% Proxy wrapper for A3 screening.
% Maps legacy A3 design vector into eval_engine(des, op, opts).

% ---- Operating points (SLS / climb / cruise) ----
op = operating_points();

% ---- Design struct for eval_engine ----
des = struct;
des.BPR = x.BPR;
des.FPR = x.FPR;
des.LPR = x.LPR;
des.HPR = x.HPR;
des.Tt4 = x.Tt4;
des.nozzle = x.nozzle;
des.cooling = x.cooling;
des.porosity = x.porosity;
des.tTBC = x.tTBC;
des.kSt = x.kSt;

% ---- Fixed assumptions / coupled-model constants ----
opts = engine_defaults();
opts.SMc = x.SMc;
if isfield(x,'Cd'); opts.Cd = x.Cd; end
if isfield(x,'alpha'); opts.alpha = x.alpha; end

y = eval_engine(des, op, opts);
end

function op = operating_points()
% A3 mission anchors with reduced required inputs for each operating point:
% altitude [ft], Mach number, required net thrust [N].
op(1).name = "SLS_takeoff"; op(1).alt_ft = 0;      op(1).M0 = 0.00; op(1).Freq = 120e3;
op(2).name = "climb";       op(2).alt_ft = 19685;  op(2).M0 = 0.40; op(2).Freq = 60e3;
op(3).name = "cruise";      op(3).alt_ft = 36089;  op(3).M0 = 0.78; op(3).Freq = 25e3;
end
