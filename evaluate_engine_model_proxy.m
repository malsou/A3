function y = evaluate_engine_model_proxy(x)
% Proxy wrapper for A3 screening.
% Thin wrapper around eval_engine(des, op, opts) to avoid duplicate physics.

op = operating_points();
opts = engine_defaults();
y = eval_engine(x, op, opts);
end

function op = operating_points()
% A3 mission anchors with reduced required inputs for each operating point:
% altitude [ft], Mach number, required net thrust [N].
op(1).name = "SLS_takeoff"; op(1).alt_ft = 0;      op(1).M0 = 0.00; op(1).Freq = 120e3;
op(2).name = "climb";       op(2).alt_ft = 19685;  op(2).M0 = 0.40; op(2).Freq = 60e3;
op(3).name = "cruise";      op(3).alt_ft = 36089;  op(3).M0 = 0.78; op(3).Freq = 25e3;
end
