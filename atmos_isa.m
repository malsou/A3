function atm = atmos_isa(alt_ft)
%ATMOS_ISA Standard-day ISA atmosphere (troposphere + tropopause).
% Valid for altitudes up to ~11 km / 36 kft.
%
% Inputs:
%   alt_ft : geometric altitude [ft]
% Outputs struct atm:
%   T [K], P [Pa], rho [kg/m^3], a [m/s], h_m [m]

% Constants
T0 = 288.15;           % K
P0 = 101325;           % Pa
g0 = 9.80665;          % m/s^2
R  = 287.05287;        % J/(kg-K)
gamma = 1.4;
L  = -0.0065;          % K/m (troposphere lapse rate)
h_trop = 11000;        % m

h = max(0, alt_ft * 0.3048);

if h <= h_trop
    T = T0 + L*h;
    P = P0 * (T/T0)^(-g0/(L*R));
else
    T11 = T0 + L*h_trop;
    P11 = P0 * (T11/T0)^(-g0/(L*R));
    T = T11;
    P = P11 * exp(-g0*(h-h_trop)/(R*T));
end

rho = P/(R*T);
a = sqrt(gamma*R*T);

atm = struct('T',T,'P',P,'rho',rho,'a',a,'h_m',h);
end
