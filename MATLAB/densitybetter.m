% h = (0:1400) * 1e3;
% rho = getDensity(h);
% semilogx(rho, h, LineWidth=2)
% grid on


M_CH4 = 12 + 4; 
M_Ar = 40;
M_N2 = 28;
M = 0.03*M_CH4 + 0.02*M_Ar + 0.95*M_N2;

NA = 6.02214076e23; % /mol

getDensity(0)

function rho_vec = getDensity(h_vec)
Z = h_vec;
kappa = 0.625;
A1 = 0.240;
A2 = 0.006;
A3 = 0.020;
Zh = 1050; % km
RT = 2574.73; % km
x = 1.76e5 * (Z-Zh) ./ ( (RT + Zh) * (RT+ Z) );
f_ch4 = A1 * (1+exp((1-kappa)*x))^(3 / (7 * (1- kappa))) + A2;
f_ar  = A3 * (1+exp((1-kappa)*x))^(-0.3/(1-kappa));

disp(f_ch4)
disp(f_ar)

return;
end