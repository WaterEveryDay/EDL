h = (0:1400) * 1e3;
T = getTemp(h);
semilogx(T, h, LineWidth=2)
grid on


function T_vec = getTemp(h_vec)
T0 = 93.28;
T1 = 70.66;
T2 = 151.21;
T3 = 178.19;
T4 = 135.49;
T5 = 175.00;

h0 = 0;
h1 = 40e3;
h2 = 120e3;
h3 = 300e3;
h4 = 550e3;
h5 = 900e3;

T_vec = zeros(size(h_vec));
for i = 1:length(h_vec)
    T = 175.00;
    h = h_vec(i);
    if (h < h1)
        Lh0 = (T1-T0)/(h1-h0);
        T = T0 + Lh0 * (h - h0);
    elseif (h < h2)
        Lh1 = (T2-T1)/(h2-h1);
        T = T1 + Lh1 * (h - h1);
    elseif (h < h3)
        Lh2 = (T3-T2)/(h3-h2);
        T = T2+ Lh2 * (h - h2);
    elseif (h < h4)
        Lh3 = (T4-T3)/(h4-h3);
        T = T3+ Lh3 * (h - h3);
    elseif (h < h5)
        Lh4 = (T5-T4)/(h5-h4);
        T = T4+ Lh4 * (h - h4);
    end
    T_vec(i) = T;
end
return;
end