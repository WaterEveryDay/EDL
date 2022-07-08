run_c

entry = table2array(readtable("../output/entry1.csv"));
close all 
figure
plot(entry(:,2)/1000, entry(:,4)/1000, LineWidth=2)
xlabel("velocity (km/s)", "FontSize",15)
ylabel("altitude (km)", "FontSize",15)
