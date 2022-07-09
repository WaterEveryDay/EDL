run_c

entry = table2array(readtable("../output/entry1.csv"));
entry_rec = table2array(readtable("../output/reconstructed_data.csv", 'VariableNamingRule','preserve'));

close all
figure
plot(entry(:,2)/1000, entry(:,4)/1000, LineWidth=2)
hold on
plot(entry_rec(:,2)/1000, entry_rec(:,3), "*", LineWidth=2)
hold off
xlabel("velocity (km/s)", "FontSize",15)
ylabel("altitude (km)", "FontSize",15)
legend("integration", "reconstruction (paper)", Location="best")

figure
plot(entry(:,1), entry(:,4)/1000, LineWidth=2)
hold on
plot(entry_rec(:,1), entry_rec(:,3), "*", LineWidth=2)
hold off
xlabel("time (s)", "FontSize",15)
ylabel("altitude (km)", "FontSize",15)
legend("integration", "reconstruction (paper)", Location="best")
xlim([-268.48, 0])
