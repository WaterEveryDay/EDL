%run_c

entryYelle = table2array(readtable("../output/entryYelle.csv"));
entryExp = table2array(readtable("../output/entryExp.csv"));
entryTrue = table2array(readtable("../output/reconstructed_data.csv", 'VariableNamingRule','preserve'));

modelYelle = table2array(readtable("../output/atm_model_entryYelle.csv"));
modelExp = table2array(readtable("../output/atm_model_entryExp.csv"));

close all
%%
figure
tiledlayout(1, 2, "TileSpacing","compact", Padding="compact")
nexttile
semilogx(modelYelle(:,2), modelYelle(:,1)/1000, LineWidth=2)
hold on
semilogx(modelExp(:,2), modelExp(:,1)/1000, '--', LineWidth=2)
hold off
xlabel("density (kg/m^3)", "FontSize",15)
ylabel("altitude (km)", "FontSize",15)
legend("Yelle Fit", "Exponential", Location="best", FontSize=15)

nexttile
plot(modelYelle(:,3), modelYelle(:,1)/1000, LineWidth=2)
hold on
plot(modelExp(:,3), modelExp(:,1)/1000, '--', LineWidth=2)
hold off
xlabel("T (K)", "FontSize",15)
ylabel("altitude (km)", "FontSize",15)
legend("Yelle Fit", "Exponential", Location="best", FontSize=15)

%%
figure
tiledlayout(1, 2, "TileSpacing","compact", Padding="compact")
nexttile
semilogx(modelYelle(:,5), modelYelle(:,1)/1000, LineWidth=2)
hold on
semilogx(modelExp(:,5), modelExp(:,1)/1000, '--', LineWidth=2)
plot([1e-3, 1e-3], [0, 1300], 'k--')
plot([10, 10], [0, 1300], 'k--')
hold off
xlabel("Kn", "FontSize",15)
ylabel("altitude (km)", "FontSize",15)
legend("Yelle Fit", "Exponential", Location="best", FontSize=15)
ylim([0,1300])
nexttile
plot(modelYelle(:,6), modelYelle(:,1)/1000, LineWidth=2)
hold on
plot(modelExp(:,6), modelExp(:,1)/1000, '--', LineWidth=2)
plot([1.4,2.2], [386.4815, 386.4815], 'k--')
plot([1.4,2.2], [862.476, 862.476], 'k--')
hold off
xlabel("C_D", "FontSize",15)
ylabel("altitude (km)", "FontSize",15)
legend("Yelle Fit", "Exponential", Location="best", FontSize=15)
ylim([0,1300])
xlim([1.4,2.2])

%%
figure
plot(entryYelle(:,2)/1000, entryYelle(:,4)/1000, LineWidth=2)
hold on
plot(entryExp(:,2)/1000, entryExp(:,4)/1000, '--', LineWidth=2)
plot(entryTrue(:,2)/1000, entryTrue(:,3), "*", color='#77AC30', LineWidth=3)
hold off
xlabel("velocity (km/s)", "FontSize",15)
ylabel("altitude (km)", "FontSize",15)
legend("integration with Yelle Fit", "integration with Exp Fit", "reconstruction (paper)", Location="best", fontsize=12)

figure
plot(entryYelle(:,1), entryYelle(:,4)/1000, LineWidth=2)
hold on
plot(entryExp(:,1), entryExp(:,4)/1000, '--', LineWidth=2)
plot(entryTrue(:,1), entryTrue(:,3), "*", color='#77AC30', LineWidth=3)
hold off
xlabel("time (s)", "FontSize",15)
ylabel("altitude (km)", "FontSize",15)
legend("integration with Yelle Fit", "integration with Exp Fit", "reconstruction (paper)", Location="best", fontsize=12)
xlim([-268.48, 0])

figure
plot(entryYelle(:,1), entryYelle(:,2)/1000, LineWidth=2)
hold on
plot(entryExp(:,1), entryExp(:,2)/1000, '--', LineWidth=2)
plot(entryTrue(:,1), entryTrue(:,2)/1000, "*", color='#77AC30', LineWidth=3)
hold off
xlabel("time (s)", "FontSize",15)
ylabel("vel (km/s)", "FontSize",15)
legend("integration with Yelle Fit", "integration with Exp Fit", "reconstruction (paper)", Location="best", fontsize=12)
xlim([-268.48, 0])


disp("altitude at t=0 Exp =" + 340972/1000)
disp("altitude at t=0 Yelle =" + 157123/1000)
disp("altitude at t=0 truth =" + 155.84331)

range_yelle = sum(entryYelle(:,2)/1000 .* cos(entryYelle(:,3)))*0.1;
disp("range yelle = " + range_yelle + "km")
