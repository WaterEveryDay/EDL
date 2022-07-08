entry = table2array(readtable("../output/entry_exp.csv"));
figure
plot(entry(:,2), entry(:,4)/1000)
