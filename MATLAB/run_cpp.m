%%
outputDir = "../output/ ";
outputFileName = "entryYelle.csv";
eigenpath = "-I /Users/sebastienhenry/Documents/GTSAM/gtsam/gtsam/3rdparty/Eigen/";
runc(outputDir, outputFileName, eigenpath)

%%
function runc(outputDir, outputFileName, eigenpath)
currentFolder = pwd;
[status, cmdout] = system("cd " + pwd);
disp(cmdout)
ccodepath = "../ccode/ccode/";
[status, cmdout] = system("clang++ -std=c++11 -stdlib=libc++ " + ...
    eigenpath + " " + ...
    ccodepath+"main.cpp " + ...
    ccodepath+"atmospheric_model.cpp " + ...
    ccodepath+"EOM.cpp " + ...
    ccodepath+"integrator.cpp " + ...
    ccodepath+"vec_functions.cpp " + ...
    ccodepath+"vehicle.cpp " + ...
    "-o ../build/to_exec");
disp(cmdout)
[status, cmdout] = system("../build/to_exec " + outputDir + " " + outputFileName);
disp(cmdout)
end
