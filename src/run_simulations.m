close all;
clearvars;

run_GaAs = true;
run_Si = true;
run_Ge = true;

if run_GaAs
    run('main_GaAs.m')
end

if run_Si
    run('main_Si.m')
end

if run_Ge
    run('main_Ge.m')
end
