%% startup_Recycle.m
%
% startup file for the HyBR recycling codes

directory = pwd;
path(directory, path)

path([directory, '/AIRToolsII-master'], path)
AIRToolsII_setup

path([directory, '/IRtools-master'], path)
IRtools_setup

path([directory, '/HyBR'], path)
path([directory, '/HyBRrecycle'], path)
clear directory




