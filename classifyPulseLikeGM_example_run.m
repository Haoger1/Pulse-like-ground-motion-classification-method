clear; clc;

load('sample_SWresults.mat');  % This file should contain SWresults.

results = classifyPulseLikeGM(SWresults, true);

disp(struct2table(results));