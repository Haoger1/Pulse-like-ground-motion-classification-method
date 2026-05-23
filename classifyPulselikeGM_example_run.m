clear; clc;

load('SWresults_sample.mat');  % This file should contain SWresults.

results = classifyPulselikeGM(SWresults_sample, true);

disp(struct2table(results));
