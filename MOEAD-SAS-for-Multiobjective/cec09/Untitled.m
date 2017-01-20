clc
% load the PF* data
PFStar  = importdata('pf_data/CF10.dat');
PFStar  = PFStar';

% randomly generate 1000 points inside the search space
PF  = importdata('b.dat');
% PF  = PF';
str  = sprintf('%.5f = IGD(PFStar, PF, C)', IGD(PFStar, PF));
disp(str);

disp('----------------------------------------');
