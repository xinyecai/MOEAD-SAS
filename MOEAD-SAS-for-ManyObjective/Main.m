%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  The  Many objective verstion of the MOEA/D-SAS
%%
%%  See the details of MOEA/D-SAS in the following paper:
%%
%%  Xinye Cai , Zhixiang Yang, Zhun Fan, Qingfu Zhang
%%
%%  Decomposition-based-Sorting and Angle-based-Selection for Evolutionary Multiobjective Optimization
%%  IEEE Transactions on Cybernetics, 2016
%%
%%  The source code of MOEAD-SAS is implemented by Zhixiang Yang 
%%
%%  If you have any questions about the code, please contact: 
%%
%%  Zhixiang Yang at xiang052@163.com
%%  Xinye Cai at xinye@nuaa.edu.cn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc
 problems = {'dtlz1','dtlz2','dtlz3','dtlz4','dtlz5','dtlz6','dtlz7'};
plength = length(problems);
%total test run.
totalrun = 1;
% path('cec09',path);
% path('moalg',path);

global numrun nobj;
seed = 10;

%%%%%
nobj = 5;
%%%%%

for i=1:7
    prob = problems{i};
    fprintf('Running on %s\n', prob);
    for j = 1:totalrun
        numrun=j;
        mop = test_MaOP(prob,0,nobj);
        pareto = moead(mop,'seed',j+seed);
        objpareto=[pareto.objective];
        fprintf('Run %d\n', j);
    end
end
