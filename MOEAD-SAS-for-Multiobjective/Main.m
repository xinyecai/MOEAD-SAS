%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  The  Multiobjective verstion of the MOEA/D-SAS
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

problems = {'uf1','uf2','uf3','uf4','uf5','uf6','uf7','uf8','uf9','uf10'};
%total test run.
totalrun = 1;

path('cec09',path);
global numrun;
seed = 10;

for i=1:10
    prob = problems{i};
    fprintf('Running on %s\n', prob);
    for j = 1:totalrun
        numrun=j;
        if i>10
            dim=10;
        else
            dim=30;
        end
        mop = testmop(prob,dim);
        pareto = moead(mop,'seed',j+seed);
        objpareto=[pareto.objective];
    end
end
