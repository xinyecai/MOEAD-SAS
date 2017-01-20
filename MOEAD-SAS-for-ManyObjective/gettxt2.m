
clear
clc
path('moalg',path);
problems = {'uf1','uf2','uf3','uf4','uf5','uf6','uf7','uf8','uf9','uf10',...
    'dtlz1','dtlz2','dtlz3','dtlz4','dtlz5','dtlz6','dtlz7'};
for k=1:17
    prob=problems{k};
    disp(prob);
    objs=[];
    objs2=[];
    t=[];
    if k<=7
        numobj=2;
    else
        numobj=3;
    end
    tmp=zeros(1,numobj);

    for i=1:30
        load(sprintf('data/%s/sub_run%d.mat',prob,i));
        
%         for j=1:length(cfsave)
            objpareto=cfsave{end};
            if k==14
                objpareto=max(objpareto,1e-30);
            end
            
            %%%%%%%%%%%
            S = size(objpareto);
            xtemp = ones(S);
            [objpareto,x] = ParetoFilter(objpareto,xtemp);
            %%%%%%%%%%%
            len=size(objpareto,2);
            tmp(end)=len;
            objs(end+1,:)=tmp;
            objs(end+1:end+len,:)=objpareto';
            objs2(end+1:end+len,:)=objpareto';
            
%             xtemp = ones(size(objpareto));
%             [CF,x] = ParetoFilter(objpareto,xtemp);
%             dlmwrite(sprintf('data/%s.txt',prob), objpareto', '-append', 'delimiter', ' ');    
%         end
    end
max(objs2)
    dlmwrite(sprintf('data/%s.txt',prob), objs, 'delimiter', ' '); 
%     save(sprintf('data/%s.txt',prob),'objs');
end
