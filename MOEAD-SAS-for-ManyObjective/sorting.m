function sorting(pop, L)
%%% 相应于子问题目标值小的子问题优先
    global subproblems idealpoint counters EP ind2subp Lev xxx;
    Indivs=[pop.objective];
    lenpop=size(Indivs,2);
    lensub=length(subproblems);
    [subproblems.indexes]=deal([]);
    counters=zeros(1,lensub);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Association
    ind2subp0=zeros(1,lenpop);
    divers_cell=cell(1,lensub);
    nums=zeros(1,lensub);
    perdist_cell=cell(1,lensub);
    weights=[subproblems.weight];
    mweights=sqrt(sum(weights.^2));
    for i=1:lenpop
        ind=Indivs(:,i)-idealpoint;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        mind=sqrt(sum(ind.^2));
        costhita=ind'*weights./mweights/mind;
        [perdist0, index_i]=min(acos(costhita));

%         dist2=ind'*weights./(mweights).^2;
%         pervct=ind(:,ones(1,lensub))-dist2(ones(1,length(ind)),:).*weights;
%         perdist=sqrt(sum(pervct.^2));
%         [~, index_i]=min(perdist);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        divers_cell{index_i}(end+1)=i;
        perdist_cell{index_i}(end+1)=perdist0;
        nums(index_i)=nums(index_i)+1;
        ind2subp0(i)=index_i;
    end
    
    for i=1:lensub
        [~,index_dist]=sort(perdist_cell{i});
        divers_cell{i}=divers_cell{i}(index_dist);
    end
    
    %%% Sorting
    allindexes=zeros(L,lensub);
    for i=1:lensub
        %%% get L solutions
        neighbour0=subproblems(i).neighbour;
        L1=min(L,length(subproblems));
        s_nums=sum(nums(neighbour0(1:L1)));
        while s_nums<L
            L1=L1+1;
            s_nums=s_nums+nums(neighbour0(L1));
        end
        neighbour0=neighbour0(1:L1);
        inx_tmp=cell2mat(divers_cell(neighbour0));
        inx_tmp=inx_tmp(1:L);
        Indivs1=Indivs(:,inx_tmp);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        weight0=subproblems(i).weight;
        weight0(weight0==0)=1e-6;
        weight01=1./weight0;
        weight1=weight01/sum(weight01);
        
        %%%
        weight1=weight0;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        weights=weight1(:,ones(1,L));
        objs=subobjective(weights, Indivs1, idealpoint, 'bi');
        [~, inx_sorted]=sort(objs);
        allindexes(:,i)=inx_tmp(inx_sorted);
    end
    
    %%% Selection
    indexes=[];
    for i=1:L
        rindexes0=allindexes(i,:);
        rindexes = setdiff(rindexes0,indexes);
        if length(indexes)+length(rindexes)<lensub
            indexes=[indexes, rindexes];
        elseif length(indexes)+length(rindexes)>lensub
            inds_ridx=Indivs(:,rindexes)-idealpoint(:,ones(1,length(rindexes)));
            inds_idx=Indivs(:,indexes)-idealpoint(:,ones(1,length(indexes)));
            thetas=getthita(inds_ridx, inds_idx);
            while length(indexes)<lensub
                [~,index_tmp]=max(thetas);
                index_r=rindexes(index_tmp);
                indexes(end+1)=index_r;
                rindexes(index_tmp)=[];
                %%%
                thetas(index_tmp)=[];
                inds_i=inds_ridx(:,index_tmp);
                inds_ridx(:,index_tmp)=[];
                thetas0=getthita(inds_ridx, inds_i);
                thetas=min(thetas,thetas0);
            end
            break;
        else
            indexes=[indexes, rindexes];
            break;
        end
    end
    
    %%%
    Lev=i;
    EP=pop(indexes);
    ind2subp=ind2subp0(indexes);
    
    %%%
    for i=1:lensub
        subproblems(i).indexes=find(ind2subp==i);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
end

function [thetas,indexes]=getthita(inds,weights)
    leninds=size(inds,2);
    thetas=zeros(1,leninds);
    indexes=zeros(1,leninds);
    mweights=sqrt(sum(weights.^2));
    for i=1:leninds
        ind=inds(:,i);
        mind=sqrt(sum(ind.^2));
        costhita=ind'*weights./mweights/mind;
        [thetas(i),indexes(i)]=min(acos(costhita));
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%
% %distance
% function [thetas,indexes]=getthita(inds,weights)
%     leninds=size(inds,2);
%     thetas=zeros(1,leninds);
%     indexes=zeros(1,leninds);
%     mweights=sqrt(sum(weights.^2));
%     for i=1:leninds
%         indsm=inds(:,i*ones(1,size(weights,2)))-weights;
%         dist=sqrt(sum(indsm.^2,1));
%         
%         [thetas(i),indexes(i)]=min(dist);
%     end
% end



