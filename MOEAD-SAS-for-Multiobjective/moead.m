function pareto = moead( mop, varargin)
%MOEAD run moea/d algorithms for the given mop.
% MOP could be obtained by function 'testmop.m'.
% the controlling parameter can be passed in by varargin, the folliwing
% parameters are defined in here. More other parameters can be passed by
% modify loadparams.m problem by problem.
%   seed: the random seed.
%   popsize: The subproblem's size.
%   niche: the neighboursize, must less then the popsize.
%   evaluation: the total evaluation of the moead algorithms before finish.
%   dynamic: whether to use dynamic resource allocation.
%   selportion: the selection portion for the dynamic resource allocation

    %global variable definition.
    global subproblems params itrCounter evalCounter rnduni EP numrun;
    %global idealpoint objDim parDim evalCounter;
    
    %load the parameters.
    params=loadparams(mop, varargin);
    
    %%%%%
    cfsave=[];
    %%%%%
    
    %set the random generator.
%     rs = RandStream.create('mt19937ar', 'Seed', params.seed);
%     RandStream.setDefaultStream(rs);
    rand('seed', sum(100*clock));
    %seed = rem((params.seed+23),1377);
    %rnduni = -seed;
    
    %the counters.
    evalCounter = 0;
    itrCounter = 0;
    
    %and Initialize the algorithm.
    init(mop);
    
    while ~terminate()
        evolve(mop); % one generation of evaluation.
        itrCounter=itrCounter+1;
        %%%%%%%%%
        obs=[EP.objective];
        cfsave=[cfsave;{obs}];
        %%%%%%%%%
        
%         if (rem(itrCounter,50)==0) % updating of the utility.
%             util_update();
%         end
        if (rem(itrCounter,200)==0) % updating of the utility.
            status(mop);
        end
    end
    
    %display the result.subproblems(i).optimal=newobj(i);
    pareto=EP;
    %%%%%%
    save(sprintf('data/%s/sub_run%d.mat',mop.name,numrun),'cfsave');
    %%%%%%
end

% The evoluation setp in MOEA/D
function evolve(mop)
	global subproblems idealpoint params rnduni EP selectionSize Lev;

    % select the subproblem according to its utility.
    % if params.dynamic is not true, then no selection is used.
%     if (params.dynamic)
% %         selindex = util_select();
%         selSize = ceil(params.popsize/params.selportion);
%         selindex = randperm(params.popsize,selSize);
%     else
%         selindex = 1:length(subproblems);
%     end
  
    selindex = 1:length(EP);
    selectionSize = length(selindex);
    
    temp_ind = get_structure('individual');
    ind = repmat(temp_ind, 1, length(selindex));

    for i=1:length(selindex)
        index = selindex(i);

        %[r, rnduni]=crandom(rnduni);
        updateneighbour = rand < params.updateprob;
        
        %new point generation using genetic operations, and evaluate it.
        ind(i) = genetic_op(index, updateneighbour, mop.domain);
        [obj,ind(i)] = evaluate(mop, ind(i));
        
        %update the idealpoint.
        idealpoint = min(idealpoint, obj);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % update external population
    pop = [EP ind];
    L=min(max(Lev+params.niche, params.niche), length(pop));
    sorting(pop, L);
end

% update the index's neighbour with the given individual.
% index is the subproblem's index in the main population.
% ind is the individual structure.
% updatenieghbour is a bool determine whether the neighbourhood of index, or the whole population should be updated.
% this procedure is also governed by a parameter from params: params.updatenb, which determine how many subproblem
% should be updated at most by this new individual.
function update(index, ind, updateneighbour)
	global subproblems idealpoint params newold newdomold EP enterEP;
  
	% collect the updation index
    if (updateneighbour)
        updateindex = subproblems(index).neighbour;
    else
        updateindex = 1:length(subproblems);
    end
      
	updateindex = random_shuffle(updateindex);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	len11=length(updateindex);
    points11=[subproblems(updateindex).curpoint];
	updateweights = [subproblems(updateindex).weight];
    newobj=subobjective(updateweights, repmat(ind.objective,1,len11),idealpoint, 'te');
    oldobj=subobjective(updateweights,[points11.objective], idealpoint, 'te');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    upindex=newobj<oldobj;
%     newold=[newold; sum(upindex)];
    upindex=find(upindex==1);
    if length(upindex)>params.updatenb
        upindex=upindex(1:params.updatenb);
    end
    index11=updateindex(upindex);
    [subproblems(index11).curpoint]=deal(ind);
    %%%%%%%%%%%
  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 	news=ind.objective;
%     olds=point11.objective;
%     if ~ismember(0,news<=olds)&&ismember(1,news<olds)
%         newdomold=[newdomold;1];
%     elseif ~ismember(0,olds<=news)&&ismember(1,olds<news)
%         newdomold=[newdomold;-1];
%     else
%         newdomold=[newdomold;0];
%     end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   for i=1:length(updateindex)
%       idx = updateindex(i);
%       updateweight = subproblems(idx).weight;
%       
%       newobj=subobjective(updateweight, ind.objective,  idealpoint, 'ws');
%       old=subobjective(updateweight, subproblems(idx).curpoint.objective,  idealpoint, 'ws');
%       
%       if (newobj<old)
%          subproblems(idx).curpoint=ind;
%          %disp(sprintf('UP -- ID: %d, Type: %d, T: %d, UI: %d -- %f', index-1, updateneighbour, time, idx-1, ind.parameter));
%          time = time+1; 
%       end
%       if (time>=params.updatenb)
%           return;
%       end
%   end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     enterEP=[enterEP;0];
%   %update EP
%   ind1=get_structure('individual');
%   for i=1:size(ind.objective,2)
%       ind1.parameter=ind.parameter(:,i);
%       ind1.objective=ind.objective(:,i);
% %       ind1.obj_cons=ind.obj_cons(:,i);
% %       if ind1.obj_cons==0
% %       if ~any((ind1.obj_cons<=1e-6)-1)    %%%%%test
%           len=length(EP);
%           if len==0
%               EP=[EP ind1];
%           else
%               %remove the vectors be dominated by ind1
%               obj=[EP.objective];
%               tmp1=repmat(ind1.objective,1,len)<=obj;
%               index1=~any(tmp1-1);
%               EP(index1)=[];
%               obj(:,index1)=[];
%               %insert into EP
%               tmp2=obj<=repmat(ind1.objective,1,length(EP));
%               index2=any(tmp2-1);
%               if ~ismember(0,index2)
%                   EP=[EP ind1];
%                   enterEP(end)=1;
%               end
%           end
% %       end
%   end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % do the comparision.
%   updateweight = [subproblems(updateindex).weight];
%   oops = [subproblems(updateindex).curpoint];
%   newobj=subobjective(updateweight, ind.objective,  idealpoint, 'te');
%   oldobj=subobjective(updateweight, [oops.objective], idealpoint, 'te' );
%   C = newobj < oldobj;
%   
%   % find the one need to be updated in betterIndex.
%   betterIndex = find(C);
%   updateNumber = params.updatenb;
%   if (length(betterIndex)>updateNumber)
%       randp = randperm(length(betterIndex));
%       randp = randp(1:updateNumber);
%       betterIndex = betterIndex(randp);
%   end
%   
%   % do the updation.
%   [subproblems(updateindex(betterIndex)).curpoint]= deal(ind);
%   clear C newobj oops oldobj;
end

function y =terminate()
    global params evalCounter;
    y = evalCounter>params.evaluation;
end

function status(mop)
    global evalCounter itrCounter selectionSize subproblems params EP idealpoint;
    %average utility
    averageutil = mean([subproblems.utility]);
    fprintf('Itr:%d\tSel:%d\tEval:%d\tUtilM:%1.4f\n', ...
        itrCounter, selectionSize, evalCounter, averageutil);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     open(sprintf('PF_figure/%s.fig',mop.name));
%     hold on
%     length(EP)
%     if ~isempty(EP)
% %         pot=[subproblems.curpoint];
%         p=[EP.objective];
% %         cons=sum([pot.obj_cons],1);
%         if mop.od==2
%             scatter(p(1,:),p(2,:),'r');
%             axis([0 max([p(1,:) 1.2]) 0 max([p(2,:) 1.2])]);
%         else
%             scatter3(p(1,:),p(2,:),p(3,:),'r');
%             axis([0 max([p(1,:) 1.2]) 0 max([p(2,:) 1.2]) 0 max([p(3,:) 1.2])]);
%             view([-45 10]);
%             grid on
%         end
%     end
% %     title(mop.name);
%     legend('Real PF','MoeadID');
    
% %     set(h,'visible','off');
%     str=sprintf('figure/%s-%d',mop.name,itrCounter/50);
%     saveas(gcf,str);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  此部分可去除
end
