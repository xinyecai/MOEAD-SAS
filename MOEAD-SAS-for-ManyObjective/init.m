function init(mop)
%Set up the initial setting for the MOEA/D.
%loading the params.

    global subproblems params idealpoint parDim objDim Lev;
    idealpoint=ones(mop.od,1)*inf;
    parDim = mop.pd;
    objDim = mop.od;
    
    %%%
    pop=[];  %in order to compare easily
    %%%
    
    subproblems = init_weights(params.popsize, params.niche, mop.od);
    params.popsize = length(subproblems);
    
    %initial the subproblem's initital state.
    for i=1:params.popsize
        ind = randompoint(mop);
        [value, ind] = evaluate( mop, ind );
        idealpoint = min(idealpoint, value);
        pop=[pop,ind];
    end
    sorting(pop, length(pop));
    Lev=2*length(subproblems);
end

function subp=init_weights(popsize, niche, objDim)
% init_weights function initialize a pupulation of subproblems structure
% with the generated decomposition weight and the neighbourhood
% relationship.
    
    subp = [];
    
    %to load the parameter from the given file.
    switch objDim
        case 4
            popsize = 220;
        case 5
            popsize = 210;
        case 6
            popsize = 182;
        case 8
            popsize = 156;
        case 10
            popsize = 275;
        case 15
            popsize = 135;
    end 
    
    filename = sprintf('weight/W%uD_%u.dat',objDim, popsize);
    
    if (exist(filename, 'file'))
        % copy from the exist file. 
        p=get_structure('subproblem');
        subp=repmat(p, popsize,1);
        allws = importdata(filename);
        for i=1:popsize
            subp(i).weight = allws(i,:)';
        end
    else
        %generate on the fly.
        for (i=1:popsize)
            if (objDim==2)
                p=get_structure('subproblem');
                weight=zeros(2,1);
                weight(1)=i/popsize;
                weight(2)=(popsize-i)/popsize;
                p.weight=weight;
                subp=[subp p];
            elseif objDim==3
            %TODO
            end
        end
    end

    %Set up the neighbourhood.
    leng=length(subp);
    distanceMatrix=zeros(leng, leng);
    for i=1:leng
        for j=i+1:leng
            A=subp(i).weight;B=subp(j).weight;
            distanceMatrix(i,j)=(A-B)'*(A-B);
            distanceMatrix(j,i)=distanceMatrix(i,j);
        end
        [s,sindex]=sort(distanceMatrix(i,:));
        subp(i).neighbour=sindex';
    end
end