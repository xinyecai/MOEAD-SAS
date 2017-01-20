function [obj,cons] = subobjective(weight, ind, idealpoint, method)%,cons
%SUBOBJECTIVE function evaluate a point's objective with a given method of
%decomposition. 

%   Two method are implemented by far is Weighted-Sum and Tchebesheff.
%   weight: is the decomposition weight.(column wise vector).
%   ind: is the individual point(column wise vector).
%   idealpoint: the idealpoint for Tchebesheff decomposition.
%   method: is the decomposition method, the default is 'te' when is
%   omitted.
%   
%   weight and ind can also be matrix. in which have two scenairos:
%   When weight is a matrix, then it's treated as a column wise set of
%   weights. in that case, if ind is a size 1 column vector, then the
%   subobjective is computed with every weight and the ind; if ind is also
%   a matrix of the same size as weight, then the subobjective is computed
%   in a column-to-column, with each column of weight computed against the
%   corresponding column of ind. 
%   A row vector of subobjective is return in both case.

%     if (nargin==2)
%         obj = ws(weight, ind);
%     elseif (nargin==3)
%         obj = te(weight, ind, idealpoint);
%     else
    global subproblems;

        if strcmp(method, 'ws')
            obj=ws(weight, ind);
        elseif strcmp(method, 'te')
            obj=te(weight, ind, idealpoint);
        elseif strcmp(method,'bi')
            obj= bi(weight, ind, idealpoint);
        end
        %%%%%%%%%%%%%%%%%%%%%%%
%         points=[subproblems.curpoints];
%         subVs=sum([points.obj_cons]);
%         Vmin=min(subVs);
%         Vmax=max(subVs);
%         T=Vmin+0.3*(Vmax-Vmin);
%         Vx=sum(cons);
%         tag=Vx<T;
%         s1=0.01;
%         s2=20;
%         obj(tag)=obj(tag)+s1*Vx(tag).^2;
%         obj(~tag)=obj(~tag)+s1*T^2+s2*(Vx(~tag)-T);
        %%%%%%%%%%%%%%%%%%%%%%%
%     end
end

% weighted sum scalarization Function
function obj = ws(weight, ind)
    obj=sum(weight.*ind);
end

% bi Scalarization Function
function d1 = bi(weight, ind, idealpoint)
    si=0.5;
    tmp1=sum((ind-repmat(idealpoint,1,size(ind,2))).*weight,1);
    tmp2=sqrt(sum(tmp1.^2,1));
    tmp3=sqrt(sum(weight.^2,1));
    d1=tmp2./tmp3;
    tmp4=ind-repmat(idealpoint,1,size(ind,2))-repmat(d1,size(weight,1),1).*weight;
%     tmp4=ind-repmat(idealpoint,1,size(ind,2))-repmat(d1./mweight,size(weight,1),1).*weight;
    d2=sqrt(sum(tmp4.^2,1));
    obj=d1+si*d2;
end

% function [d1,d2] = bi(weight, ind)
% global subproblems;
%     points=[subproblems.curpoints];
%     idealpoint1=max([points.objective],[],2);
% 
%     si=1;
%     tmp1=sum((repmat(idealpoint1,1,size(ind,2))-ind).*weight,1);
%     tmp2=sqrt(sum(tmp1.^2,1));
%     tmp3=sqrt(sum(weight.^2,1));
%     d1=tmp2./tmp3;
%     d1=100-d1;
%     tmp4=ind-repmat(idealpoint1,1,size(ind,2))+repmat(d1,size(weight,1),1).*weight;
%     d2=sqrt(sum(tmp4.^2,1));
% %     obj=d1-si*d2;
% end

% Techbycheff Scalarization Function
function obj = te(weight, ind, idealpoint)
    s = size(weight, 2);
    indsize = size(ind,2);
    
    weight((weight == 0))=0.0001;
    
    if indsize==s 
        part2 = abs(ind-idealpoint(:,ones(1, indsize)));
        obj = max(weight.*part2);
    elseif indsize ==1
        part2 = abs(ind-idealpoint);
        obj = max(weight.*part2(:,ones(1, s)));         
    else
        error('individual size must be same as weight size, or equals 1');
    end
end