function mop = test_MaOP(testname,pardim,odim,k)
%Get test multi-objective problems from a given name. 
%   The method get testing or benchmark problems for Multi-Objective
%   Optimization. The test problem will be encapsulated in a structure,
%   which can be obtained by function get_structure('testmop'). 
%   User get the corresponding test problem, which is an instance of class
%   mop, by passing the problem name and optional dimension parameters.

% Copyright (c) NUAA 
% Create by Zhixiang Yang
% WFG instances are altered edition of that are created by Luping Wang
%
% History:
% 9/10/2015
    mop=get_structure('testmop');
    
    switch lower(testname)
        case 'dtlz1'
            k=5;
            pardim=odim+k-1;
            mop=dtlz1(mop,pardim,odim);
        case 'dtlz2'
            k=10;
            pardim=odim+k-1;
            mop=dtlz2(mop,pardim,odim);
        case 'dtlz3'
            k=10;
            pardim=odim+k-1;
            mop=dtlz3(mop,pardim,odim);
        case 'dtlz4'
            k=10;
            pardim=odim+k-1;
            mop=dtlz4(mop,pardim,odim);
        case 'dtlz5'
            k=10;
            pardim=odim+k-1;
            mop=dtlz5(mop,pardim,odim);
        case 'dtlz6'
            k=10;
            pardim=odim+k-1;
            mop=dtlz6(mop,pardim,odim);
        case 'dtlz7'
            k=20;
            pardim=odim+k-1;
            mop=dtlz7(mop,pardim,odim);
        otherwise 
            error('Undefined test problem name');                
    end 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%DTLZ1 function generator
function p=dtlz1(p,pardim,odim)
p.name='DTLZ1';
p.pd=pardim;
p.od=odim;
p.domain=[zeros(pardim,1) ones(pardim,1)];
p.func=@evaluate;

%DTLZ1 evaluation function.
    function y=evaluate(x,odim)
        [dim,num] = size(x);
        u         = x(odim:dim,:)-0.5;
        tmp       = 100*(dim-odim+1+sum(u.^2-cos(20*pi*u)));
        tmp       = (1+tmp)/2;
        y(1,:)    = tmp.*prod(x(1:odim-1,:),1);
        for i=2:odim-1
            y(i,:)= tmp.*prod(x(1:odim-i,:),1).*(1-x(odim-i+1,:));
        end
        y(odim,:) = tmp.*(1-x(1,:));
    end
end

%DTLZ2 function generator
function p=dtlz2(p,pardim,odim)
p.name='DTLZ2';
p.pd=pardim;
p.od=odim;
p.domain=[zeros(pardim,1) ones(pardim,1)];
p.func=@evaluate;

%DTLZ2 evaluation function.
    function y=evaluate(x,odim)
        [dim,num] = size(x);
        tmp       = sum((x(odim:dim,:)-0.5).^2);
        tmp       = 1+tmp;
        x(1:odim-1,:)=x(1:odim-1,:)*pi/2;
        y(1,:)    = tmp.*prod(cos(x(1:odim-1,:)),1);
        for i=2:odim-1
            y(i,:)= tmp.*prod(cos(x(1:odim-i,:)),1).*sin(x(odim-i+1,:));
        end
        y(odim,:) = tmp.*sin(x(1,:));
    end
end

%%% DTLZ3
function p=dtlz3(p,pardim,odim)
p.name='DTLZ3';
p.pd=pardim;
p.od=odim;
p.domain=[zeros(pardim,1) ones(pardim,1)];
p.func=@evaluate;

%DTLZ3 evaluation function.
    function y=evaluate(x,odim)
        [dim,num] = size(x);
        u         = x(odim:dim,:)-0.5;
        tmp       = 100*(dim-odim+1+sum(u.^2-cos(20*pi*u)));
        tmp       = 1+tmp;
        x(1:odim-1,:)=x(1:odim-1,:)*pi/2;
        y(1,:)    = tmp.*prod(cos(x(1:odim-1,:)),1);
        for i=2:odim-1
            y(i,:)= tmp.*prod(cos(x(1:odim-i,:)),1).*sin(x(odim-i+1,:));
        end
        y(odim,:) = tmp.*sin(x(1,:));
    end
end

%%% DTLZ4
function p=dtlz4(p,pardim,odim)
p.name='DTLZ4';
p.pd=pardim;
p.od=odim;
p.domain=[zeros(pardim,1) ones(pardim,1)];
p.func=@evaluate;

%DTLZ4 evaluation function.
    function y=evaluate(x,odim)
        alpha0=100;
        [dim,num] = size(x);
        tmp       = sum((x(odim:dim,:)-0.5).^2);
        tmp       = 1+tmp;
        x(1:odim-1,:)=x(1:odim-1,:).^alpha0*pi/2;
        y(1,:)    = tmp.*prod(cos(x(1:odim-1,:)),1);
        for i=2:odim-1
            y(i,:)= tmp.*prod(cos(x(1:odim-i,:)),1).*sin(x(odim-i+1,:));
        end
        y(odim,:) = tmp.*sin(x(1,:));
    end
end

%%% DTLZ5
function p=dtlz5(p,pardim,odim)
p.name='DTLZ5';
p.pd=pardim;
p.od=odim;
p.domain=[zeros(pardim,1) ones(pardim,1)];
p.func=@evaluate;

%DTLZ5 evaluation function.
    function y=evaluate(x,odim)
        [dim,num] = size(x);
        tmp       = sum((x(odim:dim,:)-0.5).^2);
        tmps      = tmp(ones(odim-2,1),:);
        x(1,:)    = x(1,:)*pi/2;
        x(2:odim-1,:)= pi/4*(1+2*tmps.*x(2:odim-1,:))./(1+tmps);
        tmp       = 1+tmp;
        y(1,:)    = tmp.*prod(cos(x(1:odim-1,:)),1);
        for i=2:odim-1
            y(i,:)= tmp.*prod(cos(x(1:odim-i,:)),1).*sin(x(odim-i+1,:));
        end
        y(odim,:) = tmp.*sin(x(1,:));
    end
end

%%% DTLZ6
function p=dtlz6(p,pardim,odim)
p.name='DTLZ6';
p.pd=pardim;
p.od=odim;
p.domain=[zeros(pardim,1) ones(pardim,1)];
p.func=@evaluate;

%DTLZ6 evaluation function.
    function y=evaluate(x,odim)
        [dim,num] = size(x);
        tmp       = sum(x(odim:dim,:).^0.1);
        tmps      = tmp(ones(odim-2,1),:);
        x(1,:)    = x(1,:)*pi/2;
        x(2:odim-1,:)= pi/4*(1+2*tmps.*x(2:odim-1,:))./(1+tmps);
        tmp       = 1+tmp;
        y(1,:)    = tmp.*prod(cos(x(1:odim-1,:)),1);
        for i=2:odim-1
            y(i,:)= tmp.*prod(cos(x(1:odim-i,:)),1).*sin(x(odim-i+1,:));
        end
        y(odim,:) = tmp.*sin(x(1,:));
    end
end

%%% DTLZ7
function p=dtlz7(p,pardim,odim)
p.name='DTLZ7';
p.pd=pardim;
p.od=odim;
p.domain=[zeros(pardim,1) ones(pardim,1)];
p.func=@evaluate;

%DTLZ7 evaluation function.
    function y=evaluate(x,odim)
        [dim,num] = size(x);
        tmp       = 1+9*mean(x(odim:dim,:));
        tmps      = tmp(ones(odim-1,1),:);
        y(1:odim-1,:)  = x(1:odim-1,:);
        tmp1      = y(1:odim-1,:).*(1+sin(3*pi*y(1:odim-1,:)))./(1+tmps);
        y(odim,:)    = (1+tmp).*(odim-sum(tmp1));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The following functions are used in WFG functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
