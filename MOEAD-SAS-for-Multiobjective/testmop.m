function mop = testmop( testname, dimension )
%Get test multi-objective problems from a given name. 
%   The method get testing or benchmark problems for Multi-Objective
%   Optimization. The test problem will be encapsulated in a structure,
%   which can be obtained by function get_structure('testmop'). 
%   User get the corresponding test problem, which is an instance of class
%   mop, by passing the problem name and optional dimension parameters.

    mop=get_structure('testmop');
    
    switch lower(testname)
        case {'uf1', 'uf2','uf3','uf4','uf5','uf6','uf7'}
            mop=cecproblems(mop, testname, dimension);
            mop.od=2;
        case {'uf8','uf9','uf10'}
            mop=cecproblems(mop, testname, dimension);
            mop.od=3;
        otherwise 
            error('Undefined test problem name');                
    end 
end
%cec09 UF1 - UF10
function p=cecproblems(p, testname,dim)
 p.name=upper(testname);
 p.pd=dim;
 
 p.domain=xboundary(upper(testname),dim);
 %p.domain = [zeros(dim,1),ones(dim,1)];
 p.func=cec09(upper(testname));
end

%cec09 UF11 - UF13
function p=cecproblems2(p, testname,dim)
 p.name=upper(testname);
 p.pd=dim;
 p.od=2;
 
 p.domain=xboundary(upper(testname),dim);
 %p.domain = [zeros(dim,1),ones(dim,1)];
 p.func=cec09m(upper(testname));
end




