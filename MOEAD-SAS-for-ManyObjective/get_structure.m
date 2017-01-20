function str = get_structure( name )
%STRUCTURE Summary of this function goes here
% 
% Structure used in this toolbox.
% 
% individual structure:
% parameter: the parameter space point of the individual. it's a column-wise
% vector.
% objective: the objective space point of the individual. it's column-wise
% vector. It only have value after evaluate function is called upon the
% individual.
%
% subproblem structure:
% weight: the decomposition weight for the subproblem.
% neighbour: the index of the neighbour of this subproblem.
% optimal: the current optimal value of the current structure. 
% curpoiont: the current individual of the subproblem.
% oldpoint: the backup individual of the subproblem used to caculate the
% utility.
% utility: the searching utility of this subproblem.
%
% testmop structure:
% name: the name of the test problem.
% od: the number of objective.
% pd: the number of variable.
% domain: the domain, which is a pd*2 matrix, with domain(:,1) and
% domain(:,2) to specificy the lower and upper limit on every variable.
%
% parameter strucutre.
% the parameter setting structure is explained in loadpapams.m
%
%

switch name
    case 'individual' 
        str = struct('parameter',[],'objective',[]);
    case 'statis'
        str = struct('i',[],'c',[],'d',[]);
    case 'subproblem' 
        str = struct('weight',[],'neighbour',[],'optimal', Inf, ...
        'curpoint', [], 'utility', 1, 'oldpoint',[],'indexes',[]);
    case 'testmop'
        str = struct('name',[],'od',[],'pd',[],'domain',[],'func',[]);
    case 'parameter'
        str = struct('popsize',300,'niche',20,'dmethod', 'ts', ...
        'iteration', 2500,'evaluation', 300000, 'F',0.5, 'CR', 1, 'seed', 0, ...
        'selportion', 5, 'updatenb', 2, 'updateprob',0.9, 'dynamic',1, 'decayrate', 0.95);
    otherwise
        error('the structure name requried does not exist!');
end