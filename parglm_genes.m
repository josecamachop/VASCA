function parglmo = parglm_genes(X, F, interactions, center, n_perm)

% Parallel General Linear Model to obtain multivariate factor and interaction 
% matrices in a crossed experimental design and permutation test for multivariate 
% statistical significance based on ASCA-genes. 
%
% Related routines: asca, apca, parglm, parglmVS, parglmMC, create_design
%
% parglmo = parglm(X, F, interactions, prep, n_perm)   % complete call
%
%
% INPUTS:
%
% X: [NxM] billinear data set for model fitting, where each row is a
% measurement, each column a variable
%
% F: [NxF] design matrix, where columns correspond to factors and rows to
% levels.
%
% interactions: [Ix2] matrix where rows contain the factors for which
% interactions are to be calculated.
%
% prep: [1x1] preprocesing:
%       1: mean preping
%       2: autoscaling (default)
%
% n_perm: [1x1] number of permutations (1000 by default)
%
%
% OUTPUTS:
%
% parglmo (structure): structure with the factor and interaction
% matrices, p-values and explained variance 
%
%
% coded by: José Camacho (josecamacho@ugr.es)
% last modification: 18/Oct/22
%
% Copyright (C) 2022  José Camacho, Universidad de Granada
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

%% Arguments checking

% Set default values
routine=dbstack;
assert (nargin >= 2, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);
N = size(X, 1);
M = size(X, 2);
if nargin < 3 || isempty(interactions), interactions = []; end;
if nargin < 4 || isempty(center), center = 2; end;
if nargin < 5 || isempty(n_perm), n_perm = 1000; end;

% Validate dimensions of input data
assert (isequal(size(center), [1 1]), 'Dimension Error: 4th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(n_perm), [1 1]), 'Dimension Error: 5th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);


%% Main code
                  
n_interactions      = size(interactions,1);      % number of interactions
n_factors           = size(F,2);                 % number of factors
mtcc                = n_factors + n_interactions;       % correction for the number of tests
SSQ_factors         = zeros(n_perm*mtcc + 1,n_factors,1);      % sum of squares for factors
SSQ_interactions    = zeros(n_perm*mtcc + 1,n_interactions);   % sum of squares for interactions
p_factor            = zeros(n_factors,M);                % p-values factors
p_interaction       = zeros(n_interactions,M);           % p-values interactions

% In column space
parglmo.factors                = cell(n_factors,1);
parglmo.interactions           = cell(n_interactions,1);

% preprocess the data
[Xs,m,dt] = preprocess2D(X,center);
X = X./(ones(size(X,1),1)*dt);

SSQ_X = sum(sum(X.^2));

% Make structure with general 'variables'
parglmo.data           = X;
parglmo.design         = F;
parglmo.n_factors      = n_factors;
parglmo.n_interactions = n_interactions;

% Create Design Matrix
n = 1;
D = ones(size(X,1),1);

for f = 1 : n_factors
    uF = unique(F(:,f));
    paranovao.n_levels(f) = length(uF); 
    for i = 1:length(uF)-1
        D(find(F(:,f)==uF(i)),n+i) = 1;
    end
    parglmo.factors{f}.Dvars = n+(1:length(uF)-1);
    D(find(F(:,f)==uF(end)),parglmo.factors{f}.Dvars) = -1;
    n = n + length(uF) - 1;
end

for i = 1 : n_interactions
    for j = parglmo.factors{interactions(i,1)}.Dvars
        for k = parglmo.factors{interactions(i,2)}.Dvars
            D(:,end+1) = D(:,j).* D(:,k);
        end
    end
    parglmo.interactions{i}.Dvars = n+1:size(D,2);
    n = size(D,2);
end
    
% GLM model calibration with LS, only fixed factors
pD =  pinv(D'*D)*D';
B = pD*X;
X_residuals = X - D*B;
parglmo.D = D;
parglmo.B = B;

% Create Effect Matrices

parglmo.inter = D(:,1)*B(1,:);
SSQ_inter = sum(sum(parglmo.inter.^2));

for f = 1 : n_factors
    parglmo.factors{f}.matrix = D(:,parglmo.factors{f}.Dvars)*B(parglmo.factors{f}.Dvars,:);
    SSQ_factors(1,f) = sum(sum(parglmo.factors{f}.matrix.^2));
    %var_pca(parglmo.factors{f}.matrix); 
    %parglmo.factors{f}.pcs = input(sprintf('Number of PCS (1:n) for the factor %d: ',f));
    parglmo.factors{f}.pcs = 1; % 1 PC for the examples of the paper
end

% Interactions
for i = 1 : n_interactions
    parglmo.interactions{i}.matrix = D(:,parglmo.interactions{i}.Dvars)*B(parglmo.interactions{i}.Dvars,:);
    SSQ_interactions(1,i) = sum(sum(parglmo.interactions{i}.matrix.^2));
    var_pca(parglmo.interactions{i}.matrix); 
    parglmo.interactions{i}.pcs = input(sprintf('Number of PCS (1:n) for the interaction %d: ',i));
    
end

SSQ_residuals = sum(sum(X_residuals.^2));
parglmo.effects = 100*([SSQ_inter SSQ_factors(1,:) SSQ_interactions(1,:) SSQ_residuals]./SSQ_X);
parglmo.residuals = X_residuals;

for f = 1 : n_factors
    factors{f}.D = [];
end
for i = 1 : n_interactions
    interacts{i}.D = [];
end
        
for j = 1 : n_perm*mtcc
    
    perms = randperm(size(X,1)); % permuted data (permute whole data matrix)
      
    B = pD*X(perms, :);
    
    for f = 1 : n_factors
        factors{f}.matrix = D(:,parglmo.factors{f}.Dvars)*B(parglmo.factors{f}.Dvars,:);
        SSQ_factors(1 + j,f) = sum(sum(factors{f}.matrix.^2));
        T = loadings_pca(factors{f}.matrix,parglmo.factors{f}.pcs,0,0);
        factors{f}.D = [ factors{f}.D diag(T*T')];
    end
    
    % Interactions
    for i = 1 : n_interactions
        interacts{i}.matrix = D(:,parglmo.interactions{i}.Dvars)*B(parglmo.interactions{i}.Dvars,:);
        SSQ_interactions(1 + j,i) = sum(sum(interacts{i}.matrix.^2));
        T = loadings_pca(interacts{i}.matrix,parglmo.interactions{i}.pcs,0,0);
        interacts{i}.D = [interacts{i}.D diag(T*T')];
    end

end        % permutations



for f = 1 : n_factors
    parglmo.factors{f}.UCD = prctile(factors{f}.D,99);
    T = loadings_pca(parglmo.factors{f}.matrix,parglmo.factors{f}.pcs,0,0);
    parglmo.factors{f}.D = diag(T*T');
    parglmo.factors{f}.E = sum((parglmo.factors{f}.matrix - parglmo.factors{f}.matrix*T*T').^2,1);
    parglmo.factors{f}.E = diag(parglmo.factors{f}.E'*parglmo.factors{f}.E);
    m = mean(parglmo.factors{f}.E);
    v = var(parglmo.factors{f}.E);
    parglmo.factors{f}.UCQ = v/(2*m)*chi2inv(0.99,2*m^2/v);
    for m = 1 : M      
            p_factor(f,m) = (size(find( factors{f}.D(:,m) ...
                >= parglmo.factors{f}.D(m)), 1) + 1)/( n_perm*mtcc+1); 
    end

end
for i = 1 : n_interactions
    parglmo.interactions{i}.UCD = prctile(interacts{i}.D,99);
    T = loadings_pca(parglmo.interactions{i}.matrix,parglmo.interactions{i}.pcs,0,0);
    parglmo.interactions{i}.D = diag(T*T');
    parglmo.interaction{i}.E = sum((parglmo.interaction{i}.matrix - parglmo.interaction{i}.matrix*T*T').^2,1);
    parglmo.interaction{i}.E = diag(parglmo.interaction{i}.E'*parglmo.interaction{i}.E);
    m = mean(parglmo.interactions{i}.E);
    v = var(parglmo.interactions{i}.E);
    parglmo.interactions{i}.UCQ = v/(2*m)*chi2inv(0.99,2*m^2/v);
    for m = 1 : M      
            p_interaction(i,m) = (size(find( interacts{i}.D(:,m) ...
                >= parglmo.interaction{i}.D(m)), 1) + 1)/( n_perm*mtcc+1); 
    end
    

end

% Calculate univariate p-values and order variables by relevance
for factor = 1 : n_factors
    [~,ord_factors(factor,:)] = sort(p_factor(factor,:),'ascend');
end
for interaction = 1 : n_interactions
    [~,ord_interactions(interaction,:)] = sort(p_interaction(interaction,:),'ascend');
end

parglmo.ord_factors = ord_factors;
if n_interactions>0
    parglmo.ord_interactions = ord_interactions;
end

% Multiple test correction for several factors/interactions
p_factor = min(1,p_factor * mtcc); 
p_interaction = min(Inf,p_interaction * mtcc); 

parglmo.p = [p_factor' p_interaction'];

