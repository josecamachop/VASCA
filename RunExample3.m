%% Third simulation example in Variable-Selection ANOVA Simultaneous Component Analysis. Bioinformatics. 2022 
% Camacho J, Vitale R, Morales-Jimenez D. and Gómez-Llorente C. 
%
% We simulate a single factor with two levels and 40 subjects for which 400
% variables or responses (e.g., -omics features) are collected. We generate 
% an additive multivariate relationship between the first three variables 
% and the levels. Therefore, unlike the previous example, we need to 
% consider the combination of the three variables to properly differentiate 
% the aforementioned levels. 
%
% coded by: Jose Camacho (josecamacho@ugr.es)
% last modification: 30/Oct/2022
%
% Copyright (C) 2022  University of Granada, Granada
% Copyright (C) 2022  Jose Camacho Paez
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

%% Simulation

clear

n_obs = 40;       % number of individuals 
n_vars = 400;     % number of responses or variables
rep = 1000;        % number of repetitions of the simulation

close all
p1 = zeros(rep,1);
p2 = zeros(rep,n_vars);
p2b = zeros(rep,n_vars);
p3 = zeros(rep,n_vars);
p22 = zeros(rep,n_vars);
p22b = zeros(rep,n_vars);
p33 = zeros(rep,n_vars);
parfor i= 1:rep % Main loop
    
    X = [randn(n_obs,3) simuleMV(n_obs,n_vars-3,7)]; % Variables from 4 to 400 are independent to the class, and obtained with simuleMV with a medium correlation level
    class = sign(sum(X(:,1:3),2)); % The factor is generated from the sum of the three variables: if the sum is positive, the factor is level 2, it otherwise is level 1 
    
    s = rng(i);
    [~,parglmo] = parglm(X,class); % general linear model (GLM) factorization and (ASCA type) multivariate significance testing
    rng(s);
    [~,parglmoVS] = parglmVS(X,class); % GLM factorization and (VASCA-type) incremental multivariate significance testing
    rng(s);
    [~,parglmoMC] = parglmMC(X,class); % GLM factorization and Benjamini-Hochberg (BH) univariate significant testing
    
    p1(i) = parglmo.p;
    
    % Sorted p-values
    p2(i,:) = parglmoVS.p(parglmoVS.ord_factors);
    p3(i,:) = parglmoMC.q(parglmoMC.ord_factors);
    
    % Unsorted p-values
    p22(i,:) = parglmoVS.p;
    p33(i,:) = parglmoMC.q;
    
    p2bb = p22(i,:);
    ind = find(p22(i,:)<0.01); % VASCA + bootstrapping
    if ~isempty(ind) 
        [~,parglmo] = parglm(X(:,ind),class); 
        ascao = asca(parglmo);
        bpvals = pbootasca(X(:,ind), class, ascao, 1, 1000, 0); 
        p2bb(ind) = bpvals;
    end
    % Unsorted
    p22b(i,:) = p2bb;
    % Sorted
    p2b(i,:) = p2bb(parglmoVS.ord_factors);
end

save example3

%% Plot Figures: using sorted p-values 

load example3

minT = 1e-3;
maxT = 1;
h = figure; hold on

p1(find(p1(:)<minT)) = minT;
p1(find(p1(:)>maxT)) = maxT;
p2(find(p2(:)<minT)) = minT;
p2(find(p2(:)>maxT)) = maxT;
p3(find(p3(:)<minT)) = minT;
p3(find(p3(:)>maxT)) = maxT;

mp3 = mean(p3);
plot(mp3,'k--')
plot([1 n_vars],[mean(p1) mean(p1)],'b-.')
mp2 = mean(p2);
plot(mp2,'g-o');
mp2b = mean(p2b);
plot(mp2b,'c-o')

plot([0 n_vars],[0.05 0.05],'r:')
plot([0 n_vars],[0.01 0.01],'r--')
legend('FDR','ASCA','VASCA','VASCA + bootstrapping','\alpha=0.05','\alpha=0.01','Location','southeast')

xr = [];
yr = [];
for i=1:size(p3,2)
    xr = [xr;i*ones(1,2)];
    yr = [yr;[mean(p3(:,i))-std(p3(:,i)) mean(p3(:,i))+std(p3(:,i))]];
end
yr(find(yr<minT)) = minT;
yr(find(yr>maxT)) = maxT;
fill([xr(:,1);flipud(xr(:,2))],[yr(:,1);flipud(yr(:,2))],'k','FaceAlpha',0.2,'EdgeColor','none');
xr = [];
yr = [];
for i=1:size(p2,2)
    xr = [xr;i*ones(1,2)];
    yr = [yr;[mean(p2(:,i))-std(p2(:,i)) mean(p2(:,i))+std(p2(:,i))]];
end
yr(find(yr<minT)) = minT;
yr(find(yr>maxT)) = maxT;
fill([xr(:,1);flipud(xr(:,2))],[yr(:,1);flipud(yr(:,2))],'g','FaceAlpha',0.2,'EdgeColor','none');

xr = [];
yr = [];
for i=1:size(p2,2)
    xr = [xr;i*ones(1,2)];
    yr = [yr;[mean(p2b(:,i))-std(p2b(:,i)) mean(p2b(:,i))+std(p2b(:,i))]];
end
yr(find(yr<minT)) = minT;
yr(find(yr>maxT)) = maxT;
fill([xr(:,1);flipud(xr(:,2))],[yr(:,1);flipud(yr(:,2))],'c','FaceAlpha',0.2,'EdgeColor','none');

plot([1 n_vars],[mean(p1) mean(p1)],'b-.')

a=get(h,'CurrentAxes');
set(a,'FontSize',14)
set(a,'YScale','log')
ylabel('p-values','FontSize',18)
xlabel('Variables in selected order','FontSize',18)

saveas(gcf,'Fig/example3');
saveas(gcf,'Fig/example3.eps','epsc');


%% Plot zoom

set(a,'Box','on')
ylabel('')
xlabel('')
axis([1 10 1e-3 1])
legend HIDE

saveas(gcf,'Fig/example3_100zoom');
saveas(gcf,'Fig/example3_100zoom.eps','epsc');


%% Compute table with statistics: using un-sorted p-values

load example3

name={'VASCA','VASCA + bootstrap','FDR'}';

p2bb = p22(:,:)<0.01; % VASCA
p2_1var = sum(sum(p2bb(:,1:3),2)>0)/rep;
p2_2var = sum(sum(p2bb(:,1:3),2)>1)/rep;
p2_3var = sum(sum(p2bb(:,1:3),2)>2)/rep;
p2_FPR = sum(sum(p2bb(:,4:end)))/(rep*(n_vars-3));
p2_FDR = sum(sum(p2bb(:,4:end)))/(sum(sum(p2bb)));
p2 = [p2_1var p2_2var p2_3var p2_FPR p2_FDR];

p2bb = p22b(:,:)<0.01; % VASCA + bootstrap
p2b_1var = sum(sum(p2bb(:,1:3),2)>0)/rep;
p2b_2var = sum(sum(p2bb(:,1:3),2)>1)/rep;
p2b_3var = sum(sum(p2bb(:,1:3),2)>2)/rep;
p2b_FPR = sum(sum(p2bb(:,4:end)))/(rep*(n_vars-3));
p2b_FDR = sum(sum(p2bb(:,4:end)))/(sum(sum(p2bb)));
p2b = [p2b_1var p2b_2var p2b_3var p2b_FPR p2b_FDR];

p3b = p33(:,:)<0.01; % FDR
p3_1var = sum(sum(p3b(:,1:3),2)>0)/rep;
p3_2var = sum(sum(p3b(:,1:3),2)>1)/rep;
p3_3var = sum(sum(p3b(:,1:3),2)>2)/rep;
p3_FPR = sum(sum(p3b(:,4:end)))/(rep*(n_vars-3));
p3_FDR = sum(sum(p3b(:,4:end)))/sum(sum(p3b));
p3 = [p3_1var p3_2var p3_3var p3_FPR p3_FDR];

X = [p2;p2b;p3];

T = table(name, X(:,1), X(:,2), X(:,3), X(:,4), X(:,5), 'VariableNames', {'Method','SigVar1','SigVars2','SigVars3','FPR','FDR'})
