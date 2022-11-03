%% First simulation example, subexample B, in Variable-Selection ANOVA Simultaneous Component Analysis. Bioinformatics. 2022 
% Camacho J, Vitale R, Morales-Jimenez D. and Gómez-Llorente C. 
%
% We simulate two factors with two levels, with 40 subjects for which 400 
% variables or responses (e.g., -omics features) are collected. The  
% example illustrates the case where the data matrix and the levels coding 
% for the factor are unrelated. 
%
% coded by: Jose Camacho (josecamacho@ugr.es)
% last modification: 19/Oct/2022
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

n_obs = 48;       % number of individuals 
n_vars = 400;     % number of responses or variables
rep = 1000;        % number of repetitions of the simulation
simMV_engine = true; % simulation engine, true for simuleMV, false for randn 


close all
p1 = zeros(rep,2);
p2 = zeros(rep,n_vars,2);
p3 = zeros(rep,n_vars,2);
p4 = zeros(rep,n_vars,2);
parfor i= 1:rep

    reps = 4;
    levels = {[1,2,3,4],[1,2,3]};
    
    class = create_design(levels,reps);
    
    if simMV_engine
        X = simuleMV(n_obs,n_vars,7);       % Data is independent to the class, and obtained with simuleMV with a medium correlation level 
    else
        X = randn(n_obs,n_vars);
    end
    
    s = rng(i);
    [~,parglmo] = parglm(X,class); % general linear model (GLM) factorization and (ASCA type) multivariate significance testing
    rng(s);
    [~,parglmoVS] = parglmVS(X,class); % GLM factorization and (VASCA-type) incremental multivariate significance testing
    rng(s);
    [~,parglmoMC] = parglmMC(X,class); % GLM factorization and Benjamini-Hochberg (BH) univariate significant testing
    rng(s);
    parglmoG = parglm_genes(X,class); % GLM factorization and ASCA-genes
    
    p1(i,:) = parglmo.p;
    p2(i,:,:) = [parglmoVS.p(parglmoVS.ord_factors(1,:)',1) parglmoVS.p(parglmoVS.ord_factors(2,:)',2)];
    p3(i,:,:) = [parglmoMC.p(parglmoMC.ord_factors(1,:)',1) parglmoMC.p(parglmoMC.ord_factors(2,:)',2)];
    p4(i,:,:) = [parglmoG.p(parglmoG.ord_factors(1,:)',1) parglmoG.p(parglmoG.ord_factors(2,:)',2)];
end

save example1b

%% Plot Figures 

load example1b

minT = 1e-3;
maxT = 1;
h = figure; hold on

p1(find(p1(:)<minT)) = minT;
p1(find(p1(:)>maxT)) = maxT;
p2(find(p2(:)<minT)) = minT;
p2(find(p2(:)>maxT)) = maxT;
p3(find(p3(:)<minT)) = minT;
p3(find(p3(:)>maxT)) = maxT;
p4(find(p4(:)<minT)) = minT;
p4(find(p4(:)>maxT)) = maxT;

mp3 = squeeze(mean(p3));
mp2 = squeeze(mean(p2));
mp4 = squeeze(mean(p4));

for j=1:2
    
    h = figure; hold on
    
    plot(mp3(:,j),'k--')
    plot([1 n_vars],[mean(p1(:,j)) mean(p1(:,j))],'b-.')
    plot(mp2(:,j),'g-o')
    plot(mp4(:,j),'c-o')
    
    plot([0 n_vars],[0.05 0.05],'r:')
    plot([0 n_vars],[0.01 0.01],'r--')
    legend('FDR','ASCA','VASCA','ASCA-genes','\alpha=0.05','\alpha=0.01','Location','southeast')
    
    xr = [];
    yr = [];
    for i=1:size(p3,2)
        xr = [xr;i*ones(1,2)];
        yr = [yr;[mean(p3(:,i,j))-std(p3(:,i,j)) mean(p3(:,i,j))+std(p3(:,i,j))]];
    end
    yr(find(yr<minT)) = minT;
    yr(find(yr>maxT)) = maxT;
    fill([xr(:,1);flipud(xr(:,2))],[yr(:,1);flipud(yr(:,2))],'k','FaceAlpha',0.2,'EdgeColor','none');
    xr = [];
    yr = [];
    for i=1:size(p2,2)
        xr = [xr;i*ones(1,2)];
        yr = [yr;[mean(p2(:,i,j))-std(p2(:,i,j)) mean(p2(:,i,j))+std(p2(:,i,j))]];
    end
    yr(find(yr<minT)) = minT;
    yr(find(yr>maxT)) = maxT;
    fill([xr(:,1);flipud(xr(:,2))],[yr(:,1);flipud(yr(:,2))],'g','FaceAlpha',0.2,'EdgeColor','none');
    
    xr = [];
    yr = [];
    for i=1:size(p2,2)
        xr = [xr;i*ones(1,2)];
        yr = [yr;[mean(p4(:,i,j))-std(p4(:,i,j)) mean(p4(:,i,j))+std(p4(:,i,j))]];
    end
    yr(find(yr<minT)) = minT;
    yr(find(yr>maxT)) = maxT;
    fill([xr(:,1);flipud(xr(:,2))],[yr(:,1);flipud(yr(:,2))],'c','FaceAlpha',0.2,'EdgeColor','none');
    
    a=get(h,'CurrentAxes');
    set(a,'FontSize',14)
    set(a,'YScale','log')
    ylabel('p-values','FontSize',18)
    xlabel('Variables in selected order','FontSize',18)
    
    saveas(gcf,sprintf('Fig/example1b_100_%d',j));
    saveas(gcf,sprintf('Fig/example1b_100_%d.eps',j),'epsc');
    
    
%% Plot zoom

    set(a,'Box','on')
    ylabel('')
    xlabel('')
    axis([1 10 1e-3 1])
    legend HIDE
    
    saveas(gcf,sprintf('Fig/example1b_100zoom_%d',j));
    saveas(gcf,sprintf('Fig/example1b_100zoom_%d.eps',j),'epsc');
    
% QQ-plot

figure, qqplot(p1(:,j),makedist('Uniform')), title('QQ-plot ASCA')
saveas(gcf,'./Fig/example1b_QQ_ASCA_%d');
saveas(gcf,'./Fig/example1b_QQ_ASCA_%d.eps','epsc');

figure, qqplot(p2(:,1,j),makedist('Uniform')), title('QQ-plot VASCA')
saveas(gcf,sprintf('./Fig/example1b_QQ_VASCA_%d',j));
saveas(gcf,sprintf('./Fig/example1b_QQ_VASCA_%d.eps',j),'epsc');

figure, qqplot(p3(:,1,j),makedist('Uniform')), title('QQ-plot FDR')
saveas(gcf,sprintf('./Fig/example1b_QQ_FDR_%d',j));
saveas(gcf,sprintf('./Fig/example1b_QQ_FDR_%d.eps',j),'epsc');
end


%% Compute table with statistics

load example1b

name={'FDR_1','ASCA_1','VASCA_1','ASCA-genes_1','FDR_2','ASCA_2','VASCA_2','ASCA-genes_2'}';

X = [];
for j=1:2
    p3b = p3(:,:,j)<0.01; % FDR
    p3_FPR = sum(sum(p3b))/(rep*(n_vars));
    p1b = p1(:,j)<0.01; % ASCA
    p1_FPR = sum(p1b)/rep;
    p2b = p2(:,:,j)<0.01; % VASCA
    p2_FPR = sum(sum(p2b))/(rep*(n_vars));
    p4b = p4(:,:,j)<0.01; % ASCA-genes
    p4_FPR = sum(sum(p4b))/(rep*(n_vars));
    
    X = [X;p3_FPR;p1_FPR;p2_FPR;p4_FPR];
end

T = table(name, X, 'VariableNames', {'Method','FPR'})

