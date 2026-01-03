%% New code

close all
clc
clear

%% Range of radii and cell density

% rmin = 13;   % index of the first considered radius
% rmax = 17;   % index of the last considered radius
% elements = rmax - rmin + 1;  % number of considered radii in a window   

elements = 3;  % number of considered radii in a window

rho_c = 2.54e14;    % [cells/m^3], cell density
rho_water = 1000;   % [kg/m^3], physical density

%% Parameter ranges for the multiparametric fitting

alpha_min = 1; alpha_max = 1.; step_alpha = 0.05;  % normalization exponent for mass
beta_min = 1; beta_max = 1; step_beta = 0.05;     % normalization exponent for BMR

delta1_min = 0.9; delta1_max = 1.1; step_delta1 = 0.01;  % scaling exponent for mass
delta2_min = 0.5; delta2_max = 1.5; step_delta2 = 0.01;  % scaling exponent for BMR

%% Data importing and multiparametric fitting

results = zeros(9,15);
for rmin = 10:10                       % scorre le finestre
    
    rmax = rmin + elements - 1;
    
    % r and f contain radius and (negative) BMR distributions respectively
    [r,BMR] = DataImport_NumSamples(rmin,rmax); 
                                     
    volumes = (4*pi.*r.^3)./3;
    mass = rho_water.*volumes;       % mass contains mass distributions
    BMR = abs(BMR);                  % positive BMR distributions

    mass_av = mean(mass,1); BMR_av = mean(BMR,1);  % average masses and BMRs

    FITdist = [];

    % for all the possible combinations of Cicero eq. parameters
    for alpha_fit = alpha_min:step_alpha:alpha_max
        for beta_fit = beta_min:step_beta:beta_max
            for delta1_fit = delta1_min:step_delta1:delta1_max   
                for delta2_fit = delta2_min:step_delta2:delta2_max

                    % returns the three axes of the scaling function graph
                    [M, B, N, P, cM, cB, masses_hist, bmr_hist, prob_hist] = hist_collapseDist(mass, BMR,...
                    alpha_fit, beta_fit, delta1_fit, delta2_fit, elements);

                    % overall sum of Euclidean distances in log-scale
                    sum_dist = dist_calculator(N, masses_hist, bmr_hist, prob_hist, elements);

                    FITdist = [FITdist; alpha_fit beta_fit delta1_fit delta2_fit sum_dist];

                end
            end
        end
    end

    % minimization and parameter identification
    [FITdist_min, FITdist_index] = min(FITdist(:,5)); 
    bestCollapse = [FITdist(FITdist_index,1),FITdist(FITdist_index,2), ...
                   FITdist(FITdist_index,3),FITdist(FITdist_index,4)]

    % coordinates for best collapse
    bestMass = log10(M./(mean(abs(mass)).^bestCollapse(3)));
    bestBMR = log10(B./(mean(abs(mass)).^bestCollapse(4)));
    bestFscale = log10((cM.^bestCollapse(1)).*(cB.^bestCollapse(2)).*P);

    % visualization of collapsed joint distribution
    figure;
    for i = 1:elements
        histogram2(log10(mass(:,i)./(mean(abs(mass(:,i))).^bestCollapse(3))),...
        log10(BMR(:,i)./(mean(abs(mass(:,i))).^bestCollapse(4))),[N N],'Normalization','pdf'); hold on
        scatter3(log10(mass_av(i)./(mean(abs(mass(:,i))).^bestCollapse(3))),...
        log10(BMR_av(i)./(mean(abs(mass(:,i))).^bestCollapse(4))),400,'filled'); hold on
    end
    xlabel('m/<m>^c'); ylabel('MR/<m>^d'); zlabel('p(m,MR)*m^a*MR^b');

%     figure;
%     for i = 1:elements
%         scatter(log10(mass(:,i)./(mass_av(i).^bestCollapse(3))),log10(BMR(:,i)./(mass_av(i).^bestCollapse(4)))); hold on
%     end
%     xlabel('m/<m>^c'); ylabel('MR/<m>^d'); zlabel('p(m,MR)*m^a*MR^b');
% 
%     figure;
%     for k = 1:elements
% 
%         X = zeros(size(bestMass,1)*size(bestBMR,1),1);
%         Y = zeros(size(bestMass,1)*size(bestBMR,1),1);
%         Z = zeros(size(bestMass,1)*size(bestBMR,1),1);
% 
%         for i = 1:size(bestMass,1)
%             for j = 1:size(bestBMR,1)
% 
%                 X((i-1)*size(bestBMR,1)+j) = bestMass(i,k);
%                 Y((i-1)*size(bestBMR,1)+j) = bestBMR(j,k);     
%                 Z((i-1)*size(bestBMR,1)+j) = bestFscale(i,j,k);
% 
%             end
%         end
% 
%         scatter3(X,Y,Z,'filled'); hold on
% 
%     end
%     xlabel('m/<m>^c'); ylabel('MR/<m>^d'); zlabel('p(m,MR)*m^a*MR^b');

    % acceptable range for exponents and visualization of cost function and its sections
    [range1, range10] = err_plotter(FITdist, bestCollapse, delta1_min, delta1_max, ...
                        step_delta1, delta2_min, delta2_max, step_delta2)
    
    % statistical validation
%     [norm_p,anova_p] = statistics(mass,BMR,bestCollapse,elements);
                    
    % writing results in the output matrix       
    results(rmin,1:2) = bestCollapse(3:4);
    results(rmin,3) = FITdist_min;
    results(rmin,4:7) = [range1(1,:) range10(1,:)];
    results(rmin,8:11) = [range1(2,:) range10(2,:)];
%     results(rmin,12:end) = [norm_p anova_p(3)];
    
end

% x-y projections of joint distributions and deterministic curve
% figure;
% for i = 1:elements
%     scatter(log10(mass(:,i)),log10(BMR(:,i)),'filled'); hold on
% end
% plot(log10(mass_av(1:5)),log10(BMR_av(1:5)),'--k','LineWidth',2,'Marker','o','MarkerFaceColor','k'); hold on
% plot(log10(mass_av(7:11)),log10(BMR_av(7:11)),'--k','LineWidth',2,'Marker','o','MarkerFaceColor','k'); hold on
% plot(log10(mass_av(1:17)),log10(BMR_av(13:17)),'k','LineWidth',2,'Marker','o','MarkerFaceColor','k');
% xlabel('log10(m)'); ylabel('log10(MR)');
% 
figure;
for i = 1:elements
    histogram2(log10(mass(:,i)),log10(BMR(:,i)),[N N],'Normalization','pdf'); hold on
    scatter3(log10(mass_av(i)),log10(BMR_av(i)),400,'filled'); hold on
end
xlabel('log10(m)'); ylabel('log10(MR)'); zlabel('p(m,MR)');