function [coord_mass, coord_bmr, n_bins, pdfs, c_mass, c_bmr, m_collapsed, b_collapsed, Fscale] = hist_collapseDist(m, b, a1, b1, d1, d2, el)

    n_bins = 30;   % number of bins, arbitrarly chosen
     
    coord_mass = zeros(n_bins,el);
    coord_bmr = zeros(n_bins,el);
    
    c_mass = zeros(n_bins,n_bins,el);
    c_bmr = zeros(n_bins,n_bins,el);
    
    pdfs = zeros(n_bins,n_bins,el);
    
    % histogram generation for each joint distribution
    for i = 1:el
        
        M = m(:,i); B = b(:,i);
%         Mh = M(find(m(:,i) > 0)); Bh = B(find(b(:,i) > 0)); 
        
        % joint probability density function (p(m,b))
        H = histogram2(M,B,[n_bins n_bins],'Normalization','pdf'); hold on
        
        % only for deterministic case
%        histogram2(m(:,i),b(:,i),[n_bins n_bins],'Normalization','pdf'); hold on
        
        masses_bin = H.XBinEdges;   % bin edges on mass axis (1-by-n_bins+1 array)
        bmr_bin = H.YBinEdges;      % bin edges on BMR axis (1-by-n_bins+1 array)
        pdfs(:,:,i) = H.Values';    % z-axis, probability for each (m,b) couple 
                                    % (n_bins-by-n_bins matrix)
        
        % midpoint determination for each bin (i.e. discretized mass and BMR coordinates)
        for k = 1:n_bins       
            
            coord_mass(k,i) = (masses_bin(k) + (masses_bin(k+1) - masses_bin(k))/2);
            coord_bmr(k,i) = (bmr_bin(k) + (bmr_bin(k+1) - bmr_bin(k))/2);
            
        end

    end
    
    % collapse the x-axis (i.e. the masses) and the y-axis (i.e. the BMRs)
    m_mean = mean(m,1);    % average masses (1-by-el array)
                           % ogni elemento del vettore contiene la media
                           % della colonna
    
    m_mean_el = abs(m_mean).^d1;          % (1-by-el array)
                                          % ogni valore rappresenta la media elevata alla d1 
    m_collapsed = coord_mass./m_mean_el;  % collapsed discretized masses (n_bins-by-el matrix)
%     m_collapsed = log10(m_collapsed);   % collapsed discretized masses in log-scale (n_bins-by-el matrix)

    b_mean_el = abs(m_mean).^d2;          % (1-by-el array)
    b_collapsed = coord_bmr./b_mean_el;   % collapsed discretized BMRs (n_bins-by-el matrix)
%     b_collapsed = log10(b_collapsed);   % collapsed discretized BMRs in log-scale (n_bins-by-el matrix)
    
    % elevation by normalization exponents (n_bins-by-el matrices)
    m_hist = abs(coord_mass).^a1;
    b_hist = abs(coord_bmr).^b1;  
    
    % rearrangement in a 3D matrix (creation of (m,b) couples) 
    for i = 1:n_bins
        for k = 1:el
            
            c_mass(i,:,k) = m_hist(i,k);
            c_bmr (:,i,k) = b_hist(i,k);

        end           
    end

    Fscale = abs(pdfs.*c_mass.*c_bmr);  % scaling function according to Cicero eq.
%     Fscale(:,1:7,:) = zeros(n_bins,7,el); Fscale(:,n_bins-6:n_bins,:) = zeros(n_bins,7,el);
%     Fscale(1:7,:,:) = zeros(7,n_bins,el); Fscale(n_bins-6:n_bins,:,:) = zeros(7,n_bins,el);
    Fscale(~isfinite(Fscale)) = 0;      % delete numerical errors
        
end
    
    