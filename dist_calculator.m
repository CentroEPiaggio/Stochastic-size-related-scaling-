function overall = dist_calculator(Nbins, m_collapsed, b_collapsed, Fscale, el)
    
    % rearrangement in a 3D matrix (creation of (m,b) couples) 
    M = zeros(Nbins,Nbins,el); B = zeros(Nbins,Nbins,el);
    for i = 1:Nbins
        for k = 1:el
            
            M(i,:,k) = m_collapsed(i,k);
            B(:,i,k) = b_collapsed(i,k);

        end           
    end

    Nover = 0; overall = 0;
    for i = 1:el-1        % scorre le coppie di distribuzioni
        for j = i+1:el
            
            % distanza euclidea tra i punti corrispondenti
            dist = sqrt((log10(M(:,:,i)) - log10(M(:,:,j))).^2 + ...
                   (log10(B(:,:,i)) - log10(B(:,:,j))).^2 + ... 
                   (log10(Fscale(:,:,i)) - log10(Fscale(:,:,j))).^2);
     
            ok = isfinite(dist); NaNind = find(ok == 0);
            dist(NaNind) = 0;                              % scarto dei punti di non overlapping
                
            Nover = Nover + Nbins^2 - length(NaNind);      % somma del numero di punti di overlapping
            overall = overall + sum(sum(dist));            % somma delle distanze tra le distribuzioni

        end
    end

    overall = overall/Nover;    % normalizzazione
    
end