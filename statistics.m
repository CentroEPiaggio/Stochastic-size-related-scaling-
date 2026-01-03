function [normality,collapse] = statistics(M,B,expo,el)

    % collapsed raw data
    Mcoll = log10(abs(M)./((mean(abs(M))).^expo(3)));
    Bcoll = log10(abs(B)./((mean(abs(M))).^expo(4)));
    
    % multivariate (log)normality test for each collapsed distribution
    normality = zeros(1,el);
    AnovaNdata = zeros(length(Mcoll),2,el);
    for i = 1:el              % sweep on the distributions
        
        HZdata = [Mcoll(:,i) Bcoll(:,i)];
        normality(i) = HZmvntest(HZdata,0.05);   % p-values
        
        AnovaNdata(:,:,i) = HZdata;
        
    end
    
    % if the collapsed distributions sufficiently approaches a normal form...
   
    % 3 ways Anova for a proxy of collapse validation (look at the 3rd factor) 
    collapse = anovanTable(AnovaNdata,'varnames', {'mass', 'mr', 'subsets'});

end

