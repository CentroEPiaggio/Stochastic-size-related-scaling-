function [rangeA, rangeB] = err_plotter(stru, expo, c_min, c_max, c_step, d_min, d_max, d_step)
    
    % vettori contenenti i range degli esponenti
    c = c_min:c_step:c_max; d = d_min:d_step:d_max;
    
    % riarrangiamento dei valori della funzione costo in una matrice
    SE = zeros(length(c),length(d));
    for i = 1:length(c)
        for j = 1:length(d)
            
            SE(i,j) = stru((i-1)*length(d)+j,5);
            
        end
    end
    
    % intervallo di ammissibilità del collasso
    DerrA = 0.01; DerrB = 0.1;   % ampiezze dell'intervallo (1% e 10%)
    SEmin = min(min(SE));        % minimo della funzione costo
    
    % indici dell'intervallo   
    [rowIndA,colIndA] = find((SE - SEmin)/SEmin <= DerrA);  
    [rowIndB,colIndB] = find((SE - SEmin)/SEmin <= DerrB);  
    
    % esponenti corrispondenti agli estremi
    rangeA = [c(rowIndA(1)) c(rowIndA(end)); d(colIndA(1)) d(colIndA(end))];      
    rangeB = [c(rowIndB(1)) c(rowIndB(end)); d(colIndB(1)) d(colIndB(end))];
    
    % visualizzazione della funzione costo
%     figure; surf(d,c,SE); colormap('jet'); ylim([0.9 1.1]);       % superficie in 3D
%     xlabel('d'); ylabel('c'); zlabel('Cost function');

%     SEnorm = SE./max(max(SE));
%     figure; imshow(SEnorm,'Colormap',jet,'DisplayRange',[min(min(SEnorm)) max(max(SEnorm))]);
%     xlabel('delta2'); ylabel('delta1');                           % immagine
    
    % sezioni della funzione costo in corrispondenza degli esponenti migliori 
%     cost_c = SE(:,d == expo(4));
%     cost_d = SE(c == expo(3),:);
%     
    % visualzzazione delle sezioni
%     figure; plot(c,cost_c);
%     xlabel('c'); ylabel('Cost function');
%     
%     figure; plot(d,cost_d);
%     xlabel('d'); ylabel('Cost function');
    
end

