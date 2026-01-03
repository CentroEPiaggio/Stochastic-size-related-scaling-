function [all_radii, all_fluxes] = DataImport_NumSamples(r_min, r_max)

    s = open('DataStem001.mat');  % caricamento dati FEM
%     s_bis = open('Data001-2.mat'); s_tris = open('Data001-3.mat'); 
%     s_4 = open('Data001-4.mat'); s_5 = open('Data001-5.mat'); 
    
    % colonna1: raggi, colonna2: Km, colonna3: OCRs, colonna4: densità cellulari,
    % colonna5: BMRs, colonna6: moli di O2, colonna7: volume in cui cO2 > 0.04 mol/m^3, 
    % colonna8: volume totale
    
    % distribuzioni di raggio [m]
    R(:,1)=s.F1(:,1); R(:,2)=s.F2(:,1); R(:,3)=s.F3(:,1); R(:,4)=s.F4(:,1); R(:,5)=s.F5(:,1);
    R(:,6)=s.F6(:,1); R(:,7)=s.F7(:,1);
    R(:,8)=s.F8(:,1); R(:,9)=s.F9(:,1); R(:,10)=s.F10(:,1);
    R(:,11)=s.F11(:,1); R(:,12)=s.F12(:,1); R(:,13)=s.F13(:,1); 
    R(:,14)=s.F14(:,1); R(:,15)=s.F15(:,1); R(:,16)=s.F16(:,1); R(:,17)=s.F17(:,1); 
    
%     Rbis(:,1)=s_bis.F1(:,1); Rbis(:,2)=s_bis.F2(:,1); Rbis(:,3)=s_bis.F3(:,1); Rbis(:,4)=s_bis.F4(:,1);
%     Rbis(:,5)=s_bis.F5(:,1); Rbis(:,6)=s_bis.F6(:,1); Rbis(:,7)=s_bis.F7(:,1); Rbis(:,8)=s_bis.F8(:,1);
%     Rbis(:,9)=s_bis.F9(:,1); Rbis(:,10)=s_bis.F10(:,1); Rbis(:,11)=s_bis.F11(:,1); Rbis(:,12)=s_bis.F12(:,1);
%     Rbis(:,13)=s_bis.F13(:,1); Rbis(:,14)=s_bis.F14(:,1); Rbis(:,15)=s_bis.F15(:,1); Rbis(:,16)=s_bis.F16(:,1);
%     Rbis(:,17)=s_bis.F17(:,1); 
% 
%     Rtris(:,1)=s_tris.F1(:,1); Rtris(:,2)=s_tris.F2(:,1); Rtris(:,3)=s_tris.F3(:,1); Rtris(:,4)=s_tris.F4(:,1);
%     Rtris(:,5)=s_tris.F5(:,1); Rtris(:,6)=s_tris.F6(:,1); Rtris(:,7)=s_tris.F7(:,1); Rtris(:,8)=s_tris.F8(:,1);
%     Rtris(:,9)=s_tris.F9(:,1); Rtris(:,10)=s_tris.F10(:,1); Rtris(:,11)=s_tris.F11(:,1); Rtris(:,12)=s_tris.F12(:,1);
%     Rtris(:,13)=s_tris.F13(:,1); Rtris(:,14)=s_tris.F14(:,1); Rtris(:,15)=s_tris.F15(:,1); Rtris(:,16)=s_tris.F16(:,1);
%     Rtris(:,17)=s_tris.F17(:,1); 
% 
%     R4(:,1)=s_4.F1(:,1); R4(:,2)=s_4.F2(:,1); R4(:,3)=s_4.F3(:,1); R4(:,4)=s_4.F4(:,1);
%     R4(:,5)=s_4.F5(:,1); R4(:,6)=s_4.F6(:,1); R4(:,7)=s_4.F7(:,1); R4(:,8)=s_4.F8(:,1);
%     R4(:,9)=s_4.F9(:,1); R4(:,10)=s_4.F10(:,1); R4(:,11)=s_4.F11(:,1); R4(:,12)=s_4.F12(:,1);
%     R4(:,13)=s_4.F13(:,1); R4(:,14)=s_4.F14(:,1); R4(:,15)=s_4.F15(:,1); R4(:,16)=s_4.F16(:,1);
%     R4(:,17)=s_4.F17(:,1); 
% 
%     R5(:,1)=s_5.F1(:,1); R5(:,2)=s_5.F2(:,1); R5(:,3)=s_5.F3(:,1); R5(:,4)=s_5.F4(:,1);
%     R5(:,5)=s_5.F5(:,1); R5(:,6)=s_5.F6(:,1); R5(:,7)=s_5.F7(:,1); R5(:,8)=s_5.F8(:,1);
%     R5(:,9)=s_5.F9(:,1); R5(:,10)=s_5.F10(:,1); R5(:,11)=s_5.F11(:,1); R5(:,12)=s_5.F12(:,1);
%     R5(:,13)=s_5.F13(:,1); R5(:,14)=s_5.F14(:,1); R5(:,15)=s_5.F15(:,1); R5(:,16)=s_5.F16(:,1);
%     R5(:,17)=s_5.F17(:,1); 
    
    % distribuzioni di BMR [mol/s]
    F(:,1)=s.F1(:,3); F(:,2)=s.F2(:,3); F(:,3)=s.F3(:,3); F(:,4)=s.F4(:,3); F(:,5)=s.F5(:,3);
    F(:,6)=s.F6(:,3); F(:,7)=s.F7(:,3); 
    F(:,8)=s.F8(:,3); F(:,9)=s.F9(:,3); F(:,10)=s.F10(:,3);
    F(:,11)=s.F11(:,3); F(:,12)=s.F12(:,3); F(:,13)=s.F13(:,3);
    F(:,14)=s.F14(:,3); F(:,15)=s.F15(:,3); F(:,16)=s.F16(:,3); F(:,17)=s.F17(:,3);  
  
%     Fbis(:,1)=s_bis.F1(:,5); Fbis(:,2)=s_bis.F2(:,5); Fbis(:,3)=s_bis.F3(:,5); 
%     Fbis(:,4)=s_bis.F4(:,5); Fbis(:,5)=s_bis.F5(:,5); Fbis(:,6)=s_bis.F6(:,5); 
%     Fbis(:,7)=s_bis.F7(:,5); Fbis(:,8)=s_bis.F8(:,5); Fbis(:,9)=s_bis.F9(:,5); 
%     Fbis(:,10)=s_bis.F10(:,5); Fbis(:,11)=s_bis.F11(:,5); Fbis(:,12)=s_bis.F12(:,5); 
%     Fbis(:,13)=s_bis.F13(:,5); Fbis(:,14)=s_bis.F14(:,5); Fbis(:,15)=s_bis.F15(:,5); 
%     Fbis(:,16)=s_bis.F16(:,5); Fbis(:,17)=s_bis.F17(:,5); 
% 
%     Ftris(:,1)=s_tris.F1(:,5); Ftris(:,2)=s_tris.F2(:,5); Ftris(:,3)=s_tris.F3(:,5); 
%     Ftris(:,4)=s_tris.F4(:,5); Ftris(:,5)=s_tris.F5(:,5); Ftris(:,6)=s_tris.F6(:,5); 
%     Ftris(:,7)=s_tris.F7(:,5); Ftris(:,8)=s_tris.F8(:,5); Ftris(:,9)=s_tris.F9(:,5); 
%     Ftris(:,10)=s_tris.F10(:,5); Ftris(:,11)=s_tris.F11(:,5); Ftris(:,12)=s_tris.F12(:,5); 
%     Ftris(:,13)=s_tris.F13(:,5); Ftris(:,14)=s_tris.F14(:,5); Ftris(:,15)=s_tris.F15(:,5); 
%     Ftris(:,16)=s_tris.F16(:,5); Ftris(:,17)=s_tris.F17(:,5); 
%    
%     F4(:,1)=s_4.F1(:,5); F4(:,2)=s_4.F2(:,5); F4(:,3)=s_4.F3(:,5); F4(:,4)=s_4.F4(:,5); 
%     F4(:,5)=s_4.F5(:,5); F4(:,6)=s_4.F6(:,5);  F4(:,7)=s_4.F7(:,5); F4(:,8)=s_4.F8(:,5); 
%     F4(:,9)=s_4.F9(:,5); F4(:,10)=s_4.F10(:,5); F4(:,11)=s_4.F11(:,5); F4(:,12)=s_4.F12(:,5); 
%     F4(:,13)=s_4.F13(:,5); F4(:,14)=s_4.F14(:,5); F4(:,15)=s_4.F15(:,5); F4(:,16)=s_4.F16(:,5); 
%     F4(:,17)=s_4.F17(:,5); 
% 
%     F5(:,1)=s_5.F1(:,5); F5(:,2)=s_5.F2(:,5); F5(:,3)=s_5.F3(:,5); F5(:,4)=s_5.F4(:,5); 
%     F5(:,5)=s_5.F5(:,5); F5(:,6)=s_5.F6(:,5); F5(:,7)=s_5.F7(:,5); F5(:,8)=s_5.F8(:,5); 
%     F5(:,9)=s_5.F9(:,5); F5(:,10)=s_5.F10(:,5); F5(:,11)=s_5.F11(:,5); F5(:,12)=s_5.F12(:,5); 
%     F5(:,13)=s_5.F13(:,5); F5(:,14)=s_5.F14(:,5); F5(:,15)=s_5.F15(:,5); F5(:,16)=s_5.F16(:,5); 
%     F5(:,17)=s_5.F17(:,5); % .*(4*pi*(R(:,17).^2)); 
%     
%     Rtot = [R; Rbis; Rtris; R4; R5]; Ftot = [F; Fbis; Ftris; F4; F5];

    NumSam = 1500;                     % numero di campioni
%     [Rs,index] = sort(R);
    [Rs,index] = sort(R(NumSam+1:2*NumSam,:));    % raggi ordinati in ordine crescente
    
    Fs = zeros(NumSam,17);
    for j = 1:17
    
        Fs(:,j) = F(index(:,j),j);       % BMR ordinati conseguentemente ai raggi
        
    end
    
%     all_radii = zeros(NumSam,17);
%     all_fluxes = zeros(NumSam,17);
%     for i = r_min:r_max                  % per ogni distribuzione da considerare...
%         
%         H = histogram(Rs(:,i),30,'Normalization','probability');
%         He = H.BinEdges; Hv = H.Values;
%         
%         % elimino i dati che sono associati ad una frequenza relativa minore del 1% 
%         ind = find(Hv < 0.01);     
%         
%         % individuazione soglie di eliminazione
%         flag = [];
%         for j = 1:length(ind)-1
%             
%             diff = ind(j+1) - ind(j);
%             
%             if diff > 1  
%                flag = [flag; j];
%             end
%             
%         end
%         
%         thrL = (He(ind(flag(1))) + He(ind(flag(1))+1))/2;       % soglia coda inferiore 
%         thrH = (He(ind(flag(end)+1)) + He(ind(flag(end)+1)+1))/2;   % soglia coda superiore
%         
%         % eliminazione delle code sulla base delle soglie
%         RindL = find(Rs(:,i) < thrL); 
%         RindH = find(Rs(:,i) > thrH);
%         
%         if isempty(RindL)
%            all_radii(1:RindH(1)-1,i) = Rs(1:RindH(1)-1,i);
%            all_fluxes(1:RindH(1)-1,i) = Fs(1:RindH(1)-1,i);
%         else
%            all_radii(RindL(end)+1:RindH(1)-1,i) = Rs(RindL(end)+1:RindH(1)-1,i);
%            all_fluxes(RindL(end)+1:RindH(1)-1,i) = Fs(RindL(end)+1:RindH(1)-1,i);
%         end
%         
%     end
% %     
%     % la funzione restituisce solo le colonne necessarie
%     all_radii = all_radii(:,r_min:r_max);
%     all_fluxes = all_fluxes(:,r_min:r_max);
    
    % importo le distribuzioni della finestra in esame
    all_radii = Rs(25:NumSam-25,r_min:r_max);
    all_fluxes = Fs(25:NumSam-25,r_min:r_max);
    
    % NOTA: l'eliminazione delle code è da impostare volta per volta in funzione 
    % di NumSam
    
end
