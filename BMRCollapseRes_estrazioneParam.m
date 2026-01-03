clear all
close all
clc

%% Caricamento dati 

% colonna1: raggi, colonna2: Km, colonna3: OCRs, colonna4: densità cellulari,
% colonna5: BMRs, colonna6: moli di O2, colonna7: volume in cui cO2 > 0.04 mol/m^3, 
% colonna8: volume totale

s = open('Data001-1.mat');  % caricamento dei dati

% distribuzioni di raggio [m]
R(:,1)=s.F1(:,1); R(:,2)=s.F2(:,1); R(:,3)=s.F3(:,1); R(:,4)=s.F4(:,1);
R(:,5)=s.F5(:,1); R(:,6)=s.F6(:,1); R(:,7)=s.F7(:,1); R(:,8)=s.F8(:,1);
R(:,9)=s.F9(:,1); R(:,10)=s.F10(:,1); R(:,11)=s.F11(:,1); R(:,12)=s.F12(:,1);
R(:,13)=s.F13(:,1); R(:,14)=s.F14(:,1); R(:,15)=s.F15(:,1); R(:,16)=s.F16(:,1);
R(:,17)=s.F17(:,1); 

W = 1000;                % [kg/m^3], densità fisica dell'organoide
M = W*(4/3)*pi*(R.^3);   % [kg], distribuzioni di massa

% distribuzioni di BMR [mol/s]
BMR(:,1)=s.F1(:,5); BMR(:,2)=s.F2(:,5); BMR(:,3)=s.F3(:,5); BMR(:,4)=s.F4(:,5);
BMR(:,5)=s.F5(:,5); BMR(:,6)=s.F6(:,5); BMR(:,7)=s.F7(:,5); BMR(:,8)=s.F8(:,5); 
BMR(:,9)=s.F9(:,5); BMR(:,10)=s.F10(:,5); BMR(:,11)=s.F11(:,5); BMR(:,12)=s.F12(:,5); 
BMR(:,13)=s.F13(:,5); BMR(:,14)=s.F14(:,5); BMR(:,15)=s.F15(:,5); BMR(:,16)=s.F16(:,5); 
BMR(:,17)=s.F17(:,5);

% distribuzioni di volume [m^3]
V(:,1)=s.F1(:,8); V(:,2)=s.F2(:,8); V(:,3)=s.F3(:,8); V(:,4)=s.F4(:,8);
V(:,5)=s.F5(:,8); V(:,6)=s.F6(:,8); V(:,7)=s.F7(:,8); V(:,8)=s.F8(:,8);
V(:,9)=s.F9(:,8); V(:,10)=s.F10(:,8); V(:,11)=s.F11(:,8); V(:,12)=s.F12(:,8);
V(:,13)=s.F13(:,8); V(:,14)=s.F14(:,8); V(:,15)=s.F15(:,8); V(:,16)=s.F16(:,8);
V(:,17)=s.F17(:,8); 

% distribuzioni di volume vitale [m^3]
Vlive(:,1)=s.F1(:,7); Vlive(:,2)=s.F2(:,7); Vlive(:,3)=s.F3(:,7); Vlive(:,4)=s.F4(:,7);
Vlive(:,5)=s.F5(:,7); Vlive(:,6)=s.F6(:,7); Vlive(:,7)=s.F7(:,7); Vlive(:,8)=s.F8(:,7);
Vlive(:,9)=s.F9(:,7); Vlive(:,10)=s.F10(:,7); Vlive(:,11)=s.F11(:,7); Vlive(:,12)=s.F12(:,7);
Vlive(:,13)=s.F13(:,7); Vlive(:,14)=s.F14(:,7); Vlive(:,15)=s.F15(:,7); Vlive(:,16)=s.F16(:,7);
Vlive(:,17)=s.F17(:,7); 

% considero una distribuzione di 300 valori
R = R(1:300,:); BMR = BMR(1:300,:); 
V = V(1:300,:); Vlive(1:300,:);

%% Identificazione degli esponenti di scaling

% Rmin = 9; Rmax = 11;     % distribuzioni di raggio da considerare

elements = 3;              % ampiezza della finestra

results = zeros(15,10);
for Rmin = 1:3             % scorre le finestre
    
    Rmax = Rmin + elements - 1;
    
    % sweep sugli esponenti
    bmin = 1; bstep = 1; bmax = 1.1; dmin = 0.5; dstep = 0.001; dmax = 1.5;
    b = bmin:bstep:bmax; d = dmin:dstep:dmax;         

    Nbins = 10;
    BMRedge = zeros(Nbins+1,17);
    coord_BMR = zeros(Nbins,17);
    P = zeros(Nbins,17);
    SE = zeros(length(b),length(d));

    for bInd = 1:length(b)        % scorre le possibili combinazioni degli esponenti
        for dInd = 1:length(d)

            for i = Rmin:Rmax            % scorre le distribuzioni di raggio

                % istogramma dei BMRs 
                H = histogram(abs(BMR(:,i)),Nbins,'Normalization','pdf');
                BMRedge(:,i) = H.BinEdges;
                P(:,i) = H.Values;

                for k = 1:Nbins           % scorre i bins dell'istogramma

                   % punti medi dei bins
                   coord_BMR(k,i) = (BMRedge(k,i) + (BMRedge(k+1,i)-BMRedge(k,i))/2);

                end

            end

            % punti medi dei bins normalizzati
            BMRcoll = coord_BMR./((mean(abs(M))).^d(dInd)); 
            Fscale = (coord_BMR.^b(bInd)).*P;    % F di scaling
            Fscale(1,:) = zeros(1,17); Fscale(Nbins-2:Nbins,:) = zeros(3,17);

            % punti medi dei bins normalizzati, logaritmici e ordinati in un unico vettore
            BMRcomp = sort(log10([BMRcoll(:,Rmin); BMRcoll(:,Rmin+1); BMRcoll(:,Rmin+2); ...
                      BMRcoll(:,Rmax-1); BMRcoll(:,Rmax)]));

            Nover = 0;
            for i = Rmin:Rmax-1        % scorre le coppie di distribuzioni
                for j = i+1:Rmax

                    % interpolazione della distribuzione i-esima nei punti della j-esima in loglog
                    Finterp = interp1(log10(BMRcoll(:,i)),log10(Fscale(:,i)),log10(BMRcoll(:,j)),'linear');
    %                 figure; plot(log10(BMRcoll(:,i)),log10(Fscale(:,i)),'b*'); hold on
    %                 plot(log10(BMRcoll(:,j)),Finterp,'r*');

                    res = zeros(Nbins,1);
                    for k = 1:Nbins     % scorre i punti medi dei bins

                        % se il k-esimo punto è di overlapping...
                        if log10(BMRcoll(k,j)) >= min(log10(BMRcoll(:,i))) && ...
                           log10(BMRcoll(k,j)) <= max(log10(BMRcoll(:,i)))

                           % residui tra l'interpolazione e la distribuzione j-esima in loglog
                           res(k) = abs(Finterp(k) - log10(Fscale(k,j)));

                        end

                    end

                    res(~isfinite(res)) = 0;    % correzione di eventuali errori numerici 

                    Nover = Nover + Nbins - length(find(res == 0));  % somma del numero di punti di overlapping
                    SE(bInd,dInd) = SE(bInd,dInd) + sum(res);        % somma dei residui tra le distribuzioni

                end 
            end

            SE(bInd,dInd) = SE(bInd,dInd)/Nover;    % normalizzazione

        end
    end

    % minimizzazione della somma delle distanze euclidee tra le distribuzioni
    [SEmin,rowInd] = min(SE); [SEmin_bis,colInd] = min(SEmin);
    bestInd = [rowInd(colInd) colInd];                 % indici di SE minimo
    % bestCollapse = [b(bestInd(1)) d(bestInd(2))]       % valori degli esponenti corrispondenti
    bestCollapse = [b d(bestInd(1))]

    % intervallo di ammissibilità del collasso (1% e 10%)
    Derr1 = 0.01;                              % ampiezza dell'intervallo

    int1 = find((SE - SEmin)/SEmin <= Derr1);  % indici dell'intervallo   
    range1 = [d(int1(1)) d(int1(end))]         % esponenti corrispondenti agli estremi

    Derr10 = 0.1;                                % ampiezza dell'intervallo

    int10 = find((SE - SEmin)/SEmin <= Derr10);  % indici dell'intervallo   
    range10 = [d(int10(1)) d(int10(end))]        % esponenti corrispondenti agli estremi

    
    % PERCENTUALE DI VOLUME VITALE MEDIA E RELATIVA DISPERSIONE

    live = 100*(Vlive./V);   % distribuzioni di percentuale di volume vitale

    % per le singole distribuzioni
    perc_av = mean(live); perc_std = std(live); perc_med = median(live);

    % per l'intera finestra
    wind = live(:,Rmin:Rmax);
    perc_avw = mean(wind(:)); perc_stdw = std(wind(:)); perc_medw = median(wind(:));

    
    % STATISTICA
    
    % distribuzioni di BMR riscalate per il test di Anderson-Darling
    ADdata = log10(abs(BMR(:,Rmin:Rmax))./((mean(abs(M(:,Rmin:Rmax)))).^bestCollapse(2)));
    ADvec = ADdata(:); 
    ADind = [ones(300,1); 2*ones(300,1); 3*ones(300,1)]; % 4*ones(300,1); 5*ones(300,1);
    %         6*ones(300,1); 7*ones(300,1); 8*ones(300,1); 9*ones(300,1)];

    % test di Anderson-Darling (Trujillo - Ortiz)
    [p_value,adjp_value] = AnDarksamtest([ADvec ADind],0.05)  

    % se l'ipotesi nulla è accettata...

    % assi del collasso ottimo
    bestBMRcoll = coord_BMR./((mean(abs(M))).^bestCollapse(2));  % punti medi dei bins normalizzati
    bestFscale = (coord_BMR.^bestCollapse(1)).*P;        % F di scaling

    % assi del collasso ottimo in scala logaritmica
    xdata = log10(bestBMRcoll(:,Rmin:Rmax)); 
    ydata = log10(bestFscale(:,Rmin:Rmax)); 

    xAver = mean(xdata,2); yAver = mean(ydata,2);  % distribuzione media
    
    % scrittura dei risultati in una matrice     
    results(Rmin,1) = bestCollapse(2);
    results(Rmin,2) = SEmin;
    results(Rmin,3:6) = [range1 range10];
    results(Rmin,7:9) = [perc_avw perc_stdw perc_medw];
    results(Rmin,end) = p_value;
    
end

%% Visualizzazione 

% figure;
% for i = Rmin:Rmax
%     loglog(coord_BMR(:,i),P(:,i),'LineWidth',2); hold on               % non collassate
% end
% xlabel('MR'); ylabel('p(MR)');
% 
% figure;
% for i = Rmin:Rmax
%     loglog(bestBMRcoll(:,i),bestFscale(:,i),'LineWidth',2); hold on    % collassate
% end
% xlabel('MR/<m>^d'); ylabel('p(MR)*MR^b');


% funzione errore in 3D
% X = zeros(length(b)*length(d),1);
% Y = zeros(length(b)*length(d),1);
% Z = zeros(length(b)*length(d),1);
% for i = 1:length(b)
%     for j = 1:length(d)
%         
%         X((i-1)*length(d)+j) = b(i);
%         Y((i-1)*length(d)+j) = d(j);     % costruzione degli assi
%         Z((i-1)*length(d)+j) = SE(i,j);
%         
%     end
% end

% figure; surf(d,b,SE);
% xlabel('d'); ylabel('b'); zlabel('Sum of errors');

% sezioni 2D della funzione errore
% figure; plot(d,SE(rowInd(colInd),:)); 
% figure; plot(d,SE,'LineWidth',2); % ylim([10 80]);
% xlabel('d'); ylabel('f(d)');

% figure; plot(b,SE(:,colInd));                    
% xlabel('b'); ylabel('Sum of errors');

% sezione 2D della concavità
% figure; plot(d(1:end-2),diff2); 
% xlabel('d'); ylabel('Concavity');

% fitting della funzione di scaling
% figure; plot(linxAver,linyAver); hold on                          % lineare
% plot(linxAver,fit); xlabel('BMR/<M>^d'); ylabel('F(BMR/<M>^d)');

% figure; plot(xAver,yAver); hold on                                % log
% plot(xAver,fit); xlabel('BMR/<M>^d'); ylabel('F(BMR/<M>^d)');