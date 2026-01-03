clear all
close all
clc

elements = 9; n_wind = 17 - elements + 1;
rho_water = 1000;

[r,BMR] = DataImport_estrazioneParam(1,17); 

volumes = (4*pi.*r.^3)./3;
mass = rho_water.*volumes;       
BMR = abs(BMR);   

exponents = [ones(n_wind,2) zeros(n_wind,2)];

p_HZ = zeros(n_wind,elements); p_An3 = zeros(n_wind,3);
for i = 1:n_wind
    
    rmin = i; rmax = i + elements - 1;
    
    [p_HZ(i,:),p_An3(i,:)] = statistics(mass(:,rmin:rmax),BMR(:,rmin:rmax),exponents,elements);
    
end

results = [p_HZ p_An3(:,3)];