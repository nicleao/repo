function [cdv ,cdp] = ViscosidadePolar(perfil,Cl)

% O Cdv de uma se��o de asa � pego da polar pelo Cl dessa se��o

% fprintf(' %s Cl=%.3f\n',perfil,Cl)
rmtt = ['Polar_',perfil];
% A polar do perfil tem q tomar alguns(2) cuidados:
% � importante que haja �NGULOS NEGATIVOS na polar!
% N�O PODE ter �ngulos depois do ESTOL!

fileID = fopen(rmtt,'r');
formatSpec = '%f %f %f %f %f %f %f %f %f %f';
data = textscan(fileID,formatSpec,'Headerlines',11);

% Encontrando o index do alpha viscoso
diff = data{:,2} - Cl;
Cl_visc = min(abs(diff));
ind = find(abs(diff) == abs(Cl_visc));
ind = min(ind);
data = cell2mat(data);

cdv = data(ind,3);
cdp = data(ind,4);

% end