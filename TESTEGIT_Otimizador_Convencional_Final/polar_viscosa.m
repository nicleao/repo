function [clv, cdv ,cdp] = polar_viscosa(perfil,de,alpha)

rmtt = ['Polar_',perfil];

fileID = fopen(rmtt,'r');
formatSpec = '%f %f %f %f %f %f %f %f %f %f';
data = textscan(fileID,formatSpec,'Headerlines',11);


% Encontrando o index do alpha viscoso
diff = data{:,1} - alpha;
alpha_visc = min(abs(diff));
ind = find(abs(diff) == abs(alpha_visc));
data = cell2mat(data);

clv = data(ind,2);
cdv = data(ind,3);
cdp = data(ind,4);

% end