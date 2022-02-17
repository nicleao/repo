function [spar_height,chapeado_length,capstrip_length,nerv_surface,entelagem_length] = leitura_perfil(pontos,fator)
% 
% fileID = fopen(perfil,'r');
n(1) = size(pontos,1);
k = NaN(1,3); % 3 pontos são importantes: [origem  thick+Y  thick-Y]
                % thick+Y indica o valor Y no 1 quadrante e thick-Y no 4
                % quadrante
                % o programa roda de modo a determinar esses 3 pontos
                % thick+Y é diferente do negativo de thick-Y caso o perfil
                % não seja simétrico

k = 1; i = 1;
dist_thick = 0.25;
for p = 1:n  % Separa os pontos em 2 vetores com pontos em y positivo e negativo
    if pontos(p,2) >= 0
        positivo(k,1) = pontos(p,1);
        positivo(k,2) = pontos(p,2);
        k = k+1;
    else
        negativo(i,1) = pontos(p,1);
        negativo(i,2) = pontos(p,2);
        i = i+1;
    end
end

for i = 1:length(positivo)
    positivo(i,1) = abs(positivo(i,1)-dist_thick); % Calcula a distância entre cada ponto e dist_thick
end
for i = 1:length(negativo)
    negativo(i,1) = abs(negativo(i,1)-dist_thick); % Calcula a distância entre cada ponto e dist_thick
end

pos_min = min(positivo(:,1)); % Pega o ponto que tem a mínima distância
neg_min = min(negativo(:,1)); % Pega o ponto que tem a mínima distância

pos_positivo = find(positivo(:,1)==pos_min); % Localiza a posição em que está esse ponto
pos_negativo = find(negativo(:,1)==neg_min); % Localiza a posição em que está esse ponto

if size(pos_positivo,2) ~= 1 || size(pos_negativo,2) ~= 1 % Se houver 2 pontos com a mesma distância pega o primeiro ponto
    pos_positivo = pos_positivo(1,1);
    pos_negativo = pos_negativo(1,1);
end

y_pos = positivo(pos_positivo,2); % pega a coordenada y desse ponto
y_neg = negativo(pos_negativo,2);

p1 = min(pontos(:,1)); % Encontra o ponto do BA

% Pega a localização desses pontos
k_1 = find(pontos(:,1)==p1);
k_2 = find(pontos(:,2) == y_pos);
k_3 = find(pontos(:,2) == y_neg);

if size(k_1,2) ~= 1 || size(k_2,2) ~= 1 || size(k_3,2) ~= 1 % se tiver 2 pontos com a mesma altura (mesmo y) pega o primeiro
    k_1 = k_1(1);
    k_2 = k_2(1);
    k_3 = k_3(1);    
end

k(1) = k_1;
k(2) = k_2;
k(3) = k_3;

 
% distance(funcao,inicio,fim,n) -> n determina se deseja-se calcular a
% distancia entre dois pontos no eixo y ou se deseja-se obter o comprimento
% de uma curva
spar_height = fator*distance(pontos,k(2),k(3),1)-2; % -2mm chapeado
chapeado_length = fator*distance(pontos,k(2),k(3),2);
capstrip_length = fator*(distance(pontos,2,k(2),2) + distance(pontos,k(3),n,2));
nerv_surface = area(fator*pontos,n);
entelagem_length = chapeado_length + capstrip_length;

end


