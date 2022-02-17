%% Pesinho
function [peso] = peso_convencional(posb,cr,cm,ct,perfilr,perfilm,perfilp,b,offset,X,bev,beh,creh,h_emp,x_emp,ctev,crev)
%% Convertendo entradas para mm
cr = cr*1000; 
cm = cm*1000;
ct = ct*1000;
beh = beh*1000;
creh = creh*1000;
crev = crev*1000;
ctev = ctev*1000;
bev = bev*1000;
posb = posb.*1000;
b = b*1000;
offset = offset.*1000;
chord = [cr cm ct];
h_emp = h_emp*1000;
x_emp = x_emp*1000;

%% MATERIAIS
  % densidade % g/mm³
    densidade_1 = 190e-6; % balsa de 1mm
    densidade_238 = 125e-6; 
    densidade_honeycomb = 209e-6; %325e-6;
    densidade_H45 = 45e-6;
    densidade_manta = 469.9e-6; % 469.9 densidade superficial % 795.7e-6; % 1570e-6
    densidade_superficial_entelagem = 30e-6; % ?????
    densidade_4x3 = 9.46e-3; %g/mm
    % densidade_kevlar = 1021e-6; % g/m²
    % densidade_3 = 446e-6;
    % densidade_caixeta = 434.1e-6;
    % densidade ripa = 3.7g para 1 metro
    % massa H100+laminado por mm = 0.05941 g/mm
    
  % espessura 
    espessura_honeycomb = 5.5;
    espessura_H45 = 6;
    
%% Primeira Seção 
l_fus = 120; % valor estimado para largura da fuselagem/compartimento de carga
d_nerv = 110; % distância entre nervuras

n_vaz = floor((posb(2)-l_fus/2)/d_nerv);    % Número de espaços vazios entre nervuras %O floor arredonda para o menor valor inteiro
nerv_1 = n_vaz + 1;                         % Número de nervuras até a primeira seção (conta as 2 extremidades)
pos_nerv1 = 0:d_nerv:d_nerv*nerv_1;         % Posições das nervuras
if mod((posb(2)-l_fus/2)/d_nerv,1) >= 0.5   % Se resto for maior que metade de d_nerv, o numero de nervuras aumenta em 1 (um)
    nerv_1 = nerv_1 + 1;
    pos_nerv1(length(pos_nerv1)) = []; %Apaga o valor da ultima nervura
    pos_nerv1(length(pos_nerv1-1)) = []; %Apaga o valor da penultima nervura
    if size(pos_nerv1,2) == 0
        ult = posb(1);
    else
        ult = pos_nerv1(length(pos_nerv1));
    end
    pen_ner = ult + (posb(2)-ult)/2;         % altera o valor da penúltima nervura por um valor intermediário
    pos_nerv1 = [pos_nerv1 pen_ner posb(2)]; % Posições das nervuras
    
else
    pos_nerv1(length(pos_nerv1)) = []; %Apaga o valor da ultima nervura. Na linha abaixo, substitui a pos da ultima nervura pelo posb(2)
    pos_nerv1 = [pos_nerv1 posb(2)]; % Posições das nervuras
end

%% Segunda Seção
n_vaz = floor((posb(3)-posb(2))/d_nerv);    % Número de espaços vazios entre nervuras
nerv_2 = n_vaz;                             % Número de nervuras até a segunda seção (Conta apenas a última extremidade)
pos_nerv2 = (posb(2)+d_nerv):d_nerv:posb(2)+d_nerv*nerv_2; % Posições das nervuras
if mod((posb(3)-posb(2))/d_nerv,1) >= 0.5   % Se resto for maior que metade de d_nerv
    nerv_2 = nerv_2 +1;
    pos_nerv2(length(pos_nerv2)) = []; 
    if size(pos_nerv2,2) == 0
        ult = posb(2);
    else
        ult = pos_nerv2(length(pos_nerv2));
    end
    pen_ner = ult + (posb(3)-ult)/2;         % altera o valor da penúltima nervura por um valor intermediário
    pos_nerv2 = [pos_nerv2 pen_ner posb(3)]; % Posições das nervuras
    
else
    pos_nerv2(length(pos_nerv2)) = [];
    pos_nerv2 = [pos_nerv2 posb(3)]; % Posições das nervuras
end

numnerv = nerv_1 + nerv_2; % Número total de nervuras
pos_nerv = [pos_nerv1 pos_nerv2]; % Posição de todas as nervuras

%% Coordenadas dos perfis

% Perfil da Raiz 
fileID = fopen(perfilr,'r');
formatSpec = '%f %f %[^\n\r]'; % o arq vem nesse formato
data = textscan(fileID,formatSpec,'Headerlines',1);
perfilr = [data{:,1} data{:,2}];

%Perfil do meio 
fileID = fopen(perfilm,'r');
formatSpec = '%f %f %[^\n\r]'; % o arq vem nesse formato
data = textscan(fileID,formatSpec,'Headerlines',1);
perfilm = [data{:,1} data{:,2}];


fileID = fopen(perfilp,'r');
formatSpec = '%f %f %[^\n\r]'; % o arq vem nesse formato
data = textscan(fileID,formatSpec,'Headerlines',1);
perfilp = [data{:,1} data{:,2}];

%% Cálculos das cordas de cada nervura

[y1,c1,x1] = Corda(chord,posb,X,b,offset);

corda1 = NaN(1,nerv_1);
for k = 1:length(pos_nerv1)
    corda1(k) = interp1(y1,c1,pos_nerv1(k));
end

corda2 = NaN(1,nerv_2);
for k = 1:length(pos_nerv2)
    corda2(k) = interp1(y1,c1,pos_nerv2(k));
end

%% Primeira Seção

% Just prealocando
capstrip_raiz = NaN(1,nerv_1); nerv_raiz = NaN(1,nerv_1); chapeado_raiz = NaN(1,nerv_1); entelagem_raiz = NaN(1,nerv_1); spar_raiz = NaN(1,nerv_1);

% Perfil da Raiz primeira seção
for i = 1:nerv_1
corda = corda1(i);
[sp_raiz,chap_raiz,cap_raiz,n_raiz,ent_raiz] = leitura_perfil(perfilr,corda);
capstrip_raiz(i) = cap_raiz;
nerv_raiz(i) = n_raiz;
chapeado_raiz(i) = chap_raiz;
entelagem_raiz(i) = ent_raiz;
spar_raiz(i) = sp_raiz;
end

% Just prealocando
capstrip_m_1s = NaN(1,nerv_1); nerv_m_1s = NaN(1,nerv_1); chapeado_m_1s = NaN(1,nerv_1); entelagem_m_1s = NaN(1,nerv_1); spar_m_1s = NaN(1,nerv_1);

% Perfil meio = '_m' para Primeira seção = '_1s'
for i = 1:nerv_1
    corda = corda1(i); 
[spar_m11,chapeado_m11,capstrip_m11,nerv_m11,entelagem_m11] = leitura_perfil(perfilm,corda);
capstrip_m_1s(i) = capstrip_m11;
nerv_m_1s(i) = nerv_m11;
chapeado_m_1s(i) = chapeado_m11;
entelagem_m_1s(i) = entelagem_m11;
spar_m_1s(i) = spar_m11;
end

%% Segunda Seção

% Just prealocando
capstrip_m_2s = NaN(1,nerv_2); nerv_m_2s = NaN(1,nerv_2); chapeado_m_2s = NaN(1,nerv_2); entelagem_m_2s = NaN(1,nerv_2); spar_m_2s = NaN(1,nerv_2);

% Perfil meio = '_m' para Segunda seção = '_2s'
for i = 1:nerv_2
    corda = corda2(i);
[spar_m12,chapeado_m12,capstrip_m12,nerv_m12,entelagem_m12] = leitura_perfil(perfilm,corda);
capstrip_m_2s(i) = capstrip_m12;
nerv_m_2s(i) = nerv_m12;
chapeado_m_2s(i) = chapeado_m12;
entelagem_m_2s(i) = entelagem_m12;
spar_m_2s(i) = spar_m12;
end

% Just prealocando
capstrip_ponta = NaN(1,nerv_2); nerv_ponta = NaN(1,nerv_2); chapeado_ponta = NaN(1,nerv_2); entelagem_ponta = NaN(1,nerv_2); spar_ponta = NaN(1,nerv_2);

% Perfil ponta terceira seção
for i = 1:nerv_2
    corda = corda2(i);
[spar_p,chapeado_p,capstrip_p,nerv_p,entelagem_p] = leitura_perfil(perfilp,corda);
capstrip_ponta(i) = capstrip_p;
nerv_ponta(i) = nerv_p;
chapeado_ponta(i) = chapeado_p;
entelagem_ponta(i) = entelagem_p;
spar_ponta(i) = spar_p;
end

%% Area
Areas = NaN(1,length(chord));
for i = 1:length(chord)
    if i == 1
        n(1) = size(perfilr,1);
        Areas(i) = area(perfilr*chord(i),n);
    elseif i == 2
        n(1) = size(perfilm,1);
        Areas(i) = area(perfilm*chord(i),n);
    elseif i == 3
        n(1) = size(perfilp,1);
        Areas(i) = area(perfilp*chord(i),n);
    end
end

Area_nerv = NaN(1,numnerv);
for i = 1:numnerv
    if pos_nerv(i) <= posb(2)
        Area_nerv(i) = Areas(1) + (pos_nerv(i)-posb(1))*(Areas(2)-Areas(1))/(posb(2)-posb(1)); % interpolação
    elseif pos_nerv(i) > posb(2) 
        Area_nerv(i) = Areas(2) + (pos_nerv(i)-posb(2))*(Areas(3)-Areas(2))/(posb(3)-posb(2)); % interpolação
    end
end

%% Peso da Asa 
areahoneycomb = Area_nerv(1);
areacapstripasa = ((capstrip_raiz(1) + capstrip_ponta(nerv_2)) + hypot(posb(2),abs(cm+offset(2)-cr)) + hypot(posb(3)-posb(2),abs(ct+offset(3)-cm)))*10; %*10mm de largura
areanervuras = sum(Area_nerv);
% Calcula a area como se fosse um trapezio: 
areaentelagem = (entelagem_raiz(1) + entelagem_m_1s(nerv_1))*posb(2)/2 + (entelagem_m_1s(nerv_1) + entelagem_ponta(nerv_2))*(posb(3)-posb(2))/2;
areachapeado = (chapeado_raiz(1) + chapeado_m_1s(nerv_1))*posb(2)/2 + (chapeado_m_1s(nerv_1) + chapeado_ponta(nerv_2))*(posb(3)-posb(2))/2;
areaspar = (spar_raiz(1) + spar_m_1s(nerv_1))*posb(2)/2 + (spar_m_1s(nerv_1) + spar_ponta(nerv_2))*(posb(3)-posb(2))/2;

%% EV

fileID = fopen('NACA0012.dat','r');
formatSpec = '%f %f %[^\n\r]'; % o arq vem nesse formato
data = textscan(fileID,formatSpec,'Headerlines',1);
perfilv = [data{:,1} data{:,2}];

% vr = EV raiz e vp = EV ponta vsf forato

[~,chapeado_vr,capstrip_vr,nerv_vr,entelagem_vr] = leitura_perfil(perfilv,crev);

[~,chapeado_vp,capstrip_vp,nerv_vp,entelagem_vp] = leitura_perfil(perfilv,ctev);

num_nerv_ev = floor(bev/100)+1;
if mod(bev,100) > 0.5*100
    num_nerv_ev = num_nerv_ev +1;
end

nerv_ev = NaN(1,num_nerv_ev);
pos_nerv = 0;
for i = 1:num_nerv_ev
    nerv_ev(i) = nerv_vr + pos_nerv*(nerv_vp-nerv_vr)/bev;
    pos_nerv = i*100;
end

entelagem_ev = (entelagem_vr + entelagem_vp)*bev/2;
capstrip_ev = capstrip_vr + capstrip_vp + bev*10; %10 mm de largura
chapeado_ev = (chapeado_vr + chapeado_vp)*bev/2;
area_nerv_ev = sum(nerv_ev);


%% EH

fileID = fopen('NACA23012-24pt.dat','r');
formatSpec = '%f %f %[^\n\r]'; % o arq vem nesse formato
data = textscan(fileID,formatSpec,'Headerlines',1);
perfileh = [data{:,1} data{:,2}];

[spar_eh,chapeado_eh,capstrip_eh,nerv_eh,entelagem_eh] = leitura_perfil(perfileh,creh);

num_nerv_eh = floor(beh/100)+1;
if mod(beh,100) > 0.5*100
    num_nerv_eh = num_nerv_eh +1;
end

area_entelagem_eh = entelagem_eh*beh;
area_capstrip_eh = capstrip_eh*2*10; %10 mm de largura e 2 capstrips
area_chapeado_eh = chapeado_eh*beh;
area_nerv_eh = num_nerv_eh*nerv_eh;
area_spar_eh = spar_eh*beh;

%% Cauda supondo que é treliçada

length_cauda = hypot((x_emp - cr),h_emp); % #descubra


%% Densidade linear do boom, baseado na regular de 2018
densidade_boom = 110/(0.3651*1000); %[g/mm]

%% MASSA

m_238 = 26.9; % massa 238 - Servo 238
m_939 = 16.4; % massa 939 - Servo 939
link_metal = 2.2;
on_off = 8.2;
bateria = 64.5;
receptor = 14.8;
extensao = 10; % 10g/m de comprimento
comp_extensao = 5; %Para um convencional



massa = NaN(1,22);

%% Massa da asa
massa(1) = 2*areaspar*espessura_H45*densidade_H45;  %o fator 2 multiplica a semi-longarina                        
% massa de divynicell na longarina 0.7 de alivio
massa(2) = 2*2*areaspar*densidade_manta;                                                         
% laminação longarina % fator 2 considerar os dois lados laminados
massa(3) = 2*areanervuras*2.38*densidade_238*0.7;                                                  
% 0.7 fator alívio de nervuras + nervuras raiz e servo / 2.38 espessura da balsa
massa(4) = 2*areaentelagem*densidade_superficial_entelagem;                                        
massa(5) = 2*areachapeado*densidade_1;                                                            
massa(6) = 2*areacapstripasa*densidade_1;                                                         
massa(7) = 50; % CA na asa - estimativa mais ou menos                                              
massa(8) = 2*areahoneycomb*espessura_honeycomb*densidade_honeycomb;  

% Massa Outros
massa(9) = 210; %Tdp - estimativa mais ou menos                                                     
massa(10) = 186.9 + 790 - m_939 + 50; %sistema moto propulsor         
%(Td do elétrico(programa do Gui Chiang + (motor + servo do motor -> planilha Margonato) - servo do motor + estimativa hélice             
massa(11) = 50; % bequilha + nariz - estimativa zuada                                                                                              
massa(12) = extensao*comp_extensao;      

% Massa EV
massa(13) = entelagem_ev*densidade_superficial_entelagem;                                           
massa(14) = capstrip_ev*densidade_1;                                                                
massa(15) = chapeado_ev*densidade_1;                                                                
massa(16) = area_nerv_ev*2.38*densidade_238; % 2.38 espessura da balsa    

%Massa EH
massa(17) = area_entelagem_eh*densidade_superficial_entelagem;                                      
massa(18) = area_capstrip_eh*densidade_1;                                                           
massa(19) = area_chapeado_eh*densidade_1;                                                           
massa(20) = area_nerv_eh*2.38*densidade_238;                                                        
massa(21) = length_cauda*densidade_boom;                                                           
massa(22) = area_spar_eh*espessura_H45*densidade_H45;                                               
peso = sum(massa)/1000;
%end



%% Comentários da Tia Isa
% [22:16, 31/05/2021] Tia Isa Urubus: Se for convencional façam as centrais de honeycomb
% [22:16, 31/05/2021] Tia Isa Urubus: Asa voadora tbm
% [22:16, 31/05/2021] Tia Isa Urubus: Mas se for biplano façam duvinycell laminado







