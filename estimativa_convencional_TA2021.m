%% Otimizador TA 2021

% fclose all;clear all;clc

% Vetores que funcionam, MTOW e CP
x = [ 0.7950 0.5612 0.4999 0.6000 0.0020 2.1302 2.2193 5.4900 0.2928 1.0000 2.8766 ]; % Melhor vetor (vetor final, antes de chegar nesse debaixo)
% x = [0.88 0.5612 0.4999 0.6000 0.0020 2.1302 2.2193 5.4900 0.26 0.7 5.5]; % Vetor final do avião do TA

%function [result] = estimativa_convencional_TA2021(x,~)
try
% x(1) = X_emp [m] Distância entre bordo de ataque da asa e bordo de ataque das empenagens (raiz)
% x(2) = Seção reto-trapezoidal [% da semi envergadura]
% x(3) = Corda da raiz [m]
% x(4) = Afilamento (corda da ponta / corda da raiz)
% x(5) = Offset da ponta da asa [m]
% x(6) = Perfil da raiz
% x(7) = Perfil do meio
% x(8) = Perfil da ponta
% x(9) = Corda da EH [m]
% x(10) = Envergadura da EH [% da envergadura da asa]
% x(11) = Ângulo de incidência na corrida da asa [°]

xerro = x;
warning off;
clear a AD ai ail alfa alfamax alpha AR Area Areac ARev ARh ARhlim ARlim ARs ARsh ARvlim ats aux1 aws b ba bah beh bev c1 c1eh c_ail Cc CD_corrida CD_factor CD_polar cdi CDi_corrida CDi_corrida2 CDi_corrida_asa CDi_corrida_emp CDi_polar cdp CDp cdp_1 cdp_2 cdp_m1 cdp_t CDplane cdpr cdv CDv cdv_1 cdv_2 cdv_m1 cdv_t cdvr chord cordeh cl CL0w CL_corrida CL_corrida2 CL_corrida_asa CL_corrida_emp CL_polar CLa Clb Clda Cldl CLeh CLmax CLmaxasa CLmaxperf CLmaxperfm CLmaxperfp CLmaxperfr CLmaxrun Clp CLplane Clr cm Cma Cmde Cmq Cnb Cnda Cndl Cnp Cnr cr creh crev ct cteh ctev CYda CYdl D d_fixa de de_lim declive deda died Drag E emp epslon f Fat fator_pista fh file_aviao file_aviao_ground filename flag fmult gamalim grad_de_subida grad_sub grad_sublim H h_asa h_emp hlim i icl ih ime inc inc_massa inc_massa_dois inc_vel inutil ipa iTracao iw j k Ka Kh Kl ks ksh L L_mais_B lambda lambda1 lambda2 lambdaeh lambdaev lambdalim lb mac maceh mi MTOWt nada nEV nT offset offset1 offset2 offset_mac offset_maceh offset_motor offseteh offsetev p_disp p_req Pcg Pcglim perfilh perfilm perfilp perfilr perfilv plane_aviao planeAVL posb posbeh prop q restr_geom result ro S Seh Sev SgnDup_aileron Sum_Fx Sum_Fy themp Thrust_factor total_pista tow tracao v V_dec v_pot vel vento Vh Vhlim Vs Vv Vvlim X x1 x1_4 x1eh x1eh_4 x_emp Xcg_cr Xcg_mac Xeh Xhinge_a Xhinge_eh Xmac Xmaceh y1 y1eh Yc Ymac Ymaceh chordeh planeAVL;

tic
result = 0;
dbstop if error
save = [x(1) x(2) x(3) x(4) x(5) x(6) x(7) x(8) x(9) x(10) x(11)];
fprintf('_________________________________________________________________\n')

%% Variáveis não declaradas
d_fixa_motor = 0.05;
offset_motor = 0.05;

% Ambiente
AD      = 1100;                                 % Altitude Densidade
ro      = 1.2259*exp(-284.2*AD/2394432);        % Densidade do Ar
declive = 0;                                    % Declive da pista (+ uphill, - downhill)
váriavel_aleatória_pra_entender_git = 25;

% Fatores de correção
CD_factor     = 1.0;
Thrust_factor = 1;
fmult         = 1;

%% Pârametros da otimização - Asa

b    = 2;                 % Envergadura da asa (m)
iw   = x(11);             % Ângulo de incidência na corrida da asa [°]
died = 0;                 % Diedro [°]
X    = (b/2)*x(2);        % Transição reto-trapezoidal

% Cordas na envergadura
chord(1)    = x(3);     cr = chord(1);        % Corda da raiz
chord(2)    = x(3);     cm = chord(2);        % Corda intermediária
chord(3)    = cr*x(4);  ct = chord(3);        % Corda da ponta
chord       = [cr cm ct];                    % Vetor corda na envergadura da asa

% Posições das cordas na envergadura 
posb(1)   = 0;                                % Início da asa
posb(2)   = X;                                % Seção do meio
posb(3)   = b/2;                              % Final da asa
posb      = [posb(1) posb(2) posb(3)];        % Vetor de posições das cordas na envergadura da asa

% Offset das cordas
offset1   = 0;                      % Offset da seção do meio da asa [m]
offset2   = x(5);                   % Offset da ponta da asa [m]
offset = [0 offset1 offset2];       % Vetor offset na asa

% Perfis da asa
perfilr = round(x(6));   % Perfil da raiz da asa
perfilm = round(x(7));   % Perfil do meio da asa
perfilp = round(x(8));   % Perfil da ponta da asa

h_asa       = b/10;    % Altura da asa em relação ao solo (intradorso da raiz ate o solo)
H_CG        = 0;       % Altura CG a partir da altura da asa
% align_perf  = 0;     % Tipo de alinhamento dos perfis da raiz e ponta (1 se for pelo extradorso, 0 se for pela linha média)
align_asa   = 0;       % Tipo de alinhamento da asa (1 se for pelo BF, 0 se for pelo BA)
% thw         = 0.;    % Thickness do perfil

%% Cálculos iniciais da ASA

y1 =(0:0.01:1)*b/2;                                % Malha para calcular a área da asa 
[y1,c1,x1] = Corda(chord,posb,X,b,offset);         % Função para calcular a distribuição de corda e offset na envergadura da asa
if isnan(c1)                                       % Teste da malha de cordas  
    disp('ops')
    PrintParada(save,1,NaN,NaN,NaN,NaN,NaN,NaN,NaN);
    return
end

S          = 2*trapz(y1,c1);              % Cálculo da área da asa usando o método do trapézio
AR         = 2*b^2/S;                     % Razão de aspecto
mac        = 2/S*trapz(y1,c1.^2);         % Cálculo da MAC
Ymac       = 2/S*trapz(y1,c1.*y1);        % Posição na envergadura da MAC
x1_4       = x1 + (c1/4);                 % Posição do CA de cada corda
Xmac       = 2/S*trapz(y1,c1.*x1_4);      % Cálculo do offset do CA da MAC
offset_mac = Xmac - mac/4;                % Offset da MAC
lambda1    = cm/cr;                       % Afilamento da raíz da asa
lambda2    = ct/cm;                       % Afilamento da ponta da asa
lambda     = ct/cr;                       % Afilmento total da asa

%% Parâmetros de Otimização - EH 

beh       = x(10);                      % Envergadura da EH [m]
creh      = x(9);                       % Corda da raiz da EH
cteh      = creh;                       % Corda da ponta da EH
chordeh   = [creh cteh];                % Vetor corda EH
posbeh    = [0 beh];                    % Vetor posição das cordas na EH
offseteh  = 0;                          % Offset da ponta da EH
x_emp     = x(1);                       % Distância entre bordo de ataque da asa e bordo de ataque da EH (raiz)

Ptdp = 0.3; % Chute de posição do TDP (% da corda da raiz) [melhor errar pra menos que pra mais]

if tand(10)*(x_emp-cr) > (tand(20-iw))*(creh+x_emp-Ptdp*cr) - h_asa
    h_emp = tand(10)*(x_emp-cr); % Altura entra asa e EH (Bordo de ataque da asa até bordo de ataque da EH)
else
    h_emp = (tand(20-iw))*(creh+x_emp-Ptdp*cr) - h_asa;
end

ih        = 0;                          % Incidência da EH
perfilh   = 'Carbonarainvertido.dat';   % Perfil da EH    
Xeh       = 0;                          % Transição

%% Cálculos iniciais da EH

y1eh = (0:0.01:1)*beh/2;                                                                % Malha para calcular a área da eh
offseteh_so_para_funcao_de_corda = [0 offseteh];
[y1eh,c1eh,x1eh] = Corda(chordeh,posbeh,Xeh,beh,offseteh_so_para_funcao_de_corda);      % Função para calcular a distribuição de corda e offset na envergadura da EH
if isnan(c1eh)                                                                          % Teste da Malha de Cordas                  
     disp('eita')
     PrintParada(save,2,NaN,S,AR,NaN,NaN,NaN,NaN);
     return
end

x1eh_4 = x1eh + (c1eh/4);        % Posição do CA de cada corda

if creh == cteh                % EH retangular   
    Seh    = beh*creh;         % Área da EH
    maceh  = creh;             % Cálculo da MAC para EH
    Ymaceh = 0;                % Posição na envergadura da MAC para EH
    Xmaceh = 0.25*creh;        % Cálculo do offset do CA da MAC para EH
else                           % EH não-retangular        
    Seh          = 2*trapz(y1eh,c1eh);
    maceh        = 2/Seh*trapz(y1eh,c1eh.^2);
    Ymaceh       = 2/Seh*trapz(y1eh,c1eh.*y1eh);
    Xmaceh       = 2/Seh*trapz(y1eh,c1eh.*x1eh_4);
end

lb           = - Xmac + x_emp + Xmaceh;        % Distância entre CAs da asa e das empenagens
ARh          = beh^2/Seh;                      % Razão de aspecto da EH
lambdaeh     = cteh/creh;                      % Afilamento da EH
offset_maceh = Xmaceh - maceh/4;               % Offset da MAC

%% Parâmetros de Otimização - EV

bev      = 0.32;                % Envergadura na EV
crev     = 0.2;                 % Corda da raiz na EV
ctev     = 0.2;                 % Corda da ponta na EV     
offsetev = 0;                   % Offset da EV
nEV      = 1;                   % Quantidade de EVs
perfilv  = 'NACA0012.dat';      % Perfil da EV  

if crev > creh
    lb = lb - (crev - creh);
end

%% Cálculos iniciais da EV

Sev      = crev*bev;        % Área da EV
ARev     = bev^2/Sev;       % Razão de aspecto da EV
lambdaev = ctev/crev;       % Afilamento da EV
macev = crev;               % MAC da EV

%% Superfícies de Comando
% Aileron da asa
ail = 0.15*b;                     % Distância na semi-envergadura onde se inicia o aileron; é pra ser um absurdo propositalmente, pq sabemos que o AVL não sabe ver deflexão

if abs(ail - posb(2)) < 0.05
    if (ail - posb(2)) > 0
        ail = ail + 0.15;
    else
        ail = ail - 0.15;
    end
end

if 0.05 <= abs(ail - posb(2)) && abs(ail - posb(2)) < 0.1
    if (ail - posb(2)) > 0
        ail = ail + 0.1;
    else
        ail = ail - 0.1;
    end
end

if 0.1 <= abs(ail - posb(2)) && abs(ail - posb(2)) < 0.15
    if (ail - posb(2)) > 0
        ail = ail + 0.05;
    else
        ail = ail - 0.05;
    end
end

c_ail = interp1(y1,c1,ail);         % Corda da asa quando inicia o aileron

ba             = b/2 - ail;         % Envergadura do aileron
Xhinge_a       = 0.25;              % [%] da corda da hinge do aileron
SgnDup_aileron = -1;                % Sinal de deflexão para DUPLICATE do AVL (-1 para convencional, -2 para diferencial)
%ail_sup        = 0;                % Aileron na asa superior (0=Não; 1= Sim)

% Tab EH
bah            = beh;               % Envergadura do tab da EH
Xhinge_eh      = 0.20;              % [%] da corda da hinge da EH
themp          = 0.12;              % Thickness da EH

%% Seleção de perfil e CLmax de cada um 

switch perfilr
%     case 1
%         perfilr = 'GID19111.dat';
%         CLmaxperf = 1.9130;
%         thicknessr = 0.1352;
%     case 2
%         perfilr = 'GID18091.dat';
%         CLmaxperf = 1.7963;
%         thicknessr = 0.1352;
    case 1
        perfilr = 'G01.dat';
        CLmaxperf = 2.5296;
        thicknessr = 0.14;
    case 2
        perfilr = 'GM07.dat';
        CLmaxperf = 2.5166;
        thicknessr = 0.1419;
    case 3
        perfilr = 'GM39.dat';
        CLmaxperf = 2.4280;
        thicknessr = 0.1431;
    case 4
        perfilr = 'RR100.dat';
        CLmaxperf = 2.4054;
        thicknessr = 0.1360;
    case 5
        perfilr = 'AC03.dat';
        CLmaxperf = 2.3238;
        thicknessr = 0.1343;
    case 6
        perfilr = 'AC04.dat';
        CLmaxperf = 2.2783;
        thicknessr = 0.1320;
%     case 7
%         perfilr = 'AC05.dat';
%         CLmaxperf = 2.1959;
%         thicknessr = 0.1240;
end
switch perfilm
%     case 1
%         perfilm = 'GID19111.dat';
%         CLmaxperfm = 1.9130;
%     case 1
%         perfilm = 'GID18091.dat';
%         CLmaxperfm = 1.7963;
    case 1
        perfilm = 'G01.dat';
        CLmaxperfm = 2.5296;
    case 2
        perfilm = 'GM07.dat';
        CLmaxperfm = 2.5166;
    case 3
        perfilm = 'GM39.dat';
        CLmaxperfm = 2.4280;
    case 4
        perfilm = 'RR100.dat';
        CLmaxperfm = 2.4054;
%     case 6
%         perfilm = 'AC03.dat';
%         CLmaxperfm = 2.3238;
%     case 8
%         perfilm = 'AC04.dat';
%         CLmaxperfm = 2.2783;
%     case 9
%         perfilm = 'AC05.dat';
%         CLmaxperfm = 2.1959;  
end
switch perfilp
%     case 1
%         perfilp = 'GID19111.dat';
%         CLmaxperfp = 1.9130;
%     case 2
%         perfilp = 'GID18091.dat';
%         CLmaxperfp = 1.7963;
%     case 3
%         perfilp = 'G01.dat';
%         CLmaxperfp = 2.5296;
%     case 4
%         perfilp = 'GM07.dat';
%         CLmaxperfp = 2.5166;
    case 1
        perfilp = 'GM39.dat';
        CLmaxperfp = 2.4280;
    case 2
        perfilp = 'RR100.dat';
        CLmaxperfp = 2.4054;
    case 3
        perfilp = 'AC03.dat';
        CLmaxperfp = 2.3238;
    case 4
        perfilp = 'AC04.dat';
        CLmaxperfp = 2.2783;
    case 5
        perfilp = 'AC05.dat';
        CLmaxperfp = 2.1959;  
end

%% Restrições da otimização

ARlim       = 4;                    % Razão de aspecto limite
ARhlim      = [2.0 6.0];            % Razão de aspecto limite da EH 
ARvlim      = 1;                    % Razão de aspecto limite da EV
grad_sublim = 3.5;                  % Gradiente de subida mínimo (%)
Pcglim      = [0.30 0.42];          % Limites do CG [% da corda da raiz]
de_lim      = [4.6 40];             % Limites de deflexão do profundor
Vhlim       = [0.100 0.52];         % Limites de volume de cauda horizontal
Vvlim       = [0.025 0.05];         % Limites de volume de cauda vertical

%% Restrições geométricas de otimização
% if posb(2) >= b/2 
%     disp('Transição > envergadura')
%     return        
% end

% if cr < cm || cr < ct || cm < ct
%     disp('Cordas da asa bugadas')
%     return
% end

% if creh < cteh || crev < ctev 
%     disp('Cordas das empenagens bugadas')  
%     return
% end


if posb(2) < b/2 
    result = result-0.1;
%         fprintf('%.4f\n',result)
end

if cr >= cm && cr >= ct && cm >= ct
    result = result-0.1;
%         fprintf('%.4f\n',result)
end

if creh >= cteh && crev >= ctev 
    result = result-0.1;
%         fprintf('%.4f\n',result)
end

% if issorted(fliplr(chord)) ~= 1
%     result = 0;
%     disp('Cordas bugadas')
%     return
% end

% if chord(1) - chord(2) > 0.35 || chord(2) - chord(3) > 0.30
%     result = 0;
%     disp('Variação de corda muito grande')
%     return
% end

% if abs(b/2 - posb(2)) <= 0.10 || abs(posb(2)) <= 0.20
%     disp('Seções muito pequenas')
%     return
% end
if abs(posb(3)-posb(2)) <= 0.1 ||  abs(posb(2)-posb(1)) <= 0.1 || abs(ail-posb(2)) < 0.05
	result = 0;
	disp('Seção muito pequena')
    PrintParada(save,3,NaN,S,AR,NaN,NaN,NaN,NaN);
	return
end

if abs(b/2 - posb(2)) > 0.10 && abs(posb(2)) > 0.20
    result = result-0.1;
%         fprintf('%.4f\n',result)
end

if died < 0
    if h_asa - abs((b/2)*tand(gama)) < 0.045
    disp('Ponta de asa muito perto do chão')
    end
    if h_asa - abs((b/2)*tand(gama)) >= 0.045
        result = result-0.1;
%         fprintf('%.4f\n',result)
    end
end

if AR < ARlim            
    disp('Razão de aspecto baixa')
    PrintParada(save,4,NaN,S,AR,NaN,NaN,NaN,NaN);
    return
end
if AR >= ARlim 
        result = result-0.1;
%         fprintf('%.4f\n',result)
end

if ARh < ARhlim(1) || ARh > ARhlim(2)             
    disp('Razão de aspecto da EH')
    PrintParada(save,5,NaN,S,AR,NaN,NaN,NaN,NaN);
    return
end
if ARhlim(1) <= ARh <= ARhlim(2)  
    result = result-0.1;
%         fprintf('%.4f\n',result)
end

% if ARev < ARvlim             
%     disp('Razão de aspecto da EV')
%     PrintParada(save,7,NaN,S,AR,NaN,NaN,NaN,NaN);
%     return
% end
% if ARev >= ARvlim   
%     result = result-0.1;
% %         fprintf('%.4f\n',result)
% end

% if lambda < lambdalim
%     disp('Afilamento total grande')
%     PrintParada(save,8,NaN,S,AR,NaN,NaN,NaN,NaN);
%     return
% end
% if lambda >= lambdalim
%     result = result-0.1;
%         %fprintf('%.4f\n',result)
% end
% 
% if (lambda1 || lambda2) < 0.5
%     disp('Afilamento entre seções muito grande')
%     PrintParada(save,9,NaN,S,AR,NaN,NaN,NaN,NaN);
%     return
% end
% if (lambda1 && lambda2) >= 0.5
%     result = result-0.1;
%      %   fprintf('%.4f\n',result)
% end

% if S < 0.80
%    result = 0;
%    disp('Área insatisfatória')
%    return
% end

if cr*thicknessr > 0.1
    disp('Nervura com altura superior a 10cm')
    PrintParada(save,6,NaN,S,AR,NaN,NaN,NaN,NaN);
    return
end

if cr*thicknessr < 0.1
    result = result-0.1;
end

% Altura efetiva entre asa e EH
hh = h_emp + thicknessr*chord(1); 

if atand((h_asa + hh)/(x_emp - cr+creh)) < 18-iw        % Pré-verificação de tail-strike
    disp('Possível tail-strike') 
    PrintParada(save,7,NaN,S,AR,NaN,NaN,NaN,NaN);
    return
end

 %% Peso
 [peso] = peso_convencional_TA2021(posb,cr,cm,ct,perfilr,perfilm,perfilp,b,offset,X,bev,beh,creh,h_emp,x_emp,ctev,crev);
 
%% Volumes de cauda

Vh = lb*Seh/(S*mac);        % Volume de cauda horizontal
Vv = lb*Sev/(S*b);          % Volume de cauda vertical

if Vh < Vhlim(1) || Vh > Vhlim(2)
     disp('Volume de cauda horizontal')
     PrintParada(save,8,peso,NaN,S,AR,NaN,NaN,NaN);
     return
else       % Vh >= Vhlim(1) && Vh <= Vhlim(2)
     result = result-0.1;
     %fprintf('%.4f\n',result)
end

% if Vv < Vvlim(1) %|| Vv > Vvlim(2)
%     Vv
%      disp('Volume de cauda vertical')
%      PrintParada(save,9,peso,NaN,S,AR,NaN,NaN,NaN);
%      return
% else        % Vv >= Vvlim(1) && Vv <= Vvlim(2)
%     result = result-0.1;
%     fprintf('%.4f\n',result)
% end

%% Criação das partes do avião em .avl e cálculo analítico das derivadas longitudinais

alpha1 = -4;        % menor alpha a ser analisado
alpha2 = 12;        % maior alpha a ser analisado

[file_asa,file_eh] = convencionalAVLpartes(S,mac,b,posb,offset1,offset2,died,chord,perfilr,perfilm,perfilp,perfilh,beh,offseteh,creh,cteh,iw,ih,Xmac,Seh,maceh);
[CL,CD,Cm]         = coefsAVL(alpha1,alpha2,0,file_asa,file_eh);

ang = linspace(alpha1,alpha2,10);       % Range de ângulos analisados
n   = 0.97;                             % Queda de pressão dinâmica assumida

 for Pcg = Pcglim(2):-0.005:Pcglim(1)       % Cálculo da posição de CG a partir dos limites dados

    Xcg     = Pcg*cr;                          % Posição do CG a partir do BA [m]
    Xcg_mac = (Xcg - offset_mac)/mac;          % Posição do CG em porcentagem da mac
    Vh_cg   = Vh-(Seh/S)*(Xcg_mac-0.25);       % Volume de cauda em relação ao CG (Etkin)

    % Momentos principais
    Cm_M_asa = Cm.w(ang);                       % Cm da asa
    Cm_L_asa = CL.w(ang)*(Xcg_mac-0.25);        % Momento do CL da asa
    Cm_D_asa = -CD.w(ang)*(H_CG/mac);           % Momento do CD da asa
    
    Cm_M_eh = Cm.h(ang)*(Seh/S)*(maceh/mac);                % Cm da EH
    Cm_L_eh = -CL.h(ang)*n*Vh_cg;                           % Momento do CL da EH
    Cm_D_eh = CD.h(ang)*n*(Seh/S)*(h_emp-H_CG)/maceh;       % Momento do CD da EH
       
    Cm_D_ev = nEV*0.012796*n*(Sev/S)*((h_emp+bev/2)-H_CG)/macev;        % Momento do CD da EV

    % Coeficientes totais
    CLtot     = CL.w(ang)+n*(Seh/S)*CL.h(ang);                              % CL total do avião
    CDtot     = CD.w(ang)+n*(Seh/S)*CD.h(ang)+nEV*0.012796*n*(Sev/S);       % CD total do avião
    Cmtot_asa = Cm_M_asa+Cm_L_asa+Cm_D_asa;                                 % Cm total da asa
    Cmtot_eh  = Cm_M_eh+Cm_L_eh+Cm_D_eh;                                    % Cm total da EH
    Cmtot     = Cmtot_asa+Cmtot_eh+Cm_D_ev;                                 % Cm total do avião
    
    % Parâmetros principais
    H      = (-diff(Cmtot))./diff(CLtot)*100;       % Margem estática em porcentagem
    h_calc = H(~isnan(H));
    Xn     = Xcg_mac+h_calc/100;                    % Ponto neutro
    Cma    = diff(Cmtot);                           % Derivada de Cm por alpha
    
    Cmtotal_asa = @(alpha)interp1(ang,Cmtot_asa,alpha);        % Cm total da asa interpolado
    Cmtotal_eh  = @(alpha)interp1(ang,Cmtot_eh,alpha);         % Cm total da EH interpolado
    
    % Verificar se o alfa de trimagem é positivo para alguma incidência de EH viável:
%     for ih = 0:-1:-3
%         Cm0 = Cmtotal_asa(-iw) + Cmtotal_eh(-iw+ih)+ Cm_D_ev;
%         if Cm0 > 0
%             break
%         end
%     end
    
    % Verificar se a margem estática é positiva para todos os ângulos e se Cm0 é positivo, para o loop
%     if all(H > 0) && Cm0 > 0
%         break
%     end
 end
 
%% Verificações de Estabilidade Longitudinal

if any(H < 0)
        fprintf('Margem estática = %.4f\n',H(length(H)))
        PrintParada(save,9,peso,Pcg,S,AR,NaN,NaN,NaN);
        return
elseif all(H >= 0)
    result = result-0.1;
    %fprintf('%.4f\n',result)
end

% if Cm0 <= 0
%     fprintf('Cm0 não positivo\n')
%     PrintParada(save,10,peso,Pcg,S,AR,NaN,NaN,NaN);
%     return
% else       % Cm0 > 0
%     result = result-0.1;
% end

if any(Cma > 0)
    disp('Longitudinalmente instável')
    PrintParada(save,11,peso,Pcg,S,AR,NaN,NaN,NaN);
    return
elseif all(Cma <= 0)
    result = result-0.1;
end

%% Criação do avião inteiro em .avl e cálculo das derivadas látero-direcionais
[file_aviao] = convencionalAVL(S, mac, b, posb, offset, died, ail, c_ail, Pcg, chord, Xhinge_eh, Xhinge_a, SgnDup_aileron, perfilr, perfilm, perfilp, perfilh, perfilv, beh, offseteh, creh, cteh, bev, crev, ctev, offsetev, x_emp, h_emp, iw, ih, nEV, H_CG,1, h_asa);
[~,Clb,~,Cnb,~,~,~,~,~,~,~,~,~,~,~,~,~,~] = derivadasAVL(4,4,mac,file_aviao);

%% Verificações de Estabilidade Látero-Direcional

if Cnb < 0 
    disp('Látero-direcionalmente instável')
    delete(file_aviao)
    PrintParada(save,12,peso,Pcg,S,AR,NaN,NaN,NaN);
    return
else       % Cnb >=0
    result = result - 0.1;
end

if Clb > 0 
    disp('Látero-direcionalmente instável')
    delete(file_aviao)
    PrintParada(save,12,peso,Pcg,S,AR,NaN,NaN,NaN);
    return
else       % Clb <= 0 
    result = result-0.1;
end

%% File AVL (efeito solo)
[file_aviao_ground] = convencionalAVL(S,mac,b,posb,offset,died,ail,c_ail,Pcg,chord,Xhinge_eh,Xhinge_a,SgnDup_aileron,perfilr,perfilm,perfilp,perfilh,perfilv,beh,offseteh,creh,cteh,bev,crev,ctev,offsetev,x_emp,h_emp,iw,ih,nEV,H_CG,1,h_asa);
       %[file_aviao] = convencionalAVL(S,mac,b,posb,offset,died,ail,c_ail,Pcg,chord,Xhinge_eh,Xhinge_a,SgnDup_aileron,perfilr,perfilm,perfilp,perfilh,perfilv,beh,offseteh,creh,cteh,bev,crev,ctev,offsetev,x_emp,h_emp,iw,ih,nEV,H_CG,ef_solo,h_asa)
%fprintf('Tempo depois do AVL do efeito solo:%f\n',toc)

%% Efeito solo
% Parâmetros em cada posição da asa
[Yc,FsCL0,Area,~,~,cl,CDi_corrida,CL_corrida] = FSFT_dadosgroundAVL(-iw,0,ro,12,file_aviao_ground);
de = 0;

cdv = NaN(1,length(Area));
cdp = NaN(1,length(Area));

for j=1:length(Area)
    % Adição da viscosidade da primeira seção
    if abs(Yc(j)) <= posb(2)        
        [cdvr, cdpr] = ViscosidadePolar(perfilr, cl(j));
        cdv(j) = cdvr;
        cdp(j) = cdpr;
    % Adição da viscosidade da 2ª seção
    else
        prop = (posb(3) - abs(Yc(j)))/(posb(3)-posb(2));
        [cdv_1, cdp_1] = ViscosidadePolar(perfilr, cl(j));
        [cdv_2, cdp_2] = ViscosidadePolar(perfilp, cl(j));   

        cdv(j) = prop*cdv_1 + (1 - prop)*cdv_2;
        cdp(j) = prop*cdp_1 + (1 - prop)*cdp_2;
    end
end

% Arrasto na corrida de pista
CDvg = 4*sum(reshape(cdv,size(Area)).*Area)/S;              % Viscoso calculado através da análise de perfil no XFLR5
CDpg = 4*sum(reshape(cdp,size(Area)).*Area)/S;              % Parasita calculado através da análise de perfil no XFLR5, dobrado por serem duas asas
Cdg = CDvg+CDpg;
CD_corrida = CD_factor*(CDi_corrida + CDvg);                 % Total

if CL_corrida < 0.4 
    disp('CL corrida muito baixo')
    delete(file_aviao)
    delete(file_aviao_ground)
    delete(FsCL0)
    PrintParada(save,13,peso,Pcg,S,AR,CD_corrida,CL_corrida,NaN);
    return
end

if CL_corrida >= 0.4 
    result = result-0.1;
end

%% Efeito solo analítico
% Asa
f = 2*(h_asa)/b;
ks = 46.87*f/(f^2+40.12*f+10.86);     % Curva zuadas do livro
ARs = AR/ks;                          % AR corrigida devido ao efeito solo

% EH
fh = 2*(h_asa+h_emp)/b;                 % Cálculo do elemento f da formula descrita em ROSKAM, levando em conta a altura da fuselagem, a distancia da fuselagem ao solo e a distancia entre fuselagem e EH 
ksh = 46.87*fh/(fh^2+40.12*fh+10.86);   % Curva do livro
ARsh = ARh/ksh;                         % Razão de aspecto modificada. Uma EH com ARs em cruzeiro é equivalente a uma com AR sofrendo o efeito solo.

% Derivadas aerodinâmicas (Clalfa)
% Asa
aws = 2*pi*ARs/(2+sqrt(ARs^2/0.98^2+4));    % Coeficiente angular de Cl(alpha) da asa, obtida de Raymer (eq. 12.6), considerando o efeito solo

% EH
ats = 2*pi*ARsh/(2+sqrt(ARsh^2/0.98^2+4));  % Coeficiente angular de Cl(alpha) da EH, obtida de Raymer (eq. 12.6), considerando o efeito solo

% CL 0 da asa
[CL0w,~,~,~] = CL0asaAVL(FsCL0);

% Downwash
Ka = 1/AR-1/(1+AR^1.7);             % Fator de correção da razão de aspecto
Kl = (10-3*lambda)/7;               % Fator de correçao do afilamento
Kh = (1-h_emp/b)/((2*lb/b)^(1/3));  % Fator de correção EH
deda = 4.44*(Ka*Kl*Kh)^(1.19);      % Fórmula do gradiente de downwash
epslon = deda*(deg2rad(iw) + CL0w/aws); 

% CL corrida analítico da asa
CL_corrida_asa = CL0w + aws*deg2rad(iw);
CL_corrida_emp = 0 + ats*(deg2rad(ih) - epslon);
CL_corrida2 = CL_corrida_asa + 0.9*(Seh/S)*CL_corrida_emp;

CDi_corrida_asa = CL_corrida_asa^2/(pi*ARs*0.98);
CDi_corrida_emp = CL_corrida_emp^2/(pi*ARsh*0.98);
CDi_corrida2 = CDi_corrida_asa + 0.95*(Seh/S)*CDi_corrida_emp;

fprintf('\nCL corrida avl = %.4f\n',CL_corrida)
fprintf('CL corrida analítico = %.4f\n',CL_corrida2)
fprintf('CD corrida = %.4f\n',CD_corrida)
fprintf('CDi corrida avl = %.4f\n',CDi_corrida)
fprintf('CDi corrida analítico = %.4f\n',CDi_corrida2)
% fprintf('CDp analítico = %.4f\n',CDp)
%fprintf('k = %.4f\n', k_prandtl)
%fprintf('Tempo antes de iniciar inputs do CLmax:%f\n',toc)

%% CLmax convencional
% CLmax asa
flag = 0;
% estol_ponta_cima = 0;
% estol_ponta_baixo = 0;
beta=0;
V=12;
fileASA=file_aviao;
parfor i = 0:50
    alfa=i/2;  %#ok<PFTUSE>
    ID1 = strcat('comandosAVL_',int2str(i),'.dat');
    ID2 = strcat('FsCLmax_',int2str(i),'.dat');
    if exist(ID1,'file') == 2      % Evita arquivos já existentes
        delete(ID1)
        delete(ID2)
    end
    fileID = fopen(ID1,'w');
    fprintf(fileID,  'load %s \n',fileASA);
    fprintf(fileID,  'oper \n');
    fprintf(fileID,  'M \n');
    fprintf(fileID,  'V \n%.4f \n',V);
    fprintf(fileID,  'D \n%.4f \n\n',ro);
    fprintf(fileID,  'a \na \n');
    fprintf(fileID,  '%.4f\n',alfa);
    fprintf(fileID,  'b \nb \n');
    fprintf(fileID,  '%.4f\n',beta);
    fprintf(fileID,  'D1 \nD1\n 0\n');
    fprintf(fileID,  'D2 \nD2\n 0\n');
    fprintf(fileID,  'X \n');
    fprintf(fileID,  'FS \n%s \n',ID2);
    fprintf(fileID,  'O\n');
    fclose('all');
    % antigo (para caso de algum problema)
    % fprintf(fileID,  'load %s ',fileASA);
    % fprintf(fileID,  char([13,10]));
    % fprintf(fileID,  'oper \n');
    % fprintf(fileID,  char([13,10]));
    % fprintf(fileID,  'M \n');
    % fprintf(fileID,  char([13,10]));
    % fprintf(fileID,  'V') ;
    % fprintf(fileID,  char([13,10]));
    % fprintf(fileID, '%.4f' ,V);
    % fprintf(fileID,  char([13,10]));
    % fprintf(fileID,  'D ' );
    % fprintf(fileID,  char([13,10]));
    % fprintf(fileID,  '%.4f ' ,ro);
    % fprintf(fileID,  char([13,10]));
    % fprintf(fileID,  ['oper' char([13,10])]);
    % fprintf(fileID,  ['a' char([13,10]) 'a ' char([13,10])]);
    % fprintf(fileID,  '%.4f',alpha);
    % fprintf(fileID,  char([13,10]));
    % fprintf(fileID, [ 'b ' char([13,10]) 'b' char([13,10])]);
    % fprintf(fileID,  '%.4f',beta);
    % fprintf(fileID,  char([13,10]));
    % fprintf(fileID,  ['D1 ' char([13,10]) 'D1' char([13,10]) ' 0' char([13,10])]);
    % fprintf(fileID,  ['D2 ' char([13,10]) 'PM' char([13,10]) ' 0' char([13,10])]);
    % fprintf(fileID,  ['X ' char([13,10])]);
    % fprintf(fileID,  ['FS ' char([13,10]) ]);
    % fprintf(fileID, 'FsCLmax%u.dat ',i);
    % fprintf(fileID,  char([13,10]));
    % fprintf(fileID,  ['O' char([13,10]) char([13,10]) char([13,10]) char([13,10]) char([13,10]) char([13,10]) char([13,10]) char([13,10]) char([13,10]) char([13,10]) char([13,10]) ]);
    exec = strcat('avl.exe < ',ID1);
    [~,~] = dos(exec);
end

for k = 0:50
    alfa = k/2; 
    
    ID2 = strcat('FsCLmax_',int2str(k),'.dat');
   
    filename = ID2;
    fileID2 = fopen(filename,'r');
    data = textscan(fileID2,'CLsurf  =   %f     CDsurf  =   %f','Headerlines',15);
    CLmaxrun = data{1};
    %CDrunasa = data{2};
    
    formatSpec = '%f %f %f %f %f %f %f %f %f %f %f %f %f %[^\n\r]';
    data = textscan(fileID2,formatSpec,'Headerlines',4);
    cl = data{:,7};
    Yc = data{:,2};
    
%     clear data
%     data = textscan(fileID2,'CLsurf  =   %f     CDsurf  =   %f','Headerlines',38);
%     CLmaxrunsup = data{1};
    % CDrunsup = data{2};
    
%     clear data
%     data = textscan(fileID2,formatSpec,'Headerlines',4);
%     clsup = data{:,7};
    
%     CLmaxrun = (CLmaxruninf+CLmaxrunsup)/2;
    
    %CDrun = (CDruninf+CDrunsup)/2;
    
    %nada = 0;
    fclose all;

    for i = 1:length(cl)
        
        if Yc(i)>=posb(1) && Yc(i) <= posb(2)
            prop = (posb(2) - abs(Yc(i)))/(posb(2)-posb(1));
            CLmaxperf_int = prop*CLmaxperf + (1-prop)*CLmaxperfm;

        elseif Yc(i) >= posb(2) && Yc(i) <= posb(3)
             prop = (posb(3) - abs(Yc(i)))/(posb(3)-posb(2));
             CLmaxperf_int = prop*CLmaxperfm + (1-prop)*CLmaxperfp;
        end
        
        if cl(i) >= CLmaxperf_int 
            flag = 1;
            if cl(i) >= CLmaxperf_int
                [clmax,j] = max(cl);
                if Yc(j) >= 0.6*b/2
                    disp('Estol no aileron');
                    delete(file_aviao)
                    delete(file_aviao_ground)
                    delete(FsCL0)
                    PrintParada(save,14,peso,Pcg,S,AR,CD_corrida,CL_corrida,NaN);
                    return
                    %estol_ponta_cima = 1;
                end
            end
            break
        end
    end
    if flag
        break
    end
end
parfor i=0:50
    delete(strcat('FsCLmax_',int2str(i),'.dat'))
    delete(strcat('comandosAVL_',int2str(i),'.dat'))
end
%fprintf('Tempo antes do AVL do CLmax:%f\n',toc)

%% Verificação de CLmax e estol sem efeito solo
CLmaxasa = CLmaxrun;
alfamax = alfa;

% Intervalo de alpha
alpha = -2-iw:0.75:alfamax;

if alfamax <= iw
    result = 0;
    PrintParada(save,15,peso,Pcg,S,AR,CD_corrida,CL_corrida,NaN);
    disp('Incidência estola')
    return
end

fprintf('CL max asa = %.4f\n',CLmaxasa)
fprintf('Alfamax max = %.4f\n',alfamax+iw)
fprintf('Rotate = %.4f\n',atand((h_asa + h_emp)/(x_emp+creh)))
fprintf('Incidencia da asa = %.4f\n',iw)
fprintf('Incidencia da EH = %.4f\n',ih)

% if atand((h_asa + h_emp)/(x_emp+creh)) < alfamax-iw
%     disp('Sem rotate')
%     delete(file_aviao)
%     delete(file_aviao_ground)
%     delete(FsCL0)
%     PrintParada(save,16,peso,Pcg,S,AR,CD_corrida,CL_corrida,NaN);
%     return
% end
if atand((h_asa + h_emp)/(x_emp+creh)) >= alfamax-iw
    result = result-0.1;
end

%% Polar de arrasto convencional 
clear i
i = 1;
beta=0;
fileASA=file_aviao;

CDi_polar_ext = NaN(1,length(alpha));
CL_polar_ext = NaN(1,length(alpha));
CD_polar_ext = NaN(1,length(alpha));
CDv = NaN(1,length(alpha));
CDp = NaN(1,length(alpha));
Cd = NaN(1,length(alpha));


parfor k=1:length(alpha-1)
    ID1 = strcat('temp_commands_',int2str(k),'.dat');
    ID2 = strcat('Fsurface_polar_',int2str(k),'.dat');
    ID3 = strcat('Ftotal_polar_',int2str(k),'.dat');
    
    fileID = fopen(ID1,'w');
    
    fprintf(fileID,  'load %s \n',fileASA);
    fprintf(fileID,  'oper \n');
    fprintf(fileID,  'a \na \n');
    fprintf(fileID,  '%.4f\n',alpha(k));
    fprintf(fileID,  'b \nb \n');
    fprintf(fileID,  '%.4f\n',beta);
    fprintf(fileID,  ' D1\nD1\n 0\n');
    fprintf(fileID,  'D2 \nD2\n 0\n');
    fprintf(fileID,  'X \n');
    fprintf(fileID,  'FS \n%s \n',ID2);
    fprintf(fileID,  'FT \n%s \n',ID3);
    fclose('all');
    
    exec = strcat('avl.exe < ',ID1);  % Executa comandos no XFOIL
    [~,~] = dos(exec);
    
    delete(ID1)
end


for k = 1:length(alpha-1)  
    %% Forças por Strip, Arrasto Induzido Polar e Sustentação para o alpha dado
    ID2 = strcat('Fsurface_polar_',int2str(k),'.dat');
    fileID2 = fopen(ID2,'r');
    clear data
    data = textscan(fileID2,'CLsurf  =   %f     CDsurf  =   %f','Headerlines',15);
    CDi_asa = data{2};
    
    %% Forças por Strip
    formatSpec = '%f %f %f %f %f %f %f %f %f %f %f %f %f %[^\n\r]';
    data = textscan(fileID2,formatSpec,'Headerlines',4);
    Yc = data{:,2};
    Cc = data{:,3};
    Areac = data{:,4};
    ai = rad2deg(data{:,6});
    cl_norm = data{:,7};
    cl = data{:,8};
    %%
    clear data
    data = textscan(fileID2,'CLsurf  =   %f     CDsurf  =   %f','Headerlines',45);
    CDi_eh = data{2};
    clear data
    data = textscan(fileID2,'CLsurf  =   %f     CDsurf  =   %f','Headerlines',69);
    CDi_ev = data{2};
    
    %% Sustentação para o alpha dado
    ID3 = strcat('Ftotal_polar_',int2str(k),'.dat');
    fileID3 = fopen(ID3,'r');
    CL_polar = textscan(fileID3,'  CLtot = %f','Headerlines',23);
    CL_polar = CL_polar{1};
    %%
    
    
    CDi_polar = textscan(fileID3,'  CDtot = %f','Headerlines',0);
    CDi_polar = CDi_polar{1};
    CDi_polar = abs(CDi_polar);
    
    if CDi_eh < 0
        CDi_eh = abs(CDi_eh);
        fact = abs(CDi_polar-CDi_asa)/CDi_eh;
        CDi_polar = CDi_polar + 2*fact*CDi_eh;
    end
    
    if CDi_ev < 0
        CDi_ev = abs(CDi_ev);
        fact = abs(CDi_polar-CDi_asa)/CDi_ev;
        CDi_polar = CDi_polar + 2*fact*CDi_ev;
    end
    
    fclose all;
    delete(ID2);
    delete(ID3);
    
    
    CDi_polar_ext(i) = CDi_polar;
    CL_polar_ext(i) = CL_polar;
    
    for j=1:length(Yc)   
        % Adição da viscosidade da primeira seção
        if abs(Yc(j)) <= posb(2)
            [cdvr, cdpr] = ViscosidadePolar(perfilr, cl(j));

            cdv(j) = cdvr;
            cdp(j) = cdpr;
            
        % Adição da viscosidade da segunda seção    
        else
            prop = (posb(3) - abs(Yc(j)))/(posb(3)-posb(2));
            [cdv_m1, cdp_m1] = ViscosidadePolar(perfilr, cl(j));
            [cdv_t, cdp_t] = ViscosidadePolar(perfilp, cl(j));   

            cdv(j) = prop*cdv_m1 + (1 - prop)*cdv_t;
            cdp(j) = prop*cdp_m1 + (1 - prop)*cdp_t;    
        end
    end
    
    CDv(i) = 2*sum(reshape(cdv,size(Areac)).*Areac)/S;           % Viscoso calculado através da análise de perfil no XFLR5
    CDp(i) = 2*sum(reshape(cdp,size(Areac)).*Areac)/S;           % Parasita calculado através da análise de perfil no XFLR5
    Cd(i) = CDv(i) + CDp(i);
    CD_polar_ext(i) = CD_factor*(CDi_polar_ext(i) + CDv(i));             % Total

    % CL do ângulo analisado
    CL_polar_ext(i) = 1.0*CL_polar_ext(i);
    
    
    % CL_polar é o cl max DO AVIÃO INTEIRO (asa+eh+ev+fuselagem...) 
    % com valor que sai do avl, mas limitado pelo de estol dado pelo xf 
    % (definido no for pelo length(alpha-1), saído de arquivo Ft do avl
    
    i = i + 1;
end

% Análise de desempenho
CLmax = max(CL_polar_ext);
% monoCLmax = max(monoCL_polar);
% CLCDmax = max(CL_polar./CD_polar);

%fprintf('Tempo antes do calculo da corrida de pista:%f\n',toc)


%% Corrida de pista 
nT = 1;             % Número de trações sendo analisadas
Vs = NaN(1,nT);    
MTOWt = NaN(1,nT);
grad_sub = NaN(1,nT);

for iTracao = 1:nT
    tracao = (ro/1.1092)*[-0.0206 -0.5838   39.5820]*0.87; %13x4 07/19 entre 85% e 87.5%
    %tracao = Thrust_factor*(ro/1.0840)*[-0.0562 -0.5309 46.1025]; % Tmotor 28x9.2 U8Lite 15/02
    
    tow = 5;                % Chute inicial do TOW
    mi = 0.0875;            % Coeficiente de atrito
    total_pista = 90;       % Distância segura para rotate
    inc = 0.05;              % Tamanho de cada segmento a ser analisado
    T = total_pista/inc;
    inc_massa = 0.01;       % Incremento na massa a cada iteração

    
     % Pré alocando
    q = NaN(1,T);
    E = NaN(1,T);
    D = NaN(1,T);
    L = NaN(1,T);
    Fat = NaN(1,T);
    Sum_Fx = NaN(1,T);
    
    Sum_Fy = NaN(1,T);
    a = NaN(1,T);
    v = NaN(1,T+1);
    
    v(1) = 0;                                             % Velocidade inicial
    aux1 = 0;                                             % Variável auxiliar para quando chegar no limite do MTOW
    
    % Cálculo da corrida
    while aux1 == 0
        for i = 1:T
            
            q(i) = abs(0.5*ro*v(i)^2);                                     % Pressão dinamica
            V_dec = 1.1*sqrt((2*9.81*tow)/(ro*S*CLmax));                   % Velocidade (com correçao de 5%) de decolagem
            E(i) = abs(polyval(tracao, v(i)));                             % Tração
            D(i) = q(i)*CD_corrida*S;                                      % Cálculo do arrasto na decolagem
            L(i) = q(i)*CL_corrida*S;                                      % Cálculo da sustentação
            
            Sum_Fy(i) = L(i)-tow*9.81*cosd(declive);                        % Somatória das forças em y
            Fat(i) = abs(mi*(Sum_Fy(i)));                                   % Cálculo da força de atrito
            
            %Sum_Fx(i) = E(i)-D(i)-Fat(i)-9.81*tow*sind(declive);            % Somatória das forças em x
            Sum_Fx(i) = E(i)-D(i)-Fat(i);                                   % Somatória das forças em x
            
            % Verificação se o avião decolou com o peso escolhido
            if v(i) > V_dec && i*inc < total_pista
                fator_pista = i*inc/total_pista;
                inc_massa_dois = max([inc_massa/fator_pista  inc_massa]);
                tow = tow + inc_massa_dois;
                pista_percorrida = i*inc;
                break
            end
            
            % Caso o avião não decole
            if v(i) > V_dec || i*inc >= total_pista
                tow = tow - inc_massa;
                aux1 = 1;
            end
            
            a(i) = Sum_Fx(i)/tow;
            v(i+1) = sqrt(v(i)^2 + 2*a(i)*inc);
        end
    end
    %v_stall = v(end);
    
    % tow a partir daqui significa o MTOW do avião <<
    Vs(iTracao) = sqrt((2*9.81*tow)/(ro*S*CLmax));
    
    % Curva de potência disponível
    inc_vel = 0.5;
    v_pot = Vs(iTracao):inc_vel:2.3*Vs(iTracao);            % Envelope de velocidades
    
    % Pré alocando
    p_disp = NaN(1,length(v_pot));
    p_req = NaN(1,length(v_pot));
    grad_de_subida = NaN(1,length(v_pot));
    
    for i = 1:size(v_pot,2)
        p_disp(i) = v_pot(i)*(polyval(tracao, v_pot(i)));   % Cálculo da potência disponível do motor, levando em conta a tração do avião e a velocidade do avião.
    end
    
    % Curva de potência requerida
    for i = 1:size(v_pot,2)
        CLplane = (2*9.81*tow)/(ro*S*v_pot(i)^2);          % CL total necessário para a carga na velocidade v_pot(i)
        
        CDplane = interp1(CL_polar_ext,CD_polar_ext,CLplane);       % Arrasto interpolado da polar
        Drag = 0.5*ro*(v_pot(i)^2).*CDplane*S;             % Força de arrasto parasita da aeronave
        
        p_req(i) = v_pot(i)*Drag;                           % Potência requerida do motor, levando em conta a velocidade e o arrasto parasita.
    end
    
    for i = 1:size(v_pot,2)
        grad_de_subida(i) = 100*(p_disp(i) - p_req(i))/((tow*9.81)*v_pot(i)); % Gradiente de subida para cada velocidade [%]
        
        % Gradiente de subida logo após a decolagem
        
        if v_pot(i) > 1.1*Vs(iTracao)-inc_vel && v_pot(i) < 1.1*Vs(iTracao)+inc_vel
            grad_sub(iTracao) = grad_de_subida(i);            
        end
    end
    
MTOWt(iTracao) = tow;   % Vetor de MTOW para cada tração

end
%fprintf('Tempo depois do calculo da corrida de pista:%f\n',toc)


%% Verificação gradiente de subida 

if grad_sub <= 0
    disp('Gradiente de subida negativo fodase')
    delete(file_aviao)
    PrintParada(save,17,peso,Pcg,S,AR,CD_corrida,CL_corrida,max(MTOWt));
    return
elseif grad_sub < grad_sublim    
    disp('Gradiente de subida lamentável')
    delete(file_aviao)
    PrintParada(save,18,peso,Pcg,S,AR,CD_corrida,CL_corrida,max(MTOWt));
    return
end

j = find(MTOWt == max(MTOWt), 1);
CL_polardec = CL_polar(j);
alpha_dec = alpha(j);
tow = max(MTOWt);

%% Deflexão de profundor
alpha = 0:2:alfamax;
de = NaN(1,length(alpha));
i = 1;

for k = 1:length(alpha)        % Análise da deflexão do profundor pelo avl
    [de(i)] = deflexaoAVL(alpha(k),0,file_aviao);    
    i = i + 1;
end

if max(abs(de)) > de_lim(1) && max(abs(de)) < de_lim(2)
    fprintf('Máxima deflexão de profundor: %.4f\n',max(abs(de)))
    result = result - 0.1;
end

if max(abs(de)) < de_lim(1) || max(abs(de)) > de_lim(2)
    PrintParada(save,19,peso,Pcg,S,AR,CD_corrida,CL_corrida,max(MTOWt))
    disp('Deflexão de profundor')
    return
end

%% Peso

cargapaga = max(MTOWt) - peso;

result = -cargapaga;

fprintf('CL decolagem = %.4f\n',CL_polardec) %CL max q um ponto específico da asa gera
fprintf('Alfa decolagem = %.4f\n',alpha_dec)
fprintf('CL max = %.4f\n',CLmax)
fprintf('CL max asa = %.4f\n',CLmaxasa)
fprintf('AR = %.4f\n',AR)
fprintf('Área = %.4f\n',S)
fprintf('Gradiente de subida = %.4f\n',grad_sub(j))
fprintf('Vs = %.4f\n',Vs(j))
fprintf('Margem estática = %.4f\n',H(length(H)))
fprintf('Posiçao do CG [raiz] = %.4f\n',Pcg)
fprintf('Posiçao do CG [MAC] = %.4f\n',Xcg_mac)
fprintf('Tempo de uma iteração completa = %.4f\n',toc)
fprintf('Peso = %.4f <<<\n\n',peso)
fprintf('MTOW = %.4f <<<\n\n',max(MTOWt))
fprintf('Carga Paga = %.4f <<<\n\n', cargapaga)

fileIDBOM = fopen('backup.txt','a');
fprintf(fileIDBOM,'\n%s\n\nx = [ %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f ]',datetime,save(1),save(2),save(3),save(4),save(5),save(6),save(7),save(8),save(9),save(10),save(11));
fprintf(fileIDBOM,'\n\nMTOW = %.4f <<<\n\nPeso = %.4f\n\nCarga Paga = %.4f\n\nPosiçao do CG [raiz] = %.4f',max(MTOWt),peso ,cargapaga ,Pcg);
fprintf(fileIDBOM,'\n\nÁrea: %.4f | AR: %.4f | CDcorrida: %.4f | CLcorrida: %.4f | Parada: %f',S,AR,CD_corrida,CL_corrida,0);
% if estol_ponta_baixo == 1
%         fprintf(fileIDBOM,'\n\nCuidado: estol na ponta da asa inferior');
% end
% if estol_ponta_cima == 1
%         fprintf(fileIDBOM,'\n\nCuidado: estol na ponta da asa superior');  
% end
fprintf(fileIDBOM,'\n------------------------------------------------------------------------------------------------------------------------------');
delete(file_aviao_ground)
delete(file_aviao)

catch
fprintf('ERRO :algum erro aleatorio do avl\n');
result =0;
file_aviao = strcat('BiplanoAVL.avl');

if exist(file_aviao,'file') == 2
    delete(file_aviao)
    ID1 = strcat('temp_commands.dat');
    ID2 = strcat('derivadasAVL.dat');
    fclose('all');
    delete(ID1)
    delete(ID2)
end

PrintParada(save,20,peso,NaN,NaN,NaN,NaN,NaN,NaN);
return

end
%end
