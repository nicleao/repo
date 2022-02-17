%% Otimizador TA 2021

fclose all;clear all;clc

% Vetores que funcionam, MTOW e CP
%x = [1.35 0.4 0.4 0.75 0.07 1 1 1 0.18 0.85 5]; % Vetor de testes apenas
x = [0.88 0.5612 0.4999 0.6000 0.0020 2.1302 2.2193 5.4900 0.26 0.7 5.5];

% function [result] = estimativa_convencional(x,~)
tic
try
% x(1) = X_emp [m] Dist�ncia entre bordo de ataque da asa e bordo de ataque das empenagens (raiz)
% x(2) = Se��o reto-trapezoidal [% da semi envergadura]
% x(3) = Corda da raiz [m]
% x(4) = Afilamento (corda da ponta / corda da raiz)
% x(5) = Offset da ponta da asa [m]
% x(6) = Perfil da raiz
% x(7) = Perfil do meio
% x(8) = Perfil da ponta
% x(9) = Corda da EH [m]
% x(10) = Envergadura da EH [% da envergadura da asa]
% x(11) = �ngulo de incid�ncia na corrida da asa [�]

xerro = x;
warning off;
clear a AD ai ail alfa alfamax alpha AR Area Areac ARev ARh ARhlim ARlim ARs ARsh ARvlim ats aux1 aws b ba bah beh bev c1 c1eh c_ail Cc CD_corrida CD_factor CD_polar cdi CDi_corrida CDi_corrida2 CDi_corrida_asa CDi_corrida_emp CDi_polar cdp CDp cdp_1 cdp_2 cdp_m1 cdp_t CDplane cdpr cdv CDv cdv_1 cdv_2 cdv_m1 cdv_t cdvr chord cordeh cl CL0w CL_corrida CL_corrida2 CL_corrida_asa CL_corrida_emp CL_polar CLa Clb Clda Cldl CLeh CLmax CLmaxasa CLmaxperf CLmaxperfm CLmaxperfp CLmaxperfr CLmaxrun Clp CLplane Clr cm Cma Cmde Cmq Cnb Cnda Cndl Cnp Cnr cr creh crev ct cteh ctev CYda CYdl D d_fixa de de_lim declive deda died Drag E emp epslon f Fat fator_pista fh file_aviao file_aviao_ground filename flag fmult gamalim grad_de_subida grad_sub grad_sublim H h_asa h_emp hlim i icl ih ime inc inc_massa inc_massa_dois inc_vel inutil ipa iTracao iw j k Ka Kh Kl ks ksh L L_mais_B lambda lambda1 lambda2 lambdaeh lambdaev lambdalim lb mac maceh mi MTOWt nada nEV nT offset offset1 offset2 offset_mac offset_maceh offset_motor offseteh offsetev p_disp p_req Pcg Pcglim perfilh perfilm perfilp perfilr perfilv plane_aviao planeAVL posb posbeh prop q restr_geom result ro S Seh Sev SgnDup_aileron sum_Fx sum_Fy themp Thrust_factor total_pista tow tracao v V_dec v_pot vel vento Vh Vhlim Vs Vv Vvlim X x1 x1_4 x1eh x1eh_4 x_emp Xcg_cr Xcg_mac Xeh Xhinge_a Xhinge_eh Xmac Xmaceh y1 y1eh Yc Ymac Ymaceh chordeh planeAVL;
result = 0;
dbstop if error
save = [x(1) x(2) x(3) x(4) x(5) x(6) x(7) x(8) x(9) x(10) x(11)];
fprintf('_________________________________________________________________\n')

%% Vari�veis n�o declaradas
d_fixa_motor = 0.05;
offset_motor = 0.05;

% Ambiente
AD      = 1100;                                 % Altitude Densidade
ro      = 1.2259*exp(-284.2*AD/2394432);        % Densidade do Ar
declive = 0;                                    % Declive da pista (+ uphill, - downhill)

% Fatores de corre��o
CD_factor     = 1.0;
Thrust_factor = 1;
fmult         = 1;

%% P�rametros da otimiza��o - Asa

b    = 2;                 % Envergadura da asa (m)
iw   = x(11);             % �ngulo de incid�ncia na corrida da asa [�]
died = 0;                 % Diedro [�]
X    = (b/2)*x(2);        % Transi��o reto-trapezoidal
i_fixacao_m = 3;          % �ngulo de fixa��o do motor

% Cordas na envergadura
chord(1)    = x(3);     cr = chord(1);        % Corda da raiz
chord(2)    = x(3);     cm = chord(2);        % Corda intermedi�ria
chord(3)    = cr*x(4);  ct = chord(3);        % Corda da ponta
chord       = [cr cm ct];                     % Vetor corda na envergadura da asa

% Posi��es das cordas na envergadura 
posb(1)   = 0;                                % In�cio da asa
posb(2)   = X;                                % Se��o do meio
posb(3)   = b/2;                              % Final da asa
posb      = [posb(1) posb(2) posb(3)];        % Vetor de posi��es das cordas na envergadura da asa

% Offset das cordas
offset1   = 0;                      % Offset da se��o do meio da asa [m]
offset2   = x(5);                   % Offset da ponta da asa [m]
offset = [0 offset1 offset2];       % Vetor offset na asa

% Perfis da asa
perfilr = round(x(6));   % Perfil da raiz da asa
perfilm = round(x(7));   % Perfil do meio da asa
perfilp = round(x(8));   % Perfil da ponta da asa

h_asa       = b/10;    % Altura da asa em rela��o ao solo (intradorso da raiz ate o solo)
H_CG        = 0;       % Altura CG a partir da altura da asa
% align_perf  = 0;     % Tipo de alinhamento dos perfis da raiz e ponta (1 se for pelo extradorso, 0 se for pela linha m�dia)
align_asa   = 0;       % Tipo de alinhamento da asa (1 se for pelo BF, 0 se for pelo BA)
% thw         = 0.;    % Thickness do perfil

%% C�lculos iniciais da ASA

y1 =(0:0.01:1)*b/2;                                % Malha para calcular a �rea da asa 
[y1,c1,x1] = Corda(chord,posb,X,b,offset);         % Fun��o para calcular a distribui��o de corda e offset na envergadura da asa
if isnan(c1)                                       % Teste da malha de cordas  
    disp('ops')
    PrintParada(save,1,NaN,NaN,NaN,NaN,NaN,NaN,NaN);
    return
end

S          = 2*trapz(y1,c1);              % C�lculo da �rea da asa usando o m�todo do trap�zio
AR         = 2*b^2/S;                     % Raz�o de aspecto
mac        = 2/S*trapz(y1,c1.^2);         % C�lculo da MAC
Ymac       = 2/S*trapz(y1,c1.*y1);        % Posi��o na envergadura da MAC
x1_4       = x1 + (c1/4);                 % Posi��o do CA de cada corda
Xmac       = 2/S*trapz(y1,c1.*x1_4);      % C�lculo do offset do CA da MAC
offset_mac = Xmac - mac/4;                % Offset da MAC
lambda1    = cm/cr;                       % Afilamento da ra�z da asa
lambda2    = ct/cm;                       % Afilamento da ponta da asa
lambda     = ct/cr;                       % Afilmento total da asa

%% Par�metros de Otimiza��o - EH 

beh       = x(10)*b;                    % Envergadura da EH
creh      = x(9);                       % Corda da raiz da EH
cteh      = 0.18;                       % Corda da ponta da EH
chordeh   = [creh cteh];                % Vetor corda EH
posbeh    = [0 beh];                    % Vetor posi��o das cordas na EH
offseteh  = 0;                          % Offset da ponta da EH
x_emp     = x(1);                       % Dist�ncia entre bordo de ataque da asa e bordo de ataque da EH (raiz)
h_emp     = tand(10)*(x_emp-cr);        % Altura entra asa e EH (Bordo de ataque da asa at� bordo de ataque da EH)
ih        = 0;                          % Incid�ncia da EH
perfilh   = 'NACA23012-24pt.dat';       % Perfil da EH    
Xeh       = 0;                          % Transi��o

%% C�lculos iniciais da EH

y1eh = (0:0.01:1)*beh/2;                                                                % Malha para calcular a �rea da eh
offseteh_so_para_funcao_de_corda = [0 offseteh];
[y1eh,c1eh,x1eh] = Corda(chordeh,posbeh,Xeh,beh,offseteh_so_para_funcao_de_corda);      % Fun��o para calcular a distribui��o de corda e offset na envergadura da EH
if isnan(c1eh)                                                                          % Teste da Malha de Cordas                  
     disp('eita')
     PrintParada(save,2,NaN,S,AR,NaN,NaN,NaN,NaN);
     return
end

x1eh_4 = x1eh + (c1eh/4);        % Posi��o do CA de cada corda

if creh == cteh                % EH retangular   
    Seh    = beh*creh;         % �rea da EH
    maceh  = creh;             % C�lculo da MAC para EH
    Ymaceh = 0;                % Posi��o na envergadura da MAC para EH
    Xmaceh = 0.25*creh;        % C�lculo do offset do CA da MAC para EH
else                           % EH n�o-retangular        
    Seh          = 2*trapz(y1eh,c1eh);
    maceh        = 2/Seh*trapz(y1eh,c1eh.^2);
    Ymaceh       = 2/Seh*trapz(y1eh,c1eh.*y1eh);
    Xmaceh       = 2/Seh*trapz(y1eh,c1eh.*x1eh_4);
end

lb           = - Xmac + x_emp + Xmaceh;        % Dist�ncia entre CAs da asa e das empenagens
ARh          = beh^2/Seh;                      % Raz�o de aspecto da EH
lambdaeh     = cteh/creh;                      % Afilamento da EH
offset_maceh = Xmaceh - maceh/4;               % Offset da MAC

%% Par�metros de Otimiza��o - EV

bev      = 0.32;                % Envergadura na EV
crev     = 0.2;                 % Corda da raiz na EV
ctev     = 0.2;                 % Corda da ponta na EV     
offsetev = 0;                   % Offset da EV
nEV      = 1;                   % Quantidade de EVs
perfilv  = 'NACA0012.dat';      % Perfil da EV  

if crev > creh
    lb = lb - (crev - creh);
end

%% C�lculos iniciais da EV

Sev      = crev*bev;        % �rea da EV
ARev     = bev^2/Sev;       % Raz�o de aspecto da EV
lambdaev = ctev/crev;       % Afilamento da EV

%% Superf�cies de Comando
% Aileron da asa
ail = 0.15*b;                     % Dist�ncia na semi-envergadura onde se inicia o aileron; � pra ser um absurdo propositalmente, pq sabemos que o AVL n�o sabe ver deflex�o

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
SgnDup_aileron = -1;                % Sinal de deflex�o para DUPLICATE do AVL (-1 para convencional, -2 para diferencial)
%ail_sup        = 0;                 % Aileron na asa superior (0=N�o; 1= Sim)

% Tab EH
bah            = beh;               % Envergadura do tab da EH
Xhinge_eh      = 0.20;              % [%] da corda da hinge da EH
themp          = 0.12;              % Thickness da EH

%% Sele��o de perfil e CLmax de cada um 

switch perfilr
    case 1
        perfilr = 'GID19111.dat';
        CLmaxperf = 1.9130;
        thicknessr = 0.1352;
    case 2
        perfilr = 'GID18091.dat';
        CLmaxperf = 1.7963;
        thicknessr = 0.1352;
    case 3
        perfilr = 'G01.dat';
        CLmaxperf = 2.5296;
        thicknessr = 0.14;
    case 4
        perfilr = 'GM07.dat';
        CLmaxperf = 2.5166;
        thicknessr = 0.1419;
    case 5
        perfilr = 'GM39.dat';
        CLmaxperf = 2.4280;
        thicknessr = 0.1431;
    case 6
        perfilr = 'RR100.dat';
        CLmaxperf = 2.4054;
        thicknessr = 0.1360;
    case 7
        perfilr = 'AC03.dat';
        CLmaxperf = 2.3238;
        thicknessr = 0.1343;
    case 8
        perfilr = 'AC04.dat';
        CLmaxperf = 2.2783;
        thicknessr = 0.1320;
    case 9
        perfilr = 'AC05.dat';
        CLmaxperf = 2.1959;
        thicknessr = 0.1240;
end
switch perfilm
    case 1
        perfilr = 'GID19111.dat';
        CLmaxperfm = 1.9130;
    case 2
        perfilr = 'GID18091.dat';
        CLmaxperfm = 1.7963;
    case 3
        perfilr = 'G01.dat';
        CLmaxperfm = 2.5296;
    case 4
        perfilr = 'GM07.dat';
        CLmaxperfm = 2.5166;
    case 5
        perfilr = 'GM39.dat';
        CLmaxperfm = 2.4280;
    case 6
        perfilr = 'RR100.dat';
        CLmaxperfm = 2.4054;
    case 7
        perfilr = 'AC03.dat';
        CLmaxperfm = 2.3238;
    case 8
        perfilr = 'AC04.dat';
        CLmaxperfm = 2.2783;
    case 9
        perfilr = 'AC05.dat';
        CLmaxperfm = 2.1959;  
end
switch perfilp
    case 1
        perfilr = 'GID19111.dat';
        CLmaxperfp = 1.9130;
    case 2
        perfilr = 'GID18091.dat';
        CLmaxperfp = 1.7963;
    case 3
        perfilr = 'G01.dat';
        CLmaxperfp = 2.5296;
    case 4
        perfilr = 'GM07.dat';
        CLmaxperfp = 2.5166;
    case 5
        perfilr = 'GM39.dat';
        CLmaxperfp = 2.4280;
    case 6
        perfilr = 'RR100.dat';
        CLmaxperfp = 2.4054;
    case 7
        perfilr = 'AC03.dat';
        CLmaxperfp = 2.3238;
    case 8
        perfilr = 'AC04.dat';
        CLmaxperfp = 2.2783;
    case 9
        perfilr = 'AC05.dat';
        CLmaxperfp = 2.1959;  
end

%% Restri��es da otimiza��o

ARlim       = 4;                    % Raz�o de aspecto limite
ARhlim      = 2.0;                  % Raz�o de aspecto limite da EH 
ARvlim      = 1;                    % Raz�o de aspecto limite da EV
grad_sublim = 3.5;                  % Gradiente de subida m�nimo (%)
Pcglim      = [0.30 0.42];          % Limites do CG [% da corda da raiz]
de_lim      = [4.6 20];             % Limites de deflex�o do profundor
Vhlim       = [0.185 0.52];         % Limites de volume de cauda horizontal
Vvlim       = [0.025 0.05];         % Limites de volume de cauda vertical

%% Restri��es geom�tricas de otimiza��o
% if posb(2) >= b/2 
%     disp('Transi��o > envergadura')
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
%     disp('Varia��o de corda muito grande')
%     return
% end

% if abs(b/2 - posb(2)) <= 0.10 || abs(posb(2)) <= 0.20
%     disp('Se��es muito pequenas')
%     return
% end
if abs(posb(3)-posb(2)) <= 0.1 ||  abs(posb(2)-posb(1)) <= 0.1 || abs(ail-posb(2)) < 0.05
	result = 0;
	disp('Se��o muito pequena')
    PrintParada(save,3,NaN,S,AR,NaN,NaN,NaN,NaN);
	return
end

if abs(b/2 - posb(2)) > 0.10 && abs(posb(2)) > 0.20
    result = result-0.1;
%         fprintf('%.4f\n',result)
end


if died < 0
    if h_asa - abs((b/2)*tand(gama)) < 0.045
    disp('Ponta de asa muito perto do ch�o')
    end
    if h_asa - abs((b/2)*tand(gama)) >= 0.045
        result = result-0.1;
%         fprintf('%.4f\n',result)
    end
end

if AR < ARlim            
    disp('Raz�o de aspecto baixa')
    PrintParada(save,4,NaN,S,AR,NaN,NaN,NaN,NaN);
    return
end
if AR >= ARlim 
        result = result-0.1;
%         fprintf('%.4f\n',result)
end

if ARh < ARhlim             
    disp('Raz�o de aspecto da EH')
    PrintParada(save,5,NaN,S,AR,NaN,NaN,NaN,NaN);
    return
end
if ARh >= ARhlim  
    result = result-0.1;
%         fprintf('%.4f\n',result)
end

% if ARev < ARvlim             
%     disp('Raz�o de aspecto da EV')
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
%     disp('Afilamento entre se��es muito grande')
%     PrintParada(save,9,NaN,S,AR,NaN,NaN,NaN,NaN);
%     return
% end
% if (lambda1 && lambda2) >= 0.5
%     result = result-0.1;
%      %   fprintf('%.4f\n',result)
% end

% if S < 0.80
%    result = 0;
%    disp('�rea insatisfat�ria')
%    return
% end

if cr*thicknessr > 0.1
    result = result-0.1;
    disp('Nervura com altura superior a 10cm')
    PrintParada(save,2,NaN,S,AR,NaN,NaN,NaN,NaN);
    return
end

% Altura efetiva entre asa e EH
hh = h_emp + thicknessr*chord(1); 

if atand((h_asa + hh)/(x_emp - cr+creh)) < 18-iw        % Pr�-verifica��o de tail-strike
    result = result-0.1;
    disp('Poss�vel tail-strike') 
    PrintParada(save,18,NaN,S,AR,NaN,NaN,NaN,NaN);
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
     fprintf('%.4f\n',result)
end

if Vv < Vvlim(1) || Vv > Vvlim(2)
     disp('Volume de cauda vertical')
     PrintParada(save,9,peso,NaN,S,AR,NaN,NaN,NaN);
     return
else        % Vv >= Vvlim(1) && Vv <= Vvlim(2)
    result = result-0.1;
    fprintf('%.4f\n',result)
end

%% Cria��o das partes do avi�o em .avl e c�lculo anal�tico das derivadas longitudinais

alpha1 = -3;        % menor alpha a ser analisado
alpha2 = 12;        % maior alpha a ser analisado

[file_asa,file_eh] = convencionalAVLpartes(S,mac,b,posb,offset1,offset2,died,ail,chord,perfilr,perfilm,perfilp,perfilh,beh,offseteh,creh,cteh,iw,ih,dist_CA,Seh,maceh);
[CL,CD,Cm]         = coefsAVL(alpha1,alpha2,0,file_asa,file_eh);

ang = linspace(alpha1,alpha2,10);       % Range de �ngulos analisados
n   = 0.97;                             % Queda de press�o din�mica assumida

 for Pcg = Pcglim(2):-0.005:Pcglim(1)       % C�lculo da posi��o de CG a partir dos limites dados

    Xcg     = Pcg*cr;                          % Posi��o do CG a partir do BA [m]
    Xcg_mac = (Xcg - offset_mac)/mac;          % Posi��o do CG em porcentagem da mac
    Vh_cg   = Vh-(Seh/S)*(Xcg_mac-0.25);       % Volume de cauda em rela��o ao CG (Etkin)

    % Momentos principais
    Cm_M_asa = Cm.w(ang);                       % Cm da asa
    Cm_L_asa = CL.w(ang)*(Xcg_mac-0.25);        % Momento do CL da asa
    Cm_D_asa = -CD.w(ang)*(H_CG/mac);           % Momento do CD da asa
    
    Cm_M_eh = Cm.h(ang)*(Seh/S)*(maceh/mac);                % Cm da EH
    Cm_L_eh = -CL.h(ang)*n*Vh_cg;                           % Momento do CL da EH
    Cm_D_eh = CD.h(ang)*n*(Seh/S)*(h_emp-H_CG)/maceh;        % Momento do CD da EH
       
    Cm_D_ev = nEV*0.012796*n*(Sev/S)*((h_emp+bev/2)-H_CG)/macev;        % Momento do CD da EV

    % Coeficientes totais
    CLtot     = CL.w(ang)+n*(Seh/S)*CL.h(ang);                              % CL total do avi�o
    CDtot     = CD.w(ang)+n*(Seh/S)*CD.h(ang)+nEV*0.012796*n*(Sev/S);       % CD total do avi�o
    Cmtot_asa = Cm_M_asa+Cm_L_asa+Cm_D_asa;                                 % Cm total da asa
    Cmtot_eh  = Cm_M_eh+Cm_L_eh+Cm_D_eh;                                    % Cm total da EH
    Cmtot     = Cmtot_asa+Cmtot_eh+Cm_D_ev;                                 % Cm total do avi�o
    
    % Par�metros principais
    H      = (-diff(Cmtot))./diff(CLtot)*100;       % Margem est�tica em porcentagem
    h_calc = H(~isnan(H));
    Xn     = Xcg_mac+h_calc/100;                    % Ponto neutro
    Cma    = diff(Cmtot);                           % Derivada de Cm por alpha
    
    Cmtotal_asa = @(alpha)interp1(ang,Cmtot_asa,alpha);        % Cm total da asa interpolado
    Cmtotal_eh  = @(alpha)interp1(ang,Cmtot_eh,alpha);         % Cm total da EH interpolado
    
    % Verificar se o alfa de trimagem � positivo para alguma incid�ncia de EH vi�vel:
    for ih = 0:-1:-3
        Cm0 = Cmtotal_asa(-iw) + Cmtotal_eh(-iw+ih)+ Cm_D_ev;
        if Cm0 > 0
            break
        end
    end
    
    % Verificar se a margem est�tica � positiva para todos os �ngulos e se Cm0 � positivo, para o loop
    if all(H > 0) && Cm0 > 0
        break
    end
 end
 
%% Verifica��es de Estabilidade Longitudinal

if any(H < 0)
        fprintf('Margem est�tica = %.4f\n',H(length(H)))
        PrintParada(save,7,peso,Pcg,S,AR,NaN,NaN,NaN);
        return
elseif all(H >= 0)
    result = result-0.1;
    fprintf('%.4f\n',result)
end

if Cm0 <= 0
    fprintf('Cm0 n�o positivo\n')
    PrintParada(save,6,peso,Pcg,S,AR,NaN,NaN,NaN);
    return
else       % Cm0 > 0
    result = result-0.1;
end

if any(Cma > 0)
    disp('Longitudinalmente inst�vel')
    PrintParada(save,8,peso,Pcg,S,AR,NaN,NaN,NaN);
    return
elseif all(Cma <= 0)
    result = result-0.1;
end

%% Cria��o do avi�o inteiro em .avl e c�lculo das derivadas l�tero-direcionais
[file_aviao] = convencionalAVL(S, mac, b, posb, offset, died, ail, c_ail, Pcg, chord, Xhinge_eh, Xhinge_a, SgnDup_aileron, perfilr, perfilm, perfilp, perfilh, perfilv, beh, offseteh, creh, cteh, bev, crev, ctev, offsetev, x_eh, h_eh, x_ev, h_ev, iw, ih, nEV, H_CG, ef_solo, h_asa);
[~,Clb,~,Cnb,~,~,~,~,~,~,~,~,~,~,~,~,~,~] = derivadasAVL(4,4,mac,file_aviao);

%% Verifica��es de Estabilidade L�tero-Direcional

if Cnb < 0 
    disp('L�tero-direcionalmente inst�vel')
    delete(file_aviao)
    PrintParada(save,9,peso,Pcg,S,AR,NaN,NaN,NaN);
    return
else       % Cnb >=0
    result = result - 0.1;
end

if Clb > 0 
    disp('L�tero-direcionalmente inst�vel')
    delete(file_aviao)
    PrintParada(save,9,peso,Pcg,S,AR,NaN,NaN,NaN);
    return
else       % Clb <= 0 
    result = result-0.1;
end

%% File AVL (efeito solo)
[file_aviao_ground] = convencionalAVL(S,mac,b,posb,offset1,offset2,died,ail,c_ail,Pcg,chord,Xhinge_eh,Xhinge_a,SgnDup_aileron,perfilr,perfilm,perfilp,perfilh,perfilv,beh,offseteh,creh,cteh,bev,crev,ctev,offsetev,x_emp,h_emp,iw,ih,nEV,H_CG,1,h_asa);
%fprintf('Tempo depois do AVL do efeito solo:%f\n',toc)

%% Efeito solo
% Par�metros em cada posi��o da asa
[Yc,FsCL0,Area,~,~,cl,CDi_corrida,CL_corrida] = FSFT_dadosgroundAVL(-iw,0,ro,12,file_aviao_ground);
de = 0;

cdv = NaN(1,length(Area));
cdp = NaN(1,length(Area));

for j=1:length(Area)
    % Adi��o da viscosidade da primeira se��o
    if abs(Yc(j)) <= posb(2)        
        [cdvr, cdpr] = ViscosidadePolar(perfilr, cl(j));
        cdv(j) = cdvr;
        cdp(j) = cdpr;
    % Adi��o da viscosidade da 2� se��o
    else
        prop = (posb(3) - abs(Yc(j)))/(posb(3)-posb(2));
        [cdv_1, cdp_1] = ViscosidadePolar(perfilr, cl(j));
        [cdv_2, cdp_2] = ViscosidadePolar(perfilp, cl(j));   

        cdv(j) = prop*cdv_1 + (1 - prop)*cdv_2;
        cdp(j) = prop*cdp_1 + (1 - prop)*cdp_2;
    end
end

% Arrasto na corrida de pista
CDvg = 4*sum(reshape(cdv,size(Area)).*Area)/S;              % Viscoso calculado atrav�s da an�lise de perfil no XFLR5
CDpg = 4*sum(reshape(cdp,size(Area)).*Area)/S;              % Parasita calculado atrav�s da an�lise de perfil no XFLR5, dobrado por serem duas asas
Cdg = CDvg+CDpg;
CD_corrida = CD_factor*(CDi_corrida + CDvg);                 % Total

if CL_corrida < 0.4 
    disp('CL corrida muito baixo')
    delete(file_aviao)
    delete(file_aviao_ground)
    delete(FsCL0)
    PrintParada(save,10,peso,Pcg,S,AR,CD_corrida,CL_corrida,NaN);
    return
end

if CL_corrida >= 0.4 
    result = result-0.1;
end

%% Efeito solo anal�tico
% Asa
f = 2*(h_asa)/b;
ks = 46.87*f/(f^2+40.12*f+10.86);     % Curva zuadas do livro
ARs = AR/ks;                          % AR corrigida devido ao efeito solo

% EH
fh = 2*(h_asa+h_emp)/b;                 % C�lculo do elemento f da formula descrita em ROSKAM, levando em conta a altura da fuselagem, a distancia da fuselagem ao solo e a distancia entre fuselagem e EH 
ksh = 46.87*fh/(fh^2+40.12*fh+10.86);   % Curva do livro
ARsh = ARh/ksh;                         % Raz�o de aspecto modificada. Uma EH com ARs em cruzeiro � equivalente a uma com AR sofrendo o efeito solo.

% Derivadas aerodin�micas (Clalfa)
% Asa
aws = 2*pi*ARs/(2+sqrt(ARs^2/0.98^2+4));    % Coeficiente angular de Cl(alpha) da asa, obtida de Raymer (eq. 12.6), considerando o efeito solo

% EH
ats = 2*pi*ARsh/(2+sqrt(ARsh^2/0.98^2+4));  % Coeficiente angular de Cl(alpha) da EH, obtida de Raymer (eq. 12.6), considerando o efeito solo

% CL 0 da asa
[CL0w,~,~,~] = CL0asaAVL(FsCL0);

% Downwash
Ka = 1/AR-1/(1+AR^1.7);             % Fator de corre��o da raz�o de aspecto
Kl = (10-3*lambda)/7;               % Fator de corre�ao do afilamento
Kh = (1-h_emp/b)/((2*lb/b)^(1/3));  % Fator de corre��o EH
deda = 4.44*(Ka*Kl*Kh)^(1.19);      % F�rmula do gradiente de downwash
epslon = deda*(deg2rad(iw) + CL0w/aws); 

% CL corrida anal�o da asa
CL_corrida_asa = CL0w + aws*deg2rad(iw);
CL_corrida_emp = 0 + ats*(deg2rad(ih) - epslon);
CL_corrida2 = CL_corrida_asa + 0.9*(Seh/S)*CL_corrida_emp;

CDi_corrida_asa = CL_corrida_asa^2/(pi*ARs*0.98);
CDi_corrida_emp = CL_corrida_emp^2/(pi*ARsh*0.98);
CDi_corrida2 = CDi_corrida_asa + 0.95*(Seh/S)*CDi_corrida_emp;

fprintf('\nCL corrida avl = %.4f\n',CL_corrida)
fprintf('CL corrida analio = %.4f\n',CL_corrida2)
fprintf('CD corrida = %.4f\n',CD_corrida)
fprintf('CDi corrida avl = %.4f\n',CDi_corrida)
fprintf('CDi corrida anal�tico = %.4f\n',CDi_corrida2)
% fprintf('CDp anal�tico = %.4f\n',CDp)
fprintf('k = %.4f\n', k_prandtl)
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
    if exist(ID1,'file') == 2      % Evita arquivos j� existentes
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
                [~,j] = max(cl);
                if Yc(j) >= 0.6*b/2
                    disp('Estol no aileron');
                    delete(file_aviao)
                    delete(file_aviao_ground)
                    delete(FsCL0)
                    PrintParada(save,11,peso,Pcg,S,AR,CD_corrida,CL_corrida,NaN);
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

%% Verifica��o de CLmax e estol sem efeito solo
CLmaxasa = CLmaxrun;
alfamax = alfa;

% Intervalo de alpha
alpha = -2-iw:0.75:alfamax;

if alfamax <= iw
    result = 0;
    PrintParada(save,12,peso,Pcg,S,AR,CD_corrida,CL_corrida,NaN);
    disp('Incid�ncia estola')
    return
end

fprintf('CL max asa = %.4f\n',CLmaxasa)
fprintf('Alfamax max = %.4f\n',alfamax+iw)
fprintf('Rotate = %.4f\n',atand((h_asa + h_emp)/(x_emp+creh)))
fprintf('Incidencia da asa = %.4f\n',iw)
fprintf('Incidencia da EH = %.4f\n',ih)

if atand((h_asa + h_emp)/(x_emp+creh)) < alfamax-1
    disp('Sem rotate')
    delete(file_aviao)
    delete(file_aviao_ground)
    delete(FsCL0)
    PrintParada(save,13,peso,Pcg,S,AR,CD_corrida,CL_corrida,NaN);
    return
end
if atand((h_asa + h_emp)/(x_emp+creh)) >= alfamax-1
    result = result-0.1;
end

%% Polar de arrasto biplano
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
    %% For�as por Strip, Arrasto Induzido Polar e Sustenta��o para o alpha dado
    ID2 = strcat('Fsurface_polar_',int2str(k),'.dat');
    fileID2 = fopen(ID2,'r');
    clear data
    data = textscan(fileID2,'CLsurf  =   %f     CDsurf  =   %f','Headerlines',15);
    CDi_asa = data{2};
    
    %% For�as por Strip
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
    
    %% Sustenta��o para o alpha dado
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
        % Adi��o da viscosidade da primeira se��o
        if abs(Yc(j)) <= posb(2)
            [cdvr, cdpr] = ViscosidadePolar(perfilr, cl(j));

            cdv(j) = cdvr;
            cdp(j) = cdpr;
            
        % Adi��o da viscosidade da segunda se��o    
        else
            prop = (posb(3) - abs(Yc(j)))/(posb(3)-posb(2));
            [cdv_m1, cdp_m1] = ViscosidadePolar(perfilr, cl(j));
            [cdv_t, cdp_t] = ViscosidadePolar(perfilp, cl(j));   

            cdv(j) = prop*cdv_m1 + (1 - prop)*cdv_t;
            cdp(j) = prop*cdp_m1 + (1 - prop)*cdp_t;    
        end
    end
    
    CDv(i) = 2*sum(reshape(cdv,size(Areac)).*Areac)/S;           % Viscoso calculado atrav�s da an�lise de perfil no XFLR5
    CDp(i) = 2*sum(reshape(cdp,size(Areac)).*Areac)/S;           % Parasita calculado atrav�s da an�lise de perfil no XFLR5
    Cd(i) = CDv(i) + CDp(i);
    CD_polar_ext(i) = CD_factor*(CDi_polar_ext(i) + CDv(i));             % Total

    % CL do �ngulo analisado
    CL_polar_ext(i) = 1.0*CL_polar_ext(i);
    
    
    % CL_polar � o cl max DO AVI�O INTEIRO (asa+eh+ev+fuselagem...) 
    % com valor que sai do avl, mas limitado pelo de estol dado pelo xf 
    % (definido no for pelo length(alpha-1), sa�do de arquivo Ft do avl
    
    i = i + 1;
end

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%% An�lise de desempenho %%%%%%%%%%%%%%%%%%%%%% 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% Perfil Aerodin�mico da Asa < !!!
alpha1 = 20.10;          % �ngulo de estol da asa EM efeito solo [�]
alpha2 = 24.90;          % �ngulo de estol da asa SEM efeito solo [�] 
a1 = 'PolarASA_ef.txt';  % Nome do arquivo com os dados originais da asa EFEITO SOLO
a2 = 'PolarASA_TA.txt';  % Nome do arquivo com dados da asa FORA DE EFEITO SOLO

% Perfil Aerodin�mico da EH < !!!
alpha3 = 9.35;                 % �ngulo de estol da EH EM EFEITO SOLO [�]
alpha4 = 11.0;                 % �ngulo de estol da EH SEM EFEITO SOLO [�]
ehs = 'PolarEH1_ef_TA.txt';    % Nome do arquivo com os dados da EH - EM EFEITO SOLO
eh = 'PolarEH1_TA.txt';        % Nome do arquivo com os dados da EH - SEM EFEITO SOLO

% Perfil aerodin�mico da EV
ev = 'PolarEV.txt';     % Nome do arquivo com os dados da EV

% Par�metros Aerodin�micos 
alphar_ws = -5:0.05:alpha1;      % malha aceit�vel de �ngulo de ataque em solo - asa
alphar_w = -5:0.05:alpha2;       % malha aceit�vel de �ngulo de ataque - asa
alphar_ehs = -15:0.1:alpha3; 
alphar_eh = -15:0.1:alpha4;

% Asa < !!!
ASA_solo = ler_polar_VLMeLLT(a1);  % dados da asa extra�dos do XFLR5
ASA = ler_polar_VLMeLLT(a2);

CL0asa = 0.758012; %%%%%N�S N�O USAMOS ESSE PARAMETRO EM NENHUM LUGAR, POSSO RETIRAR????
CLstall_s = spline(ASA_solo(:,1),ASA_solo(:,2),alpha1);
CDstall_s = spline(ASA_solo(:,1),ASA_solo(:,3),alpha1);
CLstall = spline(ASA(:,1),ASA(:,2),alpha2);
CDstall = spline(ASA(:,1),ASA(:,3),alpha2);%%%%%N�S N�O USAMOS ESSE PARAMETRO EM NENHUM LUGAR, POSSO RETIRAR????
CDpasa = spline(ASA(:,1),ASA(:,4),i_fixacao_w);%%%%%N�S N�O USAMOS ESSE PARAMETRO EM NENHUM LUGAR, POSSO RETIRAR????

% CL0asa = CL para alpha igual a zero EM EFEITO SOLO
% Cmasa = Cm no centro aerodin�mico para angulo 0 EM EFEITO SOLO
% CLstall_s = CL da asa para estol EM EFEITO SOLO, de decolagem usado para v_stall2
% CDstall_s = CD da asa para estol EM EFEITO SOLO
% CLstall = CL de estol FORA DO EFEITO SOLO
% CDstall = CD de estol FORA DO EFEITO SOLO
% CDpasa = Arrasto parasita da asa para �ngulo de fixa��o

Cl_alpha_asa_solo = spline(ASA_solo(:,1),ASA_solo(:,2),alphar_ws);
%Cls_alpha_asa = Cl_alpha_ASA_solo;

Cd_alpha_asa_solo = spline(ASA_solo(:,1),ASA_solo(:,3),alphar_ws);
%Cds_alpha_asa = Cd_alpha_ASA_solo;

Cl_alpha_asa = spline(ASA(:,1),ASA(:,2),alphar_w);

Cd_alpha_asa = spline(ASA(:,1),ASA(:,3),alphar_w);


% EH < !!!
EH = ler_polar_VLMeLLT(eh);                % dados da eh extra�dos do XFLR5
EH_SOLO = ler_polar_VLMeLLT(ehs);

%NL= length(EH_SOLO(:,5));
%SI=sum(EH_SOLO(:,5));
CDpehs = spline(EH_SOLO(:,1),EH_SOLO(:,4),alpha3); % pra 11.5 graus
CDpeh = spline(EH(:,1),EH(:,4),alpha4); % pra 11.5 graus
CL0eh = spline(EH_SOLO(:,1),EH_SOLO(:,2),-i_fixacao_eh); %%%%%N�S N�O USAMOS ESSE PARAMETRO EM NENHUM LUGAR, POSSO RETIRAR????
Cmeh =  spline(EH_SOLO(:,1),EH_SOLO(:,5),9); % valor m�dio

% CDpehs = Arrasto parasita da EH em efeito solo (stall)
% CDpeh = Arrasto parasita da EH fora de efeito solo
% CL0eh = CL da EH quando �ngulo de ataque � igual a zero em efeito solo
% Cmeh = Coeficiente de momento da eh em efeito solo

Cls_alpha_eh = spline(EH_SOLO(:,1),EH_SOLO(:,2),alphar_ehs);  % Negativo para inverter o perfil, mas aqui o perfil � sim�trico
Cl_alpha_eh = spline(EH(:,1),EH(:,2),alphar_eh);              % Negativo para inverter o perfil, mas aqui o perfil � sim�trico
Cds_alpha_eh = spline(EH_SOLO(:,1),EH_SOLO(:,3),alphar_ehs);
Cd_alpha_eh = spline(EH(:,1),EH(:,3),alphar_eh);%%%%%N�S N�O USAMOS ESSE PARAMETRO EM NENHUM LUGAR, POSSO RETIRAR????

% EV
EV = ler_polar_VLMeLLT_ev(ev);                                 % dados da ev extra�dos do XFLR5
Cl_alpha_ev = spline(EV(:,1),EV(:,2),alphar_w);%%%%%N�S N�O USAMOS ESSE PARAMETRO EM NENHUM LUGAR, POSSO RETIRAR????
Cd_alpha_ev = spline(EV(:,1),EV(:,3),alphar_w);%%%%%N�S USAMOS ESSE PARAMETRO APENAS SE TRATANDO DO FREIO, POSSO RETIRAR????

CDev = 0.012586;       % Coeficiente de arrasto da EV
CLmax = max(CL_polar_ext);
% monoCLmax = max(monoCL_polar);
% CLCDmax = max(CL_polar./CD_polar);

%% Corrida de Pista - Rotate e MTOW
lim_tow = [10 20];  % Limites para an�lise de tow no programa de rotate
inc_tow = 0.25;     % Incremento de tow
MTOW    = 0;        % Flag que indica problema no rotate (se n�o entra no if, n�o h� rotate ao final)
j       = 1;
for tow = lim_tow(1):inc_tow:lim_tow(2)
    aux        = 1;
    i          = 1;
    P          = tow*g; % Peso do avi�o
    v(i,j)     = 0;     % Velocidade inicial
    alfa(i,j)  = iw;    % Incid�ncia da asa
    t_tot(i,j) = 0;     % Tempo total da corrida de pista, para cada itera��o
    d_tot(i,j) = 0;     % Dist�ncia total da corrida de pista, para cada itera��o
    omega(i,j) = 0;     % Velocidade angular de rotate em cada itera��o
    sumFy(i,j) = -1;    % Somat�ria de for�as em Y em cada itera��o
    
    while d_tot < D_max_rot-inc_pista || alfa(i,j) < stall  % sumFy negativo indica que h� rea��o normal nas rodas e ainda n�o decolou
        epslon(i,j) = alfa(i,j)*deda;                % Downwash
        ih_inst = alfa(i,j) - iw + ih - epslon(i,j); % �ngulo inst�ntaneo da EH
        im_inst = alfa(i,j) - iw + im;               % �ngulo inst�ntaneo do motor
        
        Lw(i,j) = 0.5*S*ro*v(i,j)^2*CL.w(alfa(i,j));         % Sustenta��o na asa
        Dw(i,j) = 0.5*S*ro*v(i,j)^2*CD.w(alfa(i,j));    % Arrasto na asa
        
        if d_tot(i,j)>(D_max-10)        % Quando faltar 10m para dist�ncia m�xima de decolagem, aciona-se o profundor)
            ih_inst = ih_inst + de*tau; % �ngulo inst�ntaneo com deflex�o de profundor          
        end
        Leh(i,j) = n*0.5*Seh*ro*v(i,j)^2*CL.h(ih_inst);    % Sustenta��o na EH
        Deh(i,j) = n*0.5*Seh*ro*v(i,j)^2*CD.h(ih_inst);	% Arrasto da EH
        
        Dev(i,j) = n*0.5*Sev*ro*v(i,j)^2*CD_ev; % Arrasto da EV
        
        Dfus(i,j) = 0.5*S*ro*v(i,j)^2*CD_fus;                               % Arrasto da fuselagem
        Tx(i,j) = T(v(i,j))*cosd(im_inst);                                 % Tra��o em X
        Ty(i,j) = T(v(i,j))*sind(im_inst);                                 % Tra��o em Y
        Fat(i,j) = abs(mi*(P-Ty(i,j)-Lw(i,j)-Leh(i,j)));                   % For�a de atrito
        sumFx(i,j) = Tx(i,j)-Fat(i,j)-Dw(i,j)-Deh(i,j)-Dev(i,j)-Dfus(i,j); % Somat�rio de for�as horizontais
        sumFy(i+1) = Ty(i,j)+Lw(i,j)+Leh(i,j)-P;                           % Somat�rio de for�as verticais

        if sumFy(i+1) > 0
            % fprintf('Decolou antes do stall\n')
            RN(i,j) = 0;                                                % Rea��o normal na bequilha
            RM(i,j) = P-T(v(i,j))*sind(im)-RN(i,j)-Leh(i,j)-Lw(i,j);    % Rea��o normal no TDP
            MTOW = tow;
            break
        end

        a(i,j) = sumFx(i,j)/tow;                                            % Acelera��o linear
        t(i,j) = (-v(i,j)+sqrt(v(i,j)^2+4*(a(i,j)/2)*inc_pista))/a(i,j);    % Tempo para percorrer um incremento
        v(i+1,j) = v(i,j)+a(i,j)*t(i,j);                                    % Velocidade inicial da pr�xima itera��o      
        
        MLw = Lw(i,j)*(Xcg-0.25)*mac;                       % Momento gerado pela sustenta��o da asa
        MDw = Dw(i,j)*(H_CG);                               % Momento gerado pelo arrasto da asa
        MMw = 0.5*S*ro*v(i,j)^2*Cmw(CLw(alfa(i,j)))*mac;    % Momento puro da asa de cima
            
        MLeh = -Leh(i,j)*braco;                             % Momento gerado pela sustenta��o da EH
        MDeh = Deh(i,j)*(h_emp-H_CG);                       % Momento gerado pelo arrasto da EH
        MMeh = 0.5*Seh*ro*v(i,j)^2*Cmh(CLh(ih_inst))*maceh; % Momento puro da EH
        
        % perguntar aerodinamica se h� perda de press�o dinamica na eh para o momento
        
        MDev = Dev(i,j)*(h_emp+bev/2-H_CG)*macev;   % Momento gerado pelo arrasto da EV
            
        Mmotx = -Tx(i,j)*(dist_CG_m)*sind(im_inst); % Momento gerado pela componente X da tra��o
        Mmoty = Ty(i,j)*(dist_CG_m)*cosd(im_inst);  % Momento gerado pela componente Y da tra��o

        Mt_asa = MLw+MDw+MMw;   % Momento total da asa
        Mt_eh = MLeh+MDeh+MMeh; % Momento total da EH
        Mt_ev = MDev;           % Momento total da EV
        Mt_mot = Mmotx+Mmoty;   % Momento total do motor
        
        RN(i,j) = ((P-T(v(i,j))*sind(im)-Leh(i,j)-Lw(i,j))*(d2+mi*H_CG_chao)-(Mt_asa+Mt_eh+Mt_ev+Mt_mot))/(d1+d2);  % Rea��o normal na bequilha (demonstra��o na pasta)
        RM(i,j) = P-T(v(i,j))*sind(im)-RN(i,j)-Leh(i,j)-Lw(i,j);                                                    % Rea��o normal no TDP 

        if RN(i,j)<=0   % RN negativo: bequilha n�o toca mais o ch�o e, por isso, RN = 0
            RN(i,j)=0;
            if RM(i,j) >= 0
                M_tdp = -RM(i,j)*(chord(1)*(Ptdp-Pcg)-(H_CG+alt_asa)*sind(alfa(i,j)-iw));
                M_at = -RM(i,j)*mi*((H_CG+alt_asa)*cosd(alfa(i,j)-iw)+(chord(1)*(Ptdp-Pcg))*sind(alfa(i,j)-iw));
            else
                M_tdp = 0;
                M_at = 0;
            end
            Rot(i,j) = Mt_asa + Mt_eh + Mt_ev + Mt_mot + M_tdp + M_at;	% Momento total de rota��o
            a_ang(i,j) = rad2deg(Rot(i,j)/Itotal);                      % Acelera��o angular
            dalfa(i+1,j) = omega(aux)*t(i,j)+a_ang(i,j)*t(i,j)^2/2;     % Varia��o de �ngulo
            omega(aux+1) = omega(aux)+a_ang(i,j)*t(i,j);                % Velocidade angular inicial da pr�xima itera��o
            alfa(i+1,j) = alfa(i,j)+dalfa(i+1,j);                       % Nova incid�ncia da asa;
            aux = aux+1;
        else
            alfa(i+1,j) = iw;
        end

        d_tot(i+1,j) = d_tot(i,j)+inc_pista;    % Dist�ncia total percorrida
        t_tot(i+1,j) = t_tot(i,j)+t(i,j);       % Tempo total
        i = i+1;
    end
    
    if all(sumFy <= 0)
        break
    end
    j=j+1;
end

% Verifica��o
if MTOW == 0 
    fprintf('Problema no rotate!!!')
    delete(file_aviao)
    PrintParada(save,12,peso,Pcg,S,AR,NaN,NaN,NaN);
    return
else
    result = result-0.1;
end
%% Voo em cruzeiro

for i = 1:length(alphar_ws) % < !!!
    alpha(i) = alphar_w(i);
    gama(i) = alphar_w(i) - i_fixacao_w;
    beta(i) = gama(i) + i_fixacao_m;
    CLasa(i) = Cl_alpha_asa(i);  % dados extra�dos do txt do XFLR5
    CDasa(i) = Cd_alpha_asa(i);  % dados extra�dos do txt do XFLR5
    CLehinicial(i) = 0.6;  % chute inicial para o valor de CL da EH; 
    CDehinicial(i) = CDpeh + (CLehinicial(i)^2)/(pi*ARh*0.98);    % chute inicial para o valor do CD da EH;

    for  j = 1:100
        if j == 1
           CLehi(j) = CLehinicial(i);
           CDehi(j) = CDehinicial(i);
        end

        % Equa��o do equil�brio de for�as verticais da aeronave

        vc(j) = sqrt((2*tow*9.81)/(ro*((Sw*CDasa(i) + n*Sh*CDehi(j) + n_EV*n*Sv*CDev )*tand(beta(i)) + (CLasa(i)*Sw + n*CLehi(j)*Sh)))); % estimando um valor pra velocidade --> dedu��o dessa equa��o est� no papel (falar com o Pet, qualquer coisa)                   

        % Equa��o do equil�brio de momentos em torno do CG -> + sentido hor�rio

        T(j) = ((ro*(vc(j)^2))*(Sw*CDasa(i) + n*Sh*CDehi(j) + n_EV*n*Sv*CDev ))/(2*cosd(beta(i))); % tra��o requerida do motor
        m_motor(j) = -T(j)*cosd(beta(i))*(Hm*cosd(gama(i))+dm*sind(gama(i)))+T(j)*sind(beta(i))*(dm*cosd(gama(i))-Hm*sind(gama(i))); % momento gerado pelo motor em torno do CG

        Lw(j) = 0.5*ro*(vc(j)^2)*Sw*CLasa(i);
        Dw(j) = 0.5*ro*(vc(j)^2)*Sw*CDasa(i);
        m_asa(j) = Lw(j)*((Xcg - 0.25)*mac*cosd(gama(i)) - Hw*sind(gama(i))) + Dw(j)*((Xcg - 0.25)*mac*sind(gama(i)) + Hw*cosd(gama(i)));               

        cm_asa(j) = -0.5*ro*(vc(j)^2)*Sw*(-Cmasa)*mac;

        Deh(j) = n*0.5*ro*(vc(j)^2)*Sh*CDehi(j); 
        m_eh(j) = Deh(j)*(-sind(gama(i))*lb + cosd(gama(i))*Hh); 

        cm_eh(j) =  -0.5*ro*n*(vc(j)^2)*Sh*(-Cmeh)*maceh;     

        Dev(j) = n*0.5*ro*(vc(j)^2)*Sv*CDev; 
        m_ev(j) = n_EV*Dev(j)*(-sind(gama(i))*lb + cosd(gama(i))*Hh); 
        
        
        Novo_CLeh(j) = (2*(m_motor(j) + m_asa(j) + cm_asa(j) + m_eh(j) + cm_eh(j) + m_ev(j) ))/(n*ro*(vc(j)^2)*Sh*(lb*cosd(gama(i)) + Hh*sind(gama(i))));
        Novo_CDeh(j) = CDpeh + (Novo_CLeh(j)^2)/(pi*ARh*0.98);

        % Crit�rio de converg�ncia - checando valores de for�as e momentos

        erro1 = 0.00000000001;
        erro2 = 0.00000000001;
        erro3 = 0.00000000001;

        Nova_Leh(j) = n*0.5*ro*(vc(j)^2)*Sh*Novo_CLeh(j); 

        diff1(j) = abs((T(j)*sind(beta(i)) + Lw(j) + Nova_Leh(j) - tow*9.81));
        diff2(j) = abs(Novo_CLeh(j) - CLehi(j));
        diff3(j) = (m_motor(j) + m_asa(j) + cm_asa(j) + m_eh(j) + cm_eh(j) + m_ev(j) - Nova_Leh(j)*(lb*cosd(gama(i)) + Hh*sind(gama(i))));            


        if diff1(j) <= erro1 && diff2(j) <= erro2 && diff3(j) <= erro3
            CLeh_final(i) = Novo_CLeh(j);
            CDeh_final(i) = Novo_CDeh(j);
            v_final(i) = vc(j);
            diff1_final(i) = diff1(j);
            diff2_final(i) = diff2(j);
            diff3_final(i) = diff3(j);
            Tracao(i) = T(j);
            break
        end

        CLehi(j+1) = Novo_CLeh(j);
        CDehi(j+1) = Novo_CDeh(j);

    end

    % Faixa de �ngulos de atua��o da EH

    for count = 1:length(alphar_eh)       % Para encontrar o index do CLeh
        
        diff(i) = Cl_alpha_eh(count) - CLeh_final(i);

        if abs(diff(i)) <= 0.001
            index_eh = count;
            alpha_eh(i) = alphar_eh(index_eh);
            break
        end
    end
    

    % Pot�ncia Requerida
    L_total(i) = 0.5*ro*(v_final(i)^2)*(n*CLeh_final(i)*Sh+Sw*CLasa(i));
    D_total(i) = 0.5*ro*(v_final(i)^2)*(Sw*CDasa(i) + n*Sh*CDeh_final(i) + n_EV*n*Sv*CDev );
    p_req(i) = D_total(i)*v_final(i);

    % Pot�ncia Dispon�vel
    t_disp(i) = polyval(tracao, v_final(i));
    p_disp(i) = v_final(i)*t_disp(i);  % C�lculo da pot�ncia dispon�vel do motor, levando em conta a tra��o do avi�o e a velocidade do avi�o
    p_disp_alpha(i) = p_disp(i)*cosd(beta(i));
    T_disp(i) = p_disp(i)/v_final(i);

    % Raz�o de subida e gradiente de subida
    inc_vel = 0.05; % considerado para o c�lculo de velocidade de decolagem
    v_stall1 = min(v_final);  % velocidade m�nima do avi�o estimada por meio dos c�lculos iterativos
    % v_stall3 = sqrt((2*(tow*9.81-Tracao(i)*sind(beta(i))))/(ro*Sw*CLmaxtot)); % velocidade de estol fora de efeito solo
    razao_de_subida(i) = (p_disp_alpha(i) - p_req(i))/(tow*9.81);
    grad_de_subida(i) = 100*razao_de_subida(i)/v_final(i);
    theta(i) = asind((t_disp(i)*cosd(beta(i))-D_total(i))/(tow*9.81));
    v_h(i) = v_final(i)*cosd(theta(i));

    if v_final(i) >= v_dec-inc_vel && v_final(i) < v_dec+inc_vel  
        grad_sub = grad_de_subida(i); % gradiente de subida logo ap�s a decolagem
        rc = razao_de_subida(i);
    end
end


%% Verifica��o gradiente de subida 

if grad_sub <= 0
    disp('Gradiente de subida negativo fodase')
    delete(file_aviao)
    PrintParada(save,14,peso,Pcg,S,AR,CD_corrida,CL_corrida,max(MTOWt));
    return
elseif grad_sub < grad_sublim    
    disp('Gradiente de subida lament�vel')
    delete(file_aviao)
    PrintParada(save,15,peso,Pcg,S,AR,CD_corrida,CL_corrida,max(MTOWt));
    return
end

j = find(MTOWt == max(MTOWt), 1);
CL_polardec = CL_polar(j);
alpha_dec = alpha(j);
tow = max(MTOWt);

%% Deflex�o de profundor
alpha = 0:2:alfamax;
de = NaN(1,length(alpha));
i = 1;

for k = 1:length(alpha)        % An�lise da deflex�o do profundor pelo avl
    [de(i)] = deflexaoAVL(alpha(k),0,file_aviao);    
    i = i + 1;
end

if max(abs(de)) > de_lim(1)
    fprintf('M�xima deflex�o de profundor: %.4f\n',max(abs(de)))
    result = result - 0.1;
    PrintParada(save,16,peso,Pcg,S,AR,CD_corrida,CL_corrida,max(MTOWt))
    return
end

%% Peso

cargapaga = max(MTOWt) - peso;

result = -cargapaga;

fprintf('CL decolagem = %.4f\n',CL_polardec) %CL max q um ponto espec�fico da asa gera
fprintf('Alfa decolagem = %.4f\n',alpha_dec)
fprintf('CL max = %.4f\n',CLmax)
fprintf('CL max asa = %.4f\n',CLmaxasa)
fprintf('AR = %.4f\n',AR)
fprintf('�rea = %.4f\n',S)
fprintf('Gradiente de subida = %.4f\n',grad_sub(j))
fprintf('Vs = %.4f\n',Vs(j))
fprintf('Margem est�a = %.4f\n',H(length(H)))
fprintf('Posi�ao do CG [raiz] = %.4f\n',Pcg)
fprintf('Posi�ao do CG [MAC] = %.4f\n',Xcg_mac)
fprintf('Tempo de uma itera��o completa = %.4f\n',toc)
fprintf('Peso = %.4f <<<\n\n',peso)
fprintf('MTOW = %.4f <<<\n\n',max(MTOWt))
fprintf('Carga Paga = %.4f <<<\n\n', cargapaga)

fileIDBOM = fopen('backup.txt','a');
fprintf(fileIDBOM,'\n%s\n\nx = [ %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f ]',datetime,save(1),save(2),save(3),save(4),save(5),save(6),save(7),save(8),save(9),save(10));
fprintf(fileIDBOM,'\n\nMTOW = %.4f <<<\n\nPeso = %.4f\n\nCarga Paga = %.4f\n\nPosi�ao do CG [raiz] = %.4f',max(MTOWt),peso ,cargapaga ,Pcg);
fprintf(fileIDBOM,'\n\n�rea: %.4f | AR: %.4f | CDcorrida: %.4f | CLcorrida: %.4f | Parada: %f',S,AR,CD_corrida,CL_corrida,0);
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

PrintParada(save,17,peso,NaN,NaN,NaN,NaN,NaN,NaN);
return

end
% end
