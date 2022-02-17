%% Observações gerais
% Função que cria um arquivo .avl com a configuração convencional.
% Arquivo de comandos para rodar dentro do AVL e somente criar o avião.
% IMPORTANTE: Referenciado no CG

% Entradas:
% S: Área de referência [m^2]
% mac: Corda de referência [m]
% b: Envergadura de referência [m]
% posb: Vetor com posição das cordas da asa [m]
% offset: Vetor offset da asa [m]
% died: Diedro [°]
% ail: Distância onde se inicia o aileron na semi-envergadura [m]
% c_ail: Corda da asa quando inicia o aileron
% Pcg: porcentagem do CG em relação a raíz [%]
% chord: vetor com cordas da asa [m]
% Xhinge_eh: [%] da corda da hinge da EH
% Xhinge_a: [%] da corda da hinge do aileron
% SgnDup_aileron: sinal de deflexão para superfícies duplicadas [-1,0,1]
% perfilr: Perfil da raiz da asa
% perfilm: Perfil o meio da asa
% perfilp: Perfil da ponta da asa
% perfilh: Perfil da EH  
% perfilv: Perfil da EV
% beh: Envergadura na EH
% offseteh: Offset da ponta da EH
% creh: Corda da raiz da EH
% cteh: Corda da ponta da EH
% bev: Envergadura na EV
% crev: Corda da raiz da EV
% ctev: Corda da ponta da EV
% offsetev: Offset na EV
% x_emp: Distância entre bordo de ataque da asa e bordo de ataque da empenagem (raiz)
% h_emp: Altura entra asa e EH (Bordo de ataque da asa até bordo de ataque da EH)
% iw: Ângulo de incidência da asa
% ih: Ângulo de incidência da EH
% nEV: Quantidade de EVs
% H_CG: Altura CG a partir da altura da asa
% ef_solo:
% h_asa: Altura da asa em relação ao solo (intradorso da raiz ate o solo)
% Zcg: altura (estimada) do CG em relação ao intradorso da corda da raiz [m]
% i_corrida_w: ângulo de corrida da asa [°]
% X: seção reto-trapezoidal [m]
% Xhinge_elevon: posição da hinge [% da corda]

%% Função para escrita da aeronave no AVL
function[file_aviao] = convencionalAVL(S, mac, b, posb, offset, died, ail, c_ail, Pcg, chord, Xhinge_eh, Xhinge_a, SgnDup_aileron, perfilr, perfilm, perfilp, perfilh, perfilv, beh, offseteh, creh, cteh, bev, crev, ctev, offsetev, x_emp, h_emp, iw, ih, nEV, H_CG, ef_solo, h_asa)

% Valores de entrada e condições
if ail <= posb(2)
      porcentagem = ail/posb(2);
      [Interpolado] = InterpolaXfoil(perfilr,perfilm,porcentagem);
else
      porcentagem = (ail-posb(2))/(posb(3)-posb(2));
      [Interpolado] = InterpolaXfoil(perfilm,perfilp,porcentagem); 
end

if ef_solo == 0
    file_aviao = 'convencionalAVL.avl';
else
    file_aviao = 'convencionalAVL_ground.avl';
end

% ARQUIVO DO AVIÃO
fileID = fopen(file_aviao,'w');

fprintf(fileID,  'Avião Convencional AVL 2021\n');
fprintf(fileID,  '#Mach \n');
fprintf(fileID,  '0.0 \n');

fprintf(fileID,  '#IYsym IZsym Zsym \n');
if ef_solo == 0
    fprintf(fileID, '%-8d%-8d%-8d\n\n',0,0,0);
else
    fprintf(fileID, '%-8d%-8d%-8d\n\n',0,1,-(H_CG+h_asa));
end

fprintf(fileID, '#Sref Cref Bref \n');
fprintf(fileID, '%.4f %.4f %.4f \n\n',S,mac,b);

fprintf(fileID, '#Xref Yref Zref \n');
fprintf(fileID, '%-8d%-8d%-8d\n\n',0,0,0);

fprintf(fileID, '#CDp \n');
fprintf(fileID,'%-8.4f\n\n',0);

%% Asa

fprintf(fileID,  '#------------------------------------------------------------- \n');
fprintf(fileID,  'SURFACE \n');
fprintf(fileID,  'Wing \n\n');

fprintf(fileID,  '#Nchordwise  Ccpace  Nspanwise  Sspace \n');
fprintf(fileID,  '20             1       20         0 \n\n');

fprintf(fileID,'COMPONENT\n');
fprintf(fileID,  '1 \n\n');

fprintf(fileID,  'YDUPLICATE \n');
fprintf(fileID,  '0.0 \n\n');

fprintf(fileID,'TRANSLATE\n');
fprintf(fileID,'%-8.4f%-8.4f%-8.4f\n\n',-Pcg*chord(1),0,-H_CG);

fprintf(fileID,  'ANGLE \n');
fprintf(fileID, '%.4f \n\n',iw);

% Seção da raiz
Xle = offset(1);
Yle = 0;
Zle = 0;

fprintf(fileID, '#------------------------------------------------------------- \n');
fprintf(fileID, 'SECTION \n');
fprintf(fileID, '#Xle      Yle       Zle      Chord    Ainc \n');
fprintf(fileID, '%.4f   %.4f    %.4f    %.4f    0  \n\n',Xle,Yle,Zle,chord(1));

fprintf(fileID, 'AFILE \n');
fprintf(fileID, '%s \n\n',perfilr);

% Seção do meio e começo do aileron
if ail < posb(2)
        % Aileron
        Xle = offset(2)*ail/posb(2);
        Yle = ail;
        Zle = ail*tand(died);
      

        fprintf(fileID, '#------------------------------------------------------------- \n');
        fprintf(fileID, 'SECTION \n');
        fprintf(fileID, '#Xle      Yle       Zle      Chord    Ainc \n');
        fprintf(fileID, '%.4f   %.4f    %.4f    %.4f    0  \n\n',Xle,Yle,Zle,c_ail);

        fprintf(fileID, '#Cname   Cgain  Xhinge  HingeVec     SgnDup \n');
        fprintf(fileID, 'CONTROL \n');
        fprintf(fileID, 'aileron  -1   %.4f    0    0   0    %.4f \n\n', (1-Xhinge_a),SgnDup_aileron);
           

        fprintf(fileID, 'AFILE \n');     
        fprintf(fileID, '%s \n\n',Interpolado);
        
        % Meio
        Xle = offset(2);
        Yle = posb(2);
        Zle = posb(2)*tand(died);

        fprintf(fileID, '#------------------------------------------------------------- \n');
        fprintf(fileID, 'SECTION \n');
        fprintf(fileID, '#Xle      Yle       Zle      Chord    Ainc \n');
        fprintf(fileID, '%.4f   %.4f    %.4f    %.4f    0  \n\n',Xle,Yle,Zle,chord(2));

        fprintf(fileID, '#Cname   Cgain  Xhinge  HingeVec     SgnDup \n');
        fprintf(fileID, 'CONTROL \n');
        fprintf(fileID, 'aileron  -1   %.4f    0    0   0    %.4f \n\n', (1-Xhinge_a),SgnDup_aileron);

        fprintf(fileID, 'AFILE \n');
        fprintf(fileID, '%s \n\n',perfilm);
        
        
else
        % Meio
        Xle = offset(2);
        Yle = posb(2);
        Zle = posb(2)*tand(died);

        fprintf(fileID, '#------------------------------------------------------------- \n');
        fprintf(fileID, 'SECTION \n');
        fprintf(fileID, '#Xle      Yle       Zle      Chord    Ainc \n');
        fprintf(fileID, '%.4f   %.4f    %.4f    %.4f    0  \n\n',Xle,Yle,Zle,chord(2));

        fprintf(fileID, 'AFILE \n');
        fprintf(fileID, '%s \n\n',perfilm);
        
        % Aileron
        Xle = (offset(3)*(ail-posb(2)))/((b/2)-posb(2)) + offset(2);
        Yle = ail;
        Zle = ail*tand(died);

        fprintf(fileID, '#------------------------------------------------------------- \n');
        fprintf(fileID, 'SECTION \n');
        fprintf(fileID, '#Xle      Yle       Zle      Chord    Ainc \n');
        fprintf(fileID, '%.4f   %.4f    %.4f    %.4f    0  \n\n',Xle,Yle,Zle,c_ail);

        fprintf(fileID, '#Cname   Cgain  Xhinge  HingeVec     SgnDup \n');
        fprintf(fileID, 'CONTROL \n');
        fprintf(fileID, 'aileron  -1   %.4f    0    0   0    %.4f \n\n', (1-Xhinge_a),SgnDup_aileron);

        fprintf(fileID, 'AFILE \n');        
        fprintf(fileID, '%s \n\n',Interpolado);
end

% Seção da ponta
Xle = offset(2) + offset(3);
Yle = b/2;
Zle = (b/2)*tand(died);

fprintf(fileID, '#------------------------------------------------------------- \n');
fprintf(fileID, 'SECTION \n');
fprintf(fileID, '#Xle      Yle       Zle      Chord    Ainc \n');
fprintf(fileID, '%.4f   %.4f    %.4f    %.4f    0  \n\n',Xle,Yle,Zle,chord(3));

fprintf(fileID, '#Cname   Cgain  Xhinge  HingeVec     SgnDup \n');
fprintf(fileID, 'CONTROL \n');
fprintf(fileID, 'aileron  -1   %.4f    0    0   0    %.4f \n\n', (1-Xhinge_a),SgnDup_aileron);

fprintf(fileID, 'AFILE \n');
fprintf(fileID, '%s \n\n',perfilp);

%% EH

fprintf(fileID,  '#------------------------------------------------------------- \n');
fprintf(fileID,  'SURFACE \n');
fprintf(fileID,  'EH \n\n');

fprintf(fileID,  '#Nchordwise  Ccpace  Nspanwise  Sspace \n');
fprintf(fileID,  '10             2       20         -2 \n\n');

fprintf(fileID,'COMPONENT\n');
fprintf(fileID,  '2 \n\n');

fprintf(fileID,  'YDUPLICATE \n');
fprintf(fileID,  '0.0 \n\n');

fprintf(fileID,'TRANSLATE\n');
fprintf(fileID,'%-8.4f%-8.4f%-8.4f\n\n',(-Pcg*chord(1) + x_emp),0,h_emp-H_CG);

fprintf(fileID,  'ANGLE \n');
fprintf(fileID, '%.4f \n\n',ih);

% Seção da raiz
Xle = 0;
Yle = 0;
Zle = 0;

fprintf(fileID, '#------------------------------------------------------------- \n');
fprintf(fileID, 'SECTION \n');
fprintf(fileID, '#Xle      Yle       Zle      Chord    Ainc \n');
fprintf(fileID, '%.4f   %.4f    %.4f    %.4f    0  \n\n',Xle,Yle,Zle,creh);

fprintf(fileID, '#Cname   Cgain  Xhinge  HingeVec     SgnDup \n');
fprintf(fileID, 'CONTROL \n');
fprintf(fileID, 'elevator  1   %.4f    0    0   0    1 \n\n', (1-Xhinge_eh));

fprintf(fileID, 'AFILE \n');
fprintf(fileID, '%s \n\n',perfilh);

% Seção da ponta
Xle = offseteh;
Yle = beh/2;
Zle = 0;

fprintf(fileID, '#------------------------------------------------------------- \n');
fprintf(fileID, 'SECTION \n');
fprintf(fileID, '#Xle      Yle       Zle      Chord    Ainc \n');
fprintf(fileID, '%.4f   %.4f    %.4f    %.3f    0  \n\n',Xle,Yle,Zle,cteh);

fprintf(fileID, '#Cname   Cgain  Xhinge  HingeVec     SgnDup \n');
fprintf(fileID, 'CONTROL \n');
fprintf(fileID, 'elevator  1   %.4f    0    0   0    1 \n\n', (1-Xhinge_eh));

fprintf(fileID, 'AFILE \n');
fprintf(fileID, '%s \n\n',perfilh);

%% EV

fprintf(fileID,  '#------------------------------------------------------------- \n');
fprintf(fileID,  'SURFACE \n');
fprintf(fileID,  'EV \n\n');

fprintf(fileID,  '#Nchordwise  Ccpace  Nspanwise  Sspace \n');
fprintf(fileID,  '8             2       10         -2 \n\n');

fprintf(fileID,'COMPONENT\n');
fprintf(fileID,  '3 \n\n');

if nEV ~= 1
    fprintf(fileID,  'YDUPLICATE \n');
    fprintf(fileID,  '0.0 \n\n');
end

fprintf(fileID,'TRANSLATE\n');
fprintf(fileID,'%-8.4f%-8.4f%-8.4f\n\n',(-Pcg*chord(1) + x_emp),0,h_emp-H_CG);

fprintf(fileID,  'ANGLE \n');
fprintf(fileID, '%.4f \n\n',0);

% Seção da raiz
Xle = 0;
if nEV == 1
    Yle = 0;
else 
    Yle = bh/2;
end
Zle = 0;

fprintf(fileID, '#------------------------------------------------------------- \n');
fprintf(fileID, 'SECTION \n');
fprintf(fileID, '#Xle      Yle       Zle      Chord    Ainc \n');
fprintf(fileID, '%.4f   %.4f    %.4f    %.4f    0  \n\n',Xle,Yle,Zle,crev);

fprintf(fileID, '#Cname   Cgain  Xhinge  HingeVec     SgnDup \n');
fprintf(fileID, 'CONTROL \n');
fprintf(fileID, 'rudder  1   %.4f    0    0   0    1 \n\n', 0);

fprintf(fileID, 'AFILE \n');
fprintf(fileID, '%s \n\n',perfilv);

% Ponta

Xle = offsetev;
if nEV == 1
    Yle = 0;
else 
    Yle = bh/2;
end
Zle = bev;

fprintf(fileID, '#------------------------------------------------------------- \n');
fprintf(fileID, 'SECTION \n');
fprintf(fileID, '#Xle      Yle       Zle      Chord    Ainc \n');
fprintf(fileID, '%.4f   %.4f    %.4f    %.4f    0  \n\n',Xle,Yle,Zle,ctev);

fprintf(fileID, '#Cname   Cgain  Xhinge  HingeVec     SgnDup \n');
fprintf(fileID, 'CONTROL \n');
fprintf(fileID, 'rudder  1   %.4f    0    0   0    1 \n\n', 0);

fprintf(fileID, 'AFILE \n');
fprintf(fileID, '%s \n\n',perfilv);

fclose('all');
end