%% Observa��es Gerais
% Essa fun��o retorna o valor da corda em cada posi��o da envergadura para
% diversas se��es da asa. Por motivos estruturais e construtivos, tentar
% manter ao m�ximo em 5 varia��es de corda.
% Alterem com cuidado e anotem a data e altera��es muito grandes abaixo em
% "Vers�es"

% Entradas:
%| chord: vetor com os valores das cordas nas se��es desejadas, deve ser necessariamente em ordem decrescente (da raiz para a ponta) [m]
%| posb: posi��o de cada varia��o de corda na envergadura, deve ser necessariamente em ordem crescente (de 0 at� b/2) [m]
%| X: se��o reto-trapezoidal (em [m] para uma semi-envergadura)
%| b: envegadura [m]
%| offset: offset de cada varia��o de corda na envergadura, deve ser necessariamente da ra�z para a ponta

% Saidas:
%| y: posi��o equivalente na envergadura
%| x: offset do BA em cada posi��o da envergadura
%| c: valor da corda em cada posi��o da envergadura

% Vers�es:
%(04/02/17) Vers�o 1.0 - Programa base inicial
%(06/02/17) Vers�o 1.1 - Vari�veis y,c antes das condi��es para evitar erro de output de fun��o (Gustavo)
%(01/03/18) Vers�o 1.2 - Inserido input de offsets na envergadura para obter a distribui��o de offset na envergadura --> obter output 'offset_mac' (Bal)

%% Fun��o corda

function [y,c,x] = Corda(chord,posb,X,b,offset)

    y = (0:0.01:1)*b/2;        
    c = NaN(1,length(y));
    x = NaN(1,length(y));
    
    if posb(1) ~= 0
        display('Erro: Primeira posi��o de varia��o de corda deve ser a da raiz (valor 0)')
        return
    end

%     if posb(size(posb,2)) ~= b/2
%         display('Erro: �ltima posi��o de varia��o de corda deve ser a da ponta (valor b/2)')
%         return
%     end

%     if posb(2) <= X
%         display('Erro: Primeira varia��o de corda deve ocorrer depois da se��o reto-trapezoidal')
%         return
%     end

    if issorted(posb) ~= 1 
        display('Erro: Posi��o das cordas na envergadura deve ser crescente')
        return
    end
    
    if issorted(fliplr(chord)) ~= 1
        display('Erro: Valores das cordas devem ser decrescentes')
        return
    end

    if size(posb,2) > 5
        display('Aten��o: N�o aconselh�vel realizar tantas varia��es de corda')
    end

    if length(posb) ~= length(chord)
        display('Erro: Vetor de posi��o e de cordas deve ser do mesmo tamanho')
        return
    end

%%
    j = 2;
    
    for i=1:length(y)        

        if y(i) > posb(j)
            j = j+1;
        end

        if y(i) <= X 
            c(i) = chord(1);
            x(i) = 0;
            
        elseif j == 2 
            c(i) = polyval(polyfit([X posb(j)],[chord(1) chord(j)],1),y(i));    
            x(i) = polyval(polyfit([X posb(j)],[0 offset(j)],1),y(i));
        
        else
            c(i) = polyval(polyfit([posb(j-1) posb(j)],[chord(j-1) chord(j)],1),y(i));
            x(i) = polyval(polyfit([posb(j-1) posb(j)],[offset(j-1) offset(j)],1),y(i));
        end
        
    end