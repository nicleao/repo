% Initialization and run of differential evolution optimizer.
% A simpler version with fewer explicit parameters is in run0.m
%
% Here for Rosenbrock's function
% Change relevant entries to adapt to your personal applications
%
% The file ofunc.m must also be changed 
% to return the objective function
%
clc
clear all+
close all
% delete cd_hist.nfo    
% save cd_hist.nfo -ascii

toint=cputime;

% VTR		"Value To Reach" (stop when ofunc < VTR)
		VTR = 1.e-4; 

% D		number of parameters of the objective function 
		D = 11; 

% XVmin,XVmax   vector of lower and bounds of initial population
%    		the algorithm seems to work well only if [XVmin,XVmax]
%    		covers the region where the global minimum is expected
%               *** note: these are o bound constraints!! ***
%        x_emp  posb  cr   lambda  off  raiz meio  ponta  ceh   beh    iw 
%           1     2    3      4     5     6    7     8     9     10    11
 XVmin =  [0.75  0.5  0.45  0.45  0.00  0.50  0.50  0.50  0.15  0.4  4.00];
 XVmax =  [1.1  0.68  0.55  0.58  0.15  6.49  4.49  5.49  0.25  1.00  6.50];

% Variáveis de otimização
% x(1) = x_emp [m] Distância entre bordo de ataque da asa e bordo de ataque das empenagens (raiz)
% x(2) = Seção reto-trapezoidal [% da semi envergadura]
% x(3) = Corda da raiz [m]
% x(4) = Afilamento (corda da ponta / corda da raiz)
% x(5) = Offset da ponta da asa [m]
% x(6) = Perfil da raiz
% x(7) = Perfil do meio
% x(8) = Perfil da ponta
% x(9) = Corda da EH [m]
% x(10) = Envergadura da EH [m]
% x(11) = Ângulo de incidência na corrida da asa [°]
		
% y		problem data vector (remains fixed during optimization)
		y=[]; 

% NP            number of population members
		NP = 110; 

% itermax       maximum number of iterations (generations)
		itermax = 200; 

% F             DE-stepsize F ex [0, 2]
		F = 0.8; 

% CR            crossover probabililty constant ex [0, 1]
		CR = 0.8; 

% strategy       1 --> DE/best/1/exp           6 --> DE/best/1/bin
%                2 --> DE/rand/1/exp           7 --> DE/rand/1/bin
%                3 --> DE/rand-to-best/1/exp   8 --> DE/rand-to-best/1/bin
%                4 --> DE/best/2/exp           9 --> DE/best/2/bin
%                5 --> DE/rand/2/exp           else  DE/rand/2/bin

		strategy = 8;
 
% refresh       intermediate output will be produced after "refresh"
%               iterations. No intermediate output will be produced
%               if refresh is < 1
		refresh = 1; 

[x,f,nf] = devec('estimativa_convencional_TA2021',VTR,D,XVmin,XVmax,y,NP,itermax,F,CR,strategy,refresh)

% options = optimoptions('ga','PlotFcn',{@gaplotscores,@gaplotbestindiv},...
%     'UseParallel',true, 'UseVectorized', false,'PopulationSize',NP,...
%     'Generations',itermax,'Display','iter','MaxStallGenerations',100,...
%     'FunctionTolerance',VTR,'SelectionFcn',{@selectiontournament,3});
% 
% [x,fval,exitflag,output,population,scores] = ga(@estimativa_biplana,D,[],[],[],...
%     [],XVmin,XVmax,[],[],options);

tempo_cpu = cputime-toint