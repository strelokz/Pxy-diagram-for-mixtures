%COMPUTE THE FUGACITY OF A LIQUID MIXTURE USING THE PENG-ROBINSON EOS
%AND VAN DER WALLS MIXING RULES
%composition, Tc, Pc and w must be provided in a row vector.
%the binary interaction parameters come in a square matrix
%the current mixture components are CO2 and n-pentane
clc,clear
components = {'CO2' 'n-pentane'};
T = 277.65; %Kelvin
P = 16.15e5; %Pascal
x = [0.32 0.68];%liquid mixture composition
Tc = [304.2 469.7]; %critical temperature(K)
Pc = [73.74 33.7]*1e5; %critical pressure(Pa)
w = [0.225 0.252]; %acentric factor
R = 8.314; %gas constant(J/molK)
kij = [0 0.12 ; 0.12 0]; %binary interaction parameters

%Fugacity of the liquid mixture
liquid = phase;
liquid.parameters(w,T,Tc,R,Pc,P);
liquid.mixtureparameters(x,kij,P,R,T);
liquid.fugacitycalc();

graph = table(x',liquid.fugacity','RowNames',components);
varnames = {'Xi','fugacity'};
graph.Properties.VariableNames = varnames;
graph



