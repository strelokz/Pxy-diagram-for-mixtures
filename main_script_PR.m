%COMPUTE THE Pxy DIAGRAM FOR A COMPONENT IN A LIQUID MIXTURE 
%USING THE PENG-ROBINSON EOS AND VAN DER WALLS MIXING RULES
%Tc, Pc and w must be provided in a row vector.
%the binary interaction parameters come in a square matrix
%the current mixture components are CO2 and n-pentane
clc,clear
components = {'CO2' 'n-pentane'};
T = 344.32; %Kelvin
P = 1e5; %Pascal
x1 = [0:0.01:1]';%first component composition
x2 = [1 - x1];%second component composition
x = [x1 x2];
y = [0.5 0.5]; %initial guess for vapor mixture composition
Tc = [304.2 469.7]; %critical temperature(K)
Pc = [73.74 33.7]*1e5; %critical pressure(Pa)
w = [0.225 0.252]; %acentric factor
R = 8.314; %gas constant(J/molK)
kij = [0 0.12 ; 0.12 0]; %binary interaction parameters

M = []; %matrix for the Pxy diagram
for i = 1:length(x)
    xi = x(i,:);
    for i = 1:1000
        %Fugacity of the liquid mixture
        liquid = phase;
        liquid.parameters(w,T,Tc,R,Pc,P);
        liquid.mixtureparameters(xi,kij,P,R,T);
        state = 'liquid';
        liquid.fugacitycalc(state);
        %Fugacity of the vapor mixture
        vapor = phase;
        vapor.parameters(w,T,Tc,R,Pc,P);
        vapor.mixtureparameters(y,kij,P,R,T);
        state = 'vapor';
        vapor.fugacitycalc(state);

        K = liquid.fugacity./vapor.fugacity;
        if abs(K*xi'-1) > 0.0001
            P = P*(K*xi');
            y = K.*xi/(K*xi');
        else
            break;
        end
    end

    M = [M; P/1e5 xi(1) y(1) liquid.fugacity(1) vapor.fugacity(1)];
    if abs(xi(1)- y(1)) <0.01 && xi(1) ~= 0
        break;
    end
end

varnames = {'Pressure','xCO2','yCO2','L_fugacity','V_fugacity'};
data = table(M(:,1),M(:,2),M(:,3),M(:,4),M(:,5));
data.Properties.VariableNames = varnames;
data
plot([M(:,2),M(:,3)],M(:,1))

legend({'Liquid line','Vapor line'},'location','northwest')
title('Pxy diagram for the binary mixture: CO_2- pentane at 344.32K')
ylabel('Presion(bar)')
xlabel('CO_2 molar fraction')



