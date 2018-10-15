classdef phase < handle 
    %CLASS FOR THE CREATION OF PHASE(LIQUID,VAPOR) OBJECTS WITH PROPERTIES
    %CALCULATED USING PENG-ROBINSON EOS AND VDW MIXING RULES
    properties %parameters and thermodynamic properties of the phase
        %pure component parameters
        a,b,A,B
        %mixture parameters
        am,bm,Am,Bm,C
        %calculated properties
        Zfactor,fugacity
    end
 
    methods 
        %pure component parameters method
        function obj = parameters(obj,w,T,Tc,R,Pc,P)
            m = 0.37464 + 1.54226*w - 0.26992*w.^2;
            Tr = T./Tc;
            a = 0.45724*R^2*Tc.^2.*(1+m.*(1-Tr.^0.5)).^2./Pc;
            b = 0.07780*R*Tc./Pc;
            A = a*P/(R*T)^2;
            B = b*P/(R*T);
            obj.a = a; obj.b = b;
            obj.A = A; obj.B =B;
        end
        %mixture parameters method
        function obj = mixtureparameters(obj,x,k,P,R,T)
            components = 1:length(x);
            bm = 0; %b mixture parameter
            for i = components
                sum = x(i)*obj.b(i);
                bm = bm + sum;
            end
            am = 0; %a mixture parameter
            for i = components
                a_m = 0;
               for j = components
                   sum = x(i)*x(j)*(1-k(i,j))*(obj.a(i)*obj.a(j))^0.5;
                   a_m = a_m + sum;
               end
               am = am + a_m;
            end
            Am = am*P/(R*T)^2; %A mixture parameter
            Bm = bm*P/(R*T); % B mixture parameter
            C = []; %C mixture parameter
            for i = components 
                C_ = 0;
                for j = components 
                   sum = x(j)*(obj.a(i)*obj.a(j)).^0.5*(1-k(i,j));
                   C_ = C_ + sum;
                end
                C = [C C_];
            end
            C = (Am/Bm)*(-(obj.b/bm)+(2/am)*C);
            obj.am = am; obj.bm = bm;
            obj.Am = Am; obj.Bm = Bm; obj.C = C;
        end
        %fugacity method
        function obj = fugacitycalc(obj,state)
            A = 1;
            B = (obj.Bm-1);
            C = (obj.Am-3*obj.Bm^2-2*obj.Bm);
            D = -obj.Am*obj.Bm+obj.Bm^3+obj.Bm^2;
            %compresibility factor
            Z = roots([A B C D]); %find the roots
            Z = Z(find(imag(Z)==0)); %find the real roots
            switch state
                case 'liquid'
                    Z = min(Z); %root for the liquid
                case 'vapor'
                    Z = max(Z);
            end
            phi1 = log((Z+(2^0.5+1)*obj.Bm)/(Z-(2^0.5-1)*obj.Bm));
            phi = exp(obj.b/obj.bm*(Z-1)-log(Z-obj.Bm)-obj.C/(2*2^0.5)*phi1);
            obj.Zfactor = Z; obj.fugacity = phi;
        end
    end
end