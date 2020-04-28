%Symbolic variables
syms A B C v1 v2 v3 v4; 

%Parameters and equations for this problem
I1_conc = logspace(-2,3,10);
I2_conc = logspace(-2,3,10);

A_conc = zeros(numel(I1_conc),numel(I2_conc));
I1_count = 1;
for i= I1_conc
    I2_count=1;
    for j = I2_conc
              
            I1 = i;
            I2 = j;
            
            Vmax1 = 5.0;
            Vmax2 = 5.0;
            Vmax3 = 1.0;
            Vmax4 = 1.0;
            Ks1 = 5.0;
            Ks2 = 5.0;
            Ks3 = 5.0;
            Ks4 = 5.0;
            KI1 = 1.0;
            KI2 = 1.0;
            S_total = 100.0;
            
            v1 = (Vmax1*A)/((1+(I1/KI1))*(Ks1+A));
            v2 = (Vmax2*A)/((1+(I2/KI2))*(Ks2+A));
            v3 = (Vmax3*B)/(Ks3+B);
            v4 = (Vmax4*C)/(Ks3+C);

            eqn1 = v1 - v3 == 0;
            eqn2 = v2 - v4 == 0;
            eqn3 = A + B + C == S_total;
            
            eqn4 = B >= 0;
            eqn5 = C >= 0;
            eqn6 = A >= 0;

            eqns = [eqn1,eqn2,eqn3,eqn4,eqn5,eqn6];

            S = solve(eqns,[A,B,C]);

            value_A = real(double(S.A));
            A_conc(I2_count,I1_count) = value_A;

            I2_count=I2_count+1;
    end
    I1_count = I1_count+1;
end

%3D-Plot for [A] vs. [I1] and [I2]
[X,Y]=meshgrid(log10(I1_conc),log10(I2_conc));
surf(X,Y,A_conc)
colormap autumn;
title('change of [A] with [I1] and [I2] for Ks = 5.0')
xlabel('[I1]')
ylabel('[I2]')
zlabel('[A]')