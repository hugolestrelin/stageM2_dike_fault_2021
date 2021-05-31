dfs=[0,1000,2000,3000,4000,5000,10000,15000,20000];

%sill depth :
%-3000 :
sign_tt_3 =[-0.0563   -0.0246   -0.0121   -0.0073   -0.0051   -0.0038   -0.0014   -0.0007   -0.0005]


tau_tt_3 =[0.0175    0.0205    0.0111    0.0059    0.0034    0.0021    0.0004    0.0001    0.0000]
%-1000 :
sign_tt_1 =[-0.1862   -0.0384   -0.0186   -0.0113   -0.0076   -0.0055   -0.0018   -0.0009   -0.0005]


tau_tt_1 =[0.1010    0.0250    0.0069    0.0027    0.0013    0.0007    0.0001    0.0000    0.0000]
%-5000 :
sign_tt_5 =[-0.0268   -0.0168   -0.0098   -0.0061   -0.0042   -0.0031   -0.0012   -0.0006   -0.0004]


tau_tt_5 =[0.0030    0.0100    0.0082    0.0057    0.0039    0.0026    0.0006    0.0002    0.0001]
%-10000 :
sign_tt_10 =[-0.0079   -0.0069   -0.0054   -0.0041   -0.0031   -0.0024   -0.0009   -0.0005   -0.0003]


tau_tt_10 =[-0.0014    0.0014    0.0025    0.0027    0.0024    0.0021    0.0008    0.0003    0.0002]

plot(dfs,sign_tt_1,'r','LineWidth',3)
hold on 
plot(dfs,sign_tt_3,'g','LineWidth',3)
hold on
plot(dfs,sign_tt_5,'b','LineWidth',3)
hold on
plot(dfs,sign_tt_10,'k','LineWidth',3)
xlabel('distance sill-faille (m)','FontSize', 24)
ylabel('taux de variation de \sigma_n par Pa de surpression dans le sill','FontSize',24)
legend('-1000','-3000','-5000','-10000','FontSize', 24)
title('variation de la contrainte normale en fonction de la distance sill-faille','fontsize',30)

plot(dfs,tau_tt_1,'r','LineWidth',3)
hold on 
plot(dfs,tau_tt_3,'g','LineWidth',3)
hold on
plot(dfs,tau_tt_5,'b','LineWidth',3)
hold on
plot(dfs,tau_tt_10,'k','LineWidth',3)
xlabel('distance sill-faille (m)','FontSize', 24)
ylabel('taux de variation de \tau par Pa de surpression dans le sill','FontSize', 24)
legend('-1000','-3000','-5000','-10000','FontSize', 24)
title('variation de la contrainte tangentielle en fonction de la distance sill-faille','fontsize',30)