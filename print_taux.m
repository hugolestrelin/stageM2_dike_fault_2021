dfs=[0,1000,2000,3000,4000,5000,10000,15000,20000];
dsn=[0.0422,0.0463,0.0307,0.0204,0.0143,0.0105,0.0035,0.0017,0.0010];
tau=[-0.0662,-0.0408,-0.0185,-0.0088,-0.0046,-0.0026,-3.0740e-04,-7.7977e-05,-2.8039e-05];

plot(dfs,dsn)
xlabel('distance sill-faille (m)')
ylabel('sigma_n (Pa)')
title('variation de la contrainte normale en fonction de la distance sill-faille','fontsize',15)

plot(dfs,tau)
xlabel('distance sill-faille (m)')
ylabel('tau (Pa)')
title('variation de la contrainte tangentielle en fonction de la distance sill-faille','fontsize',15)