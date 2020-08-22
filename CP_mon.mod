// Based on Dynare model file to calculate the Gertler and Karadi (2011) model
// Model section was originally created by Peter Karadi

//Model with financial accelerator and credit policy
// Parameters calibrated for Mongolia, datafile contains data for Mongolia

parameters rho_e_i sigma_e_i betta sigm hbt mu zetta gamm alfa eta_i epsil thetta theta_P kappa_p kappa_y rho_i rho_ksi sigma_ksi rho_a sigma_a rho_g sigma_g rho_psi sigma_psi tau knb varho chi b delta_c I_ss;
var S Y_obs It_obs C_obs e_i Y Yi K L I C G Q lambda Lambda Rs R N Ne Nn nu eta phi z x Pi w U X D F Z i delta Ic a ksi infl infls psi;
varexo e_a e_ksi e_g e_psi shock_e_i;

betta=0.99;           //Discount factor
sigm=0.83000000;       //Inverse elasticity of intertemporal substitution
hbt=0.48000000;        //Habit
mu=0.3817000 ;    //Inverse  Frisch elasticity of labor supply
zetta=7.20000000;     //Elasticity of marginal depreciation rate wrt utilization rate
gamm=0.9;            //Survival rate of bankers
alfa=0.51000000;      //Capital share
eta_i=1.72800000;     //Inverse elasticity of net investent to the price of capital
epsil=1.6900000;    //Elasticity of substitution - Retailers
thetta=0.71000000;       //Price stickiness
theta_P=0.39000000;     //Indexing price to past inflation
kappa_p=1.41000000;  //Taylor rule weight on inflation
kappa_y=0.14000000;   //Taylor rule weight on output gap
rho_i=0.500000;       //IR smoothing
rho_ksi=0.66000000;   //Persistence of capital quality shock
sigma_ksi=0.1000000;  //SE of capital quality shock
rho_a=0.66000000;     //Persistence of TFP shock
sigma_a=0.1000000;    //SE of TFP shock
rho_g=0.89000000;     //Persistence of gov.t spending shock
sigma_g=0.1000000;    //SE of gov.t spending shock
rho_e_i=0.3;          //Persistence of monetary policy shock
sigma_e_i=0.1;        //SE of monetary policy shock
rho_psi=0.42000000;   //Persistence of credit policy shock
sigma_psi=0.1200000;  //SE of credit policy shock
tau=0.0100000;        //multiplier of gov.t credit policy intermediation
knb=0.0008;         //new banks' transfer share
varho=0.1169;        //divertable fraction of funds
chi=1;                //relative utility weight of labor
b=0.015877023;        //depreciation rate utilization rate weight
delta_c=0.04;         //depreciation rate constant
I_ss=0.8279;          //SS investment



//Declaring the model equations
model;

//Households

//1. Consumption optimization 
exp(lambda)  =   (exp(C)-hbt*exp(C(-1)))^(-sigm)-hbt*betta*(exp(C(+1))-hbt*exp(C))^(-sigm);

//2. Consumption Euler
betta*exp(Lambda(+1))  =   1/exp(R);

//3. Discount factor for convenience
exp(Lambda)  =   exp(lambda)/exp(lambda(-1));

//4. Labor supply condition + labor demand
chi*exp(L)^mu    =   (1-alfa)*exp(lambda)*exp(Pi)*exp(Y)/exp(L);


//Financial Intermediaries
//5. Capital of financial intermediaries
exp(nu)     =   betta*exp(Lambda(+1))*(1-gamm)*(exp(Rs(+1))-exp(R))+gamm*betta*exp(Lambda(+1))*exp(nu(+1))*exp(x(+1));

//6. Net worth of financial intermediaries
exp(eta)    =   gamm*betta*exp(Lambda(+1))*exp(eta(+1))*exp(z(+1))+(1-gamm);

//7. Financial sector leverage
exp(phi)    =   exp(eta)/((varho-exp(nu))*(1-psi));

//8. Fin. intermediaries capital change
exp(z)      =   exp(phi(-1))*(1-psi(-1))*(exp(Rs)-exp(R(-1)))+exp(R(-1));

//9. Fin. intermediaries net worth change
exp(x)      =   ((1-psi)*exp(phi)/(exp(phi(-1))*(1-psi(-1))))*exp(z);

//10. Fin. intermediaries balance sheet
exp(S)     =   exp(N)*exp(phi);

//11. Fin, intermediaries net worth
exp(N)      =   exp(Nn)+exp(Ne);

//12. Net worth of surviving fin.intermediaries
exp(Ne)     =   gamm*exp(N(-1))*exp(z);

//13. Net worth of new fin. intermediaries
exp(Nn)    =   knb*exp(S(-1))*(1-psi(-1));


//Intermediate goods producers

//14. Loan
exp(S)=exp(Q)*exp(K);

//15. Production
exp(Yi)     =   exp(a)*exp(L)^(1-alfa)*(exp(ksi)*exp(U)*exp(K(-1)))^alfa;

//16. Utilization FOC
alfa*exp(Pi)*exp(Yi)/exp(U) = b*exp(ksi)*exp(K(-1))*(exp(U)-1)^zetta;

//17. Labor FOC
exp(w)      =   (1-alfa)*exp(Y)/exp(L)*exp(Pi);

//18. Capital FOC
exp(Rs)     =   (alfa*exp(Pi)*exp(Yi)/exp(K(-1))+(exp(Q)-exp(delta))*exp(ksi))/exp(Q(-1));


//Capital good producers

//19. Investment FOC
exp(Q)-1  =   eta_i*((Ic+I_ss)/(Ic(-1)+I_ss)-1)*(Ic+I_ss)/(Ic(-1)+I_ss)+eta_i/2*((Ic+I_ss)/(Ic(-1)+I_ss)-1)^2-eta_i*betta*exp(Lambda(+1))*((Ic(+1)+I_ss)/(Ic+I_ss)-1)*((Ic(+1)+I_ss)/(Ic+I_ss))^2;

//20. Endogen. depreciation
exp(delta) = delta_c+b*(exp(U)-1)^(1+zetta)/(1+zetta);

//21. Investment (minus refurbishing)
Ic  =   exp(I)-exp(ksi)*exp(delta)*exp(K(-1));

//22. Capital dynamics
exp(K)  =  Ic+ exp(ksi)*exp(K(-1)); 


//Retailers

//23. Retail demand
exp(Yi)    =   exp(Y)*exp(D);

//24. Price 
exp(D)    =   thetta*exp(infl)^epsil*exp(infl(-1))^(-theta_P*epsil)*exp(D(-1))+(1-thetta)*((1-thetta*exp(infl(-1))^(theta_P*(1-thetta))*exp(infl)^(thetta-1))/(1-thetta))^(-epsil/(1-thetta));

//25. Retail markup
exp(X)  =   1/exp(Pi);

//26. Price optimization condition
exp(F)       =   exp(Pi)*exp(Y)+thetta*betta*exp(Lambda(+1))*(exp(infl))^(-epsil*theta_P)*exp(infl(+1))^epsil*exp(F(+1));

//27.
exp(Z)       =   exp(Y)+thetta*betta*exp(Lambda(+1))*exp(infl)^(theta_P*(1-epsil))*exp(infl(+1))^(epsil-1)*exp(Z(+1));

//28. Price choice
exp(infls)   =  exp(infl)*exp(F)/exp(Z)*epsil/(epsil-1);

//29. Final price
(exp(infl))^(1-epsil)     =   (1-thetta)*(exp(infls))^(1-epsil)+thetta*exp(infl(-1))^(theta_P*(1-epsil));


//Monetary policy

//30. Fisher relationship
exp(i)  =   exp(infl(+1))*exp(R);

//31. Taylor rule
exp(i)      =   exp(i(-1))^rho_i*((1/betta)*exp(infl)^kappa_p*(exp(X)/(epsil/(epsil-1)))^(kappa_y))^(1-rho_i)*exp(e_i);

//32. Market clearing condition
exp(Y)   =   exp(C)+exp(I)+tau*psi*exp(S)+exp(G)+eta_i/2*((Ic+I_ss)/(Ic(-1)+I_ss)-1)^2*(Ic+I_ss);


//Shocks
//33. Government expenditure
G   =   rho_g*G(-1)+e_g;

//34. Credit policy shock
psi    =   rho_psi*psi(-1)+e_psi;

//35. TFP shock
a  =   rho_a*a(-1)+e_a;

//36. Capital quality shock
ksi=   rho_ksi*ksi(-1)+e_ksi;

//37. Monetary policy shock
e_i=rho_e_i*e_i(-1)+shock_e_i;


//Observation equations as the data is detrended

Y_obs=Y-steady_state(Y);
C_obs=C-steady_state(C);
It_obs=I-steady_state(I);

end;

initval;
  S      	=	2.309576344	;
    Y_obs   	=	0	;
    It_obs  	=	0	;
    C_obs   	=	0	;
   e_i     	=	0	;
    Y       	=	0.937038528	;
    Yi      	=	0.937038528	;
    K       	=	2.309576344	;
    L       	=	-0.64398914	;
    I       	=	-0.909299468	;
    C       	=	0.139418451	;
    G       	=	0	;
   Q       	=	0	;
   lambda  	=	-0.217696388	;
    Lambda  	=	0	;
    Rs      	=	0.012696341	;
    R       	=	0.01005164	;
    N       	=	-1.034780951	;
    Ne      	=	-1.057715411	;
    Nn      	=	-4.821322486	;
    nu      	=	-4.808319055	;
    eta     	=	1.125547214	;
   phi     	=	3.344357295	;
    z       	=	0.082425348	;
    x       	=	0.082425348	;
    Pi      	=	-0.895792048	;
    w       	=	-0.028113732	;  
   U       	=	0.146488767	;
    X       	=	0.895792048	;
    D       	=	5.51E-11	;
    F       	=	1.254940873	;
    Z       	=	2.150728248	;
    i       	=	0.010054878	;
    delta   	=	-3.218875812	;
    Ic      	=	0	;
    a       	=	0	;
    ksi     	=	0	;
    infl    	=	0	;
    infls	=	0	;
    psi     	=	0	;
 
end;

shocks;
var e_a=sigma_a^2;
var e_ksi=sigma_ksi^2;
var e_g=sigma_g^2;
var shock_e_i=sigma_e_i^2;
var e_psi=sigma_psi^2;
end;

varobs infl Y_obs psi It_obs C_obs;

estimated_params;
hbt,beta_pdf,0.5,0.15;
sigm,normal_pdf,0.85,0.15;
mu,gamma_pdf,0.4,0.1;
eta_i,normal_pdf,1.7,0.5;
zetta,normal_pdf,7.2,2;
theta_P,beta_pdf,0.5,0.15;
kappa_p,normal_pdf,1.5,0.125;
kappa_y,normal_pdf,0.125,0.05;
rho_i,beta_pdf,0.5,0.1;
rho_ksi,beta_pdf,0.5,0.2;
stderr e_ksi,inv_gamma_pdf,0.1,2;
rho_a,beta_pdf,0.5,0.2;
stderr e_a,inv_gamma_pdf,0.1,2;
rho_g,beta_pdf,0.5,0.2;
stderr e_g,inv_gamma_pdf,0.1,2;
rho_e_i,beta_pdf,0.5,0.2;
stderr shock_e_i,inv_gamma_pdf,0.1,2;
rho_psi,beta_pdf,0.5,0.2;
stderr e_psi,inv_gamma_pdf,0.1,2;
end;

//Run model diagnostics
//model_diagnostics;
stoch_simul (order=1) Y C K I infl;

//Check identification of the parameters at prior mode
//identification;

//Estimation of the posterior mode and parameter identification at posterior mode
//estimation(datafile='data_mongolia.xlsx',mode_compute=4,mh_replic=0,plot_priors=0);
//identification(parameter_set=posterior_mode);

//To run everything at once
//estimation(datafile='data_mongolia.xlsx',mode_compute=4,mh_replic=250000, mh_drop=0.2,mh_jscale=0.36,plot_priors=0,bayesian_irf) infl Y C I K U delta Rs R S N;

//See results
//estimation(datafile='data_mongolia.xlsx',mh_replic=0,mode_compute=0,mode_file=CP_mon_mode,load_mh_file,mh_drop=0.2,plot_priors=0,bayesian_irf) infl Y C I K Q U delta S N Rs R;


//Running counterfactual experiments using smoothed shocks
//with all shocks
y0=oo_.dr.ys;
dr=oo_.dr
iorder=1;
ex_=[oo_.SmoothedShocks.e_a oo_.SmoothedShocks.e_ksi oo_.SmoothedShocks.e_g oo_.SmoothedShocks.e_psi oo_.SmoothedShocks.shock_e_i]
save results ex_;
y_=simult_(y0,dr,ex_,1);

//setting e_psi to 0
epsi_ = zeros(49,1);
ex1_=[oo_.SmoothedShocks.e_a oo_.SmoothedShocks.e_ksi oo_.SmoothedShocks.e_g epsi_ oo_.SmoothedShocks.shock_e_i]
y1_=simult_(y0,dr,ex1_,1);

//Plotting inflation and output with counterfactual calculations
figure
plot(y_(strmatch('Y_obs',M_.endo_names,'exact'),:));
title('Output (observed)');
hold on;
plot(y1_(strmatch('Y_obs',M_.endo_names,'exact'),:)); 
figure
plot(y_(strmatch('infl',M_.endo_names,'exact'),:));
title('Inflation (observed)');
hold on;
plot(y1_(strmatch('infl',M_.endo_names,'exact'),:)); 