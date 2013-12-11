// Kiyotaki-Moore (2008) in logs - only A stochastic

//********************* CHOOSE SIMULATION SCHEME *********************//
// CHOICE BETWEEN SIMULATED AND THEORETICAL MOMENTS: 0 for SIMULATED
//                                                   1 for THEORETICAL
// MUST RUN THE SIMULATED VERSION FIRST TO COMPUTE THE STEADY STATE
// THEN RUN ssGEN_shockA_KM.m TO UPLOAD THE STEADY STATE PARAMETERS
// ONLY THEN RUN THE THEORETICAL VERSION
@#define moments  = 0

// CHOICE BETWEEN MOMENTS AND IRFS : 0 for moments
//                                   1 for irfs
@#define irfs  = 0
//********************************************************************//

@#if (moments==0)

// NUMBER OF PERIODS TO SIMULATE:
periods 20100;

// VARIABLES:
var lgY, lgY_L, lgI, lgK, lgC, lgC_w, lgC_e, lgL, lgw, lgR, lgr, lgq, lgq_r, lgp, lgD, A, lgN_s, lgbond, lgdbond, lgpbond,
        lgAA;
varexo e_A;

// PARAMETERS:
parameters pi, psi, beta, rho_A, om, lambda, theta, gamma, mu, alpha, phi;
    pi        = 0.25;
    psi       = 0.050;
    omega  = 0.01;
    // calibration as in ...
    phi       = 0.025;
    // calibration as in D'Mello and Farhat (2007)
    theta     = 0.25;
    // calibration as in King and Rebelo (1999) and Greenwood, Herkowitz, and Huffman (1988):
    beta      = 0.984;
    om        = 6.2;
    mu        = 0.6;
    lambda    = 0.975;
    gamma     = 0.333;
    rho_A     = 0.979;
    sigma_e_A = 0.0072;
    alpha     = gamma*((1+mu)/(gamma+mu));

// MODEL SPECIFICATION:
model;
    (1-pi)*(((exp(lgR(+1))+exp(lgq(+1))*lambda)/(exp(lgq)))+(1-pi)*(((exp(lgdbond(+1)+exp(lgpbond(+1)))/exp(lgpbond))))-(exp(lgp(+1))/exp(lgp)))
      *(1/((exp(lgR(+1))+exp(lgq(+1))*lambda)*exp(lgN_s)+exp(lgp(+1))*exp(lgD)+exp(lgpbond(+1))*exp(lgbond)))
    = pi*((exp(lgpbond(+1))/exp(lgpbond))-(1/exp(lgpbond))*(exp(lgdbond(+1)+exp(lgpbond(+1)*psi+(1-psi)*exp(lgq_r(+1)))))+(exp(lgp(+1))/exp(lgp))-(1/exp(lgq))*(exp(lgR(+1))+lambda*(exp(lgq(+1))*phi+(1-phi)*exp(lgq_r(+1))))) 
      *(1/((exp(lgR(+1))+lambda*(exp(lgq(+1))*phi+(1-phi)*exp(lgq_r(+1))))*exp(lgN_s)+exp(lgpbond(+1))*exp(lgbond)+exp(lgp(+1))*exp(lgD)));
    (1-theta*exp(lgq))*exp(lgI) = pi*(beta*((exp(lgR)+exp(lgq)*phi*lambda)*exp(lgK(-1))+exp(lgp)*exp(lgD))+beta*(exp(lgdbond)+exp(lgq)*psi*exp(lgbond(-1)))
      -(1-beta)*(1-phi)*lambda*exp(lgq_r)*exp(lgK(-1)))-(1-beta)*(exp(lgdbond))*(1-psi)*exp(lgq_r)*exp(lgbond(-1));
    exp(lgI)     = exp(lgK)-lambda*exp(lgK(-1));
    exp(lgL)     = (exp(lgw)/om)^(1/mu);
    (exp(lgw)/om)^(1/mu) = exp(lgK(-1))*(((1-gamma)*A)/exp(lgw))^(1/gamma);
    exp(lgR)     = exp(lgAA)*exp(lgK(-1))^(alpha-1);
    exp(lgr)     = exp(lgR)-(1-lambda);
    exp(lgY)     = A*(exp(lgK(-1))^(gamma))*(exp(lgL)^(1-gamma));
    exp(lgY)     = exp(lgI)+exp(lgC);
    exp(lgC)     = exp(lgC_e)+exp(lgC_w);
    exp(lgC_w)   = exp(lgw)*exp(lgL);
    exp(lgC_e)   = (1-beta)*((exp(lgR)+exp(lgq)*lambda)*exp(lgK(-1))+exp(lgp)*exp(lgD)+pi*(1-phi)*(exp(lgq_r)-exp(lgq))*lambda*exp(lgK(-1)))+pi*(1-psi)*(exp(lgq_r)-exp(lgq))*exp(lgbond(-1));
    exp(lgN_s)   = theta*exp(lgI)+phi*pi*lambda*exp(lgK(-1))+(1-pi)*lambda*exp(lgK(-1));
    log(A)       = rho_A*log(A(-1))+e_A;
    exp(lgAA)    = gamma*((1-gamma)/om)^((1-gamma)/(gamma+mu))*A^((1+mu)/(gamma+mu));
    exp(lgD)     = 1;
    exp(lgY_L)   = exp(lgY)/exp(lgL);
    exp(lgq_r)   = (1-theta*exp(lgq))/(1-theta);
    exp(lgbond)  = exp(lgI)+psi*pi*exp(lgK(-1))+(1-pi)*exp(lgK(-1));
    exp(lgdbond) = 0.01;
    exp(lgpbond) = exp(lgpbond(-1))*exp(lgdbond(+1))*exp(lgq(+1))*exp(psi-1)*exp(beta*pi);
end;

// INITIAL VALUES:
initval;
A = 1;
lgAA = -2.1628;
lgD = 0;
lgq = 0.0681;
lgq_r = -0.0125;
lgR = -3.3264;
lgK = 2.3254;
lgw = 0.7068;
lgL = -1.0131;
lgY = 0.0986;
lgY_L = 1.1117;
lgI = -1.3635;
lgp = -1.1085;
lgC_w = -0.3064;
lgC_e = -2.1916;
lgC = -0.1651;
lgN_s = 2.0258;
lgbond = 2.0;
lgdbond = 0.001;
lgpbond = 0.001;
end;
steady;
check;

// STOCHASTIC SHOCKS:
shocks;
    var e_A   = sigma_e_A^2;
end;

// SIMULATION ORDER:
stoch_simul(order=1,nocorr,nomoments,hp_filter=1600,IRF=0);


//-------------------------------------------------------------------------------------------------//
@#else
//-------------------------------------------------------------------------------------------------//


// VARIABLES:
var lgY, lgY_L, lgI, lgK, lgC, lgC_w, lgC_e, lgL, lgw, lgR, lgr, lgq, lgq_r, lgp, lgD, A, lgN_s, lgbond, lgdbond,
        lgAA, y, c, i, l, y_l, w, r, a, q, p;
varexo e_A;

// PARAMETERS:
parameters pi, psi, beta, rho_A, om, lambda, theta, gamma, mu, alpha, phi,
            lgYss, lgCss, lgIss, lgKss, lgLss, lgY_Lss, lgwss, lgrss, lgqss, lgpss, Ass;
    pi        = 0.25;
    phi       = 0.025;
    psi       = 0.050;
    theta     = 0.25;
    // calibration as in Gertler and Kiyotaki 2009:
    beta      = 0.984;
    om        = 6.2;
    mu        = 0.6;
    lambda    = 0.975;
    gamma     = 0.333;
    rho_A     = 0.979;
    sigma_e_A = 0.0072;
    alpha     = gamma*((1+mu)/(gamma+mu));
    load shockA_KM_RBCss;
    set_param_value('lgYss',lgYss);
    set_param_value('lgCss',lgCss);
    set_param_value('lgIss',lgIss);
    set_param_value('lgKss',lgKss);
    set_param_value('lgLss',lgLss);
    set_param_value('lgY_Lss',lgY_Lss);
    set_param_value('lgwss',lgwss);
    set_param_value('lgrss',lgrss);
    set_param_value('lgqss',lgqss);
    set_param_value('lgpss',lgpss);
    set_param_value('Ass',Ass);

// MODEL SPECIFICATION:
model;
    (1-pi)*(((exp(lgR(+1))+exp(lgq(+1))*lambda)/(exp(lgq)))-(exp(lgp(+1))/exp(lgp)))
      *(1/((exp(lgR(+1))+exp(lgq(+1))*lambda)*exp(lgN_s)+exp(lgp(+1))*exp(lgD)))
    = pi*((exp(lgp(+1))/exp(lgp))-(1/exp(lgq))*(exp(lgR(+1))+lambda*(exp(lgq(+1))*phi+(1-phi)*exp(lgq_r(+1))))) 
      *(1/((exp(lgR(+1))+lambda*(exp(lgq(+1))*phi+(1-phi)*exp(lgq_r(+1))))*exp(lgN_s)+exp(lgp(+1))*exp(lgD)));
    (1-theta*exp(lgq))*exp(lgI) = pi*(beta*((exp(lgR)+exp(lgq)*phi*lambda)*exp(lgK(-1))+exp(lgp)*exp(lgD))
      -(1-beta)*(1-phi)*lambda*exp(lgq_r)*exp(lgK(-1)));
    exp(lgI)     = exp(lgK)-lambda*exp(lgK(-1));
    exp(lgL)     = (exp(lgw)/om)^(1/mu);
    (exp(lgw)/om)^(1/mu) = exp(lgK(-1))*((1-gamma)*A/exp(lgw))^(1/gamma);
    exp(lgR)     = exp(lgAA)*exp(lgK(-1))^(alpha-1);
    exp(lgr)     = exp(lgR)-(1-lambda);
    exp(lgY)     = A*(exp(lgK(-1))^(gamma))*(exp(lgL)^(1-gamma));
    exp(lgY)     = exp(lgI)+exp(lgC);
    exp(lgC)     = exp(lgC_e)+exp(lgC_w);
    exp(lgC_w)   = exp(lgw)*exp(lgL);
    exp(lgC_e)   = (1-beta)*((exp(lgR)+exp(lgq)*lambda)*exp(lgK(-1))+exp(lgp)*exp(lgD)+pi*(1-phi)*(exp(lgq_r)-exp(lgq))*lambda*exp(lgK(-1)));
    exp(lgN_s)   = theta*exp(lgI)+phi*pi*lambda*exp(lgK(-1))+(1-pi)*lambda*exp(lgK(-1));
    log(A)       = rho_A*log(A(-1))+e_A;
    exp(lgAA)    = gamma*((1-gamma)/om)^((1-gamma)/(gamma+mu))*A^((1+mu)/(gamma+mu));
    exp(lgD)     = 1;
    exp(lgY_L)   = exp(lgY)/exp(lgL);
    exp(lgq_r)   = (1-theta*exp(lgq))/(1-theta);

@#if (irfs==0)

    r   = (exp(lgr))*100;
    y   = (lgY-lgYss)*100;
    c   = (lgC-lgCss)*100;
    i   = (lgI-lgIss)*100;
    l   = (lgL-lgLss)*100;
    y_l = (lgY_L-lgY_Lss)*100;
    a   = (log(A/Ass))*100;
    w   = (lgw-lgwss)*100;
    q   = (lgq-lgqss)*100;
    p   = (lgp-lgpss)*100;

@#else

    r   = (exp(lgr));
    y   = (lgY-lgYss);
    c   = (lgC-lgCss);
    i   = (lgI-lgIss);
    l   = (lgL-lgLss);
    y_l = (lgY_L-lgY_Lss);
    a   = (log(A/Ass));
    w   = (lgw-lgwss);
    q   = (lgq-lgqss);
    p   = (lgp-lgpss);

@#endif

end;

// INITIAL VALUES:
initval;
lgY=lgYss;
lgY_L=lgY_Lss;
lgI=lgIss;
lgK=lgKss;
lgC=lgCss;
lgL=lgLss;
lgw=lgwss;
lgr=lgrss;
lgq=lgqss;
lgp=lgpss;
A=Ass;
y=0;
c=0;
i=0;
l=0;
y_l=0;
w=0;
r=0;
a=0;
end;
steady;
check;

// STOCHASTIC SHOCKS:
shocks;

@#if (irfs==0)

    var e_A; stderr sigma_e_A;

@#else

    var e_A; stderr 0.01;

@#endif

end;

// SIMULATION ORDER:

@#if (irfs==0)

    stoch_simul(order=1,hp_filter=1600,IRF=0)
    y c i l y_l w r a q p;

@#else

    stoch_simul(order=1,nocorr,nomoments,hp_filter=1600,IRF=200)
    y c i l y_l w r a q p;

@#endif

@#endif
