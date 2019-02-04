TITLE I-h channel 

COMMENT

Start from  Magee 1998 for distal dendrites
Adaped from Li & Ascoli 2006,2008 
Adapted by A. Hanuschkin 2011
Adapted by M. Beining 2016 (added slow component and cAMP dependence with data  from Chen et al 2001 Journal of General Physiology)
ENDCOMMENT




NEURON {
	SUFFIX HCN
	NONSPECIFIC_CURRENT i
    RANGE gbar, g, taul, linf
    GLOBAL vhalfc, e, vhalfl, kl, vhalft, at, bt, q10, qfact, cAMP, qt
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
    (S)  = (siemens)
	(molar) = (1/liter)
	(uM) = (micromolar)		
}

PARAMETER {
	v 		(mV)
    e = -41.9	(mV)    : exp data		(see manuscript)
	gbar=4e-06 	(S/cm2) : dummy default

				: Boltzman fit to steady state currents
        vhalfl=-75.3	(mV)	: fitted Boltzmann	(see manuscript) 
	kl=8			: fitted Boltzmann      (see manuscript) 
        
: Fit to activation/deactivation time constance
: own fit
: E               = 0.000392656      +/- 0.0001422    (36.22%)
: F               = 0.395751         +/- 0.1126       (28.45%)
: G               = 29.385           +/- 2.715        (9.238%)
: manuscript fast, a = 0.52 +/- 0.21 ms, b = 215.1 +/- 65.7 ms, V0 = 30.4 +/- 3.4 mV, 
: E               = 0.00052
: F               = 0.2151        
: G               = 30.4      
 
    vhalft=30.4	 (mV)    : fitted 		(see manuscript)
    at=0.00052	 (/ms)   : fitted 		(see manuscript)
	bt=0.2151	 (/ms)	 :fitted               (see manuscript)

				: Temperature dependence
    celsius         (degC)  : unused
	q10=1.			: no correction for Temperature via q10 to save computation time (uncomment in rate to make T dep simulations... 
	qfact = 1 
	
	Ac = 14 (mV)
	kc 	= -0.35 (/uM)
	chalf 	= -0.8 (uM)
	cAMP = 0 (uM)
	
}

ASSIGNED {
	i (mA/cm2)
        linf      
        taul
        g (S/cm2)
		vhalfc (mV)  : voltage shift by cAMP (for HCN1+HCN2 heteromers)
		qt
}

INITIAL {
	qt = q10^((celsius-33)/10)
	rate(v)
	l1=linf
	l2=linf
}

STATE { l1 l2 }

BREAKPOINT {
	SOLVE states METHOD cnexp
	g = gbar*(0.8*l1 + 0.2*l2)  : slow component has contribution of ~ 20 % to total current
	i = g*(v-e)

}

DERIVATIVE states {     : exact when v held constant; integrates over dt step
        rate(v)
        l1' =  (linf - l1)/taul
		l2' = (linf - l2)/(taul*6.4)  : slow component has same time course as fast one but is slower by factor of ~6.4 (see Stegen&Hanuschkin 2012 Figure 3)
}

PROCEDURE rate(v (mV)) { :callable from hoc
        :LOCAL qt

	:qt = 1 	: saves computational time...
	
	if (cAMP <= 0) {
		vhalfc = 0
	}else{
		vhalfc = Ac/(1 + exp((log10(cAMP)-chalf)/kc))  : from Chen et al 2001 Journal of General Physiology
	}
    linf = 1/(1 + exp((v-vhalfl-vhalfc)/kl))
 	taul = 1/(qt * qfact * (at*exp(-v/vhalft) + bt*exp(v/vhalft) ))
}
