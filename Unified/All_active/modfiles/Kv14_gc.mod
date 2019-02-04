TITLE Kv14.mod  - modified potassium channels
 
COMMENT

This is a Hodgkin-Huxley model of inactivating K channels 
modified from hh.mod
P. Jonas, 12 September 2001
temperature dependence is given up - Q10 scaling is not used any longer
Reference: Hines and Carnevale, Expanding NEURON's repertoire of mechanisms with NMODL, Neural Computation 12, 839-851, 2000
23 October 2001, new version for Faklers Kv1.4 wild type 
18 November 2001 - final version, as published in Wissmann et al., JBC 278, 16142-16150

ENDCOMMENT
 
UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
}
 
NEURON {
    THREADSAFE
    SUFFIX Kv14                                                                                           :   KIn   
    USEION k READ ek WRITE ik                                                            :   local K may affect ion flux, and ion fluc may affect local K
    RANGE gkbar, gk, scale_a, scale_i, ik                                                                                  :   functions of position  
    GLOBAL vshift
}
 
PARAMETER {
    gkbar = 0.001 (mho/cm2)	 <0,1e9>   :   so these parameter are viewed and can be changed in GUI 
    scale_a = 1.0
    scale_i = 1.0
	vshift = 0 (mV)
}
 
STATE {
    n h                                :  n for activation, h for inactivation 
}
 
ASSIGNED {                                   
    v (mV)                          :  variables given values outside the mod file
    celsius (degC)
    ek (mV)
    
    gk (mho/cm2)               :  variables that appear on left hand side of assigment statements within the mod file
    ik (milliamp/cm2)
    ninf
    ntau (ms)
    hinf
    htau (ms) 
}
 
 
BREAKPOINT {                                      : this block is responsible for making all variables consistent at time t 
    SOLVE states METHOD cnexp
    gk = gkbar*n*n*n*n*h
    ik = gk*(v - ek)      
}
 
 
INITIAL {
    rates(v)
    n = ninf
    h = hinf
}


DERIVATIVE states {  
    rates(v)
    n' = (ninf-n)/ntau
    h' = (hinf-h)/htau
}
 
PROCEDURE rates(v(mV)) {  :Computes rate and other constants at current v.
                                                    :Call once from HOC to initialize inf at resting v.
    LOCAL  alpha, beta, sum

UNITSOFF
    alpha = scale_a*.01*vtrap((-55+vshift-v),10) 
    beta = scale_a*.125*exp((-65+vshift-v)/80)
    sum = alpha + beta
    ntau = 1/sum
    ninf = alpha/sum
    
    alpha = scale_i*0.0000256077*exp((vshift-v)/45.4217)
    beta = scale_i*0.0330402/(exp((-45.6599+vshift-v)/2.30235) + 1) :Recombinant Kv1.4
    sum = alpha + beta
    htau = 1/sum
    hinf = alpha/sum
}
 
FUNCTION vtrap(x,y) {            :Traps for 0 in denominator of rate eqns., based on three terms of infinite series expansion of exp
    if (fabs(x/y) < 1e-6) {
        vtrap = y*(1 - x/y/2)
    } else {
        vtrap = x/(exp(x/y) - 1)
    }
}
 
UNITSON
