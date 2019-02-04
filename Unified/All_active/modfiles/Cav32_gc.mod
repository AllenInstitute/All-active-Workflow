:Model from Burgess et al., 2002

NEURON {
    SUFFIX Cav32
    USEION ca READ eca WRITE ica
	USEION tca WRITE itca VALENCE 0
    RANGE g, gbar
	RANGE     kc1c2, kc2c1 ,  kc2c3  ,  kc3c2,  kc3o, koc3 ,k1i2,ki2i1,ki3i2
}

UNITS { 
	(mV) = (millivolt) 
	(S) = (siemens)
	FARADAY = (faraday) (kilocoulombs)
	R = (k-mole) (joule/degC)
}

: initialize parameters

PARAMETER {
    gbar = 2.2e-05 (S/cm2)

    kci = 0.0006  (/ms)
	kic = 0.0002  (/ms)
	
    kci2 = 0.0034  (/ms)
    kic2 = 0.0007  (/ms)
	
	kci3 = 0.058  (/ms)
	kic3 = 0.00008  (/ms)
    
}

ASSIGNED {
    v    (mV)
    eca  (mV)
    g    (S/cm2)
    ica  (milliamp/cm2)
    kc1c2   (/ms)
    kc2c1   (/ms)
    kc2c3   (/ms)
    kc3c2   (/ms)
    kc3o   (/ms)
    koc3   (/ms)
	ki1i2    (/ms)
	ki2i1   (/ms)
	ki2i3    (/ms)
	ki3i2   (/ms)
	celsius (degC)
	itca (milliamp/cm2)
}

STATE { c1 c2 c3 i1 i2 i3 io o }

BREAKPOINT {
    SOLVE kin METHOD sparse
    g = gbar*o
    ica = g*(v - eca)
	itca = ica
}

INITIAL { SOLVE kin STEADYSTATE sparse }

KINETIC kin {
    rates(v)
    ~ c1 <-> c2 (kc1c2, kc2c1)
    ~ c2 <-> c3 (kc2c3, kc3c2)
    ~ c3 <-> o (kc3o, koc3)
    ~ i1 <-> i2 (ki1i2, ki2i1)
    ~ i2 <-> i3 (ki2i3, ki3i2)
    ~ i3 <-> io (kc3o, koc3)
    ~ i1 <-> c1 (kic, kci)
    ~ i2 <-> c2 (kic2, kci2)
    ~ i3 <-> c3 (kic3, kci3)
    ~ io <-> o  (kic3, kci3)
    CONSERVE c1 + c2 + c3 + i1 + i2 + i3 + io + o = 1
}

PROCEDURE rates(v(millivolt)) {

    kc1c2 = alpha(1.6,0.72,1.82,v)
    kc2c1 = beta(0.032,0.72,1.82,v)
    
    kc2c3 = alpha(41,0.31,4.5,v)
    kc3c2 = beta(0.027,0.31,4.5,v)
    
    kc3o = 0.42
    koc3 = beta(0.015,0,0.7,v)  
    
	ki1i2 = alpha(0.032,0.72,1.82,v)
	ki2i1 = beta(0.0004,0.72,1.82,v)
	
	ki2i3 = alpha(0.25,0.31,4.5,v)
	ki3i2 = beta(0.0000011,0.31,4.5,v)
}

FUNCTION alpha(A (/ms), d, q ,  v (mV)) (/ms) {
	
	alpha = A*exp(d*q*v*FARADAY/(R*(celsius+273.15)))
}

FUNCTION beta(A (/ms), d, q ,  v (mV)) (/ms) {
	beta = A*exp(-(1-d)*q*v*FARADAY/(R*(celsius+273.15)))
}
