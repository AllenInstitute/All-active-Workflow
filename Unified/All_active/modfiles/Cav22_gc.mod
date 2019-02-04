: Ca channels (N-type) from Aradi and Holmes 1999, transferred from GENESIS to NEURON
: increased inactivation time constant and corrected calcium handling by Beining et al (2016), "A novel comprehensive and consistent electrophysiologcal model of dentate granule cells"


NEURON {
	SUFFIX Cav22
	USEION ca READ  eca WRITE ica 
	USEION nca WRITE inca VALENCE 0
	RANGE gbar, g
	GLOBAL vshift, hTau  :tau, depth, 
}

UNITS {
	(mV) = (millivolt)
	(mA) = (milliamp)
	(S) = (siemens)
}

ASSIGNED {
	inca		(mA/cm2) : instantaneous calcium current of n-type calcium channel
	v			(mV)
	ica		(mA/cm2)
	g 	(S/cm2)
	eca 		(mV)
	diam		(um)
	mInf  (1)
	hInf  (1)
	mTau (ms)
}

PARAMETER {
	gbar = 0.0003	(S/cm2)
	hTau = 80 (ms)
	vshift = 10 		(mV)  : recorrection of Jaffe  1994 compared to Fox 1987, as voltage-dependent activation curve should not depend on ion concentrations or type
}

STATE {m h} 

BREAKPOINT {
	rates()
	SOLVE state METHOD cnexp
	g = gbar*m*m*h
	ica = g*(v - eca)
	inca = ica
	
}

DERIVATIVE state {	: exact when v held constant; integrates over dt step
	m' = (mInf-m)/mTau
	h' = (hInf-h) / hTau
}

INITIAL {
	m = mInf
	h = hInf
}

PROCEDURE rates() { LOCAL alpha,beta
	alpha = f(1.9,0.1,v,19.88 + vshift) 
	beta = exponential(0.046,-0.048239,v, vshift)
	mTau = 1/(alpha+beta)
	mInf = alpha*mTau
	hInf = 1/(1+exp((v+40)/12.5)) 
}


FUNCTION f(A, k (/mV), v (mV), D (mV)) (/ms) {
	LOCAL x
	x = k*(v-D)
	if (fabs(x) > 1e-6) {
		f = A*x/(1-exp(-x))
	}else{
		f = A/(1-0.5*x)
	}
}

FUNCTION exponential(A, k (/mV), v (mV), D (mV)) (/ms) {
	exponential = A*exp(k*(v-D))
}
