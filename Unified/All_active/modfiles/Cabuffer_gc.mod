: Calcium buffer shell model with instant Calcium handling for clustered channels by Beining et al (2016), "A novel comprehensive and consistent electrophysiologcal model of dentate granule cells"

NEURON {
	SUFFIX Cabuffer
	USEION ca READ ica,cai,cao WRITE cai, cao
	USEION nca READ inca WRITE ncai VALENCE 0
	USEION lca READ ilca WRITE lcai VALENCE 0
	USEION tca READ itca WRITE tcai VALENCE 0
	GLOBAL depth,cao0,cai0
	RANGE cai,cao,ncai,lcai,brat, tau
}

UNITS {
	(molar) = (1/liter)
	(mM) = (millimolar)
	(mV) = (millivolt)
	(mA) = (milliamp)
	(S) = (siemens)
	(um) = (micrometer)
}

PARAMETER {
	tau = 240			(ms)
	depth = 0.05 		(um)
	cai0 = 4.8e-05 		(mM)	
	cao0 = 2  			(mM)	
	Fa = 96485.3365 (coulomb)
	brat = 50  : binding ratio by buffer
}

ASSIGNED {
	ica		(mA/cm2)
	diam	(um)
	VSR (um)
	ncai		(mM)
	inca		(mA/cm2) : instantaneous calcium current of n-type calcium channel
	lcai		(mM)
	ilca		(mA/cm2) : instantaneous calcium current of l-type calcium channel
	tcai		(mM)
	itca		(mA/cm2) : instantaneous calcium current of t-type calcium channel
	B 			(mM*cm2/mA)
}

STATE { 
cai (mM) 		<1e-5> 
cao (mM)}

BREAKPOINT {
	SOLVE state METHOD cnexp
}

DERIVATIVE state {	: exact when v held constant; integrates over dt step
	ncai = - inca * B  : instantaneous calcium concentration of N-type Ca channels for BK activation times sensitivity factor of BK
	lcai = - ilca * B : instantaneous calcium concentration of N-type Ca channels for BK activation times sensitivity factor of BK
	tcai = - itca * B  : instantaneous calcium concentration of N-type Ca channels for BK activation times sensitivity factor of BK
	cai' = -ica * B / brat -(cai-cai0)/tau	:(1e4)/(2*Fa*0.25*diam)	 *1e4 is for correction from um to cm to result in /cm3.. *1000 for /liter solves with /1000 for /ms ..  ja 1e4 es stimmt wirklich..june 2016
	cao' = 0
}

INITIAL {
	if (2*depth >= diam) {
		VSR = 0.25*diam : if diameter gets less than double the depth, the surface to volume ratio (here volume to surface ratio VSR) cannot be less than 1/(0.25diam) (instead of diam/(d*(diam-d)) )
	}else{
		VSR = depth*(1-depth/diam)
	}
	B = (1e4)/(2*Fa*VSR)
	cao0 	= 		cao
	cai0 	= 		cai
}