: model from Evans et al 2013, transferred from GENESIS to NEURON by Beining et al (2016), "A novel comprehensive and consistent electrophysiologcal model of dentate granule cells"
: also added Calcium dependent inactivation

NEURON {
	SUFFIX Cav13
	USEION ca READ cai, eca WRITE ica   :,cai,cao...., cai, cao
	USEION lca WRITE ilca VALENCE 0
	RANGE gbar, g
	GLOBAL kf, h2Tau, VDI
}

UNITS {
	(molar) = (1/liter)
	(mM) = (millimolar)
	(mV) = (millivolt)
	(mA) = (milliamp)
	(S) = (siemens)
	(um) = (micrometer)
}

ASSIGNED {
	ilca		(mA/cm2) : instantaneous calcium current of l-type calcium channel
	v			(mV)
	ica		(mA/cm2)
	g		(S/cm2)
	eca 		(mV)
	diam		(um)
	cai 		(mM)
	mInf  (1)
	hInf  (1)
	h2Inf (1)
	mTau (ms)
}

PARAMETER {
	hTau 	= 44.3 		(ms)
	h2Tau = 0.5 		(ms)
	gbar = 4e-06	 (S/cm2)
		vshift = 0 		(mV)
		
		:parameters for calcium-dep inactivation (CDI) 
			:f= (0.001/(0.001+[Ca]))Poirazi CA1  2003
			:f= (0.0005/(0.0005+[Ca])) Rhodes and Llinas 2001 Cort Pyr
	kf		=			0.0005 (mM)  : factor in inactivation, the higher the less sensitive. others uses 0.0002.. standen and stanfield use 0.001mM in original paper	
	VDI = 1
}

STATE {m h h2}  :a b  :cai (mM) cao (mM)

INITIAL {
	rates()
	m = mInf
	h = hInf
	h2 = h2Inf
}

BREAKPOINT {
	rates()
	SOLVE state METHOD cnexp
	g = gbar*m*h*h2 : h2 calcium dependent inactivation is taken from santhakumar 05.. tjos assumes instantaneous calcium inactivation
	ica = (g)*(v - eca) : 
	ilca = ica
	
}

DERIVATIVE state {	: exact when v held constant integrates over dt step
	m' = (mInf-m) / mTau
	h' = (hInf-h) / hTau
	h2' = (h2Inf-h2)/h2Tau
}

PROCEDURE rates(){
		LOCAL mA,mB
		mA = (39800*( v + 67.24))/( exp ( (v + 67.24)/15.005) - 1.0)
		mB = 3500* exp(v/31.4) 
		mTau = (1/(mA + mB))
		
		mInf = 1.0/((exp ( (v - (-40.0))/(-5))) + 1.0)

		hInf = VDI/( (exp ( (v - (-37))/(5))) + 1.0) + (1-VDI)
		:h2 = caIn(cai)
		h2Inf = kf/(kf+cai)
}
