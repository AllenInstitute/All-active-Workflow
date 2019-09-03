: model from Evans et al 2013, transferred from GENESIS to NEURON by Beining et al (2016), "A novel comprehensive and consistent electrophysiologcal model of dentate granule cells"
: also added Calcium dependent inactivation

NEURON {
	SUFFIX Cav12
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
	hTau 	= 44.3  (ms)
	h2Tau = 0.5 	(ms)
	gbar = 4e-05	(S/cm2)
		vshift = 0 	(mV)
		
		:parameters for calcium-dep inactivation (CDI) 
			:f= (0.001/(0.001+[Ca]))Poirazi CA1  2003
			:f= (0.0005/(0.0005+[Ca])) Rhodes and Llinas 2001 Cort Pyr
	kf		=			0.0005 (mM)  : factor in inactivation, the higher the less sensitive. others uses 0.0002.. standen and stanfield use 0.001mM in original paper	
	VDI = 0.17
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
		mA = 39800*(v + 8.124)/(exp((v + 8.124)/9.005) - 1)
		mB = 990*exp(v/31.4)
		mTau = 1/(mA + mB) 

		mInf   = 1/(1 + exp((v + 8.9)/(-6.7)))
		
		hInf   = VDI/(1 + exp((v +55)/8)) + (1-VDI)
		h2Inf = kf/(kf+cai)
}
