
:Comment :
: from Beining et al (2016), "A novel comprehensive and consistent electrophysiologcal model of dentate granule cells"

NEURON	{
	SUFFIX Kv34
	USEION k READ ek WRITE ik
	RANGE gkbar, gk, ik
	GLOBAL scale_a, Rinact, ksl, vshift,ak,ad
}

UNITS	{
	(S) = (siemens)
	(mV) = (millivolt)
	(mA) = (milliamp)
}

PARAMETER	{
	gkbar = 0.0307813 (S/cm2) 
	vshift = 0 (mV)
	vshifttau1 = 0 (mV)
	vshifttau2 = 0 (mV)
	Rinact = 0.1
	kmg = 16
	scale_a = 4
	ksl = 0.5
	ak =  -9.7 (mV)
	ad =  14 (mV)
}

ASSIGNED	{
	v	(mV)
	ek	(mV)
	ik	(mA/cm2)
	gk	(S/cm2)
	mInf
	mTau (ms)
	hInf
	hTau (ms)
	am (/ms)
	bm (/ms)
}

STATE	{ 
	m
	h
}

BREAKPOINT	{
	SOLVE states METHOD cnexp
	gk = gkbar*m*h
	ik = gk*(v-ek)
}

DERIVATIVE states	{
	rates()
	m' = (mInf-m)/mTau
	h' = (hInf-h)/hTau
}

INITIAL{
	rates()
	m = mInf
	h = hInf
}

PROCEDURE rates(){
	UNITSOFF 
	:mInf =	1/(1+exp(((v -(10 ))/(-5))) +  exp(((v -(45+vshift ))/(kmg))))	  :this is taken from Schroeter (+ mg block)
	
	mInf =	1/(1+exp(((v -(ad+vshift))/(ak)))) : This is taken from Schroeter.. is also perfect in between Kv3.3 and 3.4 (Rudy Review)
	
	::mInf =	1/(1+exp(((v -(5))/(-5))) +  exp(((v -(55))/(kmg))))	:  this is taken for the model
	
	
	: mTau from Rudy 
	am = scale_a* 1/16 * exp(0.1*ksl*(v-38))
	bm = scale_a *1/16 * exp(-0.1*ksl*(v+45))  : -0.16...v+27
	mTau = 1/(am + bm )	:


	
	hInf = Rinact + (1-Rinact)/(1+exp(((v -(-29.7 ))/(12.2))))   : leicht abgewandelt von Schroeter
	hTau = 250/(1+exp(((v -(-10))/(17))))+ 8	: hTau Rudy3

	UNITSON
}

