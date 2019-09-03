: Kv2.1
: from Beining et al (2016), "A novel comprehensive and consistent electrophysiologcal model of dentate granule cells"

NEURON	{
	SUFFIX Kv21
	USEION k READ ek WRITE ik
	RANGE gkbar, gk, ik
}

UNITS	{
	(S) = (siemens)
	(mV) = (millivolt)
	(mA) = (milliamp)
}

PARAMETER	{
	gkbar = 0.00001 (S/cm2) 
}

ASSIGNED	{
	v	(mV)
	ek	(mV)
	ik	(mA/cm2)
	gk	(S/cm2)
	alpham (ms-1)
	betam (ms-1)
	alphah (ms-1)
	betah (ms-1)
	hInf
	hTau
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
	m' = alpham*(1-m)-betam*m
	h' = alphah*(1-h)-betah*h
}

INITIAL{
	rates()
	m = alpham/(alpham+betam)
	h = alphah/(alphah+betah)
}

PROCEDURE rates(){ LOCAL x
	UNITSOFF 
		

			x = 0.1326 *(v-(11.1945))
			if (fabs(x) > 1e-6) {
				alpham = 0.0324*x/(1-exp(-x))
			}else{
				alpham = 0.0324/(1-0.5*x)
			}
			betam = 0.1302 * exp(-0.0302*(v -(-86.5604)))
		
		alphah = 0.0000151 * exp(-0.0537*(v -(-16.3280 )))
		betah = 0.0002/(1+exp((-0.2037*(v -(-22.5445))))) 
	UNITSON
}