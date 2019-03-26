: Ca-dependent K channels (BK) - alphabeta4 and alpha
: Bin Wang, Robert Brenner, and David Jaffe - Originally written May 27, 2010
: 
: June 1, 2010 - added double exponential function for voltage-dependent activation 
:
: July 3, 2010 - changed voltage-dependence for the two channels based on revised data
:
: April 2, 2011 - adjusted parameters based on updated Bin data

: Mar 2016 - added instant N-type Calcium concentration (Beining et al (2016), "A novel comprehensive and consistent electrophysiologcal model of dentate granule cells")

NEURON {
	SUFFIX BK_gc
	USEION ca READ cai
	USEION k READ ek WRITE ik
	USEION nca READ ncai VALENCE 0
	RANGE gakbar, gabkbar, gak, gabk, atau, ainf, a, ab, abinf, abtau,  ik, diff, acai
	GLOBAL base
}

UNITS {
	(molar) = (1/liter)
	(mM) = (millimolar)
	(mV) = (millivolt)
	(mA) = (milliamp)
	(S) = (siemens)
}

PARAMETER {
	diff = 100000 (1)
	gakbar = 0.0624	(S/cm2)	: maximum permeability - alpha
	gabkbar = 0.0156	(S/cm2)	: maximum permeability - alphabeta
	base = 4  	(mV)	: alphabeta4 base time constant
	
}

ASSIGNED {
	v				(mV)
	ek			(mV)
	
	ik				(mA/cm2)
      gak			(S/cm2)
      gabk		(S/cm2)
      ainf	
      atau		(ms)
      abinf	
      abtau		(ms)
	  cai 			(mM)
	  ncai		(mM)
	  acai		(mM)

}

STATE {ab a} :a

BREAKPOINT {
	SOLVE state METHOD cnexp
	gak = gakbar*a
	gabk = gabkbar*ab
	ik = (gabk+gak)*(v - ek)                               :+ gak
}

DERIVATIVE state {
	: exact when v held constant; integrates over dt step
	
	acai = (ncai) / diff : instantaneous calcium concentration of N-type Ca channels for BK activation modified by diffusion factor for nanodomains

	if (acai < cai)
		{acai = cai}
		
	rates(v,acai)				      
	a' = (ainf-a)/atau	
	ab' = (abinf-ab)/abtau
}

INITIAL {
	rates(v, cai)
	a = ainf
	ab = abinf
}

: alpha channel properties
FUNCTION shifta(ca (mM))  {
	shifta = 25 - 50.3 + (107.5*exp(-.12*ca*1e3))
}


FUNCTION peaka(ca (mM))  {
	peaka = 2.9 + (6.3*exp(-.36*ca*1e3))
}

: alpha-beta4 channel properties


FUNCTION shiftab(ca (mM))  {
	shiftab = 25 - 55.7 + 136.9*exp(-.28*ca*1e3)
}


FUNCTION peakab(ca (mM))  {
	peakab = 13.7 + 234*exp(-.72*ca*1e3)
}

: Double sigmoid function for tau voltage-dependence


FUNCTION taufunc(v (mV)) {
	 taufunc = 1 / (          (10*(exp(-v/63.6) + exp (-(150-v)/63.6)))  - 5.2                  )
	 if (taufunc <= 0.2) {	  : stop the function between 0.2 and 1
	    taufunc = 0.2
	 }
}

PROCEDURE rates(v (mV), c (mM)) { : nc (mM), 
	  LOCAL range, vv,ashift, bshift

	  : alpha model

	  ashift =  -32 + (59.2*exp(-.09*c*1e3)) + (96.7*exp(-.47*c*1e3))
	  ainf = 1/(1+exp((ashift-v)/(25/1.6)))

	  vv = v + 100 - shifta(c)
	  atau = taufunc(vv)
	  range = peaka(c)-1
	  atau = (range*((atau-.2)/.8)) + 1

	  : alpha-beta4 model

	  bshift = -56.449 + 104.52*exp(-.22964*c*1e3) + 295.68*exp(-2.1571*c*1e3)

	  abinf = 1/(1+exp((bshift-v)/(25/1.6)))

	  vv = v + 100 - shiftab(c)
	  abtau = taufunc(vv)
	  range = peakab(c)-base
	  abtau = (range*((abtau-.2)/.8)) + base		

}
