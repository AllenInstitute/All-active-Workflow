: M conductance
: from Mateos-Aparicio et al (2014)

NEURON {
	SUFFIX Kv723
	USEION k READ ek WRITE ik
	RANGE gkbar, ik, minf, tau1, tau2, i, gk, m1, m2, tadjtau
	GLOBAL Vhalf, Vshift, k, v0erev, kV, gamma
	GLOBAL Dtaumult1, Dtaumult2, tau0mult, taudiv, q10tau
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(S) = (siemens)
	(um) = (micron)
} 

PARAMETER {
	gkbar = 0.0067 			(S/cm2)
	k = 9           		(mV)
	Vhalf = -50             (mV)  :for minf(V)
	Vshift = 0              (mV)	:for g(V) and minf(V)     
	v0erev = 65             (mV)     :50-80
	kV = 40                 (mV)     
	gamma = 0.5                      :0.5,1

	temptau = 22	          (degC) :tau reference temperature 	
	q10tau  = 5
	taudiv = 1
	Dtaumult1 = 6  
	Dtaumult2 = 6 
	tau0mult = 0.2

	vmin = -100	            (mV)
	vmax = 100	            (mV)
	temp0 = 273		          (degC)
	FoverR = 11.6045039552	(degC/mV)
} 
 
ASSIGNED {
	ek		(mV)
	v 	     	(mV)
	celsius		(degC)
	Vhalf1    (mV) 
	Dtau1     (ms)
	z1               
	tau01   	(ms)	 
	Vhalf2  	(mV)	  
	Dtau2   	(ms)  
	z2               
	tau02   	(ms)	  
	alpha1				  
	beta1	  		  
	alpha2		
	beta2	
	i 	    	(mA/cm2)
	ik 	     	(mA/cm2)
	gk		      (S/cm2)
	minf
	v0        (mV)      
	tau1			(ms)
	tau2			(ms)
	tadjtau
	frt		    (/mV)
}
 
STATE { m1 m2 }

INITIAL { 
	rates(v)
	m1 = minf
	m2 = minf
}

BREAKPOINT {
  SOLVE states METHOD cnexp
	gk = gkbar*gsat(v)*(m1^2)*m2
	ik = gk*(v - ek)
	i = ik
} 

DERIVATIVE states {
	rates(v)
	m1' = (minf - m1)/tau1
	m2' = (minf - m2)/tau2
}

PROCEDURE rates(v (mV)) {
  TABLE minf, tau1, tau2
	DEPEND celsius, gamma, k, Vhalf, Vshift, taudiv, Dtaumult1, Dtaumult2, tau0mult
	FROM vmin TO vmax WITH 199
	
  IF (gamma == 0.5) {
  	z1 = 2.8
		Vhalf1 = -49.8+Vshift 	:(mV)  shifted - 20 mV (when Vshift = 0)
		tau01 = 20.7 (ms) *tau0mult	  :(ms)
		Dtau1 = 176.1 (ms) *Dtaumult1	:(ms)
		z2 = 8.9	              
		Vhalf2 = -55.5+Vshift 					:(mV)  shifted - 20 mV
		tau02 = 149 (ms) *tau0mult   					:(ms)
		Dtau2 = 1473 (ms) *Dtaumult2 	  			:(ms)
	}	
	IF (gamma == 1) {
  	z1 = 3.6
		Vhalf1 = -25.3+Vshift		:(mV)  shifted - 20 mV
		tau01 = 29.2 (ms) *tau0mult	  :(ms)
		Dtau1 = 74.6 (ms) *Dtaumult1	:(ms)
		z2 = 9.8	
		Vhalf2 = -44.7+Vshift 					:(mV)  shifted - 20 mV
		tau02 = 155 (ms) *tau0mult   					:(ms)
		Dtau2 = 549 (ms) *Dtaumult2  	  			:(ms)
	}
  tadjtau = q10tau^((celsius - temptau)/10 (degC))
	frt = FoverR/(temp0 + celsius)

  alpha1 = exp(z1*gamma*frt*(v - Vhalf1))
  beta1 = exp(-z1*(1-gamma)*frt*(v - Vhalf1))
  tau1 = (Dtau1/(alpha1 + beta1) + tau01)/(tadjtau*taudiv)
  
  alpha2 = exp(z2*gamma*frt*(v - Vhalf2))
  beta2 = exp(-z2*(1-gamma)*frt*(v - Vhalf2))
  tau2 = (Dtau2/(alpha2 + beta2) + tau02)/(tadjtau*taudiv)

  minf = 1/(1 + exp(-(v - Vhalf - Vshift)/k))
}

FUNCTION gsat (v (mV)) {
	gsat = 1
	v0 = v0erev + ek  
	IF (v > v0) {
		gsat = 1+(v0-v+kV*(1-exp(-(v-v0)/kV)))/(v-ek)
	}
}


