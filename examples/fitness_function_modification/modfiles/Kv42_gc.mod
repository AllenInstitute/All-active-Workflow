TITLE Kv4.2 with auxilliary subunits

COMMENT
This is the model of Barghaan et al. 2008: Role of N-terminal domain and accessory subunits in controlling deactivation-inactivation coupling of Kv4.2 channels.. Biophysical Journal
: implemented by M.Beining; Beining et al (2016), "A novel comprehensive and consistent electrophysiologcal model of dentate granule cells"

ENDCOMMENT

NEURON {
	SUFFIX Kv42
	USEION k READ ek WRITE ik
    RANGE  ik, gk, gkbar
	GLOBAL 	f,	a0,za,b0,zb,kco0,zco,koc0,zoc,kci,kic,koi,kio,koi2,kio2 , vshift
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
    (S)  = (siemens)
	
	(molar) = (1/liter)
	(mM) = (millimolar)
	(uM) = (micromolar)
	FARADAY = (faraday) (kilocoulombs)
	R = (k-mole) (joule/degC)
}

PARAMETER {
	v 		(mV)
	gkbar  = 0.00435                (S/cm2) : to be fitted     	
	f = 0.3  : The allosteric factor f defines the voltage-dependent coupling between activation and inactivation pathways for closed channels. deactivation
	a0 = 0.416				(/ms)
	za = 1.1				
	b0 = 0.0486				(/ms)
	zb = -1.556		
	kco0 = 0.347			(/ms)
	zco = 0
	koc0 = 1.67				(/ms)
	zoc = 0		
	kci = 0.04873			(/ms)
	kic = 0.000537			(/ms)
	koi = 0.0669			(/ms)
	kio = 0.04246			(/ms)
	koi2 = 0.0669			(/ms): is k56
	kio2 = 0.04246			(/ms): is k65
	vshift = -20			(mV)
}

STATE {
        O C0 C1 C2 C3 C4 I0 I1 I2 I3 I4  IO1 IO2 :I5 C5
}

ASSIGNED {
        ik                             (mA/cm2)
        gk                            (S/cm2)
        ek                            (mV)
		alpha   					(/ms)
		beta						(/ms)
		kco   					(/ms)
		koc						(/ms)
		celsius (degC)
}

INITIAL {
	rate(v)
	SOLVE kin STEADYSTATE sparse
}

BREAKPOINT {
	SOLVE kin METHOD sparse
	gk = gkbar * O
    ik = gk * ( v - ek )
}


KINETIC kin {
rate(v)

~ I0 <-> I1 (4*alpha/f,beta*f)
~ I1 <-> I2 (3*alpha/f,2*beta*f)
~ I2 <-> I3 (2*alpha/f,3*beta*f)
~ I3 <-> I4 (alpha/f,4*beta*f)

~ C0 <-> I0 (kci*f^4,kic/f^4)
~ C1 <-> I1 (kci*f^3,kic/f^3)
~ C2 <-> I2 (kci*f^2,kic/f^2)
~ C3 <-> I3 (kci*f,kic/f)
~ C4 <-> I4 (kci,kic)

~ C0 <-> C1 (4*alpha,beta)
~ C1 <-> C2 (3*alpha,2*beta)
~ C2 <-> C3 (2*alpha,3*beta)
~ C3 <-> C4 (alpha,4*beta)

~ O <-> IO1 (koi,kio) :kappa1,lambda1
~ IO1 <-> IO2 (koi2,kio2) :kappa2,lambda2  : is k56 and k65

~ C4 <-> O (kco,koc)

CONSERVE O + C0 + C1 + C2 + C3 + C4  + I0 + I1 + I2 + I3 + I4  + IO1 + IO2  = 1  :+ I5 + C5
}


PROCEDURE rate(v (mV)) { :callable from hoc
	alpha = exponential(a0,za,v,vshift)
	beta = exponential(b0,zb,v,vshift)
	kco = exponential(kco0,zco,v,vshift)
	koc = exponential(koc0,zoc,v,vshift)
}


FUNCTION exponential(A(/ms), z , v (mV), D (mV)) (/ms) {
	exponential = A* exp(z*(v-D)*FARADAY/(R*(celsius+273.15))) : V / (RT/F) = V * F / RT
}