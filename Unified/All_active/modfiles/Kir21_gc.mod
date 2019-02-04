TITLE Kir potassium current

COMMENT

Kir 2.1 (Mg high-affinity) model
from Beining et al (2016), "A novel comprehensive and consistent electrophysiologcal model of dentate granule cells"

based on
Yan & Ishihara (2005): Two Kir2.1 channel populations with different sensitivities to Mg(2+) and polyamine block: a model for the cardiac strong inward rectifier K(+) channel. , Journal of physiology
and Liu 2012

ENDCOMMENT

NEURON {
	SUFFIX Kir21
	USEION k READ ek WRITE ik
    RANGE  ik, gk, gkbar:, O, BS, B1, B2, B3
	GLOBAL mg_i, As, shiftmg, cas,fac, gsub, b, spm_i, vshiftbs
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
        (S)  = (siemens)
	
	(molar) = (1/liter)
	(mM) = (millimolar)
	(uM) = (micromolar)
	
}

PARAMETER {
	v 		(mV)
	gkbar  = 0.000353903                (S/cm2) :     	
	mg_i = 4 (mM)  : in Mongiat 2009
	spm_i = 1 (uM) : intracellular polyamine concentration (Yan&Ishihara 2005) Liu 2012 says physiologic is 5-10uM bzw 0.1-1uM spmd
	As = 0.2
	vshiftbs = 0 (mV)
	b= 0.105  : close to 0 makes tau big and shifts to right, b=0.1099 makes boltzmann tau to the right
	:c = -100 (mV)   : seems plausible
	
	fac = 0.005  : this influences tau a lot! make it smaller for bigger tau
	gsub = 0.25  : factor of sub state conductance 0.05-0.055 fuer spermin und 0.15-0.155 fuer spermidin
	shiftmg = 0.5 : 0 for normal 1 for shift to ek
	cas = 0.142857  
}

STATE {
        O BS B1 B2 B3 BB
}

ASSIGNED {
        : ki                              (mM)
        : ko                              (mM)
        ik                             (mA/cm2)
        gk                            (S/cm2)
        ek                            (mV)
		alpha1   					(/ms)
		beta1						(/ms)
		alphas   					(/ms)
		betas   					(/ms)
		alphas2 					(/ms)
		betas2						(/ms)
}

INITIAL {
	rate(v)
	SOLVE kin STEADYSTATE sparse
}

BREAKPOINT {
	SOLVE kin METHOD sparse
	gk = (gkbar ) * (O + 1/3 * B2 + 2/3 * B1) + (gkbar * gsub ) * BS 
    ik = gk * ( v - ek )
}


KINETIC kin {
LOCAL alpha2, alpha3, beta2, beta3
rate(v)

alpha2 = 2*alpha1
beta2 = 2 * beta1
alpha3 = 3*alpha1
beta3 = 3*beta1


~ BS <-> O (alphas,betas)
~ B1 <-> O (alpha1,beta3)
~ B2 <-> B1 (alpha2,beta2)
~ B3 <-> B2 (alpha3,beta1)
~ BB <-> BS (alphas2,betas2)

CONSERVE O + BS + BB + B1 + B2 + B3 = 1

}


PROCEDURE rate(v (mV)) { :callable from hoc
	LOCAL a,d
	
	: Mg block
	alpha1 = 12 * exp(-0.025 * (v - (shiftmg * (ek))))			: this is exactly as in paper
	beta1 = mg_i/8 * 28 * exp(0.025 * (v - (shiftmg * (ek))) ) : this is exactly as in paper  
	
	: high-affinity polyamine block
	alphas = As * 0.17 * exp(cas*-0.07 * (v - (ek) +8/8 (mV/mM) * mg_i)) :/  (1 + 0.01 * exp(0.12 * (v - (ek+vshiftbs)  +8  (mV/mM) * mg_i)))   : this is exactly as in paper, except denominator was omitted because it did not change anything in kinetics
	
	betas =  As * spm_i * 0.28 * exp(0.15 * (v - (ek) +8/8  (mV/mM) * mg_i)) :/ (1 + 0.01 * exp(0.13 * (v - (ek+vshiftbs)  + 8  (mV/mM) * mg_i)))    : this is exactly as in paper, except denominator was omitted because it did not change anything in kinetics
	
	: this is to fit two rate functions to the  kd of the paper kdd = 40 .* exp( - (v - (ek+vshiftbs)) / 9.1) 
	: b zwischen 0 und 1/9.1 ( 0.1099)

	a = - 1/9.1 + b
	
	:d = (ek  - c)/(9.1 * b) + c  : d reduces to ek if c is ek

	: low-affinity (second) polyamine block
	alphas2 = fac* 40 * exp(a*(v-(ek+vshiftbs)))  : formerly v-c..... formula is turned around compared to matlab!
	betas2 = spm_i * fac * exp(b*(v-(ek+vshiftbs))) : formerly v-d

}