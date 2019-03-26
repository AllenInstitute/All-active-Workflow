TITLE SK2 multi-state model Cerebellum Golgi Cell Model

COMMENT

Author:Sergio Solinas, Lia Forti, Egidio DAngelo
Based on data from: Hirschberg, Maylie, Adelman, Marrion J Gen Physiol 1998

Mar 2016 - added instant L-type Calcium concentration (Beining et al (2016), "A novel comprehensive and consistent electrophysiologcal model of dentate granule cells")


Published in:
             Sergio M. Solinas, Lia Forti, Elisabetta Cesana, 
             Jonathan Mapelli, Erik De Schutter and Egidio D`Angelo (2008)
             Computational reconstruction of pacemaking and intrinsic 
             electroresponsiveness in cerebellar golgi cells
             Frontiers in Cellular Neuroscience 2:2
ENDCOMMENT

NEURON{
	SUFFIX SK2
	USEION ca READ cai
	USEION k READ ek WRITE ik 
	USEION lca READ lcai VALENCE 0
	RANGE gkbar, gk, ik, acai
	GLOBAL diff, Q10, fac, invc1,invc2,invc3,invo1,invo2,diro1,diro2,dirc2,dirc3,dirc4
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(molar) = (1/liter)
	(mM) = (millimolar)
}

PARAMETER {
	gkbar = 8.33333e-05 (mho/cm2)
	Q10 = 5 (1)
	diff = 3 (1) : diffusion factor, is also 3 in the original Solinas model
	fac = 2.5 : is 1 or 0.15 in the original model , dependent on whether using high or low open probability model
	ocontrt = 0
: rates ca-independent
	invc1 = 0.32  ( /ms)
	invc2 = 0.32  ( /ms)
	invc3 = 0.09 ( /ms)

	invo1 = 1.0     ( /ms)
	invo2 = 0.1     ( /ms)
	diro1 = 0.16 ( /ms)
	diro2 = 0.1    ( /ms)

: rates ca-dependent
	dirc2 = 200 ( /ms-mM )
	dirc3 = 160 ( /ms-mM )
	dirc4 = 320  ( /ms-mM )

}

ASSIGNED{ 
	v	(mV) 
	ek	(mV) 
	cai (mM)
	celsius		(degC)
	gk	(mho/cm2) 
	ik	(mA/cm2) 
	invc1_t  ( /ms)
	invc2_t  ( /ms)
	invc3_t  ( /ms)
	invo1_t  ( /ms)
	invo2_t  ( /ms)
	diro1_t  ( /ms)
	diro2_t  ( /ms)
	dirc2_t  ( /ms-mM)
	dirc3_t  ( /ms-mM)
	dirc4_t  ( /ms-mM)
	tcorr	 (1)

	dirc2_t_ca  ( /ms)
	dirc3_t_ca  ( /ms)
	dirc4_t_ca  ( /ms)
	
	lcai		(mM)
	acai     (mM)
	:ilca		(mA/cm2) : instantaneous calcium current of l-type calcium channel
	:diam	(um)
	:VSR (um)
} 

STATE {
	c1
	c2
	c3
	c4
	o1
	o2
}

BREAKPOINT{ 
	SOLVE kin METHOD sparse 
	gk = gkbar*(o1+o2)	:(mho/cm2)
	ik = gk*(v-ek)		:(mA/cm2)
} 

INITIAL{
	rate(celsius)
	SOLVE kin STEADYSTATE sparse
	
} 

KINETIC kin{ 
	acai =  (lcai)/diff : instantaneous calcium concentration of L-type Ca channels for SK activation divided by diffusion factor (determines proximity between both channels) of SK

	if (acai < cai)
		{acai = cai}
	rates(acai) 
	
	~c1<->c2 (dirc2_t_ca, invc1_t) 
	~c2<->c3 (dirc3_t_ca, invc2_t) 
	~c3<->c4 (dirc4_t_ca, invc3_t) 
	~c3<->o1 (diro1_t, invo1_t) 
	~c4<->o2 (diro2_t, invo2_t) 
	CONSERVE c1+c2+c3+c4+o2+o1=1 
} 

FUNCTION temper (Q10, celsius (degC)) {
	temper = Q10^((celsius -23(degC)) / 10(degC)) 
}

PROCEDURE rates(c(mM)){
	dirc2_t_ca = dirc2_t*c * fac
	dirc3_t_ca = dirc3_t*c * fac
	dirc4_t_ca = dirc4_t*c * fac
} 

PROCEDURE rate (celsius(degC)) {
	tcorr = temper (Q10,celsius)
	invc1_t = invc1*tcorr  
	invc2_t = invc2*tcorr
	invc3_t = invc3*tcorr 
	invo1_t = invo1*tcorr 
	invo2_t = invo2*tcorr 
	diro1_t = diro1*tcorr 
	diro2_t = diro2*tcorr 
	dirc2_t = dirc2*tcorr
	dirc3_t = dirc3*tcorr
	dirc4_t = dirc4*tcorr
}






