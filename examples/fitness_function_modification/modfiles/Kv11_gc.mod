COMMENT
:[$URL: https://bbpteam.epfl.ch/svn/analysis/trunk/IonChannel/xmlTomod/CreateMOD.c $]
:[$Revision: 1367 $]
:[$Date: 2010-03-26 15:17:59 +0200 (Fri, 26 Mar 2010) $]
:[$Author: rajnish $]
:Comment :
:Reference :Expression of a cloned rat brain potassium channel in Xenopus oocytes. Science, 1989, 244, 221-4
ENDCOMMENT

NEURON	{
	SUFFIX Kv11
	USEION k READ ek WRITE ik
	RANGE gkbar, gk, ik 
	GLOBAL md,mk,mtA,mtd,mtk
}

UNITS	{
	(S) = (siemens)
	(mV) = (millivolt)
	(mA) = (milliamp)
}

PARAMETER	{
	gkbar = 0.00025 (S/cm2) 
	md = -30.5 (mV)
	mk = -11.3943 (mV)
	mtA = 30 (ms)
	mtd = -76.56 (mV)
	mtk = 26.1479 (mV)
}

ASSIGNED	{
	v	(mV)
	ek	(mV)
	ik	(mA/cm2)
	gk	(S/cm2)
	mInf
	mTau
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
		mInf = 1.0000/(1+ exp((v - md)/mk)) 
		mTau = mtA/(1+ exp((v - mtd)/mtk)) 
		hInf = 1.0000/(1+ exp((v - -30.0000)/27.3943)) 
		hTau = 15000.0000/(1+ exp((v - -160.5600)/-100.0000))
	UNITSON
}