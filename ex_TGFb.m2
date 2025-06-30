---------------------------------------------------------------------------------------
--  TGF-beta receptor pathway (18 species) — steady-state ideal
--  Reference: Samal S S et al., “Geometric analysis of pathways
--  dynamics: Application to versatility of TGF-β receptors”,
--  Biosystems 149 (2016) 3-14
---------------------------------------------------------------------------------------
-- This script constructs the steady-state ideal arising from
-- mass-action kinetics plus two conservation laws and then
-- computes a reduced Groebner basis with respect to a pure
-- lexicographic order.  The resulting basis can be used for
-- elimination and for detecting parameter regions that admit
-- multistationarity.
---------------------------------------------------------------------------------------

-- 1.  Coefficient field: rational functions in rate constants.
K = frac(QQ[k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13, k14, k15, k16, k17, 
	k18, k19, k20, k21, k22, k23, k24, k25, k26, k27, k28, k29, k30, k31,c1,c2]);

-- Comment out the assignments below to retain symbolic rate constants
k1 = 	 1;
k2 = 	 1;
k3 = 	 1;
k4 = 	 1;
k5 = 	 1;
k6 = 	 1;
k7 = 	 1;
k8 = 	 1;
k9 = 	 1;
k10 =	 1;
k11 =	 1;
k12 =	 1;
k13 =	 1;
k14 =	 1;
k15 =	 1;
k16 =	 1;
k17 =	 1;
k18 =	 1;
k19 =	 1;
k20 =	 1;
k21 =	 1;
k22 =	 1;
k23 =	 1;
k24 =	 1;
k25 =	 1;
k26 =	 1;
k27 =	 1;
k28 =	 1;
k29 =	 1;
k30 =	 1;
k31 =	 1;
k32 =	 1;
k33 =	 1;
k34 =	 1;
k35 =	 1;
k36 =	 1;
k37 =	 1;
c1 = 	 1;
c2 = 	 1;

-- 2.  Polynomial ring in the eighteen species, ordered
--     x18 > … > x1 (pure lex).
R = K[x18, x17, x16, x15, x14, x13, x12, x11,x10,  x9,  x8,  x7,  x6,  x5,  x4,  x3, x2, x1,
      MonomialOrder => Lex];

-- 3.  Steady-state equations f1 … f18 (d(xi)/dt = 0).
f1  =  k2*x2 - k1*x1 - k16*x1*x11;
f2  =  k1*x1 - k2*x2 + k17*k34*x6;
f3  =  k3*x4 - k3*x3 + k7*x7 + k33*k37*x18 - k6*x3*x5;
f4  =  k3*x3 - k3*x4 + k9*x8	 - k8*x4*x6;
f5  =  k5*x6 - k4*x5 + k7*x7	 + 2*k11*x9 - 2*k10*x5^2 - k6*x3*x5 + k16*x1*x11;
f6  =  k4*x5 - k5*x6 + k9*x8	 + 2*k13*x10 - 2*k12*x6^2 - k17*k34*x6  + k31*k36*x8 - k8*x4*x6;
f7  =  k6*x3*x5- (k7 + k14)*x7;
f8  =  k14*x7- (k9 + k31*k36)*x8	+ k8*x4*x6;
f9  =  k10*x5^2- (k11 + k15)*x9;
f10 =  k15*x9- k13*x10       + k12*x6^2;
f11 =  k23*x14 - k30*x11;
f12 =  k18 - (k20 + k26)*x1 + k30*x11 + k27*x15 - k22*k35*x12*x13;
f13 =  k19 - (k21 + k28)*x13 + k30*x11 + k29*x16 - k22*k35*x12*x13;
f14 =  k22*k35*x12*x13 - (k23 + k24 + k25)*x14;
f15 =  k26*x12 - k27*x15;
f16 =  k28*x13 - k29*x16;
f17 =  k31*k36*x8      - k32*x17;
f18 =  k32*x17 - k33*k37*x18;

-- 4.  Conservation relations.
g1 = x1 + x2 + x5 + x6 + x7 + x8 + 2*x9 + 2*x10       - c1;  
g2 = x3 + x4 + x7 + x8 + x17 + x18       - c2;  
-- 5.  Generate the ideal and compute a reduced Groebner basis.
I = ideal(
     f1,  f2,  f3,  f4,  f5,  f6,  f7,  f8,  f9, f10, f11, f12, f13, f14, f15, f16, f17, f18, g1,g2
    );
G = gens gb I    -- generators of the reduced Groebner basis
