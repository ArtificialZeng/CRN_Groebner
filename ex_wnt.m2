---------------------------------------------------------------------------------------
--  WNT signalling pathway (19 species)
--  Reference: Gross E. et al., “Algebraic systems biology:
--  A case study for the WNT pathway”, Bull. Math. Biol. 78 (2016) 21-51
---------------------------------------------------------------------------------------
-- This script encodes the WNT model analysed by Gross et al.,
-- imposes steady-state (mass-action) conditions together with
-- five conservation relations, and computes a reduced Gröbner
-- basis in pure lexicographic order. 
---------------------------------------------------------------------------------------

-- 1.  Coefficient field: rational functions in rate constants.
K = frac(QQ[k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13, k14, k15, k16, k17, 
	k18, k19, k20, k21, k22, k23, k24, k25, k26, k27, k28, k29, k30, k31,c1,c2,c3,c4,c5]);

-- 2. Polynomial ring in nineteen species, ordered
--    x19 > … > x1 (pure lexicographic).
R = K[
     x19, x18, x17, x16, x15, x14, x13, x12, x11, x10,
     x9,  x8,  x7,  x6,  x5,  x4,  x3,  x2,  x1,
     MonomialOrder => Lex
    ];

-- 3. Rate-constant and total-amount assignments
--   Comment out the assignments below to retain symbolic rate constants
k1 = 2;  k2 = 3;  k3 = 2;  k4 = 1;  k5 = 1;  k6 = 3;  k7 = 1;  k8 = 1;
k9 = 2;  k10 = 1; k11 = 2; k12 = 2; k13 = 1; k14 = 2; k15 = 1; k16 = 1;
k17 = 2; k18 = 1; k19 = 1; k20 = 5; k21 = 1; k22 = 1; k23 = 1; k24 = 1;
k25 = 1; k26 = 2; k27 = 2; k28 = 1; k29 = 1; k30 = 6; k31 = 1;

c1 = 5;  c2 = 16; c3 = 2;-- c4 left symbolic
c5 = 3;

-- 4. Steady-state equations  f1 … f19  (d(x)/dt = 0).
f1  = -k1*x1 + k2*x2;
f2  =  k1*x1- (k2 + k26)*x2 + k27*x3 - k3*x2*x4   + (k4 + k5)*x14;
f3  =  k26*x2        - k27*x3  - k14*x3*x6  + (k15 + k16)*x15;
f4  = -k3*x2*x4      - k9*x4*x10 + k4*x14 + k8*x16 + (k10 + k11)*x18;
f5  = -k28*x5 + k29*x7 - k6*x5*x8 + k5*x14 + k7*x16;
f6  = -k14*x3*x6 - k20*x6*x11 + k15*x15 + k19*x17 + (k21 + k22)*x19;
f7  =  k28*x5 - k29*x7 - k17*x7*x9 + k16*x15 + k18*x17;
f8  = -k6*x5*x8    + (k7 + k8)*x16;
f9  = -k17*x7*x9   + (k18 + k19)*x17;
f10 =  k12- (k13 + k30)*x10 - k9*x4*x10  + k31*x11 + k10*x18;
f11 = -k23*x11 + k30*x10 - k31*x11 - k20*x6*x11 - k24*x11*x12 + k25*x13 + k21*x19;
f12 = -k24*x11*x12  + k25*x13;
f13 =  k24*x11*x12  - k25*x13;   
f14 =  k3*x2*x4     - (k4 + k5)*x14;
f15 =  k14*x3*x6    - (k15 + k16)*x15;
f16 =  k6*x5*x8     - (k7 + k8)*x16;
f17 =  k17*x7*x9    - (k18 + k19)*x17;
f18 =  k9*x4*x10    - (k10 + k11)*x18;
f19 =  k20*x6*x11   - (k21 + k22)*x19;

-- 5. Conservation relations.
g1 =  x1 + x2 + x3 + x14 + x15 - c1;
g2 =  x4 + x5 + x6 + x7 + x14 + x15+ x16 + x17 + x18 + x19 - c2;
g3 =  x8 + x16- c3;
g4 =  x9 + x17- c4;
g5 =  x12 + x13 - c5;

-- 6. Generate the ideal and compute a reduced Groebner basis.
I = ideal(
     f1,  f2,  f3,  f4,  f5,  f6,  f7,  f8,  f9,  f10,
     f11, f12,          f14, f15, f16, f17, f18, f19,
     g1,  g2,  g3,  g4, g5
    );
time G = gens gb I	         -- lex-reduced Groebner basis generators
