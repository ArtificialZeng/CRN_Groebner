---------------------------------------------------------------------------------------
--  Minimal oscillating subnetwork of the MAPK cascade (14 species)
--  Reference: Hadač O. et al., “Minimal oscillating subnetwork in the
--  Huang–Ferrell model of the MAPK cascade”, PLOS ONE, 2017
---------------------------------------------------------------------------------------
-- The script reproduces the fourteen-species core that Hadač et al.
-- extracted from the full Huang–Ferrell MAPK model.  It imposes
-- steady-state (mass-action) conditions plus five conservation
-- relations and then computes a reduced Gröbner basis in 
-- lexicographic order.
---------------------------------------------------------------------------------------

-- 1.  Coefficient field: rational functions in rate constants.
K = frac(QQ[k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13, k14, k15, k16, k17, 
					k18,c1,c2,c3,c4,c5]);
	
-- 2. Polynomial ring in fourteen species, ordered
--    x14 > … > x1 (pure lexicographic).
R = K[
     x14, x13, x12, x11, x10, x9, x8, x7,
     x6,  x5,  x4,  x3,  x2,  x1,
     MonomialOrder => Lex
    ];

-- 3. Rate-constant assignments
--   Comment out the assignments below to retain symbolic rate constants
k1 = 1;  k2 = 1;  k3 = 1;  k4 = 1;  k5 = 1;  k6 = 1;
k7 = 1;  k8 = 1;  k9 = 1;  k10 = 1; k11 = 1; k12 = 1;
k13 = 1; k14 = 1; k15 = 1; k16 = 1; k17 = 1; k18 = 1;
c2 = 1;  c3 = 1;  c4 = 1;  c5 = 1;

-- 4. Steady-state equations  f1 … f14  (d(x) / dt = 0).
f1  =  k2*x3 + k6*x6- k1*x1*x2;
f2  =  k2*x3 + k3*x3- k1*x1*x2;
f3  =  k1*x1*x2 - k3*x3 - k2*x3;
f4  =  k3*x3 + k5*x6 + k8*x8 + k9*x8 + k14*x14 + k15*x14 - k4*x4*x5 - k7*x4*x7 - k13*x4*x9;
f5  =  k5*x6 + k6*x6- k4*x4*x5;
f6  =  k4*x4*x5 - k6*x6 - k5*x6;
f7  =  k8*x8 + k12*x11 - k7*x4*x7;
f8  =  k7*x4*x7 - k9*x8 - k8*x8;
f9  =  k9*x8 + k11*x11 + k14*x14 + k18*x13 - k13*x4*x9- k10*x9*x10;
f10 =  k11*x11 + k12*x11 + k17*x13 + k18*x13- k10*x9*x10 - k16*x10*x12;
f11 =  k10*x9*x10 - k12*x11 - k11*x11;
f12 =  k15*x14 + k17*x13 - k16*x10*x12;
f13 =  k16*x10*x12 - k18*x13 - k17*x13;
f14 =  k13*x4*x9 - k15*x14 - k14*x14;

-- 5. Conservation relations.
g1 = x2 + x3 - c1; 
g2 = x5 + x6 - c2;  
g3 = x1 + x3 + x4 + x6 + x8 + x14 - c3; 
g4 = x10 + x11 + x13 - c4; 
g5 = x7 + x8 + x9 + x11 + x12 + x13 + x14 - c5;

-- 6. Generate the ideal and compute a reduced Groebner basis.
I = ideal(
     f1,  f2,  f3,  f4,  f5,  f6,  f7,  f8,  f9,  f10, f11, f12, f13, f14,
     g1,  g2,  g3,  g4,  g5
    );
G = gens gb I   -- lex-reduced Groebner basis generators
