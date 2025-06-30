---------------------------------------------------------------------------------------
--  ERK network (9 species)  steady-state ideal
--  Reference: Loman T E et al., “Catalyst: Fast and flexible modeling
--  of reaction networks”, PLoS Comput Biol 2023.
---------------------------------------------------------------------------------------
-- The script reproduces a 9-species extracellular signal-
-- regulated kinase (ERK).  It builds the steady-state ideal
-- (mass-action kinetics plus three conservation relations) and then
-- computes a reduced Gröbner basis using a pure lexicographic
-- order.  The resulting basis can be used for
-- elimination and for detecting parameter regions that admit
-- multistationarity.
---------------------------------------------------------------------------------------

-- 1. Coefficient field: rational functions in twelve rate constants
--    k1 … k12 and three conserved totals c1, c2, c3.
K = frac(QQ[k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12, c1, c2, c3]);

-- 2. Polynomial ring in the nine species, ordered
--    x9 > x8 > … > x1 (pure lex order).
R = K[x9, x8, x7, x6, x5, x3, x2, x1, x4, MonomialOrder => Lex];

-- 3. Polynomials from steady-state equations (d(xi)/dt = 0).
f1 =  k2*x6  + k12*x9                        - k1*x1*x2;          
f2 =  k2*x6  + k4*x7 + k9*x6 + k10*x7        - k1*x1*x2 - k3*x2*x3;
f3 =  k4*x7 + k9*x6 + k8*x9                  - k3*x2*x3 - k7*x3*x5;
f4 =  k6*x8 + k10*x7 + k11*x8                - k5*x4*x5; 
f5 =  k6*x8 + k8*x9 + k11*x8 + k12*x9        - k5*x4*x5 - k7*x3*x5;
f6 =  k1*x1*x2                               - k9*x6  - k2*x6; 
f7 =  k3*x2*x3                               - k10*x7 - k4*x7;
f8 =  k5*x4*x5                               - k11*x8 - k6*x8;
f9 =  k7*x3*x5                               - k12*x9 - k8*x9;

-- 4. Conservation relations.
g1 = x2 + x6 + x7                         - c1;
g2 = x5 + x8 + x9                         - c2;
g3 = x1 + x3 + x4 + x6 + x7 + x8 + x9     - c3; 

-- 5. Build the ideal and compute a reduced Groebner basis.
I = ideal(f1, f2, f3, f4, f5, f6, f7, f8, f9, g1, g2, g3);
G = gens gb I   -- reduced Gröbner basis generators
