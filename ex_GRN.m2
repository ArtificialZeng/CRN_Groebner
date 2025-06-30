---------------------------------------------------------------------------------------
-- Steady-state analysis of a gene-regulatory CRN
-- (Conradi et al., PLoS Comput. Biol. 13 (10) 2017)
---------------------------------------------------------------------------------------
-- This script constructs the steady-state ideal arising from
-- mass-action kinetics plus two conservation laws and then
-- computes a reduced Groebner basis with respect to a pure
-- lexicographic order.  The resulting basis can be used for
-- elimination and for detecting parameter regions that admit
-- multistationarity.
---------------------------------------------------------------------------------------

-- 1.  Coefficient field:
--     a rational-function field in the ten rate constants k1…k10
--     and two conserved totals c1, c2.
K = frac(QQ[k1,k2,k3,k4,k5,k6,k7,k8,k9,k10, c1, c2]);

-- 2.  Polynomial ring in the seven species concentrations
--     ordered x7 ≻ x6 ≻ ⋯ ≻ x1 (pure lexicographic).
R = K[x7, x6, x5, x4, x3, x2, x1, MonomialOrder => Lex];

-- 3.  Polynomials from steady-state equations (d(xi)/dt = 0).
f1 =  k10*x7           - k9*x1*x6;   
f2 =  k6*x5            - k5*x2*x3;    
f3 =  k1*x1 - k3*x3    + k6*x5 - k5*x2*x3; 
f4 = -2*k7*x4^2 - k4*x4 + k2*x2 + 2*k8*x6;  
f5 =  k5*x2*x3         - k6*x5;  
f6 =  k7*x4^2 - k8*x6  + k10*x7 - k9*x1*x6;
f7 =  k9*x1*x6         - k10*x7;   

-- 4.  Conservation laws g1, g2 (total amounts are parameters c1, c2).
g1 =  x2 + x5 - c1;    -- total amount of species 2 + 5
g2 =  x1 + x7 - c2;    -- total amount of species 1 + 7

-- 5.  Generate the ideal and compute a reduced Groebner basis.
I = ideal(f1, f2, f3, f4, f5, f6, f7, g1, g2);
G = gens gb I         -- generators of the reduced Groebner basis
