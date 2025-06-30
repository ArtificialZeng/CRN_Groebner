---------------------------------------------------------------------------------------
--  Toy network (6 species)
--  Reference: Berra S. et al., “Combined Newton-Gradient Method for
--  Constrained Root-Finding in Chemical Reaction Networks”
--  (J. Optim. Theory Appl. 200 : 404-427, 2024) - Supplementary Info
---------------------------------------------------------------------------------------
-- The script encodes the six-species toy example used by Berra et al.
-- to benchmark their Newton-Gradient solver.  It builds the steady-state
-- ideal (mass-action kinetics plus three conservation laws) and then
-- computes a reduced Gröbner basis with respect to a pure
-- lexicographic order. 
---------------------------------------------------------------------------------------

-- 1. Coefficient field: rational functions in six rate constants
--    k1 … k6 and three conserved totals c1 … c3 (all kept symbolic).
K = frac(QQ[k1, k2, k3, k4, k5, k6, c1, c2, c3]);

-- 2. Polynomial ring in the nine species, ordered
--    x6 > … > x1 (pure lexicographic).
R = K[ x6, x5, x4, x3, x2, x1, MonomialOrder => Lex];

-- 3. Steady-state equations  f1 … f6  (d(x)/dt = 0).
f1 = -2*k1*x1^2           + 2*k2*x4;
f2 = -k3*x2*x4            + k4*x5;
f3 = -k5*x3*x5            + k6*x6;
f4 =  k1*x1^2             - k2*x4 - k3*x2*x4 + k4*x5;
f5 =  k3*x2*x4            - k4*x5 - k5*x3*x5 + k6*x6;
f6 =  k5*x3*x5            - k6*x6;

-- 4. Conservation relations.
g1 =  x1 + 2*x4 + 2*x5 + 2*x6          - c1;
g2 =  x2 + x5 + x6                     - c2;
g3 =  x3 + x6                          - c3;

-- 5. Generate the ideal and compute a reduced Groebner basis.
I = ideal(f1, f2, f3, f4, f5, f6, g1, g2, g3);
G = gens gb I          -- lex-reduced Groebner basis generators






​
