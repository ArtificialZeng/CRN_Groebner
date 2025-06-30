---------------------------------------------------------------------------------------
--  Toy reaction network (2 species, 1 conservation law)
--  Reference: Berra S. et al., “Combined Newton-Gradient Method for
--  Constrained Root-Finding in Chemical Reaction Networks”
--  (J. Optim. Theory Appl. 200 : 404-427, 2024) - Supplementary Info
---------------------------------------------------------------------------------------
-- The script forms the steady-state ideal (mass-action kinetics plus
-- one conservation law) and computes a reduced Gröbner basis in a
-- pure lexicographic order.
---------------------------------------------------------------------------------------

-- 1. Coefficient field: rational-function field in rate constant k and 
-- conserved total c.
K = frac(QQ[k, c]);

-- 2. Polynomial ring in two species, ordered  x2 > x1 (lex).
R = K[x2, x1, MonomialOrder => Lex];

-- 3. Non-redundant steady-state equations 
f1 = x1^2- (k-1)*x1*x2;

-- 4. Conservation relation.
g1 =  x1 + x2 - c;

-- 5. Generate the ideal and compute a reduced Gröbner basis.
I = ideal(f1, g1);
G = gens gb I          -- lex-reduced Gröbner basis generators
