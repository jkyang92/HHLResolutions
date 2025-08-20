
TEST ///
X = toricProjectiveSpace 3
Y = X ** X;
phi = diagonalToricMap X
S = ring Y;
C := makeHHLResolution(Y,matrix phi);
assert (prune HH_0 C == comodule ideal phi);
assert (prune HH_1 C == 0);
assert (prune HH_2 C == 0);
assert (prune HH_3 C == 0);
assert (length C == 3);

assert(rank C_0 == 1);
assert(rank C_1 == 6);
assert(rank C_2 == 8);
assert(rank C_3 == 3);
///


TEST ///
X = hirzebruchSurface 2
Y = X ** X;
phi = diagonalToricMap X
S = ring Y;
C := makeHHLResolution(Y,matrix phi);
--TODO test HH_0
assert (prune HH_1 C == 0);
assert (prune HH_2 C == 0);
assert (length C == 2);

assert(rank C_0 == 2);
assert(rank C_1 == 5);
assert(rank C_2 == 3);
///


TEST ///
X = hirzebruchSurface 2;
Y = X ** X;
phi = diagonalToricMap X;
G := hhlModuleGens(phi);
assert(#G == 2);
assert(entries G_0 == {0,0,0,0,0,0,0,0});
assert(entries G_1 == {-1,0,0,1,1,1,0,0});


///


-- test the example from Anderson's paper
-- https://arxiv.org/pdf/2403.09653
TEST ///
X = normalToricVariety ({{1,0},{1,1},{0,1},{-1,1},{0,-1}},{{0,1},{1,2},{2,3},{3,4},{4,0}})
Y = X**X
S = ring Y
C = andersonDiagonalResolution X

assert(set degrees C_0 == set {degree 1_S,degree ((x_1*x_9)/(x_4*x_6))})
assert(set degrees C_1 == set ({x_0*x_1*x_8,
                               (x_0*x_1*x_8*x_9)/x_6,
                               x_1*x_9,
                               x_1*x_2*x_3*x_9,
                               (x_0*x_1^2*x_2*x_9)/(x_4*x_6),
                               (x_1*x_2*x_3*x_9)/(x_4*x_6)
                               }/degree))
assert(set degrees C_2 == set ({x_0*x_1*x_8*x_9,
                                (x_0*x_1^2*x_2*x_8*x_9)/x_6,
                                (x_0*x_1^2*x_2*x_3*x_9)/(x_4*x_6),
                                x_1*x_2*x_3*x_9
                                }/degree))

assert(rank C_0 == 2);
assert(rank C_1 == 6);
assert(rank C_2 == 4);
assert(C.dd^2==0);
assert(HH_1 C == 0);
assert(HH_2 C == 0);
///


TEST///
X = hirzebruchSurface 2;
Y = X**X;
diag = diagonalToricMap X;
C = makeHHLResolution diag;
Mpsi = hhlLaurentModule diag;
--TODO worry about whether this will always work, what we want to check is isomorphism, but what ensures the order of generators is the same
assert(HH_0 C == prune Mpsi);
///

TEST ///
  -- Bruns-Gubeladze example of a triangulation of P^2 which fails projective normality
  P = convexHull transpose matrix{{0,0,0}, {1,0,0}, {0,1,0}, {1,1,2}}
  -- FIXME: fails because can't find a fundamental domain with unit volume
  -- not projectively normal, so HH_0 must be not the coordinate ring itself
  C = makeHHLResolution diagonalToricMap normalToricVariety P
  -- this one is projectively normal, so HH_0 should be isomorphic to the coordinate ring
  C = makeHHLResolution diagonalToricMap normalToricVariety(2 * P)
///
