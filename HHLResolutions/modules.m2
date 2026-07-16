--code to generate modules more directly from the subdivision

toricSemigroupGens = method();
toricSemigroupGens(ToricMap) := (phi) -> (
    --TODO this is a LOT of computation to get to the semigroup generators, there are certainly easier ways to do this
    X := target phi;
    n := #rays X;
    raysMatrix := matrix rays X;
    (cells,raysMatrix',L, fundamentalRays) := hhlPolytopes(target phi, matrix phi);
    ML := ZZ^n/image (raysMatrix * L);
    apply(n,i -> ML_i)
    )

--return the elements of ZZ^n/L generating the module
hhlModuleGens = method();
hhlModuleGens(ToricMap) := (phi) -> (
    X := target phi;
    n := #rays X;
    raysMatrix := matrix rays X;
    (cells,raysMatrix',L, fundamentalRays) := hhlPolytopes(target phi, matrix phi);
    ML := ZZ^n/(image (raysMatrix * L));
    pointToFineDegree := p -> (transpose matrix {apply(entries (raysMatrix' * p), ceiling)})_0;
    verts := unique flatten apply(cells,p -> (
            V:= vertices p;
            entries transpose V
            ));
    --apply(verts, v -> (vector v,vector(ML,entries pointToFineDegree vector v)))
    unique apply(verts, v -> vector(ML,entries pointToFineDegree vector v))
    )


andersonDiagonalModuleVertices = method()
andersonDiagonalModuleVertices(NormalToricVariety) := (X) -> (
    n := #rays X;
    d := dim X;
    raysMatrix := matrix rays X;
    (cells,raysMatrix',L, fundamentalRays) := hhlPolytopes(X, matrix toList (d:{}));
    verts := unique flatten apply(cells,p -> (
            V:= vertices p;
            entries transpose V
            ));
    apply(verts,vector)
    )

andersonVertexToExponent = method()
andersonVertexToExponent(NormalToricVariety,Vector) := (X,p) -> (
    raysMatrix := matrix rays X;
    v := (transpose matrix {apply(entries (raysMatrix * p), floor)})_0;
    w := v || -v;
    w
    )

andersonDiagonalModuleGens = method();
andersonDiagonalModuleGens(NormalToricVariety) := (X) -> (
    n := #rays X;
    d := dim X;
    raysMatrix := matrix rays X;
    (cells,raysMatrix',L, fundamentalRays) := hhlPolytopes(X, matrix toList (d:{}));
    pointToFineDegree := p -> (
        v := (transpose matrix {apply(entries (raysMatrix' * p), floor)})_0;
        v || -v);
    M := ZZ^n;
    ML := ZZ^(2*n)/image (raysMatrix || -raysMatrix);
    verts := unique flatten apply(cells,p -> (
            V:= vertices p;
            entries transpose V
            ));
    unique apply(verts, v -> vector(ML,entries pointToFineDegree vector v))
    )

gensToToricModule = method()
gensToLaurentModule = gensToToricModule
gensToToricModule(List,ToricMap) := (genExponents,phi) -> (
    imageIdeal := ideal phi;
    S := ring imageIdeal;
    degreeMatrix := transpose matrix degrees S;
    --get the maximum shift
    shift := apply(transpose genExponents,min);
    shiftDegree := entries (degreeMatrix * vector shift);
    shiftedExponents := apply(genExponents, e -> e-shift);
    print shiftedExponents;
    shiftModule := S^{-shiftDegree}/imageIdeal;
    --create the submodule
    M := sum apply(shiftedExponents,e -> S_e*shiftModule);
    M
    )

andersonModule = method();
andersonModule(NormalToricVariety) := (X) -> (
    gensToLaurentModule(andersonDiagonalModuleGens(X) / entries, diagonalToricMap X)
    )
andersonLaurentModule = andersonModule;

hhlModule = method();
hhlModule(ToricMap) := (phi) -> (
    gensToLaurentModule(hhlModuleGens(phi) / entries, phi)
    )
hhlLaurentModule = hhlModule;

--TODO get the favero huang module from the stratification...
