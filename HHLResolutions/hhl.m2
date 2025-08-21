
-- borrowed from helpers.m2
importFrom_Core {"concatCols", "concatRows", "concatBlocks"}
rows = m -> apply(numrows m, i -> m^{i})
cols = m -> apply(numcols m, j -> m_{j})
gale = m -> gens ker (if coker m == 0 then identity else transpose) m

getMinMaxRange = v -> (v = first entries v; floor min v ..< ceiling max v)

-- takes a polyhedron, and a matrix whose rows are the rays
-- returns a list of polyhedron which are the maximal cells
-- after slicing by the hyperplanes with integer constant term
-- and normals given by the rows of A
sliceByHyperplanes = method();
sliceByHyperplanes(Polyhedron,Matrix) := (P, A) -> (
    cells := {P};
    verts := vertices P;
    for rho in rows A do cells = (
	slices := hashTable for i in getMinMaxRange(rho * verts) list (
	    i => polyhedronFromHData(-rho || rho, matrix {{-i}, {i+1}}));
	--
	flatten for cell in cells list (
	    for i in getMinMaxRange(rho * vertices cell) list (
		cell' := intersection(cell, slices#i); -- ~40% of this computation, ~10s
		-- only include the top dimension cells
		if dim cell' == dim P then cell' else continue))
	);
    cells)

--takes a matrix and constructs the hyperplane for the kernel as a polyhedron
kernelPolyhedron = phi -> (
    L := transpose phi;
    n := numRows L;
    polyhedronFromHData(
	L, transpose matrix {toList (n:0)},
	L, transpose matrix {toList (n:0)})
    )

--given the maximal cells, get all polyhedra in all codimensions as a hash table.
--assumes that the cells have the correct intersection structure and all cells are compact
--returns a matrix of all the vertces and a hash table containing the cells
--represented by a tuple of indices into the vertex table.
toFacesByDimension = (cells) -> (
    verts := unique cols concatCols apply(cells, vertices);
    labels := hashTable toList(reverse \ pairs verts);
    -- Note: computing "faces cell" is very slow
    labelFaces := (verts, cell) -> applyValues(faces cell, faceList ->
	apply(faceList, face -> apply(face#0, i -> labels#(verts_{i})) ));
    facesList := apply(cells, c -> labelFaces(vertices c, c)); -- ~97% of this computation, ~315s
    faceTable := fold(facesList, (h1, h2) -> merge(h1, h2, join));
    (concatCols verts, applyValues(faceTable, unique)))

pointToRepresentative = (pt, L) ->(
    L' := L**QQ;
    c := solve(L',pt);
--    assert(all(entries c, x -> x>=0));
    c' := vector apply(entries c, x -> (x-floor x));
    L'*c'
    )

polytopeToRepresentative = (P, L) -> (
    Lcone := coneFromVData L;
    for i from 0 to numColumns L - 1 do (
        P' := P + convexHull(-L_{i});
        if contains(Lcone, P') then P = P');
    P)

partitionVertices = (pointsTable, L, polys) -> (
    values partition(p -> pointToRepresentative(pointsTable#p, L), polys))

--polyhedra on polyhedralComplexes is broken by an off by 1, so this is my horrible hack.
complexToPolytopes = PC -> applyKeys(faces PC, k -> dim PC - k - 1)

--takes a ring, a complex, and a lattice basis, the lattice basis should generate the fundemental parallelapiped
--used in the construction of the complex
--returns a hashTable containing enough information to make a resolution.
makeResolutionTable = (S,raysMatrix,cells,L) -> (
    --verts := vertices PC;
    elapsedTime (verts, faces) := toFacesByDimension cells; -- ~22% of the computation, ~325s
    d := max keys faces - 1;

    --polytopesByDimension := complexToPolytopes(PC);
    --polytopesByDimension = hashTable apply(select(keys polytopesByDimension, i -> i!=-1), k -> (k,polytopesByDimension#k));
    polytopesByDimension := applyKeys(faces,i -> d-i);
    polytopesByDimension = hashTable apply(select(keys polytopesByDimension, i -> i!=-1), k -> (k,polytopesByDimension#k));
    degreeMatrix := transpose matrix degrees S;
    pointToFineDegree := p -> (transpose matrix {apply(entries (raysMatrix * p), ceiling)})_0;
    pointToDegree := p -> (degreeMatrix * pointToFineDegree p);
    allPolyhedra := select(flatten values polytopesByDimension, p -> p!={});
    elapsedTime hulls := hashTable apply(allPolyhedra, p -> (p, convexHull verts_p)); -- ~14s
    elapsedTime pointsTable := hashTable apply(allPolyhedra, p -> (p, (1/#p * sum cols verts_p)_0)); -- ~27s, interiorPoint hulls#p is much slower
    elapsedTime modulesTable := applyValues(pointsTable, p -> S^{entries (- pointToDegree p)}); -- ~19s
    elapsedTime fineDegreeTable := applyValues(pointsTable, pointToFineDegree); -- ~16s
    elapsedTime polytopeClassesByDimension := applyValues(polytopesByDimension,
	polys -> partitionVertices(pointsTable, L, polys)); -- ~16s
    new HashTable from {
	"ring" => S,
	"verts" => verts,
	"dimension" => d,
	--polytope classes is a hash table indexed by dimension
        --where in each dimension there is a list of lists representing the equivalence classes of polytopes, 
        "polytopeClasses" => polytopeClassesByDimension,
        --modulesTable is a map from points to modules
	"modulesTable" => modulesTable,
        --fine degree table is a map from points to the actual fine degree for that point (not the degree after quotient by the lattice)
	"fineDegreeTable" => fineDegreeTable}
    )

subdivide = (K, V, A) -> (
    -- K: the kernel of the map of tori as a polyhedron
    -- V: the fundamental rays spanning the parallelogram
    -- A: the matrix of rays
    r := dim K;
    d := (mingens minors(r, V))_(0,0);
    -- take d copies of the the parallelogram cut out by the hyperplanes
    -- corresponding to the fundemental rays V
    hK := hyperplanes K;
    P := polyhedronFromHData(-V || V,
	transpose matrix {toList ((r:0) | (r-1:1) | (1:d))},
        hK#0,hK#1);
    -- print facets P;
    -- P  = intersection(K, P);

    -- the stratification
    cells := sliceByHyperplanes(P, A); -- ~26s
    cells)

makeHHLPolytopesRelative = method()
makeHHLPolytopesRelative(NormalToricVariety,ToricMap,Matrix) := (Y,phi,fundementalRays) -> (
    S := ring Y;
    raysMatrix := matrix rays Y;
    A := coker phi;
    psi:=map(A,target phi,1);
    Arays := psi*(transpose raysMatrix);
    --find a maximal rank submatrix of Arays that has deterimnant 1, if possible.
    r := rank Arays;
    assert(minors_r(fundementalRays)==1);
    K := kernelPolyhedron phi;
    assert(r == dim K);
    L := ker transpose phi;
    (subdivide(K,fundementalRays, raysMatrix),L)
    )

--find a subset of the columns of psi which has maximal rank
--and where the ideal of minors in ZZ has minimal non-zero generator
findMinimalMinorSubset = psi -> (
    -- find a maximal rank submatrix of psi that has deterimnant 1, if possible.
    bestSubset := null;
    bestDet := infinity;
    for s in subsets(numcols psi, r := rank psi) do (
        m := mingens minors(r, psi_s);
        d := if m==0 then infinity else m_(0,0);
	if d==1 then return s;
        if d<bestDet then (
            bestDet = d;
            bestSubset = s;
            ));
    print "warning: Fundamental parallelogram not cut out by hyperplanes from the rays";
    --otherwise return the subset that has the smallest minor
    bestSubset
)


makeHHLPolytopes = method()
makeHHLPolytopes(NormalToricVariety, Matrix) := (Y,    phi) -> makeHHLPolytopes(matrix rays Y,   kernel transpose phi)
makeHHLPolytopes(Ring, Matrix)               := (S,    phi) -> makeHHLPolytopes(gale effGenerators S, kernel transpose phi)
--this version takes a ring, a matrix who's rows parameterize the rays, and a matrix where the kernel of the dual gives L
makeHHLPolytopes(Ring, Matrix,       Module) := (S, A, L) -> makeHHLPolytopes(A,L)
--TODO remove the above version, the ring isn't actually used
makeHHLPolytopes(Matrix, Module) := (A, L) -> (
    g := mingens L;
    A' := A * g;
    n := rank L;
    K := polyhedronFromHData(map(ZZ^0,ZZ^n,0),map(ZZ^0,ZZ^1,0));
    
    -- the fundamental rays
    V := A'^(findMinimalMinorSubset transpose A');
    -- K := kernelPolyhedron phi;
    assert(rank A' == dim K);
    cells := subdivide(K,V,A');
    (cells,A',g))


--expects a toric variety, and a matrix mapping into the N-latice for Y, giving a toric inclusion.
makeHHLResolution = method()
makeHHLResolution ToricMap := phi -> makeHHLResolution(target phi, matrix phi)
makeHHLResolution(NormalToricVariety, Matrix) := (Y, phi) -> (
    elapsedTime (cells, raysMatrix, L) := makeHHLPolytopes(Y, phi); -- ~26s all in sliceByHyperplanes
    printerr("Cells Complete, " | #cells | " cells found");
    n := rank L;
    elapsedTime RT := makeResolutionTable(ring Y, raysMatrix, cells, mingens (ZZ^n)); -- ~28% of the computation here
    printerr "Labels Complete";
    elapsedTime makeResolution RT)                      -- ~70% of the computation here


--expects a toric map, returns the rays from HHL
hhlVectors = method()
hhlVectors(ToricMap) := phi -> hhlVectors(target phi, matrix phi)
hhlVectors(NormalToricVariety, Matrix) := (Y, phi) -> (
    raysMatrix := matrix rays Y;
    L := kernel transpose phi;
    g := mingens L;
    raysMatrix * g
    )


bondalThomsenStrata = method()
bondalThomsenStrata(List) := (normals) -> (
    normalsMatrix := matrix normals;
    L := source normalsMatrix;
    Lgens := mingens L;
    degreeSpace := cokernel normalsMatrix;
    degreeMap := inducedMap(degreeSpace, target normalsMatrix);
    hp := makeHHLPolytopes(normalsMatrix,source normalsMatrix);
    (verts,faces) := toFacesByDimension(hp#0);
    pointToFineDegree := p -> (transpose matrix {apply(entries (normalsMatrix * p), ceiling)})_0;
    allPolyhedra := select(flatten values faces, p -> p!={});
    hulls := hashTable apply(allPolyhedra, p -> (p, convexHull verts_p));
    pointsTable := hashTable apply(allPolyhedra, p -> (p, (1/#p * sum cols verts_p)_0));
    fineDegreeTable := applyValues(pointsTable, pointToFineDegree);
    polytopeClassesByCodimension := applyValues(faces,
        polys -> partitionVertices(pointsTable, Lgens, polys));
    strata := partition(x -> degreeMap (fineDegreeTable#(x#0)),
        flatten values polytopeClassesByCodimension);
    applyValues(strata, reps -> apply(flatten reps, p -> hulls#p))
    )
