makeAndersonResolutionTable = (X,raysMatrix,cells,L) -> (
    S := ring (X**X);
    --verts := vertices PC;
    elapsedTime (verts, faces) := toFacesByDimension cells; -- ~22% of the computation, ~325s
    d := max keys faces - 1;

    polytopesByDimension := applyKeys(faces,i -> d-i);
    faces = null;
    polytopesByDimension = hashTable apply(select(keys polytopesByDimension, i -> i!=-1), k -> (k,polytopesByDimension#k));
    degreeMatrix := transpose matrix degrees S;
    pointToFineDegree := p -> (
        v := vector apply(entries (raysMatrix * p), floor);
        v||-v);
    allPolyhedra := select(flatten values polytopesByDimension, p -> p!={});
    vertexLabels := hashTable apply(numColumns verts,i -> (verts_i,pointToFineDegree(verts_i)));
    fineDegreeTable := hashTable apply(allPolyhedra, p -> (
            labels := apply(p, i -> vertexLabels#(verts_i));
            (p,vector apply(entries matrix labels, max))
            ));
    vertexLabels = null;
    pointToFineDegree = null;
    
    modulesTable := applyValues(fineDegreeTable, d -> (
            S^{-entries (degreeMatrix*d)}
            ));
    
    pointsTable := hashTable apply(allPolyhedra, p -> (p, (1/#p * sum cols verts_p)_0));
    --remove some variables to help with memory useage
    allPolyhedra = null;
    polytopeClassesByDimension := applyValues(polytopesByDimension,
	polys -> partitionVertices(pointsTable, L, polys));
    
    new HashTable from {
	"ring" => S,
	"verts" => verts,
	"dimension" => d,
	"polytopeClasses" => polytopeClassesByDimension,
	"modulesTable" => modulesTable,
	"fineDegreeTable" => fineDegreeTable}
    )


andersonDiagonalResolution = method();
andersonDiagonalResolution(NormalToricVariety) := (X) -> (
    d := dim X;
    (cells,raysMatrix,L) := makeHHLPolytopes(X, matrix toList (d:{}));
    printerr("Cells Complete, " | #cells | " cells found");
    n := rank L;
    RT := makeAndersonResolutionTable(X, raysMatrix, cells, mingens (ZZ^n));
    cells = null;
    printerr "Labels Complete";
    makeResolution(RT)
    )
