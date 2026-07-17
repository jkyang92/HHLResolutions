--Code to construct monads of line bundles using the HHL resolutions
barycenter = (m) -> (
    n := numColumns m;
    (1/n)*(m*(vector toList (n:1)))
    )

-- much of this is copied from makeResolutionTable, see if we can generalize.
-- openHyperplanes should be a pair of a matrix and a vector.
cellsToComplex = (S, cells, raysMatrix, openHyperplanes) -> (
    (verts, cellFaces) := toFacesByDimension cells; -- ~22% of the computation, ~325s
    d := max keys cellFaces;
    polytopesByDimension := hashTable apply(select(keys cellFaces, i -> i!=-1), k -> (k,cellFaces#k));
    degreeMatrix := transpose matrix degrees S;
    pointToFineDegree := p -> (transpose matrix {apply(entries (raysMatrix * p), ceiling)})_0;
    pointToDegree := p -> (degreeMatrix * pointToFineDegree p);
    allPolyhedra := select(flatten values polytopesByDimension, p -> p!={});
    pointsTable := hashTable apply(allPolyhedra, p -> (p, (verts_p * transpose matrix { toList ((#p):(1/#p)) })_0)); -- ~27s, interiorPoint hulls#p is much slower
    --remove cells that are in the open part
    -- not removing open parts for now to keep the relevant information about the boundary maps
    -- pointsTable = selectValues(pointsTable, p -> all( entries ((openHyperplanes#0)*p - openHyperplanes#1), x -> x!=0));
    allPolyhedra = keys pointsTable;--reset the list of polyhedra to only the ones that matter
    hulls := hashTable apply(allPolyhedra, p -> (p, convexHull verts_p)); -- ~14s
    modulesTable := applyValues(pointsTable, p -> S^{entries (- pointToDegree p)}); -- ~19s
    fineDegreeTable := applyValues(pointsTable, pointToFineDegree); -- ~16s
    --remove extra polytopes
    polytopesByDimension = applyValues(polytopesByDimension, polys -> select(polys, p -> pointsTable#?p));
    rt := new HashTable from {
        "ring" => S,
        "verts" => verts,
        "dimension" => d,
        --polytope classes is a hash table indexed by dimension
        --where in each dimension there is a list of lists representing the equivalence classes of polytopes, 
        "polytopeClasses" => applyValues(polytopesByDimension, polys -> apply(polys, P -> {P})),
        --modulesTable is a map from points to modules
        "modulesTable" => modulesTable,
        --fine degree table is a map from points to the actual fine degree for that point (not the degree after quotient by the lattice)
        "fineDegreeTable" => fineDegreeTable};
    C := complex makeResolution rt;
    internalCells := set keys selectValues(pointsTable, p -> all( entries ((openHyperplanes#0)*p - openHyperplanes#1), x -> x!=0));
    -- print internalCells;
    internalIndices := for i from 0 to max(d,1) list (
        if i>d then {}
        else for p in pairs (rt#"polytopeClasses"#i) list (
            if not internalCells#?(p#1#0) then continue;
            p#0
            )
        );
    -- max(d,1) ensures the function complex has at least one map to work with
    C' := complex for i from 1 to max(d,1) list (
        submatrix(C.dd_i,internalIndices#(i-1),internalIndices#i)
        );
    -- this naiveTruncation trims off the potentially extra term in the case where d=0
    C' = naiveTruncation(C',0,d);
    --the inclusionMap
    incMap := map(C,C',for i from 0 to d list (
        submatrix(id_(C_i),,internalIndices#i)
        ));
    projMap := map(C',C,for i from 0 to d list (
        submatrix(id_(C_i),internalIndices#i,)
        ));
    hashTable {
        "verts" => verts,
        "cells" => rt#"polytopeClasses",
        "inc" => incMap,
        "proj" => projMap
        }
    )

--get the maps to and from the (pruned) homology to the complex C, may not exist in general
homologySummand = (C) -> (
    H := HH C;
    (a,b) := concentration C;
    K := kernel C.dd;
    kerInc := inducedMap(C,K);
    kerInc = prune kerInc;
    quotMap := inducedMap(H,K);
    quotMap = prune quotMap;
    inc := map(C,prune H,hashTable for i from a to b list (
        i => quotMap_i\\kerInc_i
        ));
    proj := map(prune H, C, hashTable for i from a to b list (
            i => kerInc_i \\ quotMap_i
            ));
    (inc,proj))


lineBundleBondalThomsenMonad = method();
lineBundleBondalThomsenMonad(NormalToricVariety,List) := (X, a) -> (
    S := ring X;
    n := #a;
    --this is the polytope for a single fundamental domain
    pointPolytopes := hhlPolytopes(X,map(ZZ^(dim X),ZZ^0,0));
    raysMatrix := matrix rays X;
    --the complex of the Alexander dual to the irrelevant ideal
    fullComplex := simplicialComplex (max X / (v -> product (gens S)_v));
    pointToSigns := (pt) -> (
        positions(entries ((-raysMatrix*pt)-vector a), v -> v<=0)
        );
    memoizedFibers := memoize(s -> inducedSubcomplex(fullComplex,(gens S)_s));
    baseCells := hashTable for s in subsets(#rays X) list (
        setS := set s;
        signMatrix := diagonalMatrix for i from 0 to (#rays X - 1) list ( if setS#?i then 1 else -1 );
        signedA := signMatrix * vector a;
        P := polyhedronFromHData(-signMatrix*raysMatrix,matrix {signedA});
        if isEmpty P then continue;
        (s,P)
        );
    --the return from HHLPolytopes is the normal vectors of the fundamental domain
    basisMatrix := pointPolytopes#3;
    extraShiftPolytope := polyhedronFromHData(-basisMatrix || basisMatrix,
	transpose matrix {toList (((dim X-1):1) | (1:(abs det basisMatrix)) | ((dim X):(0)))}); --We eventually use this to do a minkowski sum,  think about whether giving vertices would be faster.
    complexData := applyPairs(baseCells, (s,P) -> (
        setS := set s;
        if not isCompact P then return null;
        --get the cells
        shifts := latticePoints (P+extraShiftPolytope); --if P isn't full dimensional, is this enough?
        cells := flatten apply(shifts, p -> apply(pointPolytopes#0, P -> affineImage(P,p)));
        --intersection is needed for strata that aren't full dimensional, but containment is faster for full dimensional strata
        if dim P == dim X
        then cells = select(cells, Q -> contains(P, Q))
        else cells = select(apply(cells, Q -> intersect(P,Q)), Q -> not isEmpty Q);
        sComp := for i from 0 to (#rays X - 1) list (if setS#?i then continue else i);
        openRays := raysMatrix^sComp;
        --the vector transpose matrix works around the edge case where sComp is empty
        openA := vector transpose matrix {a_sComp};
        -- print (s,cells/vertices,openRays, openA);
        s => cellsToComplex(S, cells, raysMatrix, (openRays, -openA))
        ));
    -- print applyValues(baseCells, P -> (isCompact P , vertices P));
    -- print complexData;
    -- a map of complexes into and out of the specific cell
    -- the target complex is the specific strand of the E^1 page that the cell belongs to
    strandMaps := (signType, d, cell) -> (
        if not complexData#?signType then return null; --if we want to be careful, we could use zero maps.
        fiber := memoizedFibers(signType);
        --the map that includes the internal cells in the complex of all cells
        incMap := complexData#signType#"inc";
        --the map that projects onto the internal cells
        projMap := complexData#signType#"proj";
        projCellMap := submatrix(projMap_d,,{cell});
        incCellMap := submatrix(incMap_d,{cell},);
        ((map(target projMap,(complex source projCellMap)[-d],{projCellMap}))**(S ** prune HH fiber),
            (map((complex target incCellMap)[-d],source incMap,hashTable {d => incCellMap}))**(S ** prune HH fiber))
        );
    fiberTable := hashTable for s in keys baseCells list (
        s => homologySummand complex memoizedFibers s
        );
    -- this computes one step of the boundary,
    -- C is the HHL complex restricted to the chamber
    -- signType is the sign type of the target
    -- verts is the matrix of all vertices for the chamber
    -- cells is the list of all cells for the chamber
    -- cell indexes into the list of all cells
    -- d is the dimension of the cell
    -- gPrev is the previous step in the diagram chase
    boundaryStep := (C,signType,verts,cells,d,cell,gPrev) -> (
        fiber := memoizedFibers(signType);
        fiberComplex := S**(complex fiber);
        (homologyInc,homologyProj) := fiberTable#signType;
        degreeShift := (degrees C.dd_d)#1#cell;
        degreeShift' := if gPrev===null then degreeShift else (degrees ((target gPrev)_(max target gPrev)))#0;
        if degreeShift!=degreeShift'
        then (
            print ("Warning: mismatched degrees " | degreeShift | " and " | degreeShift');
            );
        -- because f will be part of a differential, it should be degree -1
        f := if gPrev === null then null else (S^{-degreeShift}**(S**homologyProj))*gPrev*inducedMap(source gPrev, (source gPrev)[-1],Degree=>-1);
        gInit := if gPrev === null
            then S^{-degreeShift}**(S**homologyInc)
            else (
                (c1,c2) := concentration source gPrev;
                map(S^{-degreeShift}**fiberComplex, (source gPrev)[-1], hashTable for i from c1 to c2 list (
                    (i+1) => gPrev_(i)//((S^{-degreeShift}**fiberComplex).dd_(i+1))
                    )));
        --these are indices into the cells of dimension d-1
        boundary := positions(entries(C.dd_d_cell), v -> v!=0);
        gs := for b in boundary list (
            currVerts := verts_(cells#(d-1)#b#0);
            n := numColumns currVerts;
            currSigns := pointToSigns barycenter currVerts;
            if currSigns == signType then continue;
            fiberMap := S**(complex map(memoizedFibers(currSigns),fiber,id_S));
            -- TODO aquire the sign/coefficient from C
            m := submatrix(C.dd_d,{b},{cell});
            --print map(,target gInit, C.dd_d_(b,cell));
            --degreeShift := (degrees m)#1#0;
            (b, ((m**target fiberMap)*(S^{-degreeShift}**fiberMap)*gInit))
            );
        (f,gs)
        );

    --this maps from vertices to indicies for each chamber
    baseCellsVertexTables := applyValues(complexData, dat -> (
        verts := dat#"verts";
        hashTable apply(numColumns verts, i -> (verts_i,i))
        ));
    --this maps from cells to indicies for each chamber
    baseCellsCellTables := applyValues(complexData, dat -> (
        cells := dat#"cells";
        applyValues(cells, cellClassList -> (
            hashTable flatten applyPairs(cellClassList,
                (i, cellList) -> apply(cellList, c -> (c,i)))))
    ));
    -- for each stratum, construct all of the appropriate maps
    -- the resulting hash table has keys given by the sign type of the source chamber
    -- the values will themselves be hash tables indexed by signs of the target chamber
    -- the values in the inside hash table will be maps of complexes (of degree -1)
    -- between the strands for the two chambers.
    connectingMaps := applyPairs(complexData, (signs,dat) -> (
        --intersect and find cells.
        verts := dat#"verts";
        cells := dat#"cells";
        --vertices of the boundary
        boundaryVerts := set select(numColumns verts, i -> (
                currSigns := pointToSigns(verts_i);
                -- complexData only contains the signs corresponding to compact cells
                -- we don't care about those which border a non-compact cell
                currSigns!=signs and complexData#?currSigns));
        --all the cells in dat that contain any vertex from the overlap
        boundaryCells := applyValues(cells, l ->
            positions(l, cl -> any(cl, c -> not isEmpty intersection(set c, boundaryVerts))));
        -- the set of cells that are in boundaryCells but are in the interior
        interiorBoundaryCells := applyPairs(boundaryCells, (d,l) -> (d,select(l,i -> (
                cell := cells#d#i#0;
                center := barycenter verts_cell;
                pointToSigns(center)==signs --this is a silly way to check the interior condition
                ))));
        C := target dat#"inc";
        --for each of the cells in interiorBoundaryCell, build the connecting maps
        allMaps := flatten apply(pairs interiorBoundaryCells, (d,l) -> (
            apply(l, cell -> (
                (sourceStrandInc,sourceStrandProj) := strandMaps(signs,d,cell);
                currDim := d;
                prevStepTables := hashTable {cell => null};
                newMaps := new MutableHashTable from {};
                while #prevStepTables != 0 do (
                    newSteps := flatten for currCell in keys prevStepTables list (
                        currCellVertIndices := cells#currDim#currCell#0;
                        currCellVerts := verts_currCellVertIndices;
                        currSigns := pointToSigns barycenter(verts_currCellVertIndices);
                        (f,gs) := boundaryStep(C,currSigns,verts,cells,currDim,currCell,prevStepTables#currCell);
                        if f=!=null and f!=0 then (
                            targetVertexTable := baseCellsVertexTables#currSigns;
                            currCellTargetVertIndices :=
                              sort apply(currCellVertIndices, i -> targetVertexTable#(verts_i));
                            targetCell :=
                              baseCellsCellTables#currSigns#currDim#currCellTargetVertIndices;
                            (targetStrandInc,targetStrandProj) := strandMaps(currSigns,currDim,targetCell);
                            liftedF := targetStrandInc*f*sourceStrandProj;
                            if newMaps#?currSigns
                            then newMaps#currSigns += liftedF
                            else newMaps#currSigns = liftedF;
                            );
                        gs
                        );
                    prevStepTables = hashTable(plus, newSteps);
                    currDim = currDim - 1;
                    );
                new HashTable from newMaps
                ))));
        --sum together maps to the common parts
        signs =>
            if #allMaps == 0 then hashTable {}
            else fold((maps1,maps2) -> (
                    merge(maps1,maps2,plus)),allMaps)
            ));
    -- print connectingMaps;
    -- these are the strands of the final complex that
    -- arise from each chamber of the base, they are given
    -- by the restricted HHL complex tensored agains the homology of the fibers
    strands := applyPairs(complexData, (s,dat) -> (
        s => (source dat#"inc")**(S**prune HH memoizedFibers(s))
        ));
    -- print strands;
    totalTarget := directSum values strands;
    C := complex map(target totalTarget.dd,source totalTarget.dd,
        for tgt in keys strands list (
            for src in keys strands list (
                if tgt==src
                then strands#src.dd
                else (
                    if connectingMaps#?src and connectingMaps#src#?tgt
                    then connectingMaps#src#tgt
                    else 0
                    )
                )
            )
        );
    C[dim X - 1]
    )

end

