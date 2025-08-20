--contains some common code to construct resolutions from cell complexes

--given a matrix, find the "leftmost" maximal rank submatrix. resulting matrix
--is square if the original matrix was maximal rank.
linearlyIndependentVectors := mat -> (
    columns := apply(numColumns mat, i -> mat_i);
    n := numRows mat;
    M := matrix toList (n:{});
    for c in columns do (
        M' := M | matrix c;
        if rank M' > rank M then M = M';
        );
    M
    )

frame := (verts,p) -> (
    if #p<=1  then return matrix (toList ((numRows verts):{}));
    vP := verts_(p);
    fP := vP_(toList(1..<(numColumns vP))) - matrix toList ((numColumns vP - 1):vP_0);
    linearlyIndependentVectors fP
    )

--given two polytopes q \subset p of codimension 1, extend a frame on one to a frame on the other
extendBoundaryFrame := (verts,q,p,fQ) -> (
    --extend the frame on Q to a frame on P by adding a inward pointing vector at the end,
    --i.e. any other vector in P
    v := (toList(set(p)-set(q)))_0;
    fQ | (verts_{v}-verts_{q_0})
    )

frameSign := (fP,fQ) -> (
    M := solve(fP,fQ);
    if det M>0 then 1 else -1
    )

--from a resolution tables tuple, make an actual complex
makeResolution = (RT) -> (
    polytopeClassesByDimension:=RT#"polytopeClasses";
    modulesTable:=RT#"modulesTable";
    fineDegreeTable:=RT#"fineDegreeTable";
    S := RT#"ring";
    verts := RT#"verts";
    allPolytopes := flatten flatten values polytopeClassesByDimension;
    --frameTable := elapsedTime hashTable apply(allPolytopes, p -> (p,frame(verts,p)));
    boundaries := for i from 0 to RT#"dimension" - 1 list elapsedTime (
	if debugLevel > 0 then printerr("Computing differential " | i+1);
        targetPolys := polytopeClassesByDimension#i;
        sourcePolys := polytopeClassesByDimension#(i+1);
        sourceModule := directSum apply(sourcePolys,src -> modulesTable#(src#0));
        targetModule := directSum apply(targetPolys,tgt -> modulesTable#(tgt#0));
        
        M := map(targetModule,sourceModule, for tgts in targetPolys list (
            for srcs in sourcePolys list (
                src := srcs#0;--it doesn't matter which source polytope we take
                tgtRep := tgts#0;
                sum for tgt in tgts list (
                    if isSubset(tgt,src) then (
                        tgtM := modulesTable#tgt;
                        srcM := modulesTable#src;
                        degreeDiff := fineDegreeTable#src - fineDegreeTable#tgt;
                        assert(all(entries degreeDiff, x -> x >= 0));
                        -- srcFrame := frameTable#src;
                        srcFrame := frame(verts,src);
                        -- tgtRepFrame := frameTable#tgtRep;
                        tgtRepFrame := frame(verts,tgtRep);
                        -- tgtFrame := frameTable#tgt;
                        tgtFrame := frame(verts,tgt);
                        extendedTgtFrame := extendBoundaryFrame(verts,tgt,src,tgtFrame);
                        --try to get the signs ....
                        frameSign(tgtRepFrame,tgtFrame)*frameSign(extendedTgtFrame,srcFrame)*S_(entries degreeDiff)
                        )
                    else (
                        0
                        )
                    )
                ))
            )
        );
    chainComplex boundaries
    )
