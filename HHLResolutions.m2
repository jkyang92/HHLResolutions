newPackage(
    "HHLResolutions",
    Version => "0.3",
    Date => "July 16, 2026",
    Authors => {
        { Name => "Jay Yang"
        , Email => "jay.k.yang@vanderbilt.edu"}
        },
    Headline => "Code to work with HHL and related resolutions",
    AuxiliaryFiles => true,
    PackageExports => {"NormalToricVarieties","Complexes"},
    PackageImports => {"SimplicialComplexes"},
--    DebuggingMode => true,
    HomePage => "https://github.com/jkyang92/HHLResolutions/"
    )

export {
    "sliceByHyperplanes",
    "makeResolutionTable",
    "makeHHLPolytopes",
    "makeHHLPolytopesRelative",
    "makeHHLResolution",
    "hhlPolytopes",
    "hhlResolution",
    "hhlVectors",
    "toFacesByDimension",
    "toricSemigroupGens",
    "hhlModuleGens",
    "andersonDiagonalModuleGens",
    "andersonDiagonalModuleVertices",
    "andersonVertexToExponent",
    "andersonLaurentModule",
    "andersonModule",
    "hhlLaurentModule",
    "hhlModule",
    "hyperplaneStratificationPolytopes",
    "gensToLaurentModule",
    "gensToToricModule",
    "andersonDiagonalResolution",
    "bondalThomsenStrata",
    "makeResolution",
    "lineBundleBondalThomsenMonad",
    }

--a deprecation helper for all of the renaming
--provide symbols!
addDeprecatedName = (oldName,newName) -> (
    -- compute these values immediately so that errors trigger in installPackage/needsPackage
    -- instead of when the function is run
    f := value newName;
    errStr := "The function " | toString oldName | " has been renamed, please use " | toString newName | " instead";
    oldName <- (args -> (
            printerr errStr;
            f args
            )))


load "./HHLResolutions/resolution_tools.m2"
load "./HHLResolutions/hhl.m2"
load "./HHLResolutions/anderson.m2"
load "./HHLResolutions/modules.m2"
load "./HHLResolutions/monads.m2"

load "./HHLResolutions/tests.m2"
load "./HHLResolutions/doc.m2"

end--

restart
installPackage("HHLResolutions", RemakeAllDocumentation=>  true)
installPackage"HHLResolutions"
needsPackage "HHLResolutions"
check HHLResolutions
viewHelp HHLResolutions


X = toricProjectiveSpace 3
Y = X ** X;
phi = diagonalToricMap X
S = ring Y;

hhlModuleGens(phi)
andersonDiagonalModuleGens(X)
andersonDiagonalModuleVertices(X)

X = hirzebruchSurface 2
Y = X ** X;
phi = diagonalToricMap X
S = ring Y;
G := hhlModuleGens(phi)
G' := andersonDiagonalModuleGens(X)
V := andersonDiagonalModuleVertices(X)
apply(V,v -> andersonVertexToExponent(X,v))

netList pairs bondalThomsenStrata rays X

ML := module G_0
ML_7
ML_6

X = normalToricVariety ({{1,0},{1,1},{0,1},{-1,1},{0,-1}},{{0,1},{1,2},{2,3},{3,4},{4,0}})
Y = X ** X;
phi = diagonalToricMap X
S = ring Y;
G := hhlModuleGens(phi)
G' := andersonDiagonalModuleGens(X)
ML := module G'_0

(ML_1+ML_9-(ML_4+ML_6))
G'_1

ML_0+ML_1+ML_8-(ML_5+ML_6+ML_3)

C := makeHHLResolution(Y,matrix phi);

M2 := trim HH_0 C

prune M1
prune M2
prune comodule M1
prune M2


X = normalToricVariety ({{1,0},{1,1},{0,1},{-1,1},{0,-1}},{{0,1},{1,2},{2,3},{3,4},{4,0}})
Y = X**X
S = ring Y

matrix rays Y
degreeTable := degrees S


degree (S_1*S_9) - degree (S_4*S_6)

C := andersonDiagonalResolution X

C.dd^2 == 0

degrees C_0
degrees C_1



primaryDecomposition ann HH_0 andersonDiagonalResolution X

decompose ideal Y
