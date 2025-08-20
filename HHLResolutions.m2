newPackage(
    "HHLResolutions",
    Version => "0.1",
    Date => "July 15, 2024",
    Authors => {
        },
    Headline => "Code to work with HHL and related resolutions",
    AuxiliaryFiles => true,
    PackageExports => {"NormalToricVarieties"},
    PackageImports => {"PruneComplex"},
    DebuggingMode => true
    )

export {
    "sliceByHyperplanes",
    "makeResolutionTable",
    "makeHHLPolytopes",
    "makeHHLPolytopesRelative",
    "makeHHLResolution",
    "hhlVectors",
    "toFacesByDimension",
    "toricSemigroupGens",
    "hhlModuleGens",
    "andersonDiagonalModuleGens",
    "andersonDiagonalModuleVertices",
    "andersonVertexToExponent",
    "andersonLaurentModule",
    "hhlLaurentModule",
    "gensToLaurentModule",
    "andersonDiagonalResolution"
    }


load "./HHLResolutions/resolution_tools.m2"
load "./HHLResolutions/hhl.m2"
load "./HHLResolutions/anderson.m2"
load "./HHLResolutions/modules.m2"
load "./HHLResolutions/tests.m2"

load "./HHLResolutions/doc.m2"

end--

restart
installPackage "HHLResolutions"
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
