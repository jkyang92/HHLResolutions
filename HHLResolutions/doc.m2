beginDocumentation()


doc ///
    Key
        HHLResolutions
    Headline
        A package to organize some of the code for the resolutions project
    Description
        Text
            TODO
///


doc ///
    Key
        makeHHLResolution
        (makeHHLResolution,ToricMap)
        (makeHHLResolution,NormalToricVariety,Matrix)
    Headline
        Construct the resolution given by HHL
    Usage
        makeHHLResolution f
        makeHHLResolution (X,g)
    Inputs
        f : ToricMap
        X : NormalToricVariety
            the target of the map g viewed as a map of toric varieties
        g : Matrix
            the matrix of the map on N-lattices
    Outputs
        : ChainComplex
    Description
        Text
            TODO
///

doc ///
    Key
        andersonDiagonalResolution
        (andersonDiagonalResolution,NormalToricVariety)
    Headline
        Construct the resolution of the diagonal given by Anderson
    Usage
        andersonDiagonalResolution X
    Inputs
        X : NormalToricVariety
            the variety for which to compute the resolution of the diagonal
    Outputs
        : ChainComplex
    Description
        Text
            TODO
///

doc ///
    Key
        andersonLaurentModule
        (andersonLaurentModule,NormalToricVariety)
    Headline
        Construct a module for the Anderson stratification
    Usage
        andersonLaurentModule X
    Inputs
        X : NormalToricVariety
            the variety for which to compute $M_\psi$ corresponding to the diagonal for
    Outputs
        : Module
    Description
        Text
            This function will return the $S$ module $M_\psi$ from the writeup
            where $\psi$ is the stratification function for the resolution of the diagonal Anderson gives.
        Text
            Note that this is NOT the module that Anderson's resolution resolves.
    SeeAlso
        hhlLaurentModule
        gensToLaurentModule
///

doc ///
    Key
        hhlLaurentModule
        (hhlLaurentModule,ToricMap)
    Headline
        Construct a module for the HHL stratification
    Usage
        andersonLaurentModule f
    Inputs
        f : ToricMap
            the map for which to compute the module $M_\psi$
    Outputs
        : Module
    Description
        Text
            This function will return the $S$ module $M_\psi$ from the writeup
            where $\psi$ is the stratification function for the HHL resolution for a toric morphism $f$
        Example
            X = toricProjectiveSpace 3;
            Y = toricProjectiveSpace 1;
            phi = map(X,Y,matrix {{1},{2},{3}})
            C = makeHHLResolution phi;
            Mpsi = hhlLaurentModule phi
            assert(prune HH_0 C == prune Mpsi);
    SeeAlso
        andersonLaurentModule
        gensToLaurentModule
///

doc ///
    Key
        hhlVectors
        (hhlVectors,ToricMap)
        (hhlVectors,NormalToricVariety,Matrix)
    Headline
        Construct the vectors for the normals for the hyperplanes of HHL
    Usage
        hhlVectors f
        hhlVectors(Y,g)
    Inputs
        f : ToricMap
        X : NormalToricVariety
            the target of the map g viewed as a map of toric varieties
        g : Matrix
            the matrix of the map on N-lattices
    Outputs
        : Matrix
    Description
        Text
            This returns the vectors which are normals to the hyperplanes from HHL as the ROWS of a matrix.
            The rows are in 1-1 correspondence with the rays of Y and are in the same order.
///
