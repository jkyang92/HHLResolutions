beginDocumentation()


doc ///
    Key
        HHLResolutions
    Headline
        A package to provide and use resolutions of toric subvarieties of toric varieties
    Description
        Text
            TODO
///


doc ///
    Key
        hhlResolution
        (hhlResolution,ToricMap)
        (hhlResolution,NormalToricVariety,Matrix)
    Headline
        Construct the resolution given by HHL
    Usage
        hhlResolution f
        hhlResolution (X,g)
    Inputs
        f : ToricMap
        X : NormalToricVariety
            the target of the map g viewed as a map of toric varieties
        g : Matrix
            the matrix of the map on $N$-lattices
    Outputs
        : Complex
    Description
        Text
            Given a toric map, this constructs the HHL resolution of the image of the toric map.
            The map can either be provided as a @TO ToricMap@ or as a pair of a normal toric variety
            for the target space and the map of $N$-lattices as a matrix.
        Example
            X = toricProjectiveSpace 2;
            Y = X ** X;
            S = ring Y
            phi = diagonalToricMap(X);
            hhlResolution phi
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
        : Complex
    Description
        Text
            TODO
///

doc ///
    Key
        andersonModule
        (andersonModule,NormalToricVariety)
    Headline
        Construct a module for the Anderson stratification
    Usage
        andersonModule X
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
        hhlModule
        gensToToricModule
///

doc ///
    Key
        hhlModule
        (hhlModule,ToricMap)
    Headline
        Construct a module for the HHL stratification
    Usage
        hhlModule f
    Inputs
        f : ToricMap
            the map for which to compute the module $M_\psi\otimes_{S[L]}S$
    Outputs
        : Module
    Description
        Text
            This function will return the $S$ module $M_\psi\otimes_{S[L]}S$ from the writeup
            where $\psi$ is the stratification function for the HHL resolution for a toric morphism $\phi$
        Example
            X = toricProjectiveSpace 3;
            Y = toricProjectiveSpace 1;
            phi = map(X,Y,matrix {{1},{2},{3}})
            C = hhlResolution phi;
            Mpsi = hhlModule phi
            assert(prune HH_0 C == prune Mpsi);
        Text
            There is a old (deprecated) name @TT "hhlLaurentModule"@. If you are using the old name please switch
    SeeAlso
        andersonModule
        gensToToricModule
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
