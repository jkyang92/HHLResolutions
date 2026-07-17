beginDocumentation()


doc ///
    Key
        HHLResolutions
    Headline
        A package to provide and use resolutions of toric subvarieties of toric varieties
    Description
        Text
            The key functions which are likely to be of interest to most
            users of this package are @TO hhlResolution@ which constructs
            the HHL resolution for a map into a toric variety and
            @TO lineBundleBondalThomsenMonad@ which constructs a monad for
            a twist of the Cox ring using only the Bondal Thomsen twists.
            The Bondal Thomsen stratification can also be obtained via
            @TO bondalThomsenStrata@.
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
        Text
            There is a old (deprecated) name @TT "makeHHLResolutions"@. If you are using the old name please switch
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
            This function will return the $S$ module $M_\psi\otimes_{S[L]}S$
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

doc ///
    Key
        lineBundleBondalThomsenMonad
        (lineBundleBondalThomsenMonad,NormalToricVariety,List)
    Headline
        Construct the the Bondal Thomsen monad for a line bundle
    Usage
        lineBundleBondalThomsenMonad(X,a)
    Inputs
        X : NormalToricVariety
            the toric variety over which to work
        a : List
            the exponent vector corresponding to a monomial twist for the line bundle
    Outputs
        : Complex
    Description
        Text
            As an example, for $\mathbb{P}^2$, and the line bundle $O(1)$, the following example illustrates
            the computation. In general the monad may only be a "virtual monad" as in it may only be a monad
            at the level of sheaves
        Example
            X = toricProjectiveSpace 2;
            S = ring X;
            C = lineBundleBondalThomsenMonad(X,{1,0,0})
            prune HH C
        Text
            For Hirzebruch surfaces the situation is more complicated, and one does get monads that are not resolutions
        Example
            X = hirzebruchSurface 2;
            S = ring X;
            C = lineBundleBondalThomsenMonad(X,{0,2,0,1})
            prune HH C
            prune sheaf HH_(-1) C
            prune sheaf HH_0 C
///

doc ///
    Key
        bondalThomsenStrata
        (bondalThomsenStrata,NormalToricVariety)
    Headline
        Construct the the Bondal Thomsen strata for a toric variety
    Usage
        bondalThomsenStrata X
    Inputs
        X : NormalToricVariety
            the toric variety over which to work
    Outputs
        : HashTable
    Description
        Text
            The returned hash table contains is indexed by the exponent vector representing
            a monomial corresponding to that stratum
            Each stratum is given as a list polyhedra, the union of the interiors of which is the strata.
            The polyhedra are in some fundamental domain of the torus.
        Example
            btStrata = bondalThomsenStrata toricProjectiveSpace 2
            applyValues(btStrata, l -> apply(l, vertices))
    Caveat
        The fundamental domain chosen is somewhat arbitrary, some more care might be taken
        so that each stratum is contained in a polyhedron. Moreover, as the code
        currently uses a closed fundamental domain, it generally produces redundant polyhedra
        when viewed on the torus.
///

doc ///
    Key
        "makeHHLResolution"
    Description
        Text
            This is a deprecated name for @TO hhlResolution@
    SeeAlso
        hhlResolution
///

doc ///
    Key
        "andersonLaurentModule"
    Description
        Text
            This is a deprecated name for @TO andersonModule@
    SeeAlso
        andersonModule
///

doc ///
    Key
        "hhlLaurentModule"
    Description
        Text
            This is a deprecated name for @TO hhlModule@
    SeeAlso
        hhlModule
///
