## Release 0.4.0

### Breaking changes
- Embeddings and reconstructions now live in `CausalityToolsBase`. The function for performing embeddings still has the same name (`customembed`), but now accepts a vector of state vectors or a `DynamicalSystems.Dataset` as input and has a new syntax explicitly specifying the mapping between input variables and embedding variables and their lags. To embed vectors of timeseries, wrap them in a `Dataset` first, for example: if `x` and `y` are time series, then let `data = Dataset(x, y)`. To embed these data with state vectors ``E = \\{(y(t+1), y(t), y(t - 2), x(t)\\}`` new syntax is `customembed(data, Positions(2, 2, 2, 1), Lags(1, 0, -2, 0)`. The `Positions` instance maps the dynamical variables (columns of a `Dataset`) to the embedding variables/column, and the `Lags` instance are the lags to use for each embedding variable. 

### New functionality
- Simplices may now be subsampled either with random points (irregular sampling), or subdividing the simplex using a shape-preserving simplex splitting routine and taking the centroids of the resulting subsimplices as the sampling points (regular sampling). If `s = Simplex(rand(3, 4)` is a simplex, then `subsample(s, n = 10)` will return ten points from within the simplex. The `sample_randomly` argument controls the type of sampling (regular/irregular). 
- Simplices may now be refined into a triangulation of subsimplices. `s = Simplex(rand(4, 5)` is a simplex, then `refine(s, k= 2)` splits the simple with a splitting factor 2. The number of subsimplices generated for a `dim`-dimensional simplex is `k^dim`.
- Triangulation partitions can be constructed from a set of points using `triangulate(pts)`. The input may be a vector of vectors, a vector of `SVector`/`MVector`, a `DynamicalSystems.dataset` instance or a `CustomReconstruction` instance (the latter is the return type of `customembed`, and is just a wrapper around a `Dataset`).

### Bug fixes

- Fixed the `SVector` simplex constructor, which was not working.
- Fixed multiple imports and some duplicate method definitions.

### Improvement
- Many minor improvements to documentation.