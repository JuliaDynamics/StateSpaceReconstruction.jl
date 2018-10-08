# StateSpaceReconstruction.jl

[![Build Status](https://travis-ci.org/kahaaga/StateSpaceReconstruction.jl.svg?branch=master)](https://travis-ci.org/kahaaga/StateSpaceReconstruction.jl)

Julia package for state space reconstruction and partitioning. This package provides necessary functionality for the [PerronFrobenius.jl](https://github.com/kahaaga/PerronFrobenius.jl),
[TransferEntropy.jl](https://github.com/kahaaga/TransferEntropy.jl) and
[CausalityTools.jl](https://github.com/kahaaga/CausalityTools.jl) packages.

Another package featuring state space reconstruction capabilities is [DynamicalSystems.jl](https://github.com/JuliaDynamics/DynamicalSystems.jl).
Although [CausalityTools.jl](https://github.com/kahaaga/CausalityTools.jl)
works seamlessly with `Dataset` instances from [DynamicalSystems.jl](https://github.com/JuliaDynamics/DynamicalSystems.jl),
the [CausalityTools.jl](https://github.com/kahaaga/CausalityTools.jl) package
has more specific needs than provided in [DynamicalSystems.jl](https://github.com/JuliaDynamics/DynamicalSystems.jl). In
particular, full reconstruction flexibility, and partitioning tools for
rectangular and simplex partitions, are needed in [CausalityTools.jl](https://github.com/kahaaga/CausalityTools.jl). Hence,
the existence of this library.


## Features
1. Generic embeddings of time series
2. SSR discretization (rectangular binning, simplex triangulations)

## Examples

Let's define three time series.

```
ts = [rand(30) for i = 1:3]
```

Embedding is done using the `embed` function, which come in two versions.
The first version takes only one argument, `ts`, which is a vector of time
series vectors (`typeof(ts)` must be `Vector{Vector{Float64}}` or
`Vector{Vector{Int}}`). Be default, this gives a non-lagged embedding.

```
# Default embedding (zero embedding lag along each dimension)
E = embed(ts)
```

For custom embeddings, you can provide two additional integer vectors `ts_inds`
and `embedding_lags`. Here, `ts_inds` specifies which time series appears
along which axis.  Embedding lags is a vector
(`length(ts_inds) == length(embedding_lags)`) specifying the embedding lag
along each dimension.

Positive, zero and negative lags are possible. Negative lags are takes as
"past affects future", and positive lags  as "future affects past".

```
# 3D embedding of only the first time series
E = embed(ts, [1, 1, 1], [0, -1, -2])

# 5D embedding of time series 1 and 2
E = embed(ts, [1, 1, 2, 2, 2], [0, -1, 0, -1, -1])

# 4D embedding of all three time series, some positive lags and some negative
# lags.
E = embed(ts, [1, 2, 3, 3], [-2, -1, 1, 0])
```
