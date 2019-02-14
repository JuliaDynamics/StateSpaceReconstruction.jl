# StateSpaceReconstruction.jl


[![Build Status](https://travis-ci.org/kahaaga/StateSpaceReconstruction.jl.svg?branch=master)](https://travis-ci.org/kahaaga/StateSpaceReconstruction.jl)

[![Docs](https://img.shields.io/badge/docs-stable-blue.svg)](https://kahaaga.github.io/StateSpaceReconstruction.jl/)

Julia package for state space reconstruction and partitioning. This package provides necessary functionality for the [PerronFrobenius.jl](https://github.com/kahaaga/PerronFrobenius.jl),
[TransferEntropy.jl](https://github.com/kahaaga/TransferEntropy.jl) and
[CausalityTools.jl](https://github.com/kahaaga/CausalityTools.jl) packages, which require fully flexible state space reconstructions, and partitioning tools for rectangular and simplex partitions.

This package works seamlessly with `Dataset` instances from [DynamicalSystems.jl](https://github.com/JuliaDynamics/DynamicalSystems.jl), which also provide methods for state space reconstruction.

See the `CausalityTools` documentation for more information.

| Documentation |
| ------------- |
| [![](https://img.shields.io/badge/docs-dev-blue.svg)](https://kahaaga.github.io/CausalityTools.jl/dev)  |
