[![Build Status](https://travis-ci.org/TorkelE/BRNEquilibrate.svg?branch=master)](https://travis-ci.org/TorkelE/BRNEquilibrate)
[![Build status](https://ci.appveyor.com/api/projects/status/f72vlmuvlpux7x6p?svg=true)](https://ci.appveyor.com/project/TorkelE/BRNEquilibrate)
[![codecov](https://codecov.io/gh/TorkelE/Why.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/TorkelE/BRNEquilibrate)
[![Coverage Status](https://coveralls.io/repos/github/TorkelE/Why.jl/badge.svg)](https://coveralls.io/github/TorkelE/BRNEquilibrate)
# BRNEquilibrate
A package for finding steady states to biochemical reaction network models. Detailed docs coming soon.

The actual solving is done by homotopycontinuation, as implemented in https://github.com/JuliaHomotopyContinuation/HomotopyContinuation.jl. Hence, if you use this package it is important that you:

### Cite HomotopyContinuation.jl
If you find HomotopyContinuation.jl useful in your work, we kindly request that you cite the [following paper](https://link.springer.com/chapter/10.1007/978-3-319-96418-8_54):

```latex
@inproceedings{HomotopyContinuation.jl,
  title={HomotopyContinuation.jl: A Package for Homotopy Continuation in Julia},
  author={Breiding, Paul and Timme, Sascha},
  booktitle={International Congress on Mathematical Software},
  pages={458--465},
  year={2018},
  organization={Springer}
}
```

A preprint of this paper is [freely available](https://arxiv.org/abs/1711.10911).

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-stable-url]: https://www.juliahomotopycontinuation.org/HomotopyContinuation.jl/stable
[docs-dev-url]: https://www.juliahomotopycontinuation.org/HomotopyContinuation.jl/dev

[build-img]: https://travis-ci.org/JuliaHomotopyContinuation/HomotopyContinuation.jl.svg?branch=master
[build-url]: https://travis-ci.org/JuliaHomotopyContinuation/HomotopyContinuation.jl
[codecov-img]: https://codecov.io/gh/juliahomotopycontinuation/HomotopyContinuation.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/juliahomotopycontinuation/HomotopyContinuation.jl

