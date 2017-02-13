# TESP
Tudat Earth Satellite Propagator (TESP) is a [Tudat](https://github.com/Tudat) app that lets propagate a body in an orbit about the Earth (optionally including the Sun and the Moon in the propagation) and export the results.

TESP takes as command-line argument a (path to a) plain-text file with extension `tespin`. This file contains the settings for the propagation: body properties, initial state, perturbations to be included in the acceleration model, propagator and integrator settings, output settings, etc. Then, by default TESP generates a `tespout` file with the same name as the input file in the same directory containing the requested results.

Even though the input files can be editted with any text editor, it is recommended to use **[Atom](https://atom.io) together with the [language-tesp](http://github.com/aleixpinardell/language-tesp) package**, which provides syntax highlighting, autocompletion snippets and as-you-write documentation for `tespin` files.

The results of the propagation can be easily processed using **[TESP's package for MATLAB](http://github.com/aleixpinardell/matlab-tesp)**.
