# PERIGEE
PERIGEE is a nonlinear dynamic finite element and isogeometric analysis code for multiphysics problems. The code has been under development since 2012 with the goal of providing a single C++ framework for parallel implementation of different physics problems using different element technologies.

## Table of Contents
- [Code Style](#Code-Style)
- [Design Pattern](#Design-Pattern)
- [References](#References)

## Code Style
1. Be clean. If you don't know what clean is, read code until you do.
2. Be consistent. Match the style of surrounding code, unless its not clean.
3. Make detailed commit message if your are commiting changes to the /src and /include folders.

### Preprocessor guard for header files
Header files are protected by an include guard, which should be the name of the file converted to preprocessor symbol style. For example, the Sys_Tools.hpp has the following at the beginning and the end.
```cpp
#ifndef SYS_TOOLS_HPP
#define SYS_TOOLS_HPP
// The actual code
#endif
```

### Const
When the variables or functions does not change, use `const` to mark them. In functions, we pass data by using references. Therefore, `const` will be used if the argument is an input. Things will be a bit more complicated when passing pointers. There can be different scenarios when passing data by a pointer. Be careful with the differences between `const double * const &data` and `double * const &data`. Please read this [page](https://isocpp.org/wiki/faq/const-correctness) for references.


## Design Pattern
The code is written in C++ using the [abstract factory pattern](https://en.wikipedia.org/wiki/Factory_method_pattern). An *abstract base class* is adopted with virtual functions and the actual behavior is defined in *derived classes*. 

As an example, the interface class for material properties is framed in [IMaterialModel.hpp](include/IMaterialModel.hpp). This base class does not define any material behavior, instead it only specifies what a material class shall define.
```cpp
virtual void get_PK(const Matrix_3x3 &F, Matrix_3x3 &P, Matrix_3x3 &S) = 0;
```
The above are two pure virtual function that every derived material model class shall define. They takes F, the deformation gradient, as an input, and output the Piloa stresses. The detailed definition of a specified material can be find in its derived class. For example, in [MaterialModel_NeoHookean_M94_Mixed.hpp](include/MaterialModel_NeoHookean_M94_Mixed.hpp), the NeoHookean material model with a volumetric energy given by Miehe in 1994 is defined. 

## References
### Finite Element Method
* [The Finite Element Method: Linear Static and Dynamic Finite Element Analysis](https://www.amazon.com/Finite-Element-Method-Mechanical-Engineering/dp/0486411818/ref=sr_1_2?keywords=the+finite+element+method&qid=1566093145&s=books&sr=1-2) by Thomas J.R. Hughes

* Incompressible Flow and the Finite Element Method, Volume 1: Advection-Diffusion and Isothermal Laminar Flow by P.M. Gresho and R.L. Sani

### C++
* [www.learncpp.com](http://www.learncpp.com).

* [Google C++ Style](https://google.github.io/styleguide/cppguide.html).

## Contact
Ju Liu, liujuy@gmail.com
