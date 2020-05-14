<h1 align="center">
    Control Volume Finite Element Method
</h1>

<p align="center">
    Project presented as part of the Final Term Paper in Exact Sciences of the <a href='http://www.ufjf.br/ufjf/'>Universidade Federal de Juiz de Fora</a>
</p>

<p align="center">
    <a href="#pencil2-what-is-cvfem">CVFEM</a>&nbsp;&nbsp;|&nbsp;&nbsp;
    <a href="#pushpin-problem-presentation">Problem presentation</a>&nbsp;&nbsp;|&nbsp;&nbsp;
    <a href="#pencil-dependencies">Dependencies</a>&nbsp;&nbsp;|&nbsp;&nbsp;
    <a href="#runner-how-to-run">How to run</a>&nbsp;&nbsp;|&nbsp;&nbsp;
    <a href="#book-references">References</a>&nbsp;&nbsp;|&nbsp;&nbsp;
    <a href="#computer-technologies">Technologies</a>&nbsp;&nbsp;
</p>

## :pencil2: What is CVFEM?

Control Volume Finite Element Method is a combination of two well known numerical methods: Finite Element Method (FEM) and Finite Volume Method (FVM). These methods aims to represent and evaluate partial differential equations (PDE) in the form of algebric equations.

The problems that these methods solve are usually two or three dimensional physical problems, such as heat transfer, fluid dynamics and magnetic flux, which are represented by complex and difficult to solve PDE. Therefore, numerical methods are generally used to discretize the domain in numerous elements and solving the equations for each of them.

CVFEM uses unstructured meshes, starting with a conservation equation in integral form. The solution domain is divided in an finite number of control volumes (CV) and the conservation equation is applied to each CV.

## :pushpin: Problem presentation

In this problem, CVFEM was used to obtain the discrete form of the conservation equation for a scalar quantity. Each integral represents, respectively, conservation of linear momentum, source, diffusion and advection:

<div align="center">

![Equation](https://render.githubusercontent.com/render/math?math=$\frac{d}{dt}\int_{V}\phi%20dV-\int_{V}Q%20dV-\int_{A}\kappa\nabla\phi\cdot\boldsymbol{n}dA%2B\int_{A}\(\boldsymbol{v}\cdot%20\boldsymbol{n}\)\phi%20dA=0)

</div>

This process is employed in the discretization of this equation generating a linear system of algebric equations. Thus, it is possible to obtain the following equation:

<div align="center">

![Equation](https://render.githubusercontent.com/render/math?math=$a_{i}\phi_{i}=\sum_{j=1}^{n_i}a_{i,j}\phi_{S_{i,j}}%2Bb_i)

</div>

Where *ai* and *aij* are the coefficients of the linear system of equations for the variable *phi*, and *bi* represents all source, transient and boundary terms contributions.


The case studied is a steady state advection-diffusion problem withoud sources. The domain geometry is a quarter of a circular crown, and considering a field of velocities and diffusivity varying radially and not including sources, the problem can be solved analytically for comparison criteria:

<div style='display:flex; justify-content:space-evenly; align-items:center'>

<img src='https://render.githubusercontent.com/render/math?math=$\phi=\frac{e^r-e^2}{e-e^2}'/>

<img src='https://res.cloudinary.com/lorransutter/image/upload/v1589420342/Domain.svg' height=150/>

</div>

![Solution](https://res.cloudinary.com/lorransutter/image/upload/v1589496924/CVFEM_solution.svg)

<!-- ![Equation](https://render.githubusercontent.com/render/math?math=$\phi=\frac{e^r-e^2}{e-e^2})

![Domain](https://res.cloudinary.com/lorransutter/image/upload/v1589420342/Domain.svg) -->

<!-- \frac{\partial v_x\phi}{\partial x} + \frac{\partial v_y\phi}{\partial y} - \frac{\partial}{\partial x}\left( \kappa \frac{\partial \phi}{\partial x} \right) - \frac{\partial}{\partial y}\left( \kappa \frac{\partial \phi}{\partial y} \right) = 0 -->

<!-- \frac{\partial \phi}{\partial r} = \frac{\partial^2 \phi}{\partial r^2} -->

<!-- v = \frac{1}{r} -->

<!-- \kappa = \frac{1}{r} -->

<!-- R_{in} \leqslant R_{out} -->


<!-- ![Domain](https://res.cloudinary.com/lorransutter/image/upload/v1589420342/Domain.svg) -->

<!-- ![Control Volume](https://res.cloudinary.com/lorransutter/image/upload/v1589420190/ControlVolume.svg) -->

## :pencil: Dependencies

Besides, of course, [Python](https://www.python.org/), you will need [NumPy](https://numpy.org/) library for numerical operations and [Matplotlib](https://matplotlib.org/) library for plotting.

Also, if you want to generate your own mesh you may use the free mesh generator [Gmsh](https://gmsh.info/).

## :runner: How to run

## :book: References

- Voller, Vaughan R, *Basic control volume finite element methods for fluids and solids*
- Versteeg, Henk Kaarle and Malalasekera, Weeratunge, *An introduction to computational fluid dynamics: The Finite Volume Method*

## :Computer: Technologies

- [Python](https://www.python.org/)
- [NumPy](https://numpy.org/)
- [Matplotlib](https://matplotlib.org/)
- [Gmsh](https://gmsh.info/)