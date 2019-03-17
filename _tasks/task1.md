---
layout: narrative
title: "Task 1 - The Big G"
rights: Public Domain
publication-date: 2019
toc:
- Problem
- Solution

---

<a id="title-page" />

<p class="centered larger">TASK 1 - The Big G</p>

---

* ToC
{:toc}

---

## PROBLEM

1. Create the mesh given below.
   > NOTE: mshr documentation can be found <a href=" https://bitbucket.org/fenics-project/mshr/wiki/API">here</a>
   ![G_mesh](../../assets/img/tasks/task1/G_mesh.png){: .center-image }

2. Solve the Poisson equation

   $$
   \begin{align}
    - \nabla^2 u(x, y) &= f(x, y), \\
    u_D &= g(x, y) \qquad \text{on}\ \ \Gamma_{D},
   \end{align}
   $$

   where

   $$
   \begin{align}
     f(x, y) &= e^{-(x^2 + y^2)} \\
     g(x, y) &= 1 \qquad \forall x = 4
   \end{align}
   $$
  
   > HINT: Use function near(x, x_val, tolerance) to mark the Dirichlet boundary condition. <a href="https://fenicsproject.org/docs/dolfin/1.6.0/python/programmers-reference/cpp/function/near.html">Source</a>  
3. Plot the solution as a colourmap and a surface plot.

---

## SOLUTION

```python
from dolfin import *
from mshr import *

# Create the mesh
rect_top = Rectangle(dolfin.Point(0., 0.), dolfin.Point(4., 1.))
rect_bottom = Rectangle(dolfin.Point(0., 4.), dolfin.Point(4., 5.))
rect_left = Rectangle(dolfin.Point(0., 1.), dolfin.Point(1., 4.))
rect_right = Rectangle(dolfin.Point(3., 1.), dolfin.Point(4., 3.))
rect_inner = Rectangle(dolfin.Point(3., 3.), dolfin.Point(2., 2.))
domain = rect_top + rect_bottom + rect_left + rect_right + rect_inner
mesh = generate_mesh(domain, 10)

V = FunctionSpace(mesh, 'P', 1)

# Set the boundary condition
u_D = Constant(1.0)

def boundary(x, on_boundary):
    return near(x[0], 4, 1e-14)

bc = DirichletBC(V, u_D, boundary)

# Setup the a and L matrices
u = TrialFunction(V)
v = TestFunction(V)
f = Expression('exp(-(x[0]*x[0] + x[1]*x[1]))', degree=2)
a = dot(grad(u),  grad(v))*dx
L = f*v*dx

# Solve the problem
u = Function(V)
solve(a == L, u, bc)

# Plot a colourmap of the solution
p = plot(u, title="Solution to Poisson")
plot(mesh)
p.set_cmap("viridis")
pyplot.colorbar(p)
pyplot.show();
```

![G_color_soln](../../assets/img/tasks/task1/G_color_soln.png){: .center-image }

```python
# Plot a surface of the solution
p = plot(u, title="Solution to Poisson", mode = "warp")
p.set_cmap("viridis")
pyplot.colorbar(p)
pyplot.show();
```
![G_surface_soln](../../assets/img/tasks/task1/G_surface_soln.png){: .center-image }

---
