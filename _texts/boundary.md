---
layout: narrative
title: "Boundary Conditions"
author: Marin Lauber, Vanessa Nehruji
rights: Public Domain
publication-date: 2019
toc:
- Introduction
- Dirichlet Boundary Conditions
- Neumann Boundary Conditions
- Marking Boundaries
- Appendix

---

<a id="title-page" />

<p class="centered larger">BOUNDARY CONDITIONS</p>

---
* ToC
{:toc}

---

## INTRODUCTION

The partial differential equations we look at generally have certain constraints at the boundaries of the domain, called boundary conditions (BC). There are many different types of conditions but in this tutorial we will only focus on applying the standard Dirichlet and Neumann conditions in FEniCS at particular boundaries.

---

## DIRICHLET BC

The Dirichlet BC enforces a fixed value on the unknown function, $$u$$, on the boundary.

For example, if we consider the standard Poisson equation on a unit square domain, 

$$
\begin{align}
  -\nabla^2 u &= f \qquad x, y \in [0, 1], \\
\end{align}
$$

with a source term of
$$
\begin{equation}
  f(x, y) = 3xy
\end{equation}
$$.

We can force all the boundaries of the domain to take the value 1.

$$
\begin{align}
  u &= 2 \qquad \text{on}\ \Gamma_{D},
\end{align}
$$

In FEniCS, this is simplemented by calling the DirichletBC function with the FunctionSpace, a BC function, and the subdomain as arguments. The subdomain is a string with C code.
> In the code below, on_boundary is a boolean value indicating whether a point is on the boundary. 

```python
mesh = UnitSquareMesh(8, 8)
V = FunctionSpace(mesh, 'P', 1)
bc = DirichletBC(V, Constant(1.0), "on_boundary")
```

Solving the Poisson equation with this condition gives the following result:

```python
## Input the Variational form for Poisson
u = TrialFunction(V)
v = TestFunction(V)
f = Expression('3*x[0]*x[1]', degree=2)
a = dot(grad(u), grad(v))*dx
L = f*v*dx

## Solve Poisson
u = Function(V)
solve(a == L, u, bc)

## Plot Solution
p = plot(u, title="Solution to Poisson")
plot(mesh)
p.set_cmap("viridis")
pyplot.colorbar(p)
pyplot.show();
```

![dirichlet_all_boundaries](../../assets/img/boundary/d_all_bounds.png)

We could modify the BC such that only the RHS boundary is set to 1.

```python
bc = DirichletBC(V, Constant(1.0), "x[0]==1")
```

![dirichlet_RHS](../../assets/img/boundary/d_RHS.png)

We could pass in a Python function defining the boundary instead of C code.
> The function <a href="https://fenicsproject.org/docs/dolfin/1.6.0/python/programmers-reference/cpp/function/near.html">near(x, x0, eps)</a> checks whether x is near x0 within a tolerance of eps and returns True or False accordingly.

```python
def boundary(x, on_boundary):
    return near(x[0], 1, eps=1e-14)

bc = DirichletBC(V, Constant(1.0), boundary)
```
---

### Mutliple Dirirchlet BC

You can apply multiple Dirichlet conditions by creating a list of DirichletBC and passing this to the solve function.

For example, if we wanted to the RHS boundary to have a value of 1 only when $$y < 0.5$$ and a value of 2 when $$y > 0.5$$, then,

```python
bcs = [DirichletBC(V, Constant(1.0), "x[0] == 1 and x[1] < 0.5"),
       DirichletBC(V, Constant(2.0), "x[0] == 1 and x[1] > 0.5")]
       
## Setup Poisson
## ...

solve(a == L, u, bcs)
```

![dirichlet_mutliple](../../assets/img/boundary/d_multiple.png)

---

## NEUMANN BC

The Neumann boundary condition forces the derivative of the unknown function to a speicific value at the boundary. This is also called the 'natural' BC since it automatically appears in the weak formulation.

For example, the strong form of the Poisson equation is

$$
\begin{equation}
  -\nabla^2 u = f \qquad \text{on}\ \Omega
\end{equation}
$$

The weak formulation from this gives

$$
\begin{equation}
  \int_{\Omega} \nabla u \cdot \nabla v \text{d}\Omega = \int_{\Omega} f v \text{d}\Omega + \int_{\Gamma_{N}} \frac{\partial u}{\partial n} v \text{d} \Gamma
\end{equation}
$$

where the flux is $$g = \frac{\partial u}{\partial n}$$

$$
\begin{equation}
  \int_{\Omega} \nabla u \cdot \nabla v \text{d}\Omega = \int_{\Omega} f v \text{d}\Omega - \int_{\Gamma_{N}} g v \text{d} \Gamma
\end{equation}
$$

In the previous examples, we have always set $$g = 0$$, (i.e. a Neumann condition of 0) so we did not include the second term on the RHS of the equation above. To have RHS boundary with $$u = 1$$, and set a Neumann BC of $$g = 5$$ on all other boundaries, we have

```python
mesh = UnitSquareMesh(8, 8)
V = FunctionSpace(mesh, 'P', 1)
bc = DirichletBC(V, Constant(1.0), "x[0] == 1")

## Input the Variational form for Poisson
u = TrialFunction(V)
v = TestFunction(V)
f = Expression('3*x[0]*x[1]', degree=2)
g = Constant(5.0)
a = dot(grad(u), grad(v))*dx
L = f*v*dx - g*v*ds	# Add the Neumann BC

## Solve Poisson
u = Function(V)
solve(a == L, u, bc)
```

![neumann_on_all](../../assets/img/boundary/n_all.png)

---

## MARKING BOUNDARIES


---

## APPENDIX

Further Reading:
+ <a href="https://fenicsproject.org/pub/tutorial/sphinx1/._ftut1005.html">FEniCS Tutorial</a>


---
