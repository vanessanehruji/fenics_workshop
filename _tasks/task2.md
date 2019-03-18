---
layout: narrative
title: "Task 2 - Coffee Cup"
rights: Public Domain
publication-date: 2019
toc:
- Problem
- Solution

---

<a id="title-page" />

<p class="centered larger">TASK 2 - The Coffee Cup</p>

---

* ToC
{:toc}

---

## PROBLEM

1. Make a "coffee cup" mesh similar to the one below. 
   ![coffee_mesh](../../assets/img/tasks/task2/cup_mesh.png){: .center-image }
2. The variational problem has been changed such that it includes Neumann boundaries

   $$
   \begin{align}
   - \nabla^2 u(x, y, z) &= f(x, y, z) \\
   u_{D} &= g(x, y, z) &\text{on}\ \ \Gamma_{D} \\
   \nabla u \cdot n &= h(x, y, z) &\text{on}\ \ \Gamma_{N}
   \end{align}
   $$
   
   where
   
   $$
   \begin{align}
    f(x, y, z) &= -e^{3z^2} \\
    g(x, y, z) &= 30 \qquad \forall z = 0 \\
    h(x, y, z) &= 1 \qquad \text{on outer boundary of cup body}
   \end{align}
   $$
   
   + TODO: Understand marking boundaries and make tutorial
3. Create a scatter plot of the solution (automatically done by the function plot())

---

## SOLUTION

```python
# Create the mesh
outer_body = Cylinder(dolfin.Point(0., 0., 0.),
	              dolfin.Point(0., 0., 1.), 0.4, 0.4)
inner_body = Cylinder(dolfin.Point(0., 0.0, 0.1),
	     	      dolfin.Point(0., 0.0, 1.), 0.3, 0.3)
outer_handle = Cylinder(dolfin.Point(0.4, -0.05, 0.5),
	                dolfin.Point(0.4, 0.05, 0.5), 0.35, 0.35)
inner_handle = Cylinder(dolfin.Point(0.4, -0.05, 0.5),
	                dolfin.Point(0.4, 0.05, 0.5), 0.25, 0.25)

domain = (outer_body + (outer_handle - inner_handle)) - inner_body 
mesh = generate_mesh(domain, 25)

V = FunctionSpace(mesh, 'P', 1)

u_D = Constant(30.0)

def boundary(x, on_boundary):
    # add dirichlet to cup inside
    return near(x[2], 0, 1e-14)

bc = DirichletBC(V, u_D, boundary)

u = TrialFunction(V)
v = TestFunction(V)
f = Expression('-exp(3*x[2]*x[2])', degree=2)
h = Constant(0.00005) # Need to apply only on outer cup boundaries
a = dot(grad(u),  grad(v))*dx
L = f*v*dx + h*v*ds

u = Function(V)
solve(a == L, u, bc)

## Plot scatter plot solution
p = plot(u, title="Solution to Poisson")
p.set_cmap("viridis")
pyplot.colorbar(p)
pyplot.show();
```

---