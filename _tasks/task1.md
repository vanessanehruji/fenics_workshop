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
   ![G_mesh](../../assets/img/tasks/task1/G_mesh.png)

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

---
