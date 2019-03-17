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

The Dirichlet BC enforces a fixed value on the unknown function, $u$, on the boundary.

For example, if we consider the standard Poisson equation on a unit square domain, we can force the RHS boundary, i.e. when $$x = 1$$, to take the value 2.

$$
\begin{equation}
  -\nabla^2 u = f \qquad \text{on} \Gamma
  u = 2 \qquad \text{when} x = 1
\end{equation}
$$

---

## NEUMANN BC


---

## MARKING BOUNDARIES


---

## APPENDIX


---
