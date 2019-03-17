---
layout: narrative
title: "The Variational Form"
author: Marin Lauber, Vanessa Nehruji
rights: Public Domain
publication-date: 2019
toc:
- Title Page
- Poisson Equation
- Appendix

---

<a id="title-page" />

<p class="centered larger">VARIATIONAL FORM</p>

---
* ToC
{:toc}

---

## POISSON EQUATION

The Poisson equation is given by:

$$
\nabla^2 u + f = 0
$$

which can be explicitly written out in 3D Cartesian form as:

$$
\left( \frac{\partial^2}{\partial x^2} + \frac{\partial^2}{\partial y^2} + \frac{\partial^2}{\partial z^2} \right) u(x, y, z) + f(x, y, z) = 0
$$

We want to transform this into the weak form. This is done by multiplying the LHS by a test function $$v$$ and integrating it over the whole domain.

$$
\int_{\Omega} (\nabla^2u + f) \,v\,\text{d}\Omega = 0
$$

```python
V = FunctionSpace(mesh, 'P', 1)

u = TrialFunction(V)
v = TestFunction(V)
```

---

## APPENDIX


---
