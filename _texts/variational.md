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


~~~ bash
$ git clone 
$ cd ed
$ gem install bundler
$ bundle install
~~~

```python
V = FunctionSpace(mesh, 'P', 1)

u = TrialFunction(V)
v = TestFunction(V)
```

---

## APPENDIX


---

[Footnotes by Frederick Douglass]


<!-- Make sure to use &#x21a9;&#xfe0e; to generate ↩︎ manually -->

<sup id="fn1">*</sup> This is the same man who gave me the roots to prevent my being whipped by Mr. Covey. He was "a clever soul." We used frequently to talk about the fight with Covey, and as often as we did so, he would claim my success as the result of the roots which he gave me. This superstition is very common among the more ignorant slaves. A slave seldom dies but that his death is attributed to trickery. [&#x21a9;&#xfe0e;](#ref1)

<sup id="fn2">*</sup> She was free. [&#x21a9;&#xfe0e;](#ref2)

<sup id="fn3">*</sup> I had changed my name from Frederick Bailey to that of Johnson. [&#x21a9;&#xfe0e;](#ref3)

<sup id="fn4">*</sup> I am told that colored persons can now get employment at calking in New Bedford—a result of anti-slavery effort. [&#x21a9;&#xfe0e;](#ref4)
