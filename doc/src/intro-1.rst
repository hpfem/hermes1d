============
Introduction
============

Mathematical Background
----------------------

When one speaks about the numerical solution of ODEs, one usually has in mind
initial value problems for equations of the form


.. math::

     {\d u_1\over\d x}=g_1(u_1, u_2, \dots, u_m, x),


.. math::
    :label: one

      \vdots


.. math::

     {\d u_m\over\d x}=g_m(u_1, u_2, \dots, u_m, x).

These are solved in a finite time interval $(0,T)$ using various time-stepping
methods. There are tons of those and some are quite sophisticated (meaning
multistep, higher-order, adaptive, etc.). But all of them have the following
common shortcomings:

* We would like to prescribe the initial value at $t = 0$ for some solution components and the end-time values at $t = T$ for others. Standard time stepping methods do not allow this.
* Global error control is problematic. One only can regulate the time step size locally -- this is something like "forward mesh refinement''. But one cannot do "backward mesh refinement'' or coarsening easily.
* We would like to prescribe a tolerance for the global error and then have the problem solved adaptively until this error tolerance is reached, without underresolving or overresolving too much. This is virtually impossible with adaptive time stepping methods.
* Standard time integration methods cannot change their order during the computation. For example, an adaptive RK4 method remains 4-order all the time. This is an analogy for $h$-refinement in FEM, and obviously it is highly inefficient. Correctly, the method should either do small low-order steps or large high-order steps to be efficient. We would like to see such an analogy of $hp$-refinement in ODE methods.
* We would like to solve more general ODEs than :eq:`one`.

This is why we decided to apply the $hp$-FEM methodology to ODEs and see what happens.

Equations
---------


We implemented the first version of Hermes1D during one day while returning
from the 2009 SIAM CSE conference. First we considered the form :eq:`one` but
then we realized that with no extra work we can actually assume a much more
general implicit form


.. math::

     f_1(u_1, u_2, \ldots, u_m, u'_1, u'_2, \ldots, u'_m, x) = 0,


.. math::
    :label: two

      \vdots


.. math::

     f_m(u_1, u_2, \ldots, u_m, u'_1, u'_2, \ldots, u'_m, x) = 0.

Note that :eq:`two` contains :eq:`one` as a special case.
In fact, :eq:`two` can be written shortly as

.. math::
    :label: qqq

      \bfF(\bfU, \bfU', x) = 0

where ${\bfU} = (u_1,\dots,u_m)$ and ${\bfF} = (f_1,\dots,f_m)$.

Boundary conditions
~~~~~~~~~~~~~~~~~~~


So far, we have considered Dirichlet boundary conditions only, which can be
imposed either at the initial time $t = 0$ or the end-time $t = T$. Exactly one
condition per solution component has to be defined.


hp-FEM discretization
---------------------


As always, the finite element discretization starts from a weak formulation.
With :eq:`two`, the situation is easy and we have


.. math::

     R_1(\bfY) = \int_0^T f_1(u_1, u_2, \ldots, u_m, u'_1, u'_2, \ldots, u'_m, x)v_1 \, \d t = 0,


.. math::
    :label: three

      \vdots


.. math::

     R_N(\bfY) = \int_0^T f_m(u_1, u_2, \ldots, u_m, u'_1, u'_2, \ldots, u'_m, x)v_N \, \d t = 0.

Here $v_1, v_2, \ldots, v_N$ are all basis functions for all solution
components (we can describe this more accurately if needed).  In the standard
sense, all basis functions corresponding to the solution component $u_i$ are
zero where $u_i$ has a Dirichlet boundary condition.  The vector $\bfY = (y_1,
y_2, \ldots, y_N)$ comprises all unknown coefficients of the finite element
basis functions for all solution components. The meshes for the solution
components $u_1, u_2, \ldots, u_m$ could (more precisely: *should*) be
different but for now we assume that they are the same.

Newton's method
---------------


We will drive the residual vector $\bfR = (R_1, R_2, \ldots, R_N)$ to zero
using the Newton's method. For that, we need the Jacobi matrix
$D\bfR/D\bfY$.

Let $1 \le i, j \le N$.
It is easy to calculate that

.. math::

     \frac{\partial R_i}{\partial y_j} = \int_0^T \frac{\partial f_{m(i)}}{\partial u_{n(j)}}(u_1, u_2, \ldots, u_m, u'_1, u'_2, \ldots, u'_m, x)v_jv_i


.. math::
    :label: newt1

      + \frac{\partial f_{m(i)}}{\partial u'_{n(j)}}(u_1, u_2, \ldots, u_m, u'_1, u'_2, \ldots, u'_m, x)v'_jv_i \, \d t = 0.

Here, the function $m(i)$ takes a global index $1 \le i \le N$ and returns the
index of the function $f_{m(i)}$ which is associated with $R_i$. Analogously,
$n(j)$ takes a global index $1 \le j \le N$ and returns the index of the
solution component $u_{n(i)}$ where the basis function $v_j$ belongs to.

The integral in :eq:`newt1` has two parts because the functions $u_s$ and
$u'_s$ depend on the same solution coefficients.  Do not be confused by the
derivatives with respect to $u'_{n(j)}$ in :eq:`newt1`.  The functions $u_s$
and $u'_s$ are used as independent variables for the differentiation.


