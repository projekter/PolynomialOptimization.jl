# Background on polynomial optimization

Given an optimization problem
```math
    \min_{x \in \mathbb R^n} \bigl\{ f(x) : h_i(x) = 0, g_j(x) \geq 0, G_k(x) \succeq 0 \ \forall i, j, k \bigr\}
```
with
```math
    \begin{align*}
        f\colon \mathbb R^n \to \mathbb R, f(x) & = \sum_{\lvert\alpha\rvert \leq d} f_\alpha x^\alpha \\
        h_i\colon \mathbb R^n \to \mathbb R, h_i(x) & = \sum_{\lvert\alpha\rvert \leq d_i^{\mathrm{eq}}} h_{i, \alpha} x^\alpha \\
        g_j\colon \mathbb R^n \to \mathbb R, g_j(x) & = \sum_{\lvert\alpha\rvert \leq d_j^{\mathrm{ineq}}} g_{j, \alpha} x^\alpha \\
        G_k\colon \mathbb R^n \to \mathbb S^{d_k}, G_k(x) & = \sum_{\lvert\alpha\rvert \leq d_j^{\mathrm{psd}}} G_{k, \alpha} x^\alpha \\
    \end{align*}
```
where ``f_\alpha, h_{i, \alpha}, g_{j, \alpha} \in \R`` and ``G_{k, \alpha} \in \mathbb S(\mathbb R^{d_k \times d_k})``.

## SOS formulation
The sums-of-squares formulation of the problem is given by
```math
    \max \biggl\{
        \bigl(f - \sigma_0 - \sum_i h_i \psi_i - \sum_j g_j \sigma_j - \sum_k \braket{G_k, S_k}\bigr)(0) :
        \sigma_0, \sigma_i \in \Sigma[x],
        \psi_i \in \mathbb R[x],
        S_k \in \Sigma^{d_k}[x]
    \biggr\}
```
where ``\Sigma[x]`` is the set of sums-of-squares polynomials in ``x``, ``\mathbb R[x]`` are all polynomials in ``x`` and
``\Sigma^{d_k}[x]`` are sums-of-squares matrices in ``x``.

This can be turned into an infinite-dimensional semidefinite program by introducing the Schauder basis of monomials indexed by
``\alpha``, leading to
```math
    \begin{align*}
        \max \biggl\{
            & f_0 - \sigma_{0, 0} - \sum_i h_{i, 0} \psi_{i, 0} - \sum_j g_{j, 0} \sigma_{j, 0} -
            \sum_k \braket{G_{k, 0}, S_{k, 0}} : \\
            & \quad f_\alpha = \sigma_{0, \alpha} + \sum_{\beta, \gamma : \beta + \gamma = \alpha} \biggl[
                \sum_i h_{i, \beta} \psi_{i, \gamma} + \sum_j g_{j, \beta} \sigma_{j, \gamma} +
                \sum_k \braket{G_{k, \beta}, S_{k, \gamma}}
            \biggr] \ \forall \alpha \neq 0, \\
            & \quad \sigma_0, \sigma_j, S_k \succeq 0,
            \psi_i \in \mathbb R^\infty
        \biggr\}
    \end{align*}
```

## Moment formulation
The SOS formulation has the dual
```math
    \min_{y \in \mathbb R^m} \biggl\{
        L_y(f) :
        y_0 = 1,
        M(y) \succeq 0,
        M(h_i y) = 0,
        M(g_j y) \succeq 0,
        M(G_k y) \succeq 0
    \biggr\}
```
where ``L_y(f)`` is the Riesz functional of ``f``, ``M(y)`` the moment matrix, and ``M(\theta y)`` are localizing matrices for
``\theta``.