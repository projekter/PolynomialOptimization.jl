export supports_rotated_quadratic, supports_quadratic, supports_lnorm, supports_lnorm_complex, supports_psd_complex,
    supports_dd, supports_dd_complex

@doc raw"""
    supports_rotated_quadratic(::AbstractSolver)

Indicates the solver support for rotated quadratic cones: if `true`, the rotated second-order cone
``2x_1x_2 \geq \sum_{i \geq 3} x_i^2`` is supported.
The default implementation returns `false`.
"""
supports_rotated_quadratic(::AbstractSolver) = false

@doc raw"""
    supports_quadratic(::AbstractSolver)

Indicates the solver support for the quadratic cone: if `true`, the second-order cone ``x_1^2 \geq \\sum_{i \geq 2} x_i^2``
is supported.
The default implementation returns the same value as [`supports_rotated_quadratic`](@ref).
"""
supports_quadratic(state::AbstractSolver) = supports_rotated_quadratic(state)

"""
    supports_psd_complex(::AbstractSolver)

This function indicates whether the solver natively supports a complex-valued PSD cone. If it returns `false` (default), the
complex-valued PSD constraints will be rewritten into real-valued PSD constraints; this is completely transparent for the
solver. If the function returns `true`, the solver must additionally implement [`add_var_psd_complex!`](@ref) and
[`add_constr_psd_complex!`](@ref).
"""
supports_psd_complex(::AbstractSolver) = false

@doc raw"""
    supports_dd(::AbstractSolver)

This function indicates whether the solver natively supports a diagonally-dominant cone (or its dual for the moment case).
If it returns `false` (default), the constraint will be rewritten in terms of multiple ``\ell_\infty``/``\ell_1`` norm
constraints (if supported, see [`supports_lnorm`](@ref)), together with slack variables and equality constraints. If these
``\ell``-norm constraints are also not supported, linear constraints will be used.
"""
supports_dd(::AbstractSolver) = false

@doc raw"""
    supports_dd_complex(::AbstractSolver)

This function indicates whether the solver natively supports a complex-valued diagonally-dominant cone (or its dual for the
moment case). If it returns `false` (default), the constraint will be rewritten in terms of quadratic constraints (if
supported, see [`supports_quadratic`](@ref)) or multiple ``\ell_\infty``/``\ell_1`` norm constraints (if supported, see
[`supports_lnorm_complex`](@ref)).
"""
supports_dd_complex(::AbstractSolver) = false

@doc raw"""
    supports_lnorm(::AbstractSolver)

Indicates the solver support for ``\ell_\infty`` (in the moment case) and ``\ell_1`` (in the SOS case) norm cones: if `true`,
the cone ``x_1 \geq \max_{i \geq 2} \lvert x_i\rvert`` or ``x_1 \geq \sum_{i \geq 2} \lvert x_i\rvert`` is supported.
The default implementation returns `false`.
"""
supports_lnorm(::AbstractSolver) = false

@doc raw"""
    supports_lnorm_complex(::AbstractSolver)

Indicates the solver support for complex-valued ``\ell_\infty`` (in the moment case) and ``\ell_1`` (in the SOS case) norm
cones: if `true`, the cone ``x_1 \geq \max_{i \geq 2} \lvert\operatorname{Re} x_i + \mathrm i \operatorname{Im} x_i\rvert`` or
``x_1 \geq \sum_{i \geq 2} \lvert\operatorname{Re} x_i + \mathrm i \operatorname{Im} x_i\rvert`` is supported.
The default implementation returns `false`.
"""
supports_lnorm_complex(::AbstractSolver) = false

@doc raw"""
    supports_sdd(::AbstractSolver)

This function indicates whether the solver natively supports a scaled diagonally-dominant cone (or its dual for the moment
case). If it returns `false` (default), the constraints will be rewritten in terms of multiple rotated quadratic or quadratic
constraints, one of which must be supported (see [`supports_rotated_quadratic`](@ref) and [`supports_quadratic`](@ref)).
"""
supports_sdd(::AbstractSolver) = false

@doc raw"""
    supports_sdd_complex(::AbstractSolver)

This function indicates whether the solver natively supports a complex-valued scaled diagonally-dominant cone (or its dual for
the moment case). If it returns `false` (default), the constraints will be rewritten in terms of multiple rotated quadratic or
quadratic constraints, one of which must be supported (see [`supports_rotated_quadratic`](@ref) and
[`supports_quadratic`](@ref)).
"""
supports_sdd_complex(::AbstractSolver) = false