export add_constr_nonnegative!, add_constr_rotated_quadratic!, add_constr_quadratic!, add_constr_linf!,
    add_constr_linf_complex!, add_constr_psd!, add_constr_psd_complex!, add_constr_dddual!, add_constr_dddual_complex!,
    add_constr_fix_prepare!, add_constr_fix!, add_constr_fix_finalize!, negate_fix, prepend_fix, fix_objective!, add_var_slack!

function add_constr_nonnegative! end

"""
    add_constr_nonnegative!(state::AbstractSolver{T,V}, indvals::Indvals{T,V}) where {T,V}

Add a nonnegative constraint to the solver that contains the decision variables (columns in the linear constraint matrix)
indexed according to `indvals`.
Falls back to the vector-valued version if not implemented.

See also [`Indvals`](@ref).
"""
add_constr_nonnegative!(state::AbstractSolver{T,V}, indvals::Indvals{T,V}) where {T,V} =
    add_constr_nonnegative!(state, IndvalsIterator(unsafe, indvals.indices, indvals.values, StackVec(length(indvals))))

"""
    add_constr_nonnegative!(state::AbstractSolver{T,V}, indvals::IndvalsIterator{T,V}) where {T,V}

Adds multiple nonnegative constraints to the solver that contain the decision variables (columns in the linear constraint
matrix) indices according to the entries in `indvals`.
Falls back to calling the scalar-valued version multiple times if not implemented.

See also [`IndvalsIterator`](@ref).
"""
function add_constr_nonnegative!(state::AbstractSolver{T,V}, iv::IndvalsIterator{T,V}) where {T,V}
    for indvals in iv
        add_constr_nonnegative!(state, indvals)
    end
    return
end

function add_constr_quadratic! end

@doc raw"""
    add_constr_quadratic!(state::AbstractSolver{T,V}, indvals::IndvalsIterator{T,V}) where {T,V}

Adds a quadratic constraint to the `N = length(indvals)` linear combinations of decision variables (columns in the conic
constraint matrix) indexed according to the `indvals`. This will read (where ``X_i`` is
``\mathit{indvals}_i.\mathit{values} \cdot x_{\mathit{indvals}_i.\mathit{indices}}``) ``X_1 \geq 0``,
``X_1^2 \geq \sum_{i = 2}^N X_i^2``.

See also [`Indvals`](@ref), [`IndvalsIterator`](@ref).

!!! note "Number of parameters"
    In the real-valued case, `indvals` is always of length three, in the complex case, it is of length four. If the scaled
    diagonally dominant representation is requested, `indvals` can have any length.

!!! warning
    This function will only be called if [`supports_quadratic`](@ref) returns `true` for the given state.
    If (rotated) quadratic constraints are unsupported, a fallback to a 2x2 PSD constraint is used.
"""
add_constr_quadratic!(::AbstractSolver{T,V}, ::IndvalsIterator{T,V}) where {T,V}

function add_constr_rotated_quadratic! end

@doc raw"""
    add_constr_rotated_quadratic!(state::AbstractSolver{T,V}, indvals::IndvalsIterator{T,V}) where {T,V}

Adds a rotated quadratic constraint to the `N = length(indvals)` linear combinations of decision variables (columns in the
conic constraint matrix) indexed according to the `indvals`. This will read (where ``X_i`` is
``\mathit{indvals}_i.\mathit{values} \cdot x_{\mathit{indvals}_i.\mathit{indices}}``) ``X_1, X_2 \geq 0``,
``2X_1 X_2 \geq \sum_{i = 3}^N X_i^2``.

See also [`Indvals`](@ref), [`IndvalsIterator`](@ref).

!!! note "Number of parameters"
    In the real-valued case, `indvals` is always of length three, in the complex case, it is of length four. If the scaled
    diagonally dominant representation is requested, `indvals` can have any length.

!!! warning
    This function will only be called if [`supports_rotated_quadratic`](@ref) returns `true` for the given state.
    If (rotated) quadratic constraints are unsupported, a fallback to a 2x2 PSD constraint is used.
"""
add_constr_rotated_quadratic!(::AbstractSolver{T,V}, ::IndvalsIterator{T,V}) where {T,V}

function add_constr_psd! end

"""
    add_constr_psd!(state::AbstractSolver{T,V}, dim::Integer, data::PSDMatrixCartesian{T,V}) where {T,V}

Add a PSD constraint of side dimension `dim` ≥ 3 to the solver. Its requested triangle is indexed according to the return value
of [`psd_indextype`](@ref)); these elements make up a linear matrix inequality with variables given by the keys when iterating
through `data`, which are of the type `T`.
Note that if [`add_constr_quadratic!`](@ref) is not implemented, `dim` may also be `2`.
This method is called if [`psd_indextype`](@ref) returns a [`PSDIndextypeMatrixCartesian`](@ref).

!!! hint "Complex-valued PSD variables"
    Note that this function will also be called for complex-valued PSD cones if [`supports_psd_complex`](@ref) returns `false`.
    The data will have been rewritten in terms of a real-valued PSD cone, which doubles the dimension.
    If the solver natively supports complex-valued PSD cones, [`add_constr_psd_complex!`](@ref) must be implemented.
"""
add_constr_psd!(::AbstractSolver{T,V}, ::Int, ::PSDMatrixCartesian{T,V}) where {T,V}

"""
    add_constr_psd!(state::AbstractSolver{T,V}, dim::Integer, data::IndvalsIterator{T,V}) where {T,V}

Add a PSD constraint of side dimension `dim` ≥ 3 to the solver. `data` is an iterable through the elements of the PSD matrix
one-by-one, in the order specified by [`psd_indextype`](@ref). The individual entries are [`Indvals`](@ref).
This method is called if [`psd_indextype`](@ref) returns a [`PSDIndextypeVector`](@ref).

!!! hint "Complex-valued PSD variables"
    Note that this function will also be called for complex-valued PSD cones if [`supports_psd_complex`](@ref) returns `false`.
    The data will have been rewritten in terms of a real-valued PSD cone, which doubles the dimension.
    If the solver natively supports complex-valued PSD cones, [`add_constr_psd_complex!`](@ref) must be implemented.
"""
add_constr_psd!(::AbstractSolver{T,V}, ::Int, ::IndvalsIterator{T,V}) where {T,V}

function add_constr_psd_complex! end

"""
    add_constr_psd_complex!(state::AbstractSolver{T,V}, dim::Int,
        data::PSDMatrixCartesian{T,Complex{V}}) where {T,V}

Add a Hermitian PSD constraint of side dimension `dim` ≥ 3 to the solver. Its requested triangle is indexed according to the
return value of [`psd_indextype`](@ref)); these elements make up a linear matrix inequality with variables given by the keys
when iterating through `data`, which are of the type `T`. The real part of any coefficient corresponds to the coefficient in
front of the real part of the matrix entry, the imaginary part is the coefficient for the imaginary part of the matrix entry.
Note that if [`add_constr_quadratic!`](@ref) is not implemented, `dim` may also be `2`.
This method is called if [`psd_indextype`](@ref) returns a [`PSDIndextypeMatrixCartesian`](@ref).

!!! warning
    This function will only be called if [`supports_psd_complex`](@ref) is defined to return `true` for the given state.
"""
add_constr_psd_complex!(::AbstractSolver{T,V}, ::Int, ::PSDMatrixCartesian{T,Complex{V}}) where {T,V}

"""
    add_constr_psd_complex!(state::AbstractSolver{T,V}, dim::Int, data::IndvalsIterator{T,V}) where {T,V}

Add a Hermitian PSD constraint of side dimension `dim` ≥ 3 to the solver. `data` is an iterable through the elements of the PSD
matrix one-by-one, in the order specified by [`psd_indextype`](@ref). The individual entries are [`Indvals`](@ref).
This method is called if [`psd_indextype`](@ref) returns a [`PSDIndextypeVector`](@ref).
Regardless of the travelling order, for diagonal elements, there will be exactly one entry, which is the real part. For
off-diagonal elements, the real part will be followed by the imaginary part. Therefore, the coefficients are real-valued.

!!! warning
    This function will only be called if [`supports_psd_complex`](@ref) is defined to return `true` for the given state.
"""
add_constr_psd_complex!(::AbstractSolver{T,V}, ::Int, ::IndvalsIterator{T,V}) where {T,V}

function add_constr_dddual! end

@doc raw"""
    add_constr_dddual!(state::AbstractSolver{T,V}, dim::Integer, data::IndvalsIterator{T,V}, u) where {T,V}

Add a constraint for membership in the dual cone to diagonally dominant matrices to the solver. `data` is an iterator through
the scaled lower triangle of the matrix. A basis change is induced by `u`, with the meaning for the primal cone that
`M ∈ DD(u) ⇔ M = uᵀ Q u` with `Q ∈ DD`.

!!! warning
    This function will only be called if [`supports_dd`](@ref) returns `true` for the given state. If diagonally dominant cones
    are not supported directly, a fallback to a columnwise representation in terms of ``\ell_\infty`` norms will be used (or
    the fallbacks if this norm is not supported).
"""
add_constr_dddual!(state::AbstractSolver{T,V}, dim::Integer, data::IndvalsIterator{T,V}, u) where {T,V}

function add_constr_dddual_complex! end

@doc raw"""
    add_constr_dddual_complex!(state::AbstractSolver{T,V}, dim::Integer, data::IndvalsIterator{T,V}, u) where {T,V}

Add a constraint for membership in the dual cone to complex-valued diagonally dominant matrices to the solver. `data` is an
iterator through the scaled lower triangle of the matrix. A basis change is induced by `u`, with the meaning for the primal
cone that `M ∈ DD(u) ⇔ M = u† Q u` with `Q ∈ DD`.
For diagonal elements, there will be exactly one entry, which is the real part. For off-diagonal elements, the real part will
be followed by the imaginary part. Therefore, the coefficients are real-valued.

!!! warning
    This function will only be called if [`supports_dd_complex`](@ref) returns `true` for the given state. If complex-valued
    diagonally dominant cones are not supported directly, a fallback to quadratic constraints on the complex-valued data is
    tried first (if supported), followed by a columnwise representation in terms of ``\ell_\infty`` norms or their fallback on
    the realification of the matrix data if not.
"""
add_constr_dddual_complex!(state::AbstractSolver{T,V}, dim::Integer, data::IndvalsIterator{T,V}, u) where {T,V}

function add_constr_linf! end

@doc raw"""
    add_constr_linf!(state::AbstractSolver{T,V}, indvals::IndvalsIterator{T,V}) where {T,V}

Adds an ``\ell_\infty`` norm constraint to the `N = length(indvals)` linear combinations of decision variables (columns in the
conic constraint matrix) indexed according to the `indvals`. This will read (where ``X_i`` is
``\mathit{indvals}_i.\mathit{values} \cdot x_{\mathit{indvals}_i.\mathit{indices}}``)
``X_1 \geq \max_{i > 2} \lvert X_i\rvert``.

See also [`Indvals`](@ref), [`IndvalsIterator`](@ref).

!!! warning
    This function will only be called if [`supports_lnorm`](@ref) returns `true` for the given state.
    If ``\ell_\infty`` norm constraints are unsupported, a fallback to multiple linear constraints will be used.
"""
add_constr_linf!(::AbstractSolver{T,V}, ::IndvalsIterator{T,V}) where {T,V}

function add_constr_linf_complex! end

@doc raw"""
    add_constr_linf_complex!(state, indvals::IndvalsIterator{T,V}) where {T,V<:Real}

Same as [`add_constr_linf!`](@ref), but now two successive items in `indvals` (starting from the second) are interpreted as
determining the real and imaginary part of a component of the ``\ell_\infty`` norm cone.

!!! warning
    This function will only be called if [`supports_lnorm_complex`](@ref) returns `true` for the given state.
    If complex-valued ``\ell_\infty`` norm constraints are unsupported, a fallback to multiple linear constraints and quadratic
    cones will be used. If [`supports_quadratic`](@ref) is not `true`, complex-valued DD cones cannot be used.
"""
add_constr_linf_complex!(::AbstractSolver{T,V}, ::IndvalsIterator{T,V}) where {T,V}

function add_constr_sdddual! end

@doc raw"""
    add_constr_sdddual!(state::AbstractSolver{T,V}, dim::Integer, data::IndvalsIterator{T,V}, u) where {T,V}

Add a constraint for membership in the dual cone to scaled diagonally dominant matrices to the solver. `data` is an iterator
through the (unscaled) lower triangle of the matrix. A basis change is induced by `u`, with the meaning for the primal cone
that `M ∈ SDD(u) ⇔ M = uᵀ Q u` with `Q ∈ SDD`.

!!! warning
    This function will only be called if [`supports_sdd`](@ref) returns `true` for the given state. If scaled diagonally
    dominant cones are not supported directly, a fallback to (rotated) quadratic cones will be used.
"""
add_constr_sdddual!(state::AbstractSolver{T,V}, dim::Integer, data::IndvalsIterator{T,V}, u) where {T,V}

function add_constr_sdddual_complex! end

@doc raw"""
    add_constr_sdddual_complex!(state::AbstractSolver{T,V}, dim::Integer, data::IndvalsIterator{T,V}, u) where {T,V}

Add a constraint for membership in the dual cone to complex-valued scaled diagonally dominant matrices to the solver. `data` is
an iterator through the (unscaled) lower triangle of the matrix. A basis change is induced by `u`, with the meaning for the
primal cone that `M ∈ SDD(u) ⇔ M = u† Q u` with `Q ∈ SDD`.
For diagonal elements, there will be exactly one entry, which is the real part. For off-diagonal elements, the real part will
be followed by the imaginary part. Therefore, the coefficients are real-valued.

!!! warning
    This function will only be called if [`supports_sdd_complex`](@ref) returns `true` for the given state. If complex-valued
    sclaed diagonally dominant cones are not supported directly, a fallback to quadratic constraints is automatically
    performed.
"""
add_constr_sdddual_complex!(state::AbstractSolver{T,V}, dim::Integer, data::IndvalsIterator{T,V}, u) where {T,V}

"""
    add_constr_fix_prepare!(state::AbstractSolver, num::Int)

Prepares to add exactly `num` constraints that are fixed to a certain value; the actual data is then put into the solver by
subsequent calls of [`add_constr_fix!`](@ref) and the whole transaction is completed by [`add_constr_fix_finalize!`](@ref).
The return value of this function is passed on as `constrstate` to [`add_constr_fix!`](@ref).
The default implementation does nothing.
"""
add_constr_fix_prepare!(::AbstractSolver, _) = nothing

"""
    add_constr_fix!(state::AbstractSolver{T,V}, constrstate, indvals::Indvals{T,V}, rhs::V) where {T,V}

Add a constraint fixed to `rhs` to the solver that is composed of all variables (columns in the linear constraint matrix)
indexed according to `indvals`.
The parameter `constrstate` is, upon first call, the value returned by [`add_constr_fix_prepare!`](@ref); and on all further
calls, it will be the return value of the previous call.
Note that `rhs` will almost always be zero, so if the right-hand side is represented by a sparse vector, it is worth checking
for this value (the compiler will be able to remove the check).

See also [`Indvals`](@ref).
"""
function add_constr_fix! end

"""
    add_constr_fix_finalize!(state::AbstractSolver, constrstate)

Finishes the addition of fixed constraints to `state`; the value of `constrstate` is the return value of the last call to
[`add_constr_fix!`](@ref).
The default implementation does nothing.
"""
add_constr_fix_finalize!(::AbstractSolver, _) = nothing

"""
    negate_fix(::AbstractSolver)

Depending on the exact definition of equality constraints (where to put the minus), the dual solutions, i.e., the SOS
decomposition, may yield wrong values; then define this function to return `true`, which will pass the equality constraints
with opposite sign.
"""
negate_fix(::AbstractSolver) = false

"""
    prepend_fix(::AbstractSolver)

If this method yields `true` (default), fixed constraints will be created before all others. If it is `false`, PSD (or
nonnegative) constraints will be created first.

!!! info
    Note that this does not imply a global order; if DD or SDD cones are used, fixed constraints may still appear at an
    arbitrary position. However, this method guarantees the order with respect to cones that will use the same monomial
    indices.
"""
prepend_fix(::AbstractSolver) = true

"""
    fix_objective!(state, indvals::Indvals)

Puts the variables indexed according to `indvals` into the objective (that is to be minimized).
This function will be called exactly once by [`moment_setup!`](@ref) after all variables and constraints have been set up.

See also [`Indvals`](@ref).
"""
function fix_objective! end

"""
    add_var_slack!(state::AbstractSolver{T}, num::Int)

Creates `num` slack variables in the problem. Slack variables must be free. The result should be an abstract vector (typically
a unit range) that contains the indices of all created slack variables. The indices should be of type `T`.
"""
function add_var_slack! end