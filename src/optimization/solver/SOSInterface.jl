export add_var_nonnegative!, add_var_rotated_quadratic!, add_var_quadratic!, add_var_l1!, add_var_l1_complex!, add_var_psd!,
    add_var_psd_complex!, add_var_dd!, add_var_free_prepare!, add_var_free!, add_var_free_finalize!, negate_free,
    fix_constraints!, add_constr_slack!

function add_var_nonnegative! end

"""
    add_var_nonnegative!(state::AbstractSolver{T,V}, indvals::Indvals{T,V}) where {T,V}

Add a nonnegative decision variable to the solver and put its value into the linear constraints (rows in the linear constraint
matrix) indexed according to `indvals`.
Falls back to the vector-valued version if not implemented.

See also [`Indvals`](@ref).
"""
add_var_nonnegative!(state::AbstractSolver{T,V}, indvals::Indvals{T,V}) where {T,V} =
    add_var_nonnegative!(state, IndvalsIterator(unsafe, indvals.indices, indvals.values, StackVec(length(indvals))))

"""
    add_var_nonnegative!(state::AbstractSolver{T,V}, indvals::IndvalsIterator{T,V}) where {T,V}

Add multiple nonnegative decision variables to the solver and put their values into the linear constraints (rows in the linear
constraint matrix) indexed according to the entries in `indvals`.
Falls back to calling the scalar-valued version multiple times if not implemented.

See also [`IndvalsIterator`](@ref).
"""
function add_var_nonnegative!(state::AbstractSolver{T,V}, iv::IndvalsIterator{T,V}) where {T,V}
    for indvals in iv
        add_var_nonnegative!(state, indvals)
    end
    return
end

function add_var_quadratic! end

@doc raw"""
    add_var_quadratic!(state::AbstractSolver{T,V}, indvals::IndvalsIterator{T,V}) where {T,V}

Adds decision variables in a quadratic cone to the solver and put their values into the linear constraints (rows in the linear
constraint matrix), indexed according to `indvals`. The `N = length(indvals)` variables will satisfy ``x_1 \geq 0``,
``x_1^2 \geq \sum_{i = 2}^N x_i^2``.

See also [`Indvals`](@ref), [`IndvalsIterator`](@ref).

!!! note "Number of parameters"
    In the real-valued case, `indvals` is always of length three, in the complex case, it is of length four. If the scaled
    diagonally dominant representation is requested, `indvals` can have any length.

!!! warning
    This function will only be called if [`supports_quadratic`](@ref) returns `true` for the given state.
    If (rotated) quadratic constraints are unsupported, a fallback to a 2x2 PSD variable is used.
"""
add_var_quadratic!(::AbstractSolver{T,V}, ::IndvalsIterator{T,V}) where {T,V}

function add_var_rotated_quadratic! end

@doc raw"""
    add_var_rotated_quadratic!(state::AbstractSolver{T,V}, indvals::IndvalsIterator{T,V}) where {T,V}

Adds decision variables in a rotated quadratic cone to the solver and put their values into the linear constraints (rows in
the linear constraint matrix), indexed according to `indvals`. The `N = length(indvals)` variables will satisfy
``x_1, x_2 \geq 0``, ``2x_1 x_2 \geq \sum_{i = 3}^N x_i^2``.

See also [`Indvals`](@ref), [`IndvalsIterator`](@ref).

!!! note "Number of parameters"
    In the real-valued case, `indvals` is always of length three, in the complex case, it is of length four. If the scaled
    diagonally dominant representation is requested, `indvals` can have any length.

!!! warning
    This function will only be called if [`supports_quadratic`](@ref) returns `true` for the given state.
    If (rotated) quadratic constraints are unsupported, a fallback to a 2x2 PSD variable is used.
"""
add_var_rotated_quadratic!(::AbstractSolver{T,V}, ::IndvalsIterator{T,V}) where {T,V}

function add_var_psd! end

"""
    add_var_psd!(state::AbstractSolver{T,V}, dim::Int, data::PSDMatrixCartesian{T,V}) where {T,V}

Add a PSD variable of side dimension `dim` ≥ 3 to the solver. Its requested triangle is indexed according to the return value
of [`psd_indextype`](@ref)); these elements of the matrix are put into the linear constraints (rows in the linear constraint
matrix) indicated by the keys when iterating through `data`, which are of the type `T`, at positions and with coefficients
given by their values.
Note that if [`add_var_quadratic!`](@ref) is not implemented, `dim` may also be `2`.
This method is called if [`psd_indextype`](@ref) returns a [`PSDIndextypeMatrixCartesian`](@ref).

!!! hint "Complex-valued PSD variables"
    Note that this function will also be called for complex-valued PSD cones if [`supports_psd_complex`](@ref) returns `false`.
    The data will have been rewritten in terms of a real-valued PSD cone, which doubles the dimension.
    If the solver natively supports complex-valued PSD cones, [`add_var_psd_complex!`](@ref) must be implemented.
"""
add_var_psd!(::AbstractSolver{T,V}, ::Int, ::PSDMatrixCartesian{T,V}) where {T,V}

"""
    add_var_psd!(state::AbstractSolver{T,V}, dim::Int, data::IndvalsIterator{T,V}) where {T,V}

Conceptually the same as above; but now, `data` is an iterable through the elements of the PSD variable one-by-one. The
individual entries are [`Indvals`](@ref).
This method is called if [`psd_indextype`](@ref) returns a [`PSDIndextypeVector`](@ref).

!!! hint "Complex-valued PSD variables"
    Note that this function will also be called for complex-valued PSD cones if [`supports_psd_complex`](@ref) returns `false`.
    The data will have been rewritten in terms of a real-valued PSD cone, which doubles the dimension.
    If the solver natively supports complex-valued PSD cones, [`add_var_psd_complex!`](@ref) must be implemented.
"""
add_var_psd!(::AbstractSolver{T,V}, ::Int, ::IndvalsIterator{T,V}) where {T,V}

function add_var_psd_complex! end

"""
    add_var_psd_complex!(state::AbstractSolver{T,V}, dim::Int, data::PSDMatrixCartesian{T,Complex{V}}) where {T,V}

Add a Hermitian PSD variable of side dimension `dim` ≥ 3 to the solver. Its requested triangle is indexed according to the
return value of [`psd_indextype`](@ref)); these elements of the matrix are put into the linear constraints (rows in the linear
constraint matrix) indicated by the keys when iterating through `data`, which are of the type `T`, at positions and with
coefficients given by their values. The real part of the coefficient corresponds to the coefficient in front of the real part
of the matrix entry, the imaginary part is the coefficient for the imaginary part of the matrix entry.
Note that if [`add_var_quadratic!`](@ref) is not implemented, `dim` may also be `2`.
This method is called if [`psd_indextype`](@ref) returns a [`PSDIndextypeMatrixCartesian`](@ref).

!!! warning
    This function will only be called if [`supports_psd_complex`](@ref) is defined to return `true` for the given state.
"""
add_var_psd_complex!(::AbstractSolver{T,V}, ::Int, ::PSDMatrixCartesian{T,Complex{V}}) where {T,V}

"""
    add_var_psd_complex!(state::AbstractSolver{T,V}, dim::Int, data::IndvalsIterator{T,V}) where {T,V}

Conceptually the same as above; but now, `data` is an iterable through the elements of the PSD variable one-by-one. The
individual entries are [`Indvals`](@ref).
This method is called if [`psd_indextype`](@ref) returns a [`PSDIndextypeVector`](@ref).
Regardless of the travelling order, for diagonal elements, there will be exactly one entry, which is the real part. For
off-diagonal elements, the real part will be followed by the imaginary part. Therefore, the coefficients are real-valued.

!!! warning
    This function will only be called if [`supports_psd_complex`](@ref) is defined to return `true` for the given state.
"""
add_var_psd_complex!(::AbstractSolver{T,V}, ::Int, ::IndvalsIterator{T,V}) where {T,V}

function add_var_dd! end

@doc raw"""
    add_var_dd!(state::AbstractSolver{T,V}, dim::Integer, data::IndvalsIterator{T,V}, u) where {T,V}

Add a constraint for membership in the cone of diagonally dominant matrices to the solver. `data` is an iterator through the
scaled lower triangle of the matrix. A basis change is induced by `u`, with the meaning that `M ∈ DD(u) ⇔ M = uᵀ Q u` with
`Q ∈ DD`.

!!! warning
    This function will only be called if [`supports_dd`](@ref) returns `true` for the given state. If diagonally dominant cones
    are not supported directly, a fallback to a columnwise representation in terms of ``\ell_1`` norms will be used (or the
    fallbacks if this norm is not supported).
"""
add_var_dd!(::AbstractSolver{T,V}, ::Integer, ::IndvalsIterator{T,V}, u) where {T,V}

function add_var_dd_complex! end

@doc raw"""
    add_var_dd_complex!(state::AbstractSolver{T,V}, dim::Integer, data::IndvalsIterator{T,V}, u) where {T,V}

Add a constraint for membership in the cone of complex-valued diagonally dominant matrices to the solver. `data` is an iterator
hrough the scaled lower triangle of the matrix. A basis change is induced by `u`, with the meaning that
`M ∈ DD(u) ⇔ M = u† Q u` with `Q ∈ DD`.
For diagonal elements, there will be exactly one entry, which is the real part. For off-diagonal elements, the real part will
be followed by the imaginary part. Therefore, the coefficients are real-valued.

!!! warning
    This function will only be called if [`supports_dd_complex`](@ref) returns `true` for the given state. If complex-valued
    diagonally dominant cones are not supported directly, a fallback to quadratic cones on the complex-valued data is tried
    first (if supported), followed by a columnwise representation in terms of ``\ell_1`` norms or their fallback on the
    realification of the matrix data if not.
"""
add_var_dd_complex!(::AbstractSolver{T,V}, ::Integer, ::IndvalsIterator{T,V}, u) where {T,V}

function add_var_l1! end

@doc raw"""
    add_var_l1!(state::AbstractSolver{T,V}, indvals::IndvalsIterator{T,V}) where {T,V}

Adds decision variables in an ``\ell_1`` norm cone to the solver and put their values into the linear constraints (rows in
the linear constraint matrix), indexed according to the `indvals`. The `N = length(indvals)` variables will satisfy
``x_1 \geq \sum_{i = 2}^N \lvert x_i\rvert``.

See also [`Indvals`](@ref), [`IndvalsIterator`](@ref).

!!! warning
    This function will only be called if [`supports_lnorm`](@ref) returns `true` for the given state.
    If ``\ell_\infty`` norm cones are unsupported, a fallback to multiple nonnegative variables will be used.
"""
add_var_l1!(::AbstractSolver{T,V}, ::IndvalsIterator{T,V}) where {T,V}

function add_var_l1_complex! end

@doc raw"""
    add_var_l1_complex!(state::AbstractSolver{T,V}, indvals::IndvalsIterator{T,V}) where {T,V}

Same as [`add_var_l1!`](@ref), but now two successive items in `indvals` (starting from the second) are interpreted as
determining the real and imaginary part of a component of the ``\ell_1`` norm variable.

!!! warning
    This function will only be called if [`supports_lnorm_complex`](@ref) returns `true` for the given state.
    If complex-valued ``\ell_1`` norm cones are unsupported, a fallback to multiple nonnegative and quadratic variables will be
    used.
"""
add_var_l1_complex!(::AbstractSolver{T,V}, ::IndvalsIterator{T,V}) where {T,V}

function add_var_sdd! end

@doc raw"""
    add_var_sdd!(state::AbstractSolver{T,V}, dim::Integer, data::IndvalsIterator{T,V}, u) where {T,V}

Add a constraint for membership in the cone of scaled diagonally dominant matrices to the solver. `data` is an iterator through
the (unscaled) lower triangle of the matrix. A basis change is induced by `u`, with the meaning that `M ∈ SDD(u) ⇔ M = uᵀ Q u`
with `Q ∈ SDD`.

!!! warning
    This function will only be called if [`supports_sdd`](@ref) returns `true` for the given state. If scaled diagonally
    dominant cones are not supported directly, a fallback to (rotated) quadratic cones will be used.
"""
add_var_sdd!(::AbstractSolver{T,V}, ::Integer, ::IndvalsIterator{T,V}, u) where {T,V}

function add_var_sdd_complex! end

@doc raw"""
    add_var_sdd_complex!(state::AbstractSolver{T,V}, dim::Integer, data::IndvalsIterator{T,V}, u) where {T,V}

Add a constraint for membership in the cone of complex-valued scaled diagonally dominant matrices to the solver. `data` is an
iterator through the (unscaled) lower triangle of the matrix. A basis change is induced by `u`, with the meaning that
`M ∈ SDD(u) ⇔ M = u† Q u` with `Q ∈ SDD`.
For diagonal elements, there will be exactly one entry, which is the real part. For off-diagonal elements, the real part will
be followed by the imaginary part. Therefore, the coefficients are real-valued.

!!! warning
    This function will only be called if [`supports_sdd_complex`](@ref) returns `true` for the given state. If complex-valued
    scaled diagonally dominant cones are not supported directly, a fallback to quadratic cones is automatically performed.
"""
add_var_sdd_complex!(::AbstractSolver{T,V}, ::Integer, ::IndvalsIterator{T,V}, u) where {T,V}

"""
    add_var_free_prepare!(state::AbstractSolver, num::Int)

Prepares to add exactly `num` free variables that may become part of the objective; the actual data is then put into the solver
by subsequent calls of [`add_var_free!`](@ref) and the whole transaction is completed by [`add_var_free_finalize!`](@ref).
The return value of this function is passed on as `eqstate` to [`add_var_free!`](@ref).
The default implementation does nothing.
"""
add_var_free_prepare!(::AbstractSolver, _) = nothing

"""
    add_var_free!(state::AbstractSolver{T,V}, eqstate, indvals::Indvals{T,V}, obj::V) where {T,V}

Add a free variable to the solver and put its value into the linear constraints (rows in the linear constraint matrix), indexed
according to `indvals`.
The variable should also be put into the objective with coefficient `obj` (which is likely to be zero).
The parameter `eqstate` is, upon first call, the value returned by [`add_var_free_prepare!`](@ref); and on all further calls,
it will be the return value of the previous call.

See also [`Indvals`](@ref).
"""
function add_var_free! end

"""
    add_var_free_finalize!(state::AbstractSolver, eqstate)

Finishes the addition of free variables to `state`; the value of `eqstate` is the return value of the last call to
[`add_var_free!`](@ref).
The default implementation does nothing.
"""
add_var_free_finalize!(::AbstractSolver, _) = nothing

"""
    negate_free(::AbstractSolver)

Depending on the exact definition of equality constraints (where to put the minus), the dual solutions, i.e., the SOS
decomposition, may yield wrong values; then define this function to return `true`, which will flip the sign of free variables.
"""
negate_free(::AbstractSolver) = false

"""
    fix_constraints!(state::AbstractSolver{T,V}, indvals::Indvals{T,V}) where {T,V}

Ensures that all constraints in the optimization problem are fixed to the values according to `indvals`.
This function will be called exactly once by [`sos_setup!`](@ref) after all variables and constraints have been set up.

See also [`Indvals`](@ref).
"""
function fix_constraints! end

"""
    add_constr_slack!(state::AbstractSolver{T}, num::Int)

Creates `num` linear fix-to-zero slack constraints in the problem (i.e., constraints that do not correspond to moments). The
result should be an abstract vector (typically a unit range) that contains the indices of type `T` of all created slack
constraints.
"""
function add_constr_slack! end