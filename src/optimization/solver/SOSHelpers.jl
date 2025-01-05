export sos_add_matrix!, sos_add_equality!, sos_setup!

struct SOSWrapper{S<:AbstractSolver} # We don't make this an AbstractSolver, which allows to disambiguate all our substitution
                                     # methods without actually specifying all the arguments in detail. However, this means
                                     # that we have to allow the setup methods to take a Union of AbstractSolver and
                                     # SOSWrapper.
    state::S
end

const AnySolver{T,V} = Union{<:AbstractSolver{T,V},SOSWrapper{<:AbstractSolver{T,V}}}

mindex(state::SOSWrapper, args...) = mindex(state.state, args...)

for mapfn in (
    :supports_rotated_quadratic, :supports_quadratic, :supports_lnorm, :supports_lnorm_complex, :supports_psd_complex,
    :supports_dd, :supports_dd_complex, :supports_sdd, :supports_sdd_complex,
    :psd_indextype
)
    @eval $mapfn(state::SOSWrapper) = $mapfn(state.state)
end

for mapfn in (
    :nonnegative!, :rotated_quadratic!, :quadratic!, :psd!, :psd_complex!
)
    @eval $(Symbol(:add_constr_, mapfn))(state::SOSWrapper, args...) = $(Symbol(:add_var_, mapfn))(state.state, args...)
end

for (newname, oldname) in (
    (:add_var_dd!, :add_constr_dddual!), (:add_var_dd_complex!, :add_constr_dddual_complex!),
    (:add_var_l1!, :add_constr_linf!), (:add_var_l1_complex!, :add_constr_linf_complex!),
    (:add_var_sdd!, :add_constr_sdddual!), (:add_var_sdd_complex!, :add_constr_sdddual_complex!),
    (:add_var_free_prepare!, :add_constr_fix_prepare!), (:add_var_free!, :add_constr_fix!),
    (:add_var_free_finalize!, :add_constr_fix_finalize!), (:negate_free, :negate_fix),
    (:add_constr_slack!, :add_var_slack!),
    (:fix_constraints!, :fix_objective!),
    (:addtocounter!, :addtocounter!)
)
    @eval $oldname(state::SOSWrapper, args...) = $newname(state.state, args...)
end

Base.getproperty(c::Counters, f::Symbol) = getfield(c, f === :free ? :fix : f)
Base.setproperty!(c::Counters, f::Symbol, v) = setfield!(c, f === :free ? :fix : f, v)

addtocounter!(state::SOSWrapper, counters::Counters, ::Val{:fix}, args...) =
    addtocounter!(state.state, counters, Val(:free), args...)

"""
    sos_add_matrix!(state::AbstractSolver, grouping::SimpleMonomialVector,
        constraint::Union{<:SimplePolynomial,<:AbstractMatrix{<:SimplePolynomial}},
        representation::RepresentationMethod=RepresentationPSD())

Parses a SOS constraint with a basis given in `grouping` (this might also be a partial basis due to sparsity), premultiplied by
`constraint` (which may be the unit polynomial for the SOS cone membership) and calls the appropriate solver functions to set
up the problem structure according to `representation`.

To make this function work for a solver, implement the following low-level primitives:
- [`mindex`](@ref)
- [`add_var_nonnegative!`](@ref)
- [`add_var_rotated_quadratic!`](@ref) (optional, then set [`supports_rotated_quadratic`](@ref) to `true`)
- [`add_var_quadratic!`](@ref) (optional, then set [`supports_quadratic`](@ref) to `true`)
- [`add_var_psd!`](@ref)
- [`add_var_psd_complex!`](@ref) (optional, then set [`supports_psd_complex`](@ref) to `true`)
- [`add_var_dd!`](@ref) (optional, then set [`supports_dd`](@ref) to `true`)
- [`add_var_dd_complex!`](@ref) (optional, then set [`supports_dd_complex`](@ref) to `true`)
- [`add_var_l1!`](@ref) (optional, then set [`supports_lnorm`](@ref) to `true`)
- [`add_var_l1_complex!`](@ref) (optional, then set [`supports_lnorm_complex`](@ref) to `true`)
- [`add_var_sdd!`](@ref) (optional, then set [`supports_sdd`](@ref) to `true`)
- [`add_var_sdd_complex!`](@ref) (optional, then set [`supports_sdd_complex`](@ref) to `true`)
- [`psd_indextype`](@ref)
- [`add_constr_slack!`](@ref)

Usually, this function does not have to be called explicitly; use [`sos_setup!`](@ref) instead.

See also [`sos_add_equality!`](@ref).
"""
sos_add_matrix!(state::AbstractSolver, args...) = moment_add_matrix!(SOSWrapper(state), args...)

"""
    sos_add_equality!(state::AbstractSolver, grouping::SimpleMonomialVector, constraint::SimplePolynomial)

Parses a polynomial equality constraint for sums-of-squares and calls the appropriate solver functions to set up the problem
structure. `grouping` contains the basis that will be squared in the process to generate the prefactor.

To make this function work for a solver, implement the following low-level primitives:
- [`add_var_free_prepare!`](@ref) (optional)
- [`add_var_free!`](@ref) (required)
- [`add_var_free_finalize!`](@ref) (optional)

Usually, this function does not have to be called explicitly; use [`sos_setup!`](@ref) instead.

See also [`sos_add_matrix!`](@ref).
"""
sos_add_equality!(state::AbstractSolver, args...) = moment_add_equality!(SOSWrapper(state), args...)

"""
    sos_setup!(state::AbstractSolver, relaxation::AbstractRelaxation, groupings::RelaxationGroupings[; representation])

Sets up all the necessary SOS matrices, free variables, objective, and constraints of a polynomial optimization problem
`problem` according to the values given in `grouping` (where the first entry corresponds to the basis of the objective, the
second of the equality, the third of the inequality, and the fourth of the PSD constraints).
The function returns a `Vector{<:Vector{<:Tuple{Symbol,Any}}}` that contains internal information on the problem. This
information is required to obtain dual constraints and re-optimize the problem and should be stored in the `state`.

The following methods must be implemented by a solver to make this function work:
- [`mindex`](@ref)
- [`add_var_nonnegative!`](@ref)
- [`add_var_rotated_quadratic!`](@ref) (optional, then set [`supports_rotated_quadratic`](@ref) to `true`)
- [`add_var_quadratic!`](@ref) (optional, then set [`supports_quadratic`](@ref) to `true`)
- [`add_var_psd!`](@ref)
- [`add_var_psd_complex!`](@ref) (optional, then set [`supports_psd_complex`](@ref) to `true`)
- [`add_var_dd!`](@ref) (optional, then set [`supports_dd`](@ref) to `true`)
- [`add_var_dd_complex!`](@ref) (optional, then set [`supports_dd_complex`](@ref) to `true`)
- [`add_var_l1!`](@ref) (optional, then set [`supports_lnorm`](@ref) to `true`)
- [`add_var_l1_complex!`](@ref) (optional, then set [`supports_lnorm_complex`](@ref) to `true`)
- [`add_var_sdd!`](@ref) (optional, then set [`supports_sdd`](@ref) to `true`)
- [`add_var_sdd_complex!`](@ref) (optional, then set [`supports_sdd_complex`](@ref) to `true`)
- [`psd_indextype`](@ref)
- [`add_var_free_prepare!`](@ref) (optional)
- [`add_var_free!`](@ref)
- [`add_var_free_finalize!`](@ref) (optional)
- [`fix_constraints!`](@ref)
- [`add_constr_slack!`](@ref)

!!! warning "Indices"
    The constraint indices used in all solver functions directly correspond to the indices given back by [`mindex`](@ref).
    However, in a sparse problem there may be far fewer indices present; therefore, when the problem is finally given to the
    solver, care must be taken to eliminate all unused indices. The functionality provided by [`AbstractAPISolver`](@ref) and
    [`AbstractSparseMatrixSolver`](@ref) already takes care of this.

!!! info "Order"
    The individual variable types can be added in any order (including interleaved).

!!! info "Representation"
    This function may also be used to describe simplified cones such as the (scaled) diagonally dominant one. The
    `representation` parameter can be used to define a representation that is employed for the individual groupings. This may
    either be an instance of a [`RepresentationMethod`](@ref) - which requires the method to be independent of the dimension of
    the grouping - or a callable. In the latter case, it will be passed as a first parameter an identifier[^3] of the current
    conic variable, and as a second parameter the side dimension of its matrix. The method must then return a
    [`RepresentationMethod`](@ref) instance.

See also [`moment_setup!`](@ref), [`sos_add_matrix!`](@ref), [`sos_add_equality!`](@ref),
[`RepresentationMethod`](@ref).
"""
sos_setup!(state::AbstractSolver, args...; kwargs...) = moment_setup!(SOSWrapper(state), args...; kwargs...)