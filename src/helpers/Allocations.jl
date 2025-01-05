if Sys.iswindows()
    mutable struct ProcessMemoryCountersEx
        cb::Cuint
        PageFaultCount::Cuint
        PeakWorkingSetSize::Csize_t
        WorkingSetSize::Csize_t
        QuotaPeakPagedPoolUsage::Csize_t
        QuotaPagedPoolUsage::Csize_t
        QuotaPeakNonPagedPoolUsage::Csize_t
        QuotaNonPagedPoolUsage::Csize_t
        PagefileUsage::Csize_t
        PeakPagefileUsage::Csize_t
        PrivateUsage::Csize_t

        ProcessMemoryCountersEx() = new()
    end

    function get_own_memory()
        Process = ccall((:GetCurrentProcess, "kernel32"), stdcall, Ptr{Cvoid}, ())
        # proc is a pseudo handle, should be -1
        ppsmemCounters = ProcessMemoryCountersEx()
        iszero(ccall((:GetProcessMemoryInfo, "Psapi"), stdcall, Cint,
            (Ptr{Cvoid}, Ref{ProcessMemoryCountersEx}, Cuint),
            Process, ppsmemCounters, sizeof(ppsmemCounters)
        )) && error("Unable to obtain process memory info")
        return ppsmemCounters.PrivateUsage
    end
elseif Sys.isapple()
    # completely untested, I don't have a Mac
    struct TimeValue
        seconds::Cint
        microseconds::Cint
    end

    mutable struct MachTaskBasicInfo
        virtual_size::Culonglong # mach_vm_size_t
        resident_size::Culonglong # mach_vm_size_t
        resident_size_max::Culonglong # resident_size_max
        user_time::TimeValue
        system_time::TimeValue
        policy::Cint # policy_t
        suspend_count::Cint # integer_t

        TaskBasicInfo() = new()
    end

    function get_own_memory()
        task_name = ccall((:mach_task_self, "kern"), Ptr{Cvoid}, ()) # returns pointer to ipc_port structure
        t_info = MachTaskBasicInfo()
        t_info_count = Cuint(sizeof(t_info) ÷ sizeof(Cuint)) # MACH_TASK_BASIC_INFO_COUNT
        iszero(ccall((:task_info, "kern"), Cint, # kern_return_t
            (Ptr{Cvoid}, # task_name_t = task*
             Cuint, # task_flavor_t = natural_t
             Ptr{Cint}, # task_into_t = integer_t*
             Ref{Cuint} # mach_msg_type_number_t* = natural_t*
            ),
            task_name,
            20, # MACH_TASK_BASIC_INFO
            Ref(t_info),
            Ref(t_info_count)
        )) || error("Unable to obtain task info")
        return t_info.virtual_size
    end
elseif Sys.isunix()
    function get_own_memory()
        f = open("/proc/self/status", "r")
        try
            readuntil(f, "VmSize:\t")
            return parse(UInt, lstrip(readline(f))[1:end-3]) << 10 # should end in "kB"
        finally
            close(f)
        end
    end
else
    function get_own_memory()
        error("get_own_memory not implemented for the current operating system")
    end
end

"""
    @allocdiff

A macro to evaluate an expression, discarding the resulting value, instead returning the difference in the number of bytes
allocated after vs. before evaluation of the expression (which is not guaranteed to be the peak, but allows to capture
allocations done in third-party libraries that don't use Julia's GC, contrary to `@allocated`).
In order to provide consistent results, Julia's GC is disabled while the expression is evaluated.
Note that on *nix systems, the value is only accurate up to the kibibyte.
"""
macro allocdiff(ex)
    quote
        local oldgc = GC.enable(false)
        local b0 = get_own_memory()
        $(esc(ex))
        local result = reinterpret(Int, get_own_memory() - b0)
        GC.enable(oldgc)
        result
    end
end