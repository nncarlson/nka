#if defined(NDEBUG)
#define ASSERT(x) !! Assert( x )
#else
#define ASSERT(x) if (.not.(x)) call f90_assert(__FILE__,__LINE__)
#endif
