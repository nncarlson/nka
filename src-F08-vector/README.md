The other versions of `nka_type.F90` operate on mathematical vectors stored
as contiguous rank-1 arrays. This is simple, but frequently awkward to use in
practice. Applications may naturally store their vectors in multi-rank arrays
and/or may include ghost or halo elements and thus be non-contiguous. To use
those versions it requires copying vectors between this natural storage format
and contiguous rank-1 arrays. This can be seen in `src-F08/nka_example.F90`.

This directory contains an alternative implementation of `nka_type.F90` that
operates on mathematical vectors stored as polymorphic objects of the abstract
base class `vector` defined in `vector_class.F90`. The base class has methods
for the usual vector space operations, including dot product and norm reduction
operations. To use this version an application must define and use an extention
of the `vector` base class that is adapted to the needs of the application.
The file `grid_vector.F90` does exactly this for the `nka_example.F90` program.

In a parallel application where its vectors are distributed across processing
elements, all versions of `nka_type.F90` work without modification when called
collectively, passing each processing element's portion of the vector to the
update procedure. With the other versions, however, a parallel-aware dot
product procedure must be explicitly provided to NKA. But in this version this
is not necessary as the implementation of the vector base class reduction
methods will necessarily be parallel-aware.
