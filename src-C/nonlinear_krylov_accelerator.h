/* nonlinear_krylov_accelerator.h */

typedef struct nka_state *NKA;
extern NKA nka_create (int, int, double, double (*dp)(int, double *, double *));
extern void nka_destroy (NKA);
extern void nka_correction (NKA, double *);
extern void nka_restart (NKA);
extern void nka_relax (NKA);
