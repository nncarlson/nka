/* nonlinear_krylov_accelerator.h */

typedef struct nka_state *NKA;
extern NKA nka_init (int, int, double, double (*dp)(int, double *, double *));
extern void nka_delete (NKA);
extern void nka_accel_update (NKA, double *);
extern void nka_restart (NKA);
extern void nka_relax (NKA);
extern int nka_num_vec (NKA);
extern int nka_max_vec (NKA);
extern int nka_vec_len (NKA);
extern double nka_vec_tol (NKA);
