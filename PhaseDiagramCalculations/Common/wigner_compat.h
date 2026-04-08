#ifndef POLLEN_WIGNER_COMPAT_H
#define POLLEN_WIGNER_COMPAT_H

// Minimal compatibility surface for the subset of WIGXJPF used here.
void wig_table_init(int max_two_j, int wigner_type);
void wig_temp_init(int max_two_j);
void wig_thread_temp_init(int max_two_j);
double wig3jj(int two_j1, int two_j2, int two_j3, int two_m1, int two_m2, int two_m3);
void wig_temp_free();
void wig_table_free();

#endif
