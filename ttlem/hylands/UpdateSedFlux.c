/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: UpdateSedFlux.c
 *
 * MATLAB Coder version            : 4.3
 * C/C++ source code generated on  : 13-Mar-2020 23:17:35
 */

/* Include Files */
#include "UpdateSedFlux.h"
#include "UpdateSedFlux_emxutil.h"

/* Function Definitions */

/*
 * Arguments    : const emxArray_uint32_T *i
 *                const emxArray_uint32_T *k
 *                double phi
 *                const emxArray_real_T *E_reg
 *                double dx2
 *                double Ff
 *                const emxArray_real_T *E_bed
 *                const emxArray_real_T *V
 *                const emxArray_real_T *Q
 *                emxArray_real_T *Qs
 *                emxArray_real_T *Qs_in
 * Return Type  : void
 */
void UpdateSedFlux(const emxArray_uint32_T *i, const emxArray_uint32_T *k,
                   double phi, const emxArray_real_T *E_reg, double dx2, double
                   Ff, const emxArray_real_T *E_bed, const emxArray_real_T *V,
                   const emxArray_real_T *Q, emxArray_real_T *Qs,
                   emxArray_real_T *Qs_in)
{
  unsigned int unnamed_idx_0;
  unsigned int unnamed_idx_1;
  int b_i;
  int loop_ub;
  int i1;
  int i2;
  unnamed_idx_0 = (unsigned int)E_bed->size[0];
  unnamed_idx_1 = (unsigned int)E_bed->size[1];
  b_i = Qs->size[0] * Qs->size[1];
  Qs->size[0] = (int)unnamed_idx_0;
  Qs->size[1] = (int)unnamed_idx_1;
  emxEnsureCapacity_real_T(Qs, b_i);
  loop_ub = (int)unnamed_idx_0 * (int)unnamed_idx_1;
  for (b_i = 0; b_i < loop_ub; b_i++) {
    Qs->data[b_i] = 0.0;
  }

  unnamed_idx_0 = (unsigned int)E_bed->size[0];
  unnamed_idx_1 = (unsigned int)E_bed->size[1];
  b_i = Qs_in->size[0] * Qs_in->size[1];
  Qs_in->size[0] = (int)unnamed_idx_0;
  Qs_in->size[1] = (int)unnamed_idx_1;
  emxEnsureCapacity_real_T(Qs_in, b_i);
  loop_ub = (int)unnamed_idx_0 * (int)unnamed_idx_1;
  for (b_i = 0; b_i < loop_ub; b_i++) {
    Qs_in->data[b_i] = 0.0;
  }

  b_i = i->size[0];
  for (loop_ub = 0; loop_ub < b_i; loop_ub++) {
    i1 = (int)i->data[loop_ub] - 1;
    Qs->data[i1] = ((Qs_in->data[i1] + (1.0 - phi) * E_reg->data[i1] * dx2) +
                    (1.0 - Ff) * E_bed->data[i1] * dx2) / (V->data[i1] * dx2 /
      Q->data[i1] + 1.0);
    i2 = (int)k->data[loop_ub] - 1;
    Qs_in->data[i2] += Qs->data[i1];
    Qs_in->data[i1] = 0.0;
  }
}

/*
 * File trailer for UpdateSedFlux.c
 *
 * [EOF]
 */
