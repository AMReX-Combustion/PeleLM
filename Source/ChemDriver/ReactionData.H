#ifndef _ReactionData_
#define _ReactionData_

typedef enum
{
  FWD_A, FWD_BETA, FWD_EA,
  LOW_A, LOW_BETA, LOW_EA,
  REV_A, REV_BETA, REV_EA,
  TROE_A, TROE_TS, TROE_TSS, TROE_TSSS,
  SRI_A, SRI_B, SRI_C, SRI_D, SRI_E,
  THIRD_BODY
} REACTION_PARAMETER;

#endif
