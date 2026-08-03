FUNCTION DIFF_ROT2, X, P, ERR
;+
; Name: DIFF_ROT2
; Purpose: Evaluate the asymmetric fourth-order differential-rotation model
;   A + B sin(latitude) + C sin^2(latitude) + D sin^3(latitude)
;   + E sin^4(latitude).
; Calling Sequence: omega = DIFF_ROT2(latitude, parameters, error)
; Inputs:
;   X   - Latitude in degrees; scalar or array.
;   P   - Five-element coefficient vector [A, B, C, D, E].
;   ERR - Placeholder required by the MPFITFUN model-function interface.
; Return Value: Model angular rotation rate with the same shape as X.
; Side Effects: None.
;-

  IF N_ELEMENTS(P) LT 5 THEN MESSAGE, $
    'DIFF_ROT2 requires parameters [A, B, C, D, E].'

  SIN_LATITUDE = SIN(X * !DTOR)
  RETURN, P[0] + P[1] * SIN_LATITUDE + P[2] * SIN_LATITUDE^2 + $
    P[3] * SIN_LATITUDE^3 + P[4] * SIN_LATITUDE^4
END
