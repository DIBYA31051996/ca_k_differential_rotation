FUNCTION DIFF_ROT1, X, P, ERR
;+
; Name: DIFF_ROT1
; Purpose: Evaluate the symmetric differential-rotation model
;   A + B sin^2(latitude) + C sin^4(latitude).
; Calling Sequence: omega = DIFF_ROT1(latitude, parameters, error)
; Inputs:
;   X   - Latitude in degrees; scalar or array.
;   P   - Three-element coefficient vector [A, B, C].
;   ERR - Placeholder required by the MPFITFUN model-function interface.
; Return Value: Model angular rotation rate with the same shape as X.
; Side Effects: None.
;-

  IF N_ELEMENTS(P) LT 3 THEN MESSAGE, $
    'DIFF_ROT1 requires parameters [A, B, C].'

  SIN_LATITUDE = SIN(X * !DTOR)
  RETURN, P[0] + P[1] * SIN_LATITUDE^2 + P[2] * SIN_LATITUDE^4
END
