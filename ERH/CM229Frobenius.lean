import Std

/-! Finite Frobenius actions on the six embeddings of a cubic CM compositum.
Each profile lists the multiplicities for six times the slope in [0, 2, 3, 4, 6]. -/

namespace CM229Frobenius

abbrev Embedding := Nat × Bool

def embeddings : List Embedding :=
  [(0, false), (0, true), (1, false), (1, true), (2, false), (2, true)]

def frobenius (perm : List Nat) (flip : Bool) (e : Embedding) : Embedding :=
  (perm[e.1]?.getD e.1, if flip then !e.2 else e.2)

def visits (perm : List Nat) (flip : Bool) (mask : Nat) :
    Nat → Embedding → Nat
  | 0, _ => 0
  | n + 1, e =>
    (if ((mask / 2 ^ e.1 % 2 == 1) == e.2) then 1 else 0) +
      visits perm flip mask n (frobenius perm flip e)

def profile (perm : List Nat) (flip : Bool) (mask : Nat) : List Nat :=
  let slopes := embeddings.map (visits perm flip mask 6)
  [0, 2, 3, 4, 6].map fun s => slopes.count s

def profiles : List (List Nat) :=
  [[3, 0, 0, 0, 3], [0, 0, 6, 0, 0], [2, 0, 2, 0, 2],
   [1, 0, 4, 0, 1], [0, 3, 0, 3, 0]]

def histogram (perm : List Nat) (flip : Bool) : List Nat :=
  profiles.map fun s => ((List.range 8).filter fun m => profile perm flip m == s).length

def permutations : List (List Nat) :=
  [[0, 1, 2], [1, 0, 2], [2, 1, 0], [0, 2, 1], [1, 2, 0], [2, 0, 1]]

set_option maxRecDepth 32768 in
set_option maxHeartbeats 4000000 in
theorem orbit_table :
    permutations.map (fun g => (histogram g false, histogram g true)) =
      [([8, 0, 0, 0, 0], [0, 8, 0, 0, 0]),
       ([4, 0, 0, 4, 0], [0, 4, 4, 0, 0]),
       ([4, 0, 0, 4, 0], [0, 4, 4, 0, 0]),
       ([4, 0, 0, 4, 0], [0, 4, 4, 0, 0]),
       ([2, 0, 0, 0, 6], [0, 8, 0, 0, 0]),
       ([2, 0, 0, 0, 6], [0, 8, 0, 0, 0])] := by decide +kernel

/-- At 37 the cubic splits and both -11 and -263 are squares. -/
theorem common_split_witness :
    ((List.range 37).filter fun x => (x^3 + 1) % 37 == 4*x % 37) = [8, 13, 16] ∧
    (10^2 + 11) % 37 = 0 ∧ (12^2 + 263) % 37 = 0 ∧
    11 % 8 = 3 ∧ 263 % 16 = 7 ∧ 264 % 11 = 0 := by decide

/-- The fixed cubic has a single root at 7 and three roots at 347. -/
theorem local_examples :
    ((List.range 7).filter fun x => (x^3 + 1) % 7 == 4*x % 7) = [4] ∧
    ((List.range 347).filter fun x => (x^3 + 1) % 347 == 4*x % 347) = [9, 92, 246] := by
  decide +kernel

end CM229Frobenius
