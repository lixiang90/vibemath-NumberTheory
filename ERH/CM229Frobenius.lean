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

/-- A mixed CM type has one plus sign, at the indicated real embedding. -/
def inMixedType (minority : Nat) (e : Embedding) : Bool :=
  (e.1 == minority) == e.2

/-- Edges leaving the CM type in the circular-word description of A[p]. -/
def aNumber (perm : List Nat) (flip : Bool) (minority : Nat) : Nat :=
  (embeddings.filter fun e =>
    inMixedType minority e && !(inMixedType minority (frobenius perm flip e))).length

theorem mixed_a_numbers :
    permutations.map (fun g =>
      ((List.range 3).map (aNumber g false), (List.range 3).map (aNumber g true))) =
      [([0,0,0], [3,3,3]), ([2,2,0], [1,1,3]), ([2,0,2], [1,3,1]),
       ([0,2,2], [3,1,1]), ([2,2,2], [1,1,1]), ([2,2,2], [1,1,1])] := by
  decide +kernel

/-- In a CM elliptic module only the edge starting at false survives modulo p. -/
def ellipticNext (flip e : Bool) : Bool := if flip then !e else e

def ellipticANumber (flip : Bool) : Nat :=
  ([false, true].filter fun e => !e && ellipticNext flip e).length

def ellipticFrobeniusSquareRank (flip : Bool) : Nat :=
  ([false, true].filter fun e => !e && !(ellipticNext flip e)).length

/-- Each top exterior power has dimension one precisely when its a-number is three. -/
def localTensorLength (g : List Nat) (qInert dInert iInert : Bool) : Nat :=
  ((List.range 3).map fun j => if aNumber g qInert j == 3 then 1 else 0).prod *
    ellipticANumber dInert * ellipticFrobeniusSquareRank iInert

theorem local_tensor_indicator (j : Fin 6) (q d i : Bool) :
    localTensorLength (permutations[j.val]?.getD []) q d i =
      (if j == 0 && q && d && !i then 1 else 0) := by
  revert j q d i
  decide +kernel

theorem two_mixed_types_suffice :
    ∀ (j : Fin 6) (q : Bool) (a b : Fin 3), a ≠ b →
      ((aNumber (permutations[j.val]?.getD []) q a.val == 3) &&
       (aNumber (permutations[j.val]?.getD []) q b.val == 3)) = (j == 0 && q) := by
  decide +kernel

/-- Parity of the permutation on the three real embeddings, via inversions. -/
def oddPermutation (g : List Nat) : Bool :=
  let values := fun i => g[i]?.getD i
  let inversions := [(0,1), (0,2), (1,2)].filter fun (i,j) => values i > values j
  inversions.length % 2 == 1

/-- A single mixed reduction, together with the two quadratic characters. -/
def mixedRelativeSplit (g : List Nat) (iInert : Bool) (minority : Nat) : Bool :=
  oddPermutation g || (aNumber g iInert minority == (if iInert then 3 else 0))

/-- Relative Frobenius is the base residue-degree power of rational Frobenius. -/
def relativeFrobTrivial (g : List Nat) (iInert : Bool) : Bool :=
  let step := frobenius g iInert
  embeddings.all fun e =>
    (if oddPermutation g || iInert then step (step e) else step e) == e

theorem single_mixed_relative_split (j : Fin 6) (iInert : Bool) (minority : Fin 3) :
    mixedRelativeSplit (permutations[j.val]?.getD []) iInert minority.val =
      relativeFrobTrivial (permutations[j.val]?.getD []) iInert := by
  revert j iInert minority
  decide +kernel

/-- Trace of the standard representation on the three roots. -/
def standardTrace (g : List Nat) : Int :=
  Int.ofNat ((List.range 3).filter fun i => g[i]?.getD i == i).length - 1

def squarePermutation (g : List Nat) : List Nat :=
  (List.range 3).map fun i => g[g[i]?.getD i]?.getD i

/-- The first prime coefficient of the induced cubic character. -/
theorem first_carry_trace (j : Fin 6) (iInert : Bool) :
    let g := permutations[j.val]?.getD []
    (if iInert then 0 else 2) * standardTrace g =
      2 * (if oddPermutation g || iInert then 0
           else if mixedRelativeSplit g iInert 0 then 2 else -1) := by
  revert j iInert
  decide +kernel

/-- The trivial constituent at squares accounts for the square-root drift. -/
theorem square_carry_trace (j : Fin 6) :
    let g := permutations[j.val]?.getD []
    standardTrace (squarePermutation g) =
      1 - (if oddPermutation g then (-1 : Int) else 1) + standardTrace g := by
  revert j
  decide +kernel

def cubicRoots (p : Nat) : List Nat :=
  (List.range p).filter fun x => (x^3 + 1) % p == (4*x) % p

def minusNonSquare (a p : Nat) : Bool :=
  !((List.range p).any fun x => (x*x + a) % p == 0)

theorem mixed_reduction_examples :
    ([193,277,29].all fun p => p % 4 == 1 && minusNonSquare 11 p && minusNonSquare 263 p) = true ∧
    cubicRoots 193 = [6,36,151] ∧ cubicRoots 277 = [] ∧ cubicRoots 29 = [13] := by
  decide +kernel

end CM229Frobenius
