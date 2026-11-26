import Std

/-!
The natural action of A₄ on the four points of the affine plane over F₂,
its action on the six edges, and its cyclic quotient on the three directions.
The affine plane here is the two-dimensional vector space over F₂.
-/

namespace Tetra3136

abbrev Point := Bool × Bool

def vertices : List Point :=
  [(false, false), (true, false), (false, true), (true, true)]

def directions : List Point := vertices.drop 1

def edges : List (Point × Point) :=
  vertices.flatMap fun x => (vertices.drop ((vertices.idxOf x) + 1)).map (x, ·)

def rotate (j : Fin 3) (v : Point) : Point :=
  match j.val with
  | 0 => v
  | 1 => (v.2, v.1.xor v.2)
  | _ => (v.1.xor v.2, v.1)

def act (a b : Bool) (j : Fin 3) (v : Point) : Point :=
  let r := rotate j v
  (a.xor r.1, b.xor r.2)

def trace (g : Point → Point) : Int :=
  ((vertices.filter fun v => g v == v).length : Int) - 1

def fixedDirections (j : Fin 3) : Int :=
  (directions.filter fun d => rotate j d == d).length

def fixedEdges (g : Point → Point) : Int :=
  (edges.filter fun (x, y) =>
    (g x == x && g y == y) || (g x == y && g y == x)).length

def identityIndicator (a b : Bool) (j : Fin 3) : Int :=
  if !a && !b && j == 0 then 1 else 0

/-- Prime and prime-square traces in the degree-three constituent. -/
theorem prime_traces : ∀ (a b : Bool) (j : Fin 3),
    let g := act a b j
    trace g = 4 * identityIndicator a b j - (if j = 0 then 1 else 0) ∧
    trace (g ∘ g) = fixedDirections j ∧
    fixedDirections j = (if j = 0 then 3 else 0) := by
  decide

/-- Characters for the index-six field and the full twelve-point regular action. -/
theorem subfield_characters : ∀ (a b : Bool) (j : Fin 3),
    let g := act a b j
    fixedEdges g = fixedDirections j + trace g ∧
    12 * identityIndicator a b j = fixedDirections j + 3 * trace g := by
  decide

end Tetra3136
