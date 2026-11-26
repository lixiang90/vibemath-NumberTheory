import Std

/-! Cofactor support for the degree-two factor of the cubic field of discriminant 229.
The input is a list of prime powers with their Euler types. -/

namespace CM229CofactorSupport

inductive EulerType
  | split | transposition | cycle | ramified
  deriving DecidableEq

def coefficient : EulerType → Nat → Int
  | .split, e => (e : Int) + 1
  | .transposition, e => if e % 2 = 0 then 1 else 0
  | .cycle, e => if e % 3 = 0 then 1 else if e % 3 = 1 then -1 else 0
  | .ramified, _ => 1

def inverse : EulerType → Nat → Int
  | _, 0 => 1
  | .split, 1 => -2
  | .transposition, 1 => 0
  | .cycle, 1 => 1
  | .ramified, 1 => -1
  | .split, 2 | .cycle, 2 => 1
  | .transposition, 2 => -1
  | _, _ => 0

/-- A local term in u(d) a(h/d) can survive exactly at the admissible exponents. -/
theorem local_support (kind : EulerType) (e : Nat) :
    (∃ j, j ≤ e ∧ inverse kind j ≠ 0 ∧ coefficient kind (e-j) ≠ 0) ↔
    (kind = .transposition → e % 2 = 0) := by
  cases kind with
  | split =>
    constructor
    · intro _ contra
      cases contra
    · intro _
      exact ⟨0, Nat.zero_le _, by decide, by simp [coefficient]; omega⟩
  | ramified =>
    constructor
    · intro _ contra
      cases contra
    · intro _
      exact ⟨0, Nat.zero_le _, by decide, by simp [coefficient]⟩
  | cycle =>
    constructor
    · intro _ contra
      cases contra
    · intro _
      by_cases hm : e % 3 = 2
      · have hrem : (e - 1) % 3 = 1 := by omega
        exact ⟨1, by omega, by decide, by simp [coefficient, hrem]⟩
      · have hrem : e % 3 = 0 ∨ e % 3 = 1 := by omega
        rcases hrem with h | h <;> exact ⟨0, Nat.zero_le _, by decide, by simp [coefficient, h]⟩
  | transposition =>
    constructor
    · rintro ⟨j, hj, hu, ha⟩ _
      have evenRest : (e - j) % 2 = 0 := by
        by_cases hp : (e - j) % 2 = 0
        · exact hp
        · simp [coefficient, hp] at ha
      have hcases : j = 0 ∨ j = 2 := by
        cases j with
        | zero => simp
        | succ j =>
          cases j with
          | zero => contradiction
          | succ j =>
            cases j with
            | zero => simp
            | succ j => simp [inverse] at hu
      omega
    · intro h
      exact ⟨0, Nat.zero_le _, by decide, by simp [coefficient, h rfl]⟩

/-- One entry for each divisor of the integer described by the prime powers. -/
def divisorTerms : List (EulerType × Nat) → List Int
  | [] => [1]
  | (kind, e) :: rest =>
    (List.range (e+1)).flatMap fun j =>
      (divisorTerms rest).map fun v => inverse kind j * coefficient kind (e-j) * v

theorem support_on_factorization (factors : List (EulerType × Nat)) :
    (∃ v ∈ divisorTerms factors, v ≠ 0) ↔
    ∀ f ∈ factors, f.1 = .transposition → f.2 % 2 = 0 := by
  induction factors with
  | nil => simp [divisorTerms]
  | cons f fs ih =>
    rcases f with ⟨kind, e⟩
    simp only [List.forall_mem_cons]
    constructor
    · rintro ⟨v, hv, hn⟩
      simp only [divisorTerms, List.mem_flatMap, List.mem_range, List.mem_map] at hv
      rcases hv with ⟨j, hj, w, hw, rfl⟩
      refine ⟨(local_support kind e).mp ⟨j, by omega, ?_, ?_⟩, ih.mp ⟨w, hw, ?_⟩⟩
      all_goals intro hz; simp [hz] at hn
    · rintro ⟨he, hfs⟩
      obtain ⟨j, hj, hu, ha⟩ := (local_support kind e).mpr he
      obtain ⟨w, hw, hn⟩ := ih.mpr hfs
      refine ⟨_, ?_, Int.mul_ne_zero (Int.mul_ne_zero hu ha) hn⟩
      simp only [divisorTerms, List.mem_flatMap, List.mem_range, List.mem_map]
      exact ⟨j, by omega, w, hw, rfl⟩

/-- The cofactor for h > 1 is the negative prefix sum of its divisor terms. -/
def cutoff (term : Nat → Int) : Nat → Int
  | 0 => 0
  | z+1 => cutoff term z - term (z+1)

def divisorTerm (a u : Nat → Int) (h d : Nat) : Int :=
  if d ∣ h then u d * a (h/d) else 0

def residual (a u : Nat → Int) (h z : Nat) : Int :=
  (if h = 1 then 1 else 0) + cutoff (divisorTerm a u h) z

private theorem cutoff_scale (f g : Nat → Int) (b : Int) (z : Nat)
    (h : ∀ d, 0 < d → d ≤ z → f d = g d * b) :
    cutoff f z = cutoff g z * b := by
  induction z with
  | zero => simp [cutoff]
  | succ z ih =>
    have prev := ih (fun d hd hle => h d hd (by omega))
    simp only [cutoff, prev, h (z+1) (by omega) (Nat.le_refl _), Int.sub_mul]

/-- Pull the rough part t through the actual divisor cutoff. -/
theorem smooth_part_formula (a u : Nat → Int) (s t z : Nat)
    (product_ne_one : s*t ≠ 1) (cop : s.Coprime t)
    (rough : ∀ d, 0 < d → d ≤ z → d.Coprime t)
    (mul : ∀ m n, m.Coprime n → a (m*n) = a m * a n) :
    residual a u (s*t) z =
      (residual a u s z - (if s = 1 then 1 else 0)) * a t := by
  have scaled : cutoff (divisorTerm a u (s*t)) z = cutoff (divisorTerm a u s) z * a t := by
    apply cutoff_scale
    intro d hd hle
    by_cases hds : d ∣ s
    · have hdst := Nat.dvd_mul_right_of_dvd hds t
      have hdiv : s*t/d = (s/d)*t := by
        rw [Nat.mul_comm s t, Nat.mul_div_assoc t hds, Nat.mul_comm t (s/d)]
      have hcop := Nat.Coprime.coprime_dvd_left (Nat.div_dvd_of_dvd hds) cop
      simp [divisorTerm, hds, hdst, hdiv, mul (s/d) t hcop, Int.mul_assoc]
    · have hn : ¬d ∣ s*t := fun h => hds ((rough d hd hle).dvd_of_dvd_mul_right h)
      simp [divisorTerm, hds, hn]
  unfold residual
  rw [if_neg product_ne_one, Int.zero_add, scaled]
  congr 1
  omega

theorem some_cutoff_iff (term : Nat → Int) (Z : Nat) :
    (∃ z, z ≤ Z ∧ cutoff term z ≠ 0) ↔
    (∃ d, 0 < d ∧ d ≤ Z ∧ term d ≠ 0) := by
  induction Z with
  | zero =>
    simp [cutoff]
    omega
  | succ Z ih =>
    constructor
    · rintro ⟨z, hz, hc⟩
      by_cases hz' : z ≤ Z
      · obtain ⟨d, hd, hle, hn⟩ := ih.mp ⟨z, hz', hc⟩
        exact ⟨d, hd, by omega, hn⟩
      · have heq : z = Z+1 := by omega
        subst z
        by_cases ht : term (Z+1) = 0
        · have hprev : cutoff term Z ≠ 0 := by simpa [cutoff, ht] using hc
          obtain ⟨d, hd, hle, hn⟩ := ih.mp ⟨Z, Nat.le_refl _, hprev⟩
          exact ⟨d, hd, by omega, hn⟩
        · exact ⟨Z+1, by omega, Nat.le_refl _, ht⟩
    · rintro ⟨d, hd, hle, hn⟩
      by_cases hd' : d ≤ Z
      · obtain ⟨z, hz, hc⟩ := ih.mpr ⟨d, hd, hd', hn⟩
        exact ⟨z, by omega, hc⟩
      · have heq : d = Z+1 := by omega
        subst d
        by_cases hprev : cutoff term Z = 0
        · exact ⟨Z+1, Nat.le_refl _, by simpa [cutoff, hprev] using hn⟩
        · exact ⟨Z, by omega, hprev⟩

end CM229CofactorSupport
