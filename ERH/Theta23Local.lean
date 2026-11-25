import Mathlib.Algebra.Polynomial.AlgebraMap
import Mathlib.Algebra.Polynomial.Roots
import Mathlib.Algebra.Polynomial.Taylor
import Mathlib.Data.ZMod.Basic
import Mathlib.Tactic

/-! Algebraic steps in the two-weight theta content argument for discriminant -23.
The statements about harmonic theta coefficients and prime distribution are in the paper. -/

namespace Theta23

/-- The second column in units of the first weight, at a principal prime. -/
def triple {R : Type*} [CommRing R] (s p : R) : R :=
  s ^ 3 - 3 * p * s ^ 2 + 3 * p ^ 3

theorem triple_identity {R : Type*} [CommRing R] (x y : R) :
    x ^ 6 + (x * y) ^ 3 + y ^ 6 = triple (x ^ 2 + x * y + y ^ 2) (x * y) := by
  unfold triple
  ring

theorem principal_difference {R : Type*} [CommRing R] (s p : R) :
    s ^ 3 - triple s p = 3 * p * (s ^ 2 - p ^ 2) := by
  unfold triple
  ring

/-- Away from residue characteristic 3 the principal row has no common zero for p a unit. -/
theorem principal_row_nonzero {F : Type*} [Field F] (s p : F)
    (hp : p ≠ 0) (h3 : (3 : F) ≠ 0) : s ^ 3 ≠ 0 ∨ triple s p ≠ 0 := by
  by_cases hs : s = 0
  · right
    simpa [triple, hs] using mul_ne_zero h3 (pow_ne_zero 3 hp)
  · exact Or.inl (pow_ne_zero 3 hs)

/-- At 3 the principal row has valuation exactly one when s is divisible by 3 and p is a unit. -/
theorem principal_mod_nine : ∀ t p : ZMod 9, IsUnit p →
    triple (3 * t) p ≠ 0 ∧ triple (3 * t) p = 3 * p ^ 3 := by
  decide

theorem inert_difference {R : Type*} [CommRing R] (p : R) :
    (2 - p) ^ 3 - (2 - p ^ 3) = 6 * (p - 1) ^ 2 := by ring

/-- A prime common to both inert entries has residue characteristic 2 or 3. -/
theorem inert_common_zero {F : Type*} [Field F] (p : F)
    (h : (2 - p) ^ 3 = 0) (h' : 2 - p ^ 3 = 0) : (2 : F) = 0 ∨ (3 : F) = 0 := by
  have h0 : 2 - p = 0 := by simpa using h
  have hp : p = 2 := (sub_eq_zero.mp h0).symm
  rw [hp] at h'
  apply mul_eq_zero.mp
  linear_combination -h'

/-- Reading a common coefficient and putting it back on a normalized line. -/
def project {F V : Type*} [Field F] [AddCommGroup V] [Module F V]
    (read : V →ₗ[F] F) (w x : V) : V := read x • w

theorem projected_step {F V : Type*} [Field F] [AddCommGroup V] [Module F V]
    (read : V →ₗ[F] F) (w : V) (hw : read w = 1) (A : V →ₗ[F] V) (r : F) (x : V) :
    project read w (A (project read w x) + r • project read w x) =
      (read (A w) + r) • project read w x := by
  simp [project, hw, smul_smul, mul_add, mul_comm]

open Polynomial
open scoped Polynomial

/-- The first possible Taylor coefficient after factoring out h copies of X-a;
the coefficient is nonzero precisely when q(a) is. -/
theorem shifted_factor_jets {R : Type*} [CommRing R] (q : R[X]) (a : R) (h : ℕ) :
    (∀ j < h, (hasseDeriv j ((X - C a) ^ h * q)).eval a = 0) ∧
    (hasseDeriv h ((X - C a) ^ h * q)).eval a = q.eval a := by
  constructor
  · intro j hj
    rw [← taylor_coeff, taylor_mul, taylor_pow, taylor_X_sub_C, coeff_X_pow_mul']
    simp [Nat.not_le_of_lt hj]
  · rw [← taylor_coeff, taylor_mul, taylor_pow, taylor_X_sub_C, coeff_X_pow_mul']
    simp [taylor_coeff_zero]

/-- A finite set of Hasse jets detects every nonzero polynomial of the stated degree,
also in positive characteristic. The roots in the finite set are pairwise distinct. -/
theorem exists_nonzero_jet {F : Type*} [Field F] (f : F[X]) (hf : f ≠ 0)
    (points : Finset F) (m : ℕ) (hd : f.natDegree < points.card * m) :
    ∃ a ∈ points, ∃ j < m, (hasseDeriv j f).eval a ≠ 0 := by
  classical
  by_contra! hall
  have hm (a : F) (ha : a ∈ points) : m ≤ f.roots.count a := by
    rw [count_roots, le_rootMultiplicity_iff hf, X_sub_C_pow_dvd_iff, X_pow_dvd_iff]
    intro j hj
    rw [← taylor_apply, taylor_coeff]
    exact hall a ha j hj
  have hsum : points.card * m ≤ ∑ a ∈ points, f.roots.count a := by
    simpa using Finset.sum_le_sum hm
  have hext : ∑ a ∈ points, f.roots.count a ≤
      ∑ a ∈ points ∪ f.roots.toFinset, f.roots.count a :=
    Finset.sum_le_sum_of_subset_of_nonneg Finset.subset_union_left (by intros; omega)
  have htotal : (∑ a ∈ points ∪ f.roots.toFinset, f.roots.count a) = f.roots.card :=
    Multiset.sum_count_eq_card (by intro a ha; exact Finset.mem_union_right _ (by simpa using ha))
  exact (Nat.not_lt_of_ge (hsum.trans (hext.trans (htotal.le.trans f.card_roots')))) hd

end Theta23
