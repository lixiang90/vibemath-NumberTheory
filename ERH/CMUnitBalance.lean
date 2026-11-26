import Mathlib.Algebra.Order.Round
import Mathlib.Algebra.Order.BigOperators.Group.Finset
import Mathlib.Basic.Real.Basic
import Mathlib.Algebra.Order.Archimedean.Real.Basic
import Mathlib.Tactic

/-! Rounding in a real logarithmic unit lattice and its effect on weighted norms. -/

namespace CMUnitBalance

open scoped BigOperators

/-- Columns are logarithmic unit vectors. One choice of integer powers works
for every infinite place, provided the centered vector lies in their real span. -/
theorem balance {ι κ : Type*} [Fintype κ]
    (L : ι → κ → ℝ) (x : ι → ℝ) (mean : ℝ)
    (span : ∃ t : κ → ℝ, ∀ i, x i - mean = ∑ j, t j * L i j) :
    ∃ n : κ → ℤ, ∀ i,
      |x i - mean - ∑ j, (n j : ℝ) * L i j| ≤ (∑ j, |L i j|) / 2 := by
  classical
  obtain ⟨t, ht⟩ := span
  refine ⟨fun j => round (t j), ?_⟩
  intro i
  rw [ht i, ← Finset.sum_sub_distrib]
  simp_rw [← sub_mul]
  calc
    _ ≤ ∑ j, |(t j - (round (t j) : ℝ)) * L i j| := Finset.abs_sum_le_sum_abs _ _
    _ ≤ ∑ j, (1 / 2 : ℝ) * |L i j| := by
      apply Finset.sum_le_sum
      intro j _
      rw [abs_mul]
      exact mul_le_mul_of_nonneg_right (abs_sub_round (t j)) (abs_nonneg _)
    _ = _ := by rw [← Finset.mul_sum]; ring

/-- Only differences of archimedean weights contribute to a centered error. -/
theorem weighted_error {ι : Type*} [Fintype ι] (r w B : ι → ℝ) (w₀ : ℝ)
    (centered : ∑ i, r i = 0) (bounded : ∀ i, |r i| ≤ B i) :
    |∑ i, w i * r i| ≤ ∑ i, |w i - w₀| * B i := by
  classical
  have h : ∑ i, (w i - w₀) * r i = ∑ i, w i * r i := by
    simp only [sub_mul, Finset.sum_sub_distrib, ← Finset.mul_sum, centered, mul_zero, sub_zero]
  rw [← h]
  calc
    _ ≤ ∑ i, |(w i - w₀) * r i| := Finset.abs_sum_le_sum_abs _ _
    _ ≤ _ := by
      apply Finset.sum_le_sum
      intro i _
      rw [abs_mul]
      exact mul_le_mul_of_nonneg_left (bounded i) (abs_nonneg _)

end CMUnitBalance
