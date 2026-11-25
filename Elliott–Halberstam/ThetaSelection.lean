import Mathlib.Algebra.Order.BigOperators.Group.Finset
import Mathlib.Algebra.BigOperators.Field
import Mathlib.Data.Complex.Basic
import Mathlib.Tactic.Linarith

/-! The finite selection and centering steps in the translated cyclic test. -/
namespace EHThetaSelection
open scoped BigOperators

private theorem supported_le_mean {ι : Type*} [Fintype ι] (p cost : ι → ℝ)
    (hp : ∀ i, 0 ≤ p i) (hmass : ∑ i, p i = 1) {v : ℝ}
    (hmean : ∑ i, p i * cost i ≤ v) : ∃ i, 0 < p i ∧ cost i ≤ v := by
  classical
  by_contra h
  have hbad (i : ι) (hi : 0 < p i) : v < cost i :=
    lt_of_not_ge (fun hle => h ⟨i, hi, hle⟩)
  have hex : ∃ j, 0 < p j := by
    by_contra hn
    have hz (i : ι) : p i = 0 :=
      le_antisymm (le_of_not_gt (fun hi => hn ⟨i, hi⟩)) (hp i)
    simp [hz] at hmass
  obtain ⟨j, hj⟩ := hex
  have hlt : ∑ i, p i * v < ∑ i, p i * cost i := by
    refine Finset.sum_lt_sum (s := Finset.univ)
      (f := fun i : ι => p i * v) (g := fun i : ι => p i * cost i) ?_ ?_
    · intro i _
      by_cases hi : 0 < p i
      · exact mul_le_mul_of_nonneg_left (hbad i hi).le (hp i)
      · simp [le_antisymm (le_of_not_gt hi) (hp i)]
    · exact ⟨j, Finset.mem_univ _, mul_lt_mul_of_pos_left (hbad j hj) hj⟩
  rw [← Finset.sum_mul, hmass, one_mul] at hlt
  exact (not_lt.mpr hmean) hlt

/-- Select one outcome of positive probability for a whole finite family of budgets.
Zero-frequency count errors, marked energies and orbit loss may be separate entries. -/
theorem select_costs {ι κ : Type*} [Fintype ι] [Fintype κ]
    (p : ι → ℝ) (cost : κ → ι → ℝ) (budget : κ → ℝ)
    (hp : ∀ i, 0 ≤ p i) (hmass : ∑ i, p i = 1)
    (hn : ∀ j i, 0 ≤ cost j i) (hb : ∀ j, 0 < budget j)
    (he : ∑ j, (∑ i, p i * cost j i) / budget j ≤ 1) :
    ∃ i, 0 < p i ∧ ∀ j, cost j i ≤ budget j := by
  classical
  have hm : ∑ i, p i * (∑ j, cost j i / budget j) ≤ 1 := by
    simp only [Finset.mul_sum]
    rw [Finset.sum_comm]
    simpa only [Finset.sum_div, mul_div_assoc] using he
  obtain ⟨i, hi, hc⟩ := supported_le_mean p _ hp hmass hm
  refine ⟨i, hi, fun j => ?_⟩
  have h := Finset.single_le_sum (f := fun l => cost l i / budget l)
    (fun l _ => div_nonneg (hn l i) (hb l).le) (Finset.mem_univ j)
  exact (div_le_one (hb j)).mp (h.trans hc)

/-- The empirical total may fluctuate. Centering the test first removes its contribution. -/
theorem centered_mark {κ : Type*} [Fintype κ] (counts means test : κ → ℂ) (w : ℂ)
    (hcenter : ∑ j, means j * (test j - w) = 0) :
    ∑ j, (counts j - means j) * (test j - w) =
      (∑ j, counts j * test j) - (∑ j, counts j) * w := by
  simp_rw [sub_mul, Finset.sum_sub_distrib, hcenter, sub_zero, mul_sub,
    Finset.sum_sub_distrib, Finset.sum_mul]

end EHThetaSelection

#print axioms EHThetaSelection.select_costs
#print axioms EHThetaSelection.centered_mark
