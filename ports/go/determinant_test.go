package matrix

import (
	"math"
	"math/rand"
	"testing"
)

func permutationDeterminant(values []float64, size int) int64 {
	permutation, used := make([]int, size), make([]bool, size)
	var result int64
	var visit func(int)
	visit = func(row int) {
		if row == size {
			product, inversions := int64(1), 0
			for i := 0; i < size; i++ {
				product *= int64(values[i*size+permutation[i]])
				for j := i + 1; j < size; j++ {
					if permutation[i] > permutation[j] {
						inversions++
					}
				}
			}
			if inversions%2 == 0 {
				result += product
			} else {
				result -= product
			}
			return
		}
		for col := 0; col < size; col++ {
			if !used[col] {
				used[col] = true
				permutation[row] = col
				visit(row + 1)
				used[col] = false
			}
		}
	}
	visit(0)
	return result
}

func TestDeterminantAlgorithmsMatchIndependentPermutationOracle(t *testing.T) {
	rng := rand.New(rand.NewSource(781))
	for size := 0; size <= 6; size++ {
		for sample := 0; sample < 6; sample++ {
			values := make([]float64, size*size)
			for i := range values {
				values[i] = float64(rng.Intn(7) - 3)
			}
			m := mustMatrix(t, size, size, values)
			expected := float64(permutationDeterminant(values, size))
			for _, algorithm := range []DeterminantAlgorithm{DeterminantAuto, DeterminantLU} {
				got, err := m.DeterminantWith(algorithm)
				if err != nil || math.Abs(got-expected) > 1e-10*math.Max(1, math.Abs(expected)) {
					t.Fatalf("size%d/%d/%d: got%g error%v want%g", size, sample, algorithm, got, err, expected)
				}
			}
			assertValues(t, m, values)
			for row := 0; row < size; row++ {
				for col := row; col < size; col++ {
					value := float64(rng.Intn(7) - 3)
					if row == col {
						value += float64(4 * size)
					}
					values[row*size+col], values[col*size+row] = value, value
				}
			}
			m = mustMatrix(t, size, size, values)
			expected = float64(permutationDeterminant(values, size))
			got, err := m.DeterminantWith(DeterminantCholesky)
			if err != nil || math.Abs(got-expected) > 1e-10*math.Max(1, math.Abs(expected)) {
				t.Fatalf("SPD size%d/%d: got%g error%v want%g", size, sample, got, err, expected)
			}
			assertValues(t, m, values)
		}
	}
}

func TestDeterminantExtremeScalesAndUnderflowingFactors(t *testing.T) {
	cases := []struct {
		size   int
		values []float64
		want   float64
	}{
		{4, []float64{1e300, 0, 0, 0, 0, 1e300, 0, 0, 0, 0, 1e-300, 0, 0, 0, 0, 1e-300}, 1},
		{2, []float64{1e308, 1e308, 1e-308, 2e-308}, 1},
		{2, []float64{1e-200, 1e-200, 3e-124, 6e-124}, math.SmallestNonzeroFloat64},
		{3, []float64{1e308, 1e308, 1e-308, 1e-308, 2e-308, 1e308, 0, 0, 1e-308}, 1e-308},
		{3, []float64{math.MaxFloat64, 0, 0, 0, math.MaxFloat64, 0, 0, 0, 0}, 0},
	}
	for _, tc := range cases {
		m := mustMatrix(t, tc.size, tc.size, tc.values)
		for _, algorithm := range []DeterminantAlgorithm{DeterminantAuto, DeterminantLU} {
			got, err := m.DeterminantWith(algorithm)
			if err != nil || tc.want == 0 && got != 0 || tc.want != 0 && math.Abs(got/tc.want-1) > 1e-12 {
				t.Fatalf("%v algorithm%d: got%g error%v want%g", tc.values, algorithm, got, err, tc.want)
			}
		}
	}
	for _, values := range [][]float64{{4, 2, 2, 3}, {1e300, 0, 0, 1e-300}, {math.SmallestNonzeroFloat64}} {
		size := 2
		if len(values) == 1 {
			size = 1
		}
		m := mustMatrix(t, size, size, values)
		want, _ := m.Determinant()
		got, err := m.DeterminantWith(DeterminantCholesky)
		if err != nil || math.Abs(got/want-1) > 1e-12 {
			t.Fatalf("Cholesky %v got%g error%v want%g", values, got, err, want)
		}
	}
}

func TestCholeskyRejectsInvalidInputAndUnknownAlgorithms(t *testing.T) {
	for _, values := range [][]float64{{4, 2, 1, 3}, {1, 2, 2, 1}, {1, 1, 1, 1}, {-1, 0, 0, -1}} {
		if _, err := mustMatrix(t, 2, 2, values).DeterminantWith(DeterminantCholesky); err == nil {
			t.Fatalf("accepted invalid Cholesky input%v", values)
		}
	}
	m := mustMatrix(t, 1, 1, []float64{1})
	if _, err := m.DeterminantWith(DeterminantAlgorithm(99)); err == nil {
		t.Fatal("accepted invalid algorithm")
	}
	m.values[0] = math.NaN()
	for _, algorithm := range []DeterminantAlgorithm{DeterminantAuto, DeterminantLU, DeterminantCholesky} {
		if _, err := m.DeterminantWith(algorithm); err == nil {
			t.Fatal("accepted nonfinite input")
		}
	}
}
