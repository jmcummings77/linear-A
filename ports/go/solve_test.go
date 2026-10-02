package matrix

import (
	"math"
	"testing"
)

func TestReusableFactorSnapshots(t *testing.T) {
	for _, algorithm := range []string{"lu", "cholesky", "qr"} {
		a, _ := New(2, 2, []float64{4, 1, 1, 3})
		b, _ := New(2, 2, []float64{6, 5, 7, 4})
		factor, err := factorize(a, algorithm)
		if err != nil {
			t.Fatal(err)
		}
		a.values[0] = 99
		for run := 0; run < 3; run++ {
			x, err := factor.Solve(b)
			if err != nil {
				t.Fatal(err)
			}
			for i, want := range []float64{1, 1, 2, 1} {
				if math.Abs(x.values[i]-want) > 1e-12 {
					t.Fatal(x.values)
				}
			}
		}
		condition, err := factor.ReciprocalCondition()
		if err != nil || math.Abs(condition-.44) > 1e-12 {
			t.Fatal(condition, err)
		}
		if b.values[0] != 6 {
			t.Fatal("RHS changed")
		}
		bad, _ := Zeros(1, 1)
		if _, err := factor.Solve(bad); err == nil {
			t.Fatal("wrong RHS accepted")
		}
	}
}
