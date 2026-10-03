package main

import (
	"encoding/json"
	"errors"
	"fmt"
	"io"
	"math"
	"os"
	"runtime/pprof"
	"strconv"
	"strings"
	"time"

	matrix "github.com/jmcummings77/linear-A/ports/go"
)

type operation string

func parseOperation(value string) (operation, error) {
	switch value {
	case "pseudoinverse", "spectral_diagnostics", "solve_minimum_norm", "svd", "svd_one_sweep", "add", "subtract", "scale", "transpose", "multiply", "solve", "solve_cholesky", "least_squares", "rcond", "cross", "rotation2d", "rotation3d", "trace", "determinant", "determinant_lu", "determinant_spd_lu", "determinant_cholesky", "eigen_symmetric", "eigen_general", "triangular":
		return operation(value), nil
	default:
		return "", fmt.Errorf("unknown operation: %s", value)
	}
}

func (op operation) binary() bool {
	return op == "add" || op == "subtract" || op == "multiply" || op == "cross" || op == "solve" || op == "solve_cholesky" || op == "least_squares" || op == "solve_minimum_norm"
}

type outcome struct {
	matrix         *matrix.Matrix
	eigen          *matrix.SymmetricEigenDecomposition
	general        *matrix.GeneralEigenDecomposition
	value          float64
	classification bool
	upper, lower   bool
}

func execute(op operation, a, b *matrix.Matrix, scalar float64) (outcome, error) {
	var result outcome
	var err error
	switch op {
	case "add":
		result.matrix, err = a.Add(b)
	case "subtract":
		result.matrix, err = a.Subtract(b)
	case "scale":
		result.matrix, err = a.Scale(scalar)
	case "transpose":
		result.matrix = a.Transpose()
	case "multiply":
		result.matrix, err = a.Multiply(b)
	case "pseudoinverse":
		result.matrix, err = a.PseudoinverseWith(scalar)
	case "solve_minimum_norm":
		result.matrix, err = a.SolveMinimumNorm(b)
	case "spectral_diagnostics":
		var d *matrix.SpectralDiagnostics
		d, err = a.SpectralDiagnosticsWith(scalar)
		if err == nil {
			result.matrix, err = matrix.New(1, 3, []float64{float64(d.Rank), d.ReciprocalCondition, d.RetainedReciprocalCondition})
		}
	case "svd", "svd_one_sweep":
		sweeps := 100
		if op == "svd_one_sweep" {
			sweeps = 1
		}
		var r *matrix.SingularValueDecomposition
		r, err = a.SVDWith(1e-12, sweeps)
		if err == nil {
			values := append(r.U.Values(), r.Values...)
			values = append(values, r.Vt.Transpose().Values()...)
			result.matrix, err = matrix.New(a.Rows()+1+a.Cols(), len(r.Values), values)
		}
	case "solve":
		result.matrix, err = a.Solve(b)
	case "least_squares":
		result.matrix, err = a.LeastSquares(b)
	case "solve_cholesky", "rcond":
		var f *matrix.Factorization
		if op == "rcond" {
			f, err = a.FactorLU()
		} else {
			f, err = a.FactorCholesky()
		}
		if err == nil {
			if op == "rcond" {
				result.value, err = f.ReciprocalCondition()
			} else {
				result.matrix, err = f.Solve(b)
			}
		}
	case "cross":
		result.matrix, err = a.Cross(b)
	case "rotation2d":
		if a.Rows() != 0 || a.Cols() != 0 {
			return result, errors.New("rotation2d requires the empty 0x0 input shape")
		}
		result.matrix, err = matrix.Rotation2D(scalar)
	case "rotation3d":
		result.matrix, err = matrix.RotationAxisAngle(a, scalar)
	case "trace":
		result.value, err = a.Trace()
	case "determinant":
		result.value, err = a.Determinant()
	case "determinant_lu", "determinant_spd_lu":
		result.value, err = a.DeterminantWith(matrix.DeterminantLU)
	case "determinant_cholesky":
		result.value, err = a.DeterminantWith(matrix.DeterminantCholesky)
	case "eigen_general":
		result.general, err = a.EigenGeneral()
	case "eigen_symmetric":
		result.eigen, err = a.EigenSymmetric()
	case "triangular":
		result.classification = true
		result.upper, result.lower = a.Triangular()
	default:
		err = errors.New("unknown operation")
	}
	return result, err
}

func finite(value float64) bool { return !math.IsNaN(value) && !math.IsInf(value, 0) }

func (result outcome) checksum() (float64, error) {
	if result.classification {
		return 0, errors.New("triangular is not a benchmark operation")
	}
	value := result.value
	if result.general != nil {
		g := result.general
		for i, r := range g.EigenvaluesReal {
			value += float64(i+1) * (r + math.Abs(g.EigenvaluesImag[i]))
		}
		for row := 0; row < g.EigenvectorsReal.Rows(); row++ {
			for col := 0; col < g.EigenvectorsReal.Cols(); col++ {
				r, _ := g.EigenvectorsReal.At(row, col)
				im, _ := g.EigenvectorsImag.At(row, col)
				value += r*r + im*im
			}
		}
	}
	if result.eigen != nil {
		for i, eigenvalue := range result.eigen.Eigenvalues {
			value += float64(i+1) * eigenvalue
		}
		q := result.eigen.Eigenvectors
		for row := 0; row < q.Rows(); row++ {
			for col := 0; col < q.Cols(); col++ {
				entry, _ := q.At(row, col)
				value += entry * entry
			}
		}
	}
	if result.matrix != nil {
		m := result.matrix
		length := m.Rows() * m.Cols()
		value = 0
		if length > 0 {
			first, _ := m.At(0, 0)
			middle, _ := m.At((length/2)/m.Cols(), (length/2)%m.Cols())
			last, _ := m.At(m.Rows()-1, m.Cols()-1)
			value = first + middle + last
		}
	}
	if !finite(value) {
		return 0, errors.New("checksum overflowed")
	}
	return value, nil
}

func (result outcome) write() error {
	encoder := json.NewEncoder(os.Stdout)
	if result.general != nil {
		g := result.general
		payload := func(m *matrix.Matrix) any {
			return map[string]any{"rows": m.Rows(), "cols": m.Cols(), "values": m.Values()}
		}
		return encoder.Encode(map[string]any{"eigenvalues_real": g.EigenvaluesReal, "eigenvalues_imag": g.EigenvaluesImag, "eigenvectors_real": payload(g.EigenvectorsReal), "eigenvectors_imag": payload(g.EigenvectorsImag)})
	}
	if result.eigen != nil {
		q := result.eigen.Eigenvectors
		return encoder.Encode(struct {
			Eigenvalues  []float64 `json:"eigenvalues"`
			Eigenvectors struct {
				Rows   int       `json:"rows"`
				Cols   int       `json:"cols"`
				Values []float64 `json:"values"`
			} `json:"eigenvectors"`
		}{result.eigen.Eigenvalues, struct {
			Rows   int       `json:"rows"`
			Cols   int       `json:"cols"`
			Values []float64 `json:"values"`
		}{q.Rows(), q.Cols(), q.Values()}})
	}
	if result.matrix != nil {
		return encoder.Encode(struct {
			Rows   int       `json:"rows"`
			Cols   int       `json:"cols"`
			Values []float64 `json:"values"`
		}{result.matrix.Rows(), result.matrix.Cols(), result.matrix.Values()})
	}
	if result.classification {
		return encoder.Encode(struct {
			Upper bool `json:"upper"`
			Lower bool `json:"lower"`
		}{result.upper, result.lower})
	}
	return encoder.Encode(struct {
		Value float64 `json:"value"`
	}{result.value})
}

func integer(value, name string) (int, error) {
	if value == "" {
		return 0, fmt.Errorf("%s must be a nonnegative integer", name)
	}
	for _, digit := range value {
		if digit < '0' || digit > '9' {
			return 0, fmt.Errorf("%s must be a nonnegative integer", name)
		}
	}
	parsed, err := strconv.ParseInt(value, 10, strconv.IntSize)
	if err != nil {
		return 0, fmt.Errorf("%s is too large", name)
	}
	return int(parsed), nil
}

func number(value string) (float64, error) {
	parsed, err := strconv.ParseFloat(value, 64)
	if err != nil {
		return 0, fmt.Errorf("invalid number: %s", value)
	}
	if !finite(parsed) {
		return 0, errors.New("numbers must be finite")
	}
	return parsed, nil
}

func length(rows, cols int) (int, error) {
	maxElements := int(^uint(0)>>1) / 8
	if cols != 0 && rows > maxElements/cols {
		return 0, errors.New("matrix dimensions overflow storage size")
	}
	return rows * cols, nil
}

func check(args []string) error {
	if len(args) < 4 {
		return errors.New("usage: check OP ROWS COLS [BROWS BCOLS | SCALAR]")
	}
	op, err := parseOperation(args[1])
	if err != nil {
		return err
	}
	expected := 4
	if op.binary() {
		expected = 6
	} else if op == "scale" || op == "rotation2d" || op == "rotation3d" || op == "pseudoinverse" || op == "spectral_diagnostics" {
		expected = 5
	}
	if len(args) != expected {
		return errors.New("wrong number of arguments for check operation")
	}
	rows, err := integer(args[2], "rows")
	if err != nil {
		return err
	}
	cols, err := integer(args[3], "columns")
	if err != nil {
		return err
	}
	aLength, err := length(rows, cols)
	if err != nil {
		return err
	}
	bRows, bCols := 0, 0
	if op.binary() {
		bRows, err = integer(args[4], "right rows")
		if err != nil {
			return err
		}
		bCols, err = integer(args[5], "right columns")
		if err != nil {
			return err
		}
	}
	bLength, err := length(bRows, bCols)
	if err != nil {
		return err
	}
	if aLength > int(^uint(0)>>1)-bLength {
		return errors.New("input count overflow")
	}
	scalar := 1.25
	if op == "scale" || op == "rotation2d" || op == "rotation3d" || op == "pseudoinverse" || op == "spectral_diagnostics" {
		scalar, err = number(args[4])
		if err != nil {
			return err
		}
	}
	input, err := io.ReadAll(os.Stdin)
	if err != nil {
		return err
	}
	fields := strings.Fields(string(input))
	if len(fields) != aLength+bLength {
		return fmt.Errorf("expected %d input values, received %d", aLength+bLength, len(fields))
	}
	values := make([]float64, len(fields))
	for index, field := range fields {
		values[index], err = number(field)
		if err != nil {
			return err
		}
	}
	a, err := matrix.New(rows, cols, values[:aLength])
	if err != nil {
		return err
	}
	var b *matrix.Matrix
	if op.binary() {
		b, err = matrix.New(bRows, bCols, values[aLength:])
		if err != nil {
			return err
		}
	}
	result, err := execute(op, a, b, scalar)
	if err != nil {
		return err
	}
	return result.write()
}

func generated(size, seed int) (*matrix.Matrix, error) {
	result, err := matrix.Zeros(size, size)
	if err != nil {
		return nil, err
	}
	for row := 0; row < size; row++ {
		for col := 0; col < size; col++ {
			index := row*size + col
			numerator := (index%101*17+seed%101*13)%101 - 50
			if err := result.Set(row, col, float64(numerator)/16); err != nil {
				return nil, err
			}
		}
	}
	return result, nil
}

type measurement struct {
	ElapsedNS  int64   `json:"elapsed_ns"`
	Iterations int     `json:"iterations"`
	Checksum   float64 `json:"checksum"`
}

func measure(op operation, a, b *matrix.Matrix, iterations int) (measurement, error) {
	scalar := 1.25
	if op == "rotation2d" || op == "rotation3d" {
		scalar = 0.5
	}
	for i := 0; i < max(5, min(iterations, 100)); i++ {
		result, err := execute(op, a, b, scalar)
		if err != nil {
			return measurement{}, err
		}
		if _, err = result.checksum(); err != nil {
			return measurement{}, err
		}
	}
	started := time.Now()
	checksum := 0.0
	for i := 0; i < iterations; i++ {
		result, err := execute(op, a, b, scalar)
		if err != nil {
			return measurement{}, err
		}
		value, err := result.checksum()
		if err != nil {
			return measurement{}, err
		}
		checksum += value
		if !finite(checksum) {
			return measurement{}, errors.New("accumulated checksum overflowed")
		}
	}
	return measurement{time.Since(started).Nanoseconds(), iterations, checksum}, nil
}

func bench(args []string) error {
	if len(args) != 5 {
		return errors.New("usage: bench OP SIZE ITERATIONS SEED")
	}
	op, err := parseOperation(args[1])
	if err != nil {
		return err
	}
	if op == "triangular" {
		return errors.New("triangular is not a benchmark operation")
	}
	size, err := integer(args[2], "size")
	if err != nil {
		return err
	}
	iterations, err := integer(args[3], "iterations")
	if err != nil {
		return err
	}
	seed, err := integer(args[4], "seed")
	if err != nil {
		return err
	}
	if size == 0 || iterations == 0 {
		return errors.New("size and iterations must be positive")
	}
	if ((op == "cross" || op == "rotation3d") && size != 3) || (op == "rotation2d" && size != 2) {
		return errors.New("cross/rotation3d require size 3; rotation2d requires size 2")
	}
	if seed > 2147483646 {
		return errors.New("seed must be between 0 and 2147483646")
	}
	a, err := generated(size, seed)
	if err != nil {
		return err
	}
	b, err := generated(size, seed+1)
	if err != nil {
		return err
	}
	if op == "cross" {
		av, bv := a.Values(), b.Values()
		a, err = matrix.New(3, 1, av[:3])
		if err != nil {
			return err
		}
		b, err = matrix.New(3, 1, []float64{bv[0], bv[1], -bv[2]})
		if err != nil {
			return err
		}
	} else if op == "rotation2d" {
		a, err = matrix.Zeros(0, 0)
		if err != nil {
			return err
		}
	} else if op == "rotation3d" {
		a, err = matrix.New(3, 1, []float64{1, 2, 3})
		if err != nil {
			return err
		}
	}
	if op == "eigen_general" {
		for i := 0; i < size; i++ {
			for j := 0; j < size; j++ {
				k := i / 2
				a0 := 1 + float64(seed%17)/16 + float64(k)/8
				b0 := .5 + float64(k)/16
				value := 0.0
				if i == j {
					value = a0
				} else if i%2 == 0 && j == i+1 {
					value = -b0
				} else if i%2 == 1 && j == i-1 {
					value = b0
				} else if i < j {
					value = float64((i*3+j*5+seed)%11-5) / 32
				}
				if err := a.Set(i, j, value); err != nil {
					return err
				}
			}
		}
	}
	if op == "eigen_symmetric" {
		for row := 0; row < size; row++ {
			for col := 0; col < size; col++ {
				value := 0.0
				if row == col {
					value = 2 + float64(seed%17)/16
				} else if row-col == 1 || col-row == 1 {
					value = -1
				}
				if err := a.Set(row, col, value); err != nil {
					return err
				}
			}
		}
	}
	if op == "determinant_spd_lu" || op == "determinant_cholesky" {
		for row := 0; row < size; row++ {
			for col := row + 1; col < size; col++ {
				left, _ := a.At(row, col)
				right, _ := a.At(col, row)
				value := (left + right) / 2
				if err := a.Set(row, col, value); err != nil {
					return err
				}
				if err := a.Set(col, row, value); err != nil {
					return err
				}
			}
		}
	}
	if op == "determinant" || op == "determinant_lu" || op == "determinant_spd_lu" || op == "determinant_cholesky" {
		for index := 0; index < size; index++ {
			value, _ := a.At(index, index)
			if err := a.Set(index, index, value+float64(size)*4); err != nil {
				return err
			}
		}
	}
	// Profiling is opt-in and occurs on separate runs from ordinary measurements.
	// Samples cover warmup and the timed loop, excluding input preparation.
	var profile *os.File
	if path := os.Getenv("LINEAR_A_CPU_PROFILE"); path != "" {
		profile, err = os.Create(path)
		if err != nil {
			return fmt.Errorf("create CPU profile: %w", err)
		}
		if err = pprof.StartCPUProfile(profile); err != nil {
			profile.Close()
			return fmt.Errorf("start CPU profile: %w", err)
		}
	}
	result, measureErr := measure(op, a, b, iterations)
	if profile != nil {
		pprof.StopCPUProfile()
		if closeErr := profile.Close(); measureErr == nil {
			measureErr = closeErr
		}
	}
	if measureErr != nil {
		return measureErr
	}
	return json.NewEncoder(os.Stdout).Encode(result)
}

func run(args []string) error {
	if len(args) == 0 {
		return errors.New("usage: runner {check|bench} OP ...")
	}
	switch args[0] {
	case "check":
		return check(args)
	case "bench":
		return bench(args)
	}
	return errors.New("usage: runner {check|bench} OP ...")
}

func main() {
	if err := run(os.Args[1:]); err != nil {
		fmt.Fprintln(os.Stderr, "error:", err)
		os.Exit(1)
	}
}
