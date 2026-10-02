// Package matrix provides owned, row-major float64 matrices.
package matrix

import (
	"errors"
	"math"
)

// Matrix owns its storage. Its zero value is the empty 0x0 matrix. Copy produces
// an independent matrix; assigning a Matrix value is a shallow Go struct copy.
type Matrix struct {
	rows, cols int
	values     []float64
}

func count(rows, cols int) (int, error) {
	if rows < 0 || cols < 0 {
		return 0, errors.New("matrix dimensions must be nonnegative")
	}
	maxElements := int(^uint(0)>>1) / 8
	if cols != 0 && rows > maxElements/cols {
		return 0, errors.New("matrix dimensions overflow storage size")
	}
	return rows * cols, nil
}

func isFinite(value float64) bool {
	return !math.IsNaN(value) && !math.IsInf(value, 0)
}

// New copies row-major values, rejecting invalid shapes and nonfinite numbers.
func New(rows, cols int, values []float64) (*Matrix, error) {
	length, err := count(rows, cols)
	if err != nil {
		return nil, err
	}
	if len(values) != length {
		return nil, errors.New("value count does not match matrix dimensions")
	}
	for _, value := range values {
		if !isFinite(value) {
			return nil, errors.New("matrix values must be finite")
		}
	}
	storage := make([]float64, length)
	copy(storage, values)
	return &Matrix{rows: rows, cols: cols, values: storage}, nil
}

func Zeros(rows, cols int) (*Matrix, error) {
	length, err := count(rows, cols)
	if err != nil {
		return nil, err
	}
	return &Matrix{rows: rows, cols: cols, values: make([]float64, length)}, nil
}

func Identity(size int) (*Matrix, error) {
	result, err := Zeros(size, size)
	if err != nil {
		return nil, err
	}
	for index := 0; index < size; index++ {
		result.values[index*size+index] = 1
	}
	return result, nil
}

func (m *Matrix) Rows() int { return m.rows }
func (m *Matrix) Cols() int { return m.cols }

// Values returns a copy in row-major order.
func (m *Matrix) Values() []float64 {
	values := make([]float64, len(m.values))
	copy(values, m.values)
	return values
}

func (m *Matrix) Copy() *Matrix {
	return &Matrix{rows: m.rows, cols: m.cols, values: m.Values()}
}

func (m *Matrix) offset(row, col int) (int, error) {
	if row < 0 || row >= m.rows || col < 0 || col >= m.cols {
		return 0, errors.New("matrix index is out of range")
	}
	return row*m.cols + col, nil
}

func (m *Matrix) At(row, col int) (float64, error) {
	index, err := m.offset(row, col)
	if err != nil {
		return 0, err
	}
	return m.values[index], nil
}

func (m *Matrix) Set(row, col int, value float64) error {
	index, err := m.offset(row, col)
	if err != nil {
		return err
	}
	if !isFinite(value) {
		return errors.New("matrix values must be finite")
	}
	m.values[index] = value
	return nil
}

func (m *Matrix) Row(row int) ([]float64, error) {
	if row < 0 || row >= m.rows {
		return nil, errors.New("row index is out of range")
	}
	values := make([]float64, m.cols)
	copy(values, m.values[row*m.cols:(row+1)*m.cols])
	return values, nil
}

func (m *Matrix) Column(col int) ([]float64, error) {
	if col < 0 || col >= m.cols {
		return nil, errors.New("column index is out of range")
	}
	values := make([]float64, m.rows)
	for row := range values {
		values[row] = m.values[row*m.cols+col]
	}
	return values, nil
}

func (m *Matrix) sameShape(other *Matrix) error {
	if other == nil || m.rows != other.rows || m.cols != other.cols {
		return errors.New("matrix dimensions must match")
	}
	return nil
}

func (m *Matrix) Add(other *Matrix) (*Matrix, error) {
	if err := m.sameShape(other); err != nil {
		return nil, err
	}
	result, _ := Zeros(m.rows, m.cols)
	for index, value := range m.values {
		result.values[index] = value + other.values[index]
		if !isFinite(result.values[index]) {
			return nil, errors.New("matrix addition overflowed")
		}
	}
	return result, nil
}

func (m *Matrix) Subtract(other *Matrix) (*Matrix, error) {
	if err := m.sameShape(other); err != nil {
		return nil, err
	}
	result, _ := Zeros(m.rows, m.cols)
	for index, value := range m.values {
		result.values[index] = value - other.values[index]
		if !isFinite(result.values[index]) {
			return nil, errors.New("matrix subtraction overflowed")
		}
	}
	return result, nil
}

func (m *Matrix) Scale(scalar float64) (*Matrix, error) {
	if !isFinite(scalar) {
		return nil, errors.New("scale must be finite")
	}
	result, _ := Zeros(m.rows, m.cols)
	for index, value := range m.values {
		result.values[index] = value * scalar
		if !isFinite(result.values[index]) {
			return nil, errors.New("matrix scaling overflowed")
		}
	}
	return result, nil
}

func (m *Matrix) Transpose() *Matrix {
	result, _ := Zeros(m.cols, m.rows)
	for row := 0; row < m.rows; row++ {
		for col := 0; col < m.cols; col++ {
			result.values[col*m.rows+row] = m.values[row*m.cols+col]
		}
	}
	return result
}

func (m *Matrix) Multiply(other *Matrix) (*Matrix, error) {
	if other == nil || m.cols != other.rows {
		return nil, errors.New("left columns must equal right rows")
	}
	result, err := Zeros(m.rows, other.cols)
	if err != nil {
		return nil, err
	}
	for row := 0; row < m.rows; row++ {
		for inner := 0; inner < m.cols; inner++ {
			left := m.values[row*m.cols+inner]
			for col := 0; col < other.cols; col++ {
				index := row*other.cols + col
				// Keep multiplication and addition separate for reproducible float64
				// rounding across implementations (Go otherwise permits fusion).
				product := float64(left * other.values[inner*other.cols+col])
				result.values[index] += product
				if !isFinite(result.values[index]) {
					return nil, errors.New("matrix multiplication overflowed")
				}
			}
		}
	}
	return result, nil
}

func (m *Matrix) square() error {
	if m.rows != m.cols {
		return errors.New("operation requires a square matrix")
	}
	return nil
}

func (m *Matrix) Trace() (float64, error) {
	if err := m.square(); err != nil {
		return 0, err
	}
	result := 0.0
	for index := 0; index < m.rows; index++ {
		result += m.values[index*m.cols+index]
		if !isFinite(result) {
			return 0, errors.New("matrix trace overflowed")
		}
	}
	return result, nil
}

// Triangular returns (upper, lower), using exact zeros. Diagonal and 0x0
// matrices are both; rectangular matrices are neither.
func (m *Matrix) Triangular() (upper, lower bool) {
	if m.rows != m.cols {
		return false, false
	}
	upper, lower = true, true
	for row := 0; row < m.rows; row++ {
		for col := 0; col < m.cols; col++ {
			if row > col && m.values[row*m.cols+col] != 0 {
				upper = false
			}
			if row < col && m.values[row*m.cols+col] != 0 {
				lower = false
			}
		}
	}
	return upper, lower
}
