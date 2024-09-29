// src/quantum_state/mod.rs

use num_complex::Complex64;
use ndarray::{Array1, Array2};

/// Represents a pure quantum state as a complex vector
pub struct PureState {
    amplitudes: Array1<Complex64>,
}

impl PureState {
    /// Create a new pure state with given amplitudes
    pub fn new(amplitudes: Vec<Complex64>) -> Self {
        PureState {
            amplitudes: Array1::from(amplitudes),
        }
    }

    /// Get the number of qubits in the state
    pub fn num_qubits(&self) -> usize {
        (self.amplitudes.len() as f64).log2() as usize
    }

    /// Calculate the probability of measuring a specific basis state
    pub fn probability(&self, basis_state: usize) -> f64 {
        self.amplitudes[basis_state].norm_sqr()
    }
}

/// Represents a mixed quantum state as a density matrix
pub struct MixedState {
    density_matrix: Array2<Complex64>,
}

impl MixedState {
    /// Create a new mixed state with given density matrix
    pub fn new(density_matrix: Array2<Complex64>) -> Self {
        MixedState {
            density_matrix,
        }
    }

    /// Get the number of qubits in the state
    pub fn num_qubits(&self) -> usize {
        (self.density_matrix.nrows() as f64).log2() as usize
    }

    /// Calculate the probability of measuring a specific basis state
    pub fn probability(&self, basis_state: usize) -> f64 {
        self.density_matrix[[basis_state, basis_state]].re
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pure_state() {
        let state = PureState::new(vec![
            Complex64::new(1.0 / 2.0_f64.sqrt(), 0.0),
            Complex64::new(0.0, 1.0 / 2.0_f64.sqrt()),
        ]);
        assert_eq!(state.num_qubits(), 1);
        assert_eq!(state.probability(0), 0.5);
        assert_eq!(state.probability(1), 0.5);
    }

    #[test]
    fn test_mixed_state() {
        let density_matrix = Array2::from_shape_vec(
            (2, 2),
            vec![
                Complex64::new(0.5, 0.0), Complex64::new(0.0, 0.0),
                Complex64::new(0.0, 0.0), Complex64::new(0.5, 0.0),
            ],
        ).unwrap();
        let state = MixedState::new(density_matrix);
        assert_eq!(state.num_qubits(), 1);
        assert_eq!(state.probability(0), 0.5);
        assert_eq!(state.probability(1), 0.5);
    }
}
