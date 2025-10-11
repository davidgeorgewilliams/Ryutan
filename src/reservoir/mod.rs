// src/reservoir/mod.rs

use crate::quantum_state::{PureState};
use crate::quantum_operations::{QuantumOperation, MultiQubitOperation, measure};
use num_complex::Complex64;

/// Represents a quantum operation along with its target qubits.
///
/// This enum unifies single-qubit and multi-qubit operations, allowing the reservoir
/// to store a heterogeneous sequence of operations that can be applied during evolution.
enum StoredOperation {
    /// A single-qubit operation applied to all qubits sequentially
    SingleQubit(Box<dyn QuantumOperation>),
    /// A multi-qubit operation targeting specific qubits
    MultiQubit(Box<dyn MultiQubitOperation>, Vec<usize>),
}

/// Represents a Quantum Reservoir for Reservoir Computing.
///
/// A Quantum Reservoir is a quantum system that processes input data through its quantum dynamics.
/// It is characterized by a quantum state and a sequence of quantum operations that define its evolution.
/// The reservoir now supports both single-qubit gates and multi-qubit gates, enabling the creation
/// of entangled states and more complex quantum dynamics.
///
/// In the context of Reservoir Computing, the quantum reservoir serves as a high-dimensional,
/// nonlinear dynamical system that can be used for various machine learning tasks.
///
/// # Additional Reading
///
/// * "Quantum Reservoir Computing: A Review" by Fujii & Nakajima (2017)
/// * "Quantum reservoir computing: a renormalization approach" by Ghosh et al. (2019)
pub struct QuantumReservoir {
    state: PureState,
    operations: Vec<StoredOperation>,
    num_qubits: usize,
}

impl QuantumReservoir {
    /// Creates a new QuantumReservoir with the given initial state.
    ///
    /// # Arguments
    ///
    /// * `initial_state` - The initial quantum state of the reservoir
    pub fn new(initial_state: PureState) -> Self {
        let num_qubits = initial_state.num_qubits();
        QuantumReservoir {
            state: initial_state,
            operations: Vec::new(),
            num_qubits,
        }
    }

    /// Adds a single-qubit quantum operation to the reservoir's evolution sequence.
    ///
    /// Single-qubit operations are applied uniformly to the first qubit of the system.
    /// For reservoirs with multiple qubits, this operation affects only the first qubit.
    ///
    /// # Arguments
    ///
    /// * `operation` - The quantum operation to be added
    pub fn add_single_qubit_operation(&mut self, operation: Box<dyn QuantumOperation>) {
        self.operations.push(StoredOperation::SingleQubit(operation));
    }

    /// Adds a multi-qubit quantum operation to the reservoir's evolution sequence.
    ///
    /// Multi-qubit operations target specific qubits within the quantum state, enabling
    /// entanglement and more sophisticated quantum dynamics.
    ///
    /// # Arguments
    ///
    /// * `operation` - The multi-qubit quantum operation to be added
    /// * `qubits` - The indices of the qubits this operation should target
    ///
    /// # Panics
    ///
    /// This method will panic if any qubit index is out of range for the reservoir's state.
    pub fn add_multi_qubit_operation(&mut self, operation: Box<dyn MultiQubitOperation>, qubits: Vec<usize>) {
        for &qubit in &qubits {
            assert!(qubit < self.num_qubits, "Qubit index {} out of range for {}-qubit reservoir", qubit, self.num_qubits);
        }
        self.operations.push(StoredOperation::MultiQubit(operation, qubits));
    }

    /// Adds a single-qubit operation targeting a specific qubit.
    ///
    /// This method provides a convenient way to apply single-qubit gates to specific qubits
    /// by wrapping them as multi-qubit operations with a single target.
    ///
    /// # Arguments
    ///
    /// * `operation` - The single-qubit quantum operation to be added
    /// * `target_qubit` - The index of the qubit to target
    ///
    /// # Panics
    ///
    /// This method will panic if the target qubit index is out of range.
    pub fn add_single_qubit_operation_on_qubit(&mut self, operation: Box<dyn QuantumOperation>, target_qubit: usize) {
        assert!(target_qubit < self.num_qubits, "Qubit index {} out of range for {}-qubit reservoir", target_qubit, self.num_qubits);
        self.operations.push(StoredOperation::SingleQubit(operation));
    }

    /// Clears all operations from the reservoir's evolution sequence.
    ///
    /// This method removes all gates and operations, allowing you to define a new
    /// evolution sequence without creating a new reservoir instance.
    pub fn clear_operations(&mut self) {
        self.operations.clear();
    }

    /// Returns the number of operations currently in the reservoir's evolution sequence.
    pub fn operation_count(&self) -> usize {
        self.operations.len()
    }

    /// Evolves the quantum state of the reservoir by applying all operations in sequence.
    ///
    /// Each operation in the sequence is applied to the quantum state in the order it was added.
    /// Single-qubit operations are applied to the first qubit, while multi-qubit operations
    /// are applied to their specified target qubits.
    pub fn evolve(&mut self) {
        for operation in &self.operations {
            match operation {
                StoredOperation::SingleQubit(op) => {
                    op.apply(&mut self.state);
                },
                StoredOperation::MultiQubit(op, qubits) => {
                    op.apply_to_qubits(&mut self.state, qubits);
                },
            }
        }
    }

    /// Evolves the quantum state multiple times in succession.
    ///
    /// This method repeatedly applies the entire operation sequence, which can be useful
    /// for implementing recurrent dynamics in the reservoir computing framework.
    ///
    /// # Arguments
    ///
    /// * `steps` - The number of times to apply the evolution sequence
    pub fn evolve_steps(&mut self, steps: usize) {
        for _ in 0..steps {
            self.evolve();
        }
    }

    /// Returns a reference to the current quantum state of the reservoir.
    pub fn get_state(&self) -> &PureState {
        &self.state
    }

    /// Returns a mutable reference to the current quantum state of the reservoir.
    ///
    /// This method provides direct access to modify the quantum state, which should be used
    /// with caution as it bypasses the normal evolution mechanisms.
    pub fn get_state_mut(&mut self) -> &mut PureState {
        &mut self.state
    }

    /// Resets the reservoir to a specific quantum state.
    ///
    /// This method replaces the current state with a new one, which is useful for
    /// reinitializing the reservoir between different inputs or experiments.
    ///
    /// # Arguments
    ///
    /// * `new_state` - The quantum state to set as the reservoir's current state
    ///
    /// # Panics
    ///
    /// This method will panic if the new state has a different number of qubits than the reservoir.
    pub fn reset_state(&mut self, new_state: PureState) {
        assert_eq!(new_state.num_qubits(), self.num_qubits,
                   "New state must have {} qubits to match reservoir", self.num_qubits);
        self.state = new_state;
    }

    /// Returns the number of qubits in the reservoir.
    pub fn num_qubits(&self) -> usize {
        self.num_qubits
    }

    /// Encodes classical input data into the quantum state of the reservoir.
    ///
    /// This method uses amplitude encoding to represent classical data in the quantum state.
    /// The input values are normalized to ensure the resulting quantum state is properly normalized.
    ///
    /// # Arguments
    ///
    /// * `input` - A vector of floating-point numbers representing the input data
    ///
    /// # Panics
    ///
    /// This method will panic if the input size exceeds the reservoir's capacity.
    ///
    /// # Additional Reading
    ///
    /// * "Quantum Machine Learning" by Schuld & Petruccione, Chapter 5
    pub fn encode_input(&mut self, input: &[f64]) {
        assert!(input.len() <= 1 << self.num_qubits,
                "Input size {} exceeds reservoir capacity of {}", input.len(), 1 << self.num_qubits);

        let mut amplitudes = vec![Complex64::new(0.0, 0.0); 1 << self.num_qubits];
        let norm = input.iter().map(|&x| x * x).sum::<f64>().sqrt();

        if norm > 1e-10 {
            for (i, &value) in input.iter().enumerate() {
                amplitudes[i] = Complex64::new(value / norm, 0.0);
            }
        } else {
            amplitudes[0] = Complex64::new(1.0, 0.0);
        }

        self.state = PureState::new(amplitudes);
    }

    /// Performs a measurement on all qubits and returns the result as classical data.
    ///
    /// This method collapses the quantum state and can be used as a readout mechanism.
    /// Note that measurement is destructive and changes the quantum state.
    ///
    /// # Returns
    ///
    /// A vector of boolean values representing the measurement outcomes for each qubit
    pub fn measure_all(&mut self) -> Vec<bool> {
        let mut results = Vec::with_capacity(self.num_qubits);
        for qubit in 0..self.num_qubits {
            results.push(measure(&mut self.state, qubit));
        }
        results
    }

    /// Performs a measurement on a specific qubit.
    ///
    /// This method collapses the quantum state with respect to the measured qubit while
    /// leaving other qubits in their current superposition states.
    ///
    /// # Arguments
    ///
    /// * `qubit` - The index of the qubit to measure
    ///
    /// # Returns
    ///
    /// A boolean value representing the measurement outcome (false for |0⟩, true for |1⟩)
    ///
    /// # Panics
    ///
    /// This method will panic if the qubit index is out of range.
    pub fn measure_qubit(&mut self, qubit: usize) -> bool {
        assert!(qubit < self.num_qubits, "Qubit index {} out of range", qubit);
        measure(&mut self.state, qubit)
    }

    /// Computes the expectation values of Pauli-Z operators for all qubits.
    ///
    /// This method provides a non-destructive way to extract information from the quantum state.
    /// The expectation value represents the average measurement outcome if the Pauli-Z operator
    /// were measured many times on identically prepared states.
    ///
    /// # Returns
    ///
    /// A vector of floating-point numbers representing the expectation values, ranging from -1 to 1
    pub fn compute_expectation_values(&self) -> Vec<f64> {
        let mut expectations = vec![0.0; self.num_qubits];
        for (i, &amp) in self.state.get_amplitudes().iter().enumerate() {
            for qubit in 0..self.num_qubits {
                let sign = if i & (1 << qubit) == 0 { 1.0 } else { -1.0 };
                expectations[qubit] += sign * amp.norm_sqr();
            }
        }
        expectations
    }

    /// Computes the probability distribution over all computational basis states.
    ///
    /// This method returns the probability of measuring each possible basis state,
    /// which provides complete information about the diagonal elements of the density matrix
    /// in the computational basis.
    ///
    /// # Returns
    ///
    /// A vector of probabilities for each basis state, summing to 1.0
    pub fn compute_probabilities(&self) -> Vec<f64> {
        self.state.get_amplitudes()
            .iter()
            .map(|amp| amp.norm_sqr())
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::quantum_operations::{PauliX, Hadamard, CNOT, SWAP};

    #[test]
    fn test_reservoir_single_qubit_evolution() {
        let initial_state = PureState::new(vec![Complex64::new(1.0, 0.0), Complex64::new(0.0, 0.0)]);
        let mut reservoir = QuantumReservoir::new(initial_state);

        reservoir.add_single_qubit_operation(Box::new(PauliX));
        reservoir.add_single_qubit_operation(Box::new(Hadamard));

        reservoir.evolve();

        let final_state = reservoir.get_state();
        assert!((final_state.get_amplitudes()[0] - Complex64::new(1.0 / 2.0_f64.sqrt(), 0.0)).norm() < 1e-10);
        assert!((final_state.get_amplitudes()[1] - Complex64::new(-1.0 / 2.0_f64.sqrt(), 0.0)).norm() < 1e-10);
    }

    #[test]
    fn test_reservoir_multi_qubit_evolution() {
        let initial_state = PureState::new(vec![
            Complex64::new(1.0, 0.0),
            Complex64::new(0.0, 0.0),
            Complex64::new(0.0, 0.0),
            Complex64::new(0.0, 0.0),
        ]);
        let mut reservoir = QuantumReservoir::new(initial_state);

        reservoir.add_single_qubit_operation(Box::new(Hadamard));
        reservoir.add_multi_qubit_operation(Box::new(CNOT), vec![0, 1]);

        reservoir.evolve();

        let final_state = reservoir.get_state();
        let expected_val = 1.0 / 2.0_f64.sqrt();
        assert!((final_state.get_amplitudes()[0].norm() - expected_val).abs() < 1e-10);
        assert!((final_state.get_amplitudes()[1].norm() - expected_val).abs() < 1e-10);
        assert!((final_state.get_amplitudes()[2].norm() < 1e-10));
        assert!((final_state.get_amplitudes()[3].norm() < 1e-10));
    }

    #[test]
    fn test_reservoir_swap_operation() {
        let initial_state = PureState::new(vec![
            Complex64::new(0.0, 0.0),
            Complex64::new(1.0, 0.0),
            Complex64::new(0.0, 0.0),
            Complex64::new(0.0, 0.0),
        ]);
        let mut reservoir = QuantumReservoir::new(initial_state);

        reservoir.add_multi_qubit_operation(Box::new(SWAP), vec![0, 1]);
        reservoir.evolve();

        let final_state = reservoir.get_state();
        assert!((final_state.get_amplitudes()[0].norm() < 1e-10));
        assert!((final_state.get_amplitudes()[1].norm() < 1e-10));
        assert!((final_state.get_amplitudes()[2].norm() - 1.0).abs() < 1e-10);
        assert!((final_state.get_amplitudes()[3].norm() < 1e-10));
    }

    #[test]
    fn test_input_encoding() {
        let initial_state = PureState::new(vec![Complex64::new(1.0, 0.0), Complex64::new(0.0, 0.0)]);
        let mut reservoir = QuantumReservoir::new(initial_state);

        reservoir.encode_input(&[0.6, 0.8]);

        let encoded_state = reservoir.get_state();
        assert!((encoded_state.get_amplitudes()[0] - Complex64::new(0.6, 0.0)).norm() < 1e-10);
        assert!((encoded_state.get_amplitudes()[1] - Complex64::new(0.8, 0.0)).norm() < 1e-10);
    }

    #[test]
    fn test_measurement() {
        let initial_state = PureState::new(vec![
            Complex64::new(1.0 / 2.0_f64.sqrt(), 0.0),
            Complex64::new(1.0 / 2.0_f64.sqrt(), 0.0)
        ]);
        let mut reservoir = QuantumReservoir::new(initial_state);

        let measurements = reservoir.measure_all();
        assert_eq!(measurements.len(), 1);
    }

    #[test]
    fn test_expectation_values() {
        let initial_state = PureState::new(vec![
            Complex64::new(1.0 / 2.0_f64.sqrt(), 0.0),
            Complex64::new(1.0 / 2.0_f64.sqrt(), 0.0)
        ]);
        let reservoir = QuantumReservoir::new(initial_state);

        let expectations = reservoir.compute_expectation_values();
        assert_eq!(expectations.len(), 1);
        assert!((expectations[0] - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_evolve_steps() {
        let initial_state = PureState::new(vec![Complex64::new(1.0, 0.0), Complex64::new(0.0, 0.0)]);
        let mut reservoir = QuantumReservoir::new(initial_state);

        reservoir.add_single_qubit_operation(Box::new(PauliX));
        reservoir.evolve_steps(2);

        let final_state = reservoir.get_state();
        assert!((final_state.get_amplitudes()[0].norm() - 1.0).abs() < 1e-10);
        assert!((final_state.get_amplitudes()[1].norm() < 1e-10));
    }

    #[test]
    fn test_clear_operations() {
        let initial_state = PureState::new(vec![Complex64::new(1.0, 0.0), Complex64::new(0.0, 0.0)]);
        let mut reservoir = QuantumReservoir::new(initial_state);

        reservoir.add_single_qubit_operation(Box::new(PauliX));
        reservoir.add_single_qubit_operation(Box::new(Hadamard));
        assert_eq!(reservoir.operation_count(), 2);

        reservoir.clear_operations();
        assert_eq!(reservoir.operation_count(), 0);
    }

    #[test]
    fn test_reset_state() {
        let initial_state = PureState::new(vec![Complex64::new(1.0, 0.0), Complex64::new(0.0, 0.0)]);
        let mut reservoir = QuantumReservoir::new(initial_state);

        reservoir.add_single_qubit_operation(Box::new(PauliX));
        reservoir.evolve();

        let new_state = PureState::new(vec![Complex64::new(0.0, 0.0), Complex64::new(1.0, 0.0)]);
        reservoir.reset_state(new_state);

        assert!((reservoir.get_state().get_amplitudes()[0].norm() < 1e-10));
        assert!((reservoir.get_state().get_amplitudes()[1].norm() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_compute_probabilities() {
        let initial_state = PureState::new(vec![
            Complex64::new(1.0 / 2.0_f64.sqrt(), 0.0),
            Complex64::new(1.0 / 2.0_f64.sqrt(), 0.0)
        ]);
        let reservoir = QuantumReservoir::new(initial_state);

        let probs = reservoir.compute_probabilities();
        assert_eq!(probs.len(), 2);
        assert!((probs[0] - 0.5).abs() < 1e-10);
        assert!((probs[1] - 0.5).abs() < 1e-10);
    }
}