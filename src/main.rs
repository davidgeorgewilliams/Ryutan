mod reservoir;
mod quantum_state;
mod quantum_operations;
mod hardware;
mod visualization;
mod training;

use quantum_state::PureState;
use quantum_operations::{Hadamard, PauliX, CNOT};
use reservoir::QuantumReservoir;
use training::{Trainer, TrainingData, TaskType, FeatureExtractionStrategy};
use num_complex::Complex64;

fn main() {
    println!("=== Quantum Reservoir Computing System ===\n");

    // Initialize a two-qubit quantum reservoir
    println!("Initializing quantum reservoir with 2 qubits...");
    let initial_state = PureState::new(vec![
        Complex64::new(1.0, 0.0),
        Complex64::new(0.0, 0.0),
        Complex64::new(0.0, 0.0),
        Complex64::new(0.0, 0.0),
    ]);
    let mut reservoir = QuantumReservoir::new(initial_state);

    // Configure quantum operations for reservoir dynamics
    reservoir.add_single_qubit_operation(Box::new(Hadamard));
    reservoir.add_multi_qubit_operation(Box::new(CNOT), vec![0, 1]);
    reservoir.add_single_qubit_operation(Box::new(PauliX));

    println!("Reservoir configured with {} quantum operations\n", reservoir.operation_count());

    // Generate training data for a nonlinear function approximation task
    println!("Generating training data...");
    let mut training_inputs = Vec::new();
    let mut training_outputs = Vec::new();

    for i in 0..20 {
        let x = i as f64 / 20.0;
        let y = (x * std::f64::consts::PI).sin();
        training_inputs.push(vec![x, 1.0 - x]);
        training_outputs.push(vec![y]);
    }

    let training_data = TrainingData::new(training_inputs, training_outputs);
    println!("Created training dataset with {} samples\n", training_data.len());

    // Initialize trainer with ridge regression
    println!("Initializing trainer...");
    let mut trainer = Trainer::new(
        TaskType::Regression,
        FeatureExtractionStrategy::Combined,
        0.1,
    );
    println!("Task type: Regression");
    println!("Feature strategy: Combined (expectations + probabilities)");
    println!("Ridge parameter: 0.1\n");

    // Train the model
    println!("Training quantum reservoir computing model...");
    match trainer.train(&mut reservoir, &training_data) {
        Ok(()) => println!("Training completed successfully!\n"),
        Err(e) => {
            eprintln!("Training failed: {}", e);
            return;
        }
    }

    // Generate test data
    println!("Generating test data...");
    let mut test_inputs = Vec::new();
    let mut test_outputs = Vec::new();

    for i in 0..10 {
        let x = (i as f64 + 0.5) / 10.0;
        let y = (x * std::f64::consts::PI).sin();
        test_inputs.push(vec![x, 1.0 - x]);
        test_outputs.push(vec![y]);
    }

    let test_data = TrainingData::new(test_inputs, test_outputs);
    println!("Created test dataset with {} samples\n", test_data.len());

    // Evaluate model performance
    println!("Evaluating model performance...");
    match trainer.evaluate(&mut reservoir, &test_data) {
        Ok((mse, mae)) => {
            println!("Test Set Performance:");
            println!("  Mean Squared Error: {:.6}", mse);
            println!("  Mean Absolute Error: {:.6}\n", mae);
        }
        Err(e) => {
            eprintln!("Evaluation failed: {}", e);
            return;
        }
    }

    // Demonstrate predictions on specific inputs
    println!("Sample predictions:");
    let sample_inputs = vec![
        vec![0.25, 0.75],
        vec![0.50, 0.50],
        vec![0.75, 0.25],
    ];

    for input in sample_inputs {
        match trainer.predict(&mut reservoir, &input) {
            Ok(prediction) => {
                let x = input[0];
                let true_value = (x * std::f64::consts::PI).sin();
                println!("  Input: [{:.2}, {:.2}] → Prediction: {:.4}, True: {:.4}",
                         input[0], input[1], prediction[0], true_value);
            }
            Err(e) => eprintln!("  Prediction failed: {}", e),
        }
    }

    println!("\n=== Demonstration Complete ===");
}