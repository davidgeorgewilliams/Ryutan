// src/training/mod.rs

use crate::reservoir::QuantumReservoir;
use ndarray::{Array1, Array2, s};
use ndarray_linalg::Inverse;

/// Represents the type of machine learning task.
///
/// The training module supports both regression tasks, where outputs are continuous values,
/// and classification tasks, where outputs are discrete class labels.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum TaskType {
    /// Regression task with continuous outputs
    Regression,
    /// Classification task with discrete class labels
    Classification,
}

/// Strategy for extracting features from quantum states.
///
/// Different feature extraction methods capture different aspects of the quantum state
/// and may be more suitable for different types of problems.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum FeatureExtractionStrategy {
    /// Use expectation values of Pauli-Z operators for each qubit
    ExpectationValues,
    /// Use probability distributions over computational basis states
    Probabilities,
    /// Use both expectation values and probabilities concatenated
    Combined,
}

/// Training data structure holding input-output pairs.
///
/// This structure organizes the dataset for supervised learning tasks in quantum
/// reservoir computing. Each input is processed through the quantum reservoir to
/// extract features, which are then used to train the classical readout layer.
#[derive(Clone, Debug)]
pub struct TrainingData {
    /// Input samples, where each sample is a vector of features
    pub inputs: Vec<Vec<f64>>,
    /// Output targets corresponding to each input
    pub outputs: Vec<Vec<f64>>,
}

impl TrainingData {
    /// Creates a new TrainingData structure.
    ///
    /// # Arguments
    ///
    /// * `inputs` - Vector of input samples
    /// * `outputs` - Vector of output targets
    ///
    /// # Panics
    ///
    /// Panics if the number of inputs does not match the number of outputs.
    pub fn new(inputs: Vec<Vec<f64>>, outputs: Vec<Vec<f64>>) -> Self {
        assert_eq!(
            inputs.len(),
            outputs.len(),
            "Number of inputs must match number of outputs"
        );
        TrainingData { inputs, outputs }
    }

    /// Returns the number of samples in the dataset.
    pub fn len(&self) -> usize {
        self.inputs.len()
    }

    /// Returns whether the dataset is empty.
    pub fn is_empty(&self) -> bool {
        self.inputs.is_empty()
    }

    /// Returns the dimensionality of input features.
    pub fn input_dim(&self) -> usize {
        if self.inputs.is_empty() {
            0
        } else {
            self.inputs[0].len()
        }
    }

    /// Returns the dimensionality of output targets.
    pub fn output_dim(&self) -> usize {
        if self.outputs.is_empty() {
            0
        } else {
            self.outputs[0].len()
        }
    }
}

/// Quantum Reservoir Trainer for machine learning tasks.
///
/// The Trainer manages the complete training pipeline for quantum reservoir computing,
/// including feature extraction from quantum states, optimization of readout weights,
/// and prediction on new data. The quantum reservoir itself remains fixed during training,
/// with only the classical linear readout layer being optimized.
///
/// # Additional Reading
///
/// * "Quantum Reservoir Computing: A Review" by Fujii & Nakajima (2017)
/// * "Reservoir Computing approaches to recurrent neural network training" by Lukoševičius & Jaeger (2009)
pub struct Trainer {
    /// Weights of the linear readout layer
    weights: Option<Array2<f64>>,
    /// Bias terms for the readout layer
    bias: Option<Array1<f64>>,
    /// Type of machine learning task
    task_type: TaskType,
    /// Strategy for extracting features from quantum states
    feature_strategy: FeatureExtractionStrategy,
    /// Regularization parameter for ridge regression
    ridge_parameter: f64,
}

impl Trainer {
    /// Creates a new Trainer with specified configuration.
    ///
    /// # Arguments
    ///
    /// * `task_type` - The type of machine learning task (regression or classification)
    /// * `feature_strategy` - Strategy for extracting features from quantum states
    /// * `ridge_parameter` - Regularization parameter (lambda) for ridge regression
    pub fn new(
        task_type: TaskType,
        feature_strategy: FeatureExtractionStrategy,
        ridge_parameter: f64,
    ) -> Self {
        Trainer {
            weights: None,
            bias: None,
            task_type,
            feature_strategy,
            ridge_parameter,
        }
    }

    /// Extracts features from a quantum reservoir state.
    ///
    /// This method applies the configured feature extraction strategy to obtain
    /// a classical feature vector from the quantum state. The feature vector serves
    /// as input to the classical readout layer.
    ///
    /// # Arguments
    ///
    /// * `reservoir` - Reference to the quantum reservoir
    ///
    /// # Returns
    ///
    /// A feature vector extracted from the current quantum state
    fn extract_features(&self, reservoir: &QuantumReservoir) -> Vec<f64> {
        match self.feature_strategy {
            FeatureExtractionStrategy::ExpectationValues => {
                reservoir.compute_expectation_values()
            }
            FeatureExtractionStrategy::Probabilities => {
                reservoir.compute_probabilities()
            }
            FeatureExtractionStrategy::Combined => {
                let mut features = reservoir.compute_expectation_values();
                features.extend(reservoir.compute_probabilities());
                features
            }
        }
    }

    /// Trains the readout layer using ridge regression.
    ///
    /// This method implements the training procedure for quantum reservoir computing.
    /// For each input sample, the reservoir is reset, the input is encoded, the reservoir
    /// evolves according to its quantum dynamics, and features are extracted. These features
    /// are then used to train a linear readout layer via ridge regression, which provides
    /// a closed-form solution that is both efficient and numerically stable.
    ///
    /// # Arguments
    ///
    /// * `reservoir` - Mutable reference to the quantum reservoir
    /// * `data` - Training data containing input-output pairs
    ///
    /// # Returns
    ///
    /// Result indicating success or failure of the training process
    ///
    /// # Additional Reading
    ///
    /// * "Ridge Regression: Biased Estimation for Nonorthogonal Problems" by Hoerl & Kennard (1970)
    pub fn train(
        &mut self,
        reservoir: &mut QuantumReservoir,
        data: &TrainingData,
    ) -> Result<(), String> {
        if data.is_empty() {
            return Err("Training data is empty".to_string());
        }

        let num_samples = data.len();
        let output_dim = data.output_dim();

        let mut feature_matrix = Vec::new();
        let mut target_matrix = Vec::new();

        for i in 0..num_samples {
            reservoir.encode_input(&data.inputs[i]);
            reservoir.evolve();

            let features = self.extract_features(reservoir);
            feature_matrix.push(features);
            target_matrix.push(data.outputs[i].clone());
        }

        let feature_dim = feature_matrix[0].len();

        let x = Array2::from_shape_fn((num_samples, feature_dim + 1), |(i, j)| {
            if j < feature_dim {
                feature_matrix[i][j]
            } else {
                1.0
            }
        });

        let y = Array2::from_shape_fn((num_samples, output_dim), |(i, j)| {
            target_matrix[i][j]
        });

        let xt = x.t();
        let xtx = xt.dot(&x);

        let mut regularization = Array2::eye(feature_dim + 1);
        regularization[[feature_dim, feature_dim]] = 0.0;
        let lambda_i = regularization * self.ridge_parameter;

        let xtx_reg = xtx + lambda_i;

        match xtx_reg.inv() {
            Ok(xtx_inv) => {
                let w = xtx_inv.dot(&xt).dot(&y);

                self.weights = Some(w.slice(s![0..feature_dim, ..]).to_owned());
                self.bias = Some(w.slice(s![feature_dim, ..]).to_owned());

                Ok(())
            }
            Err(_) => Err("Failed to invert matrix during training".to_string()),
        }
    }

    /// Makes a prediction on a single input sample.
    ///
    /// This method processes the input through the quantum reservoir, extracts features,
    /// and applies the trained readout layer to produce a prediction. For classification
    /// tasks, the output represents class probabilities or logits that can be converted
    /// to class labels using argmax.
    ///
    /// # Arguments
    ///
    /// * `reservoir` - Mutable reference to the quantum reservoir
    /// * `input` - Input sample to make prediction on
    ///
    /// # Returns
    ///
    /// Result containing the prediction vector or an error message
    pub fn predict(
        &self,
        reservoir: &mut QuantumReservoir,
        input: &[f64],
    ) -> Result<Vec<f64>, String> {
        if self.weights.is_none() || self.bias.is_none() {
            return Err("Model has not been trained yet".to_string());
        }

        reservoir.encode_input(input);
        reservoir.evolve();

        let features = self.extract_features(reservoir);
        let feature_array = Array1::from_vec(features);

        let weights = self.weights.as_ref().unwrap();
        let bias = self.bias.as_ref().unwrap();

        let output = weights.t().dot(&feature_array) + bias;

        Ok(output.to_vec())
    }

    /// Makes predictions on multiple input samples.
    ///
    /// This method applies the trained model to a batch of inputs, returning predictions
    /// for each sample. This is more efficient than calling predict repeatedly for large
    /// datasets.
    ///
    /// # Arguments
    ///
    /// * `reservoir` - Mutable reference to the quantum reservoir
    /// * `inputs` - Vector of input samples
    ///
    /// # Returns
    ///
    /// Result containing vector of prediction vectors or an error message
    pub fn predict_batch(
        &self,
        reservoir: &mut QuantumReservoir,
        inputs: &[Vec<f64>],
    ) -> Result<Vec<Vec<f64>>, String> {
        let mut predictions = Vec::with_capacity(inputs.len());

        for input in inputs {
            let prediction = self.predict(reservoir, input)?;
            predictions.push(prediction);
        }

        Ok(predictions)
    }

    /// Evaluates model performance on a test dataset.
    ///
    /// This method computes appropriate metrics based on the task type. For regression
    /// tasks, it calculates mean squared error and mean absolute error. For classification
    /// tasks, it computes accuracy.
    ///
    /// # Arguments
    ///
    /// * `reservoir` - Mutable reference to the quantum reservoir
    /// * `data` - Test data for evaluation
    ///
    /// # Returns
    ///
    /// Result containing a tuple of evaluation metrics or an error message
    pub fn evaluate(
        &self,
        reservoir: &mut QuantumReservoir,
        data: &TrainingData,
    ) -> Result<(f64, f64), String> {
        if data.is_empty() {
            return Err("Evaluation data is empty".to_string());
        }

        let predictions = self.predict_batch(reservoir, &data.inputs)?;

        match self.task_type {
            TaskType::Regression => {
                let (mse, mae) = self.compute_regression_metrics(&predictions, &data.outputs);
                Ok((mse, mae))
            }
            TaskType::Classification => {
                let accuracy = self.compute_classification_accuracy(&predictions, &data.outputs);
                Ok((accuracy, 0.0))
            }
        }
    }

    /// Computes regression metrics (MSE and MAE).
    fn compute_regression_metrics(
        &self,
        predictions: &[Vec<f64>],
        targets: &[Vec<f64>],
    ) -> (f64, f64) {
        let n = predictions.len() as f64;
        let mut mse = 0.0;
        let mut mae = 0.0;

        for (pred, target) in predictions.iter().zip(targets.iter()) {
            for (p, t) in pred.iter().zip(target.iter()) {
                let error = p - t;
                mse += error * error;
                mae += error.abs();
            }
        }

        mse /= n;
        mae /= n;

        (mse, mae)
    }

    /// Computes classification accuracy.
    fn compute_classification_accuracy(
        &self,
        predictions: &[Vec<f64>],
        targets: &[Vec<f64>],
    ) -> f64 {
        let mut correct = 0;
        let total = predictions.len();

        for (pred, target) in predictions.iter().zip(targets.iter()) {
            let pred_class = argmax(pred);
            let target_class = argmax(target);

            if pred_class == target_class {
                correct += 1;
            }
        }

        correct as f64 / total as f64
    }

    /// Returns whether the model has been trained.
    pub fn is_trained(&self) -> bool {
        self.weights.is_some() && self.bias.is_some()
    }

    /// Returns the current task type.
    pub fn task_type(&self) -> TaskType {
        self.task_type
    }

    /// Returns the current feature extraction strategy.
    pub fn feature_strategy(&self) -> FeatureExtractionStrategy {
        self.feature_strategy
    }
}

/// Helper function to find the index of the maximum value in a vector.
fn argmax(values: &[f64]) -> usize {
    values
        .iter()
        .enumerate()
        .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
        .map(|(index, _)| index)
        .unwrap_or(0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::quantum_state::PureState;
    use crate::quantum_operations::{Hadamard, CNOT};
    use num_complex::Complex64;

    fn create_test_reservoir() -> QuantumReservoir {
        let initial_state = PureState::new(vec![
            Complex64::new(1.0, 0.0),
            Complex64::new(0.0, 0.0),
            Complex64::new(0.0, 0.0),
            Complex64::new(0.0, 0.0),
        ]);
        let mut reservoir = QuantumReservoir::new(initial_state);
        reservoir.add_single_qubit_operation(Box::new(Hadamard));
        reservoir.add_multi_qubit_operation(Box::new(CNOT), vec![0, 1]);
        reservoir
    }

    #[test]
    fn test_training_data_creation() {
        let inputs = vec![vec![1.0, 2.0], vec![3.0, 4.0]];
        let outputs = vec![vec![0.5], vec![1.5]];
        let data = TrainingData::new(inputs, outputs);

        assert_eq!(data.len(), 2);
        assert_eq!(data.input_dim(), 2);
        assert_eq!(data.output_dim(), 1);
    }

    #[test]
    fn test_trainer_creation() {
        let trainer = Trainer::new(
            TaskType::Regression,
            FeatureExtractionStrategy::ExpectationValues,
            0.01,
        );

        assert_eq!(trainer.task_type(), TaskType::Regression);
        assert_eq!(
            trainer.feature_strategy(),
            FeatureExtractionStrategy::ExpectationValues
        );
        assert!(!trainer.is_trained());
    }

    #[test]
    fn test_feature_extraction() {
        let mut reservoir = create_test_reservoir();
        let trainer = Trainer::new(
            TaskType::Regression,
            FeatureExtractionStrategy::ExpectationValues,
            0.01,
        );

        reservoir.encode_input(&[0.6, 0.8]);
        reservoir.evolve();

        let features = trainer.extract_features(&reservoir);
        assert_eq!(features.len(), 2);
    }

    #[test]
    fn test_training_regression() {
        let mut reservoir = create_test_reservoir();
        let mut trainer = Trainer::new(
            TaskType::Regression,
            FeatureExtractionStrategy::ExpectationValues,
            0.01,
        );

        let inputs = vec![
            vec![1.0, 0.0],
            vec![0.0, 1.0],
            vec![0.5, 0.5],
            vec![0.8, 0.2],
        ];
        let outputs = vec![vec![1.0], vec![0.0], vec![0.5], vec![0.8]];
        let data = TrainingData::new(inputs, outputs);

        let result = trainer.train(&mut reservoir, &data);
        assert!(result.is_ok());
        assert!(trainer.is_trained());
    }

    #[test]
    fn test_prediction() {
        let mut reservoir = create_test_reservoir();
        let mut trainer = Trainer::new(
            TaskType::Regression,
            FeatureExtractionStrategy::ExpectationValues,
            0.01,
        );

        let inputs = vec![
            vec![1.0, 0.0],
            vec![0.0, 1.0],
            vec![0.5, 0.5],
        ];
        let outputs = vec![vec![1.0], vec![0.0], vec![0.5]];
        let data = TrainingData::new(inputs, outputs);

        trainer.train(&mut reservoir, &data).unwrap();

        let prediction = trainer.predict(&mut reservoir, &[0.7, 0.3]);
        assert!(prediction.is_ok());
        assert_eq!(prediction.unwrap().len(), 1);
    }

    #[test]
    fn test_batch_prediction() {
        let mut reservoir = create_test_reservoir();
        let mut trainer = Trainer::new(
            TaskType::Regression,
            FeatureExtractionStrategy::ExpectationValues,
            0.01,
        );

        let inputs = vec![
            vec![1.0, 0.0],
            vec![0.0, 1.0],
        ];
        let outputs = vec![vec![1.0], vec![0.0]];
        let data = TrainingData::new(inputs, outputs);

        trainer.train(&mut reservoir, &data).unwrap();

        let test_inputs = vec![vec![0.8, 0.2], vec![0.3, 0.7]];
        let predictions = trainer.predict_batch(&mut reservoir, &test_inputs);

        assert!(predictions.is_ok());
        assert_eq!(predictions.unwrap().len(), 2);
    }

    #[test]
    fn test_evaluation() {
        let mut reservoir = create_test_reservoir();
        let mut trainer = Trainer::new(
            TaskType::Regression,
            FeatureExtractionStrategy::ExpectationValues,
            0.01,
        );

        let inputs = vec![
            vec![1.0, 0.0],
            vec![0.0, 1.0],
            vec![0.5, 0.5],
        ];
        let outputs = vec![vec![1.0], vec![0.0], vec![0.5]];
        let data = TrainingData::new(inputs.clone(), outputs.clone());

        trainer.train(&mut reservoir, &data).unwrap();

        let test_data = TrainingData::new(inputs, outputs);
        let result = trainer.evaluate(&mut reservoir, &test_data);

        assert!(result.is_ok());
        let (mse, _mae) = result.unwrap();
        assert!(mse >= 0.0);
    }

    #[test]
    fn test_classification_task() {
        let mut reservoir = create_test_reservoir();
        let mut trainer = Trainer::new(
            TaskType::Classification,
            FeatureExtractionStrategy::Probabilities,
            0.01,
        );

        let inputs = vec![
            vec![1.0, 0.0],
            vec![0.0, 1.0],
            vec![0.9, 0.1],
            vec![0.1, 0.9],
        ];
        let outputs = vec![
            vec![1.0, 0.0],
            vec![0.0, 1.0],
            vec![1.0, 0.0],
            vec![0.0, 1.0],
        ];
        let data = TrainingData::new(inputs, outputs);

        let result = trainer.train(&mut reservoir, &data);
        assert!(result.is_ok());
    }

    #[test]
    fn test_combined_feature_extraction() {
        let mut reservoir = create_test_reservoir();
        let trainer = Trainer::new(
            TaskType::Regression,
            FeatureExtractionStrategy::Combined,
            0.01,
        );

        reservoir.encode_input(&[0.6, 0.8]);
        reservoir.evolve();

        let features = trainer.extract_features(&reservoir);
        assert!(features.len() > 2);
    }
}