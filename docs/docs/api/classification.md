# Classification API

Functions for trajectory classification and survival analysis.

## Overview

The classification module provides:

- Trajectory-based classification
- Survival prediction
- Classifier training and evaluation
- Feature importance analysis

## Main Classification Function

### classify_trajectories

Train classifier to predict disease progression from trajectories.

::: renalprog.modeling.predict.classify_trajectories

**Example Usage:**

```python
from renalprog.modeling.predict import classify_trajectories
import pandas as pd
from pathlib import Path

# Load trajectory data
trajectories = pd.read_csv("data/processed/trajectories.csv")
labels = pd.read_csv("data/processed/progression_labels.csv")

# Train classifier
classifier, metrics = classify_trajectories(
    trajectories=trajectories.values,
    labels=labels['progressed'].values,
    output_dir=Path("models/trajectory_classifier"),
    model_type='random_forest',  # or 'logistic', 'svm', 'gradient_boosting'
    test_size=0.2,
    random_state=42
)

# Print performance
print(f"Accuracy: {metrics['accuracy']:.3f}")
print(f"AUC-ROC: {metrics['auc_roc']:.3f}")
print(f"Precision: {metrics['precision']:.3f}")
print(f"Recall: {metrics['recall']:.3f}")
print(f"F1 Score: {metrics['f1']:.3f}")
```

## Classification Models

The `classify_trajectories` function supports multiple model types:

### Random Forest

Default choice for interpretability and feature importance:

```python
classifier, metrics = classify_trajectories(
    trajectories, labels,
    model_type='random_forest',
    n_estimators=100,
    max_depth=None,
    min_samples_split=2
)

# Feature importance
importances = classifier.feature_importances_
```

### Logistic Regression

For linear decision boundaries:

```python
classifier, metrics = classify_trajectories(
    trajectories, labels,
    model_type='logistic',
    C=1.0,
    penalty='l2',
    max_iter=1000
)
```

### Support Vector Machine

For complex decision boundaries:

```python
classifier, metrics = classify_trajectories(
    trajectories, labels,
    model_type='svm',
    kernel='rbf',
    C=1.0,
    gamma='scale'
)
```

### Gradient Boosting

For maximum performance:

```python
classifier, metrics = classify_trajectories(
    trajectories, labels,
    model_type='gradient_boosting',
    n_estimators=100,
    learning_rate=0.1,
    max_depth=3
)
```

## Evaluation Metrics

The classification function returns comprehensive metrics:

| Metric | Description |
|--------|-------------|
| `accuracy` | Overall classification accuracy |
| `auc_roc` | Area under ROC curve |
| `precision` | Positive predictive value |
| `recall` | Sensitivity/True positive rate |
| `f1` | Harmonic mean of precision and recall |
| `confusion_matrix` | 2×2 confusion matrix |
| `classification_report` | Detailed per-class metrics |

## Visualization

### plot_confusion_matrix

Visualize classifier performance:

::: renalprog.plots.plot_confusion_matrix

**Example:**

```python
from renalprog.plots import plot_confusion_matrix
from pathlib import Path

# Plot confusion matrix
plot_confusion_matrix(
    confusion_matrix=metrics['confusion_matrix'],
    class_names=['Non-progressing', 'Progressing'],
    output_path=Path("reports/figures/confusion_matrix.png"),
    title="Trajectory Classification"
)
```

## Complete Classification Workflow

```python
import torch
import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.preprocessing import StandardScaler
from renalprog.modeling.train import VAE
from renalprog.modeling.predict import (
    apply_vae,
    generate_trajectories,
    classify_trajectories
)
from renalprog.plots import plot_confusion_matrix

# 1. Load model and data
model = VAE(input_dim=20000, mid_dim=1024, features=128)
model.load_state_dict(torch.load("models/my_vae/best_model.pt"))

expr = pd.read_csv("data/interim/split/test_expression.tsv", sep="\t", index_col=0)
clinical = pd.read_csv("data/interim/split/test_clinical.tsv", sep="\t", index_col=0)

# 2. Generate trajectories for each patient
early_mask = clinical['stage'] == 'early'
late_mask = clinical['stage'] == 'late'

trajectories = generate_trajectories(
    model=model,
    start_data=expr.values[early_mask],
    end_data=expr.values[late_mask],
    n_steps=50,
    interpolation='spherical',
    device='cuda'
)

# 3. Extract trajectory features
# Option A: Use trajectory statistics (mean, std, slope)
trajectory_features = []
for traj in trajectories:
    features = np.concatenate([
        traj.mean(axis=0),  # Mean expression
        traj.std(axis=0),   # Variance
        (traj[-1] - traj[0])  # Net change
    ])
    trajectory_features.append(features)
trajectory_features = np.array(trajectory_features)

# Option B: Use latent trajectory
results = apply_vae(model, expr.values, device='cuda')
latent_trajectories = results['latent']

# 4. Create labels (e.g., based on survival)
labels = clinical['progressed'].values[early_mask]

# 5. Train classifier
classifier, metrics = classify_trajectories(
    trajectories=trajectory_features,
    labels=labels,
    output_dir=Path("models/trajectory_classifier"),
    model_type='random_forest',
    test_size=0.2,
    random_state=42
)

# 6. Visualize results
plot_confusion_matrix(
    confusion_matrix=metrics['confusion_matrix'],
    class_names=['Stable', 'Progressed'],
    output_path=Path("reports/figures/classification_cm.png")
)

# 7. Feature importance (for tree-based models)
if hasattr(classifier, 'feature_importances_'):
    importances = pd.DataFrame({
        'feature': [f'feature_{i}' for i in range(len(classifier.feature_importances_))],
        'importance': classifier.feature_importances_
    }).sort_values('importance', ascending=False)
    
    print("Top 10 important features:")
    print(importances.head(10))

# 8. Save results
results_df = pd.DataFrame({
    'patient_id': clinical.index[early_mask],
    'true_label': labels,
    'predicted_label': classifier.predict(trajectory_features),
    'prediction_proba': classifier.predict_proba(trajectory_features)[:, 1]
})
results_df.to_csv("reports/classification_results.csv", index=False)

print("\nClassification Performance:")
print(f"Accuracy: {metrics['accuracy']:.3f}")
print(f"AUC-ROC: {metrics['auc_roc']:.3f}")
print(f"F1 Score: {metrics['f1']:.3f}")
```

## Cross-Validation

For robust performance estimation:

```python
from sklearn.model_selection import cross_val_score
from sklearn.ensemble import RandomForestClassifier

# Create classifier
rf = RandomForestClassifier(n_estimators=100, random_state=42)

# Cross-validation
cv_scores = cross_val_score(
    rf, trajectory_features, labels,
    cv=5,
    scoring='roc_auc'
)

print(f"CV AUC-ROC: {cv_scores.mean():.3f} ± {cv_scores.std():.3f}")
```

## Hyperparameter Tuning

```python
from sklearn.model_selection import GridSearchCV
from sklearn.ensemble import RandomForestClassifier

# Define parameter grid
param_grid = {
    'n_estimators': [50, 100, 200],
    'max_depth': [None, 10, 20, 30],
    'min_samples_split': [2, 5, 10],
    'min_samples_leaf': [1, 2, 4]
}

# Grid search
rf = RandomForestClassifier(random_state=42)
grid_search = GridSearchCV(
    rf, param_grid,
    cv=5,
    scoring='roc_auc',
    n_jobs=-1,
    verbose=1
)

grid_search.fit(trajectory_features, labels)

print(f"Best parameters: {grid_search.best_params_}")
print(f"Best CV AUC-ROC: {grid_search.best_score_:.3f}")
```

## See Also

- [Trajectories API](trajectories.md) - Generate trajectories
- [Prediction API](prediction.md) - Apply VAE models
- [Plots API](plots.md) - Visualization tools
- [Complete Pipeline Tutorial](../tutorials/complete-pipeline.md)

