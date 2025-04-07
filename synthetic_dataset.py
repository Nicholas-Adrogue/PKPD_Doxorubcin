import pandas as pd
import numpy as np
from sklearn.impute import KNNImputer

# Load original dataset
df = pd.read_csv("clinical_pharmacological_parameters.csv")

# KNN Imputation
imputer = KNNImputer(n_neighbors=5)
df_imputed = pd.DataFrame(imputer.fit_transform(df), columns=df.columns)

# Synthetic data generation function
def generate_synthetic_data(df, num_samples=100):
    synthetic_data = []
    num_existing_samples = len(df)

    for i in range(num_samples):
        idx1, idx2 = np.random.choice(num_existing_samples, 2, replace=False)
        patient1, patient2 = df.iloc[idx1], df.iloc[idx2]

        synthetic_patient = {}
        synthetic_patient["Patient number"] = f"P{1000 + i}"  # Assign new string ID like "P1000"

        for col in df.columns:
            if col == "Patient number":
                continue  # Skip interpolating this column

            if df[col].dtype in [np.float64, np.int64]:
                alpha = np.random.uniform(0.4, 0.6)
                noise_scale = max(0.05 * abs(patient1[col] - patient2[col]), 1e-6)
                noise = np.random.normal(0, noise_scale)
                synthetic_patient[col] = alpha * patient1[col] + (1 - alpha) * patient2[col] + noise
            else:
                synthetic_patient[col] = np.random.choice([patient1[col], patient2[col]])

        synthetic_data.append(synthetic_patient)

    return pd.DataFrame(synthetic_data)

# Generate and save
synthetic_df = generate_synthetic_data(df_imputed, num_samples=100)
synthetic_df.to_csv("synthetic_patient_data.csv", index=False)

print("✅ Synthetic dataset saved as 'synthetic_patient_data.csv'")


