import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# === Load synthetic patient data ===
df = pd.read_csv("synthetic_patient_data.csv")

# === AUC Calculations ===
# Prevent division by zero or missing values
df["AUC_DOX (ng·min/ml)"] = df["Maximum doxorubicin concentration (ng/ml)"] / df["Mean clearance (ml/min/m^2)"]
df["AUC_DOXOL (ng·min/ml)"] = df["Maximum doxorubicinol concentration )ng/ml)"] / df["Mean clearance (ml/min/m^2)"]

# Flag high AUC exposure (custom thresholds)
df["High AUC DOX"] = (df["AUC_DOX (ng·min/ml)"] > 0.04).astype(int)
df["High AUC DOXOL"] = (df["AUC_DOXOL (ng·min/ml)"] > 0.015).astype(int)

# === Toxicity Classification ===
def classify_toxicity(row):
    cardio = (row["Maximum doxorubicinol concentration )ng/ml)"] > 5.0) or (row["Mean clearance (ml/min/m^2)"] < 400)
    hema = row["Nadir WBC  (x10^3 μl"] < 2.0
    stomatitis = row["Maximum Stomatitis (grade)"] >= 2
    auc_dox = row["High AUC DOX"]
    auc_doxol = row["High AUC DOXOL"]
    any_tox = cardio or hema or stomatitis or auc_dox or auc_doxol

    return pd.Series({
        "Cardiotoxicity Risk": int(cardio),
        "Hematological Toxicity Risk": int(hema),
        "Severe Stomatitis": int(stomatitis),
        "Any Toxicity": int(any_tox)
    })

tox_df = df.apply(classify_toxicity, axis=1)
final_df = pd.concat([df, tox_df], axis=1)

# === Save final data ===
final_df.to_csv("synthetic_toxicity_data.csv", index=False)
print("✅ Data saved as 'synthetic_toxicity_data.csv'")

# === Visualization 1: Toxicity Bar Plot ===
plt.figure(figsize=(8, 5))
tox_df.sum().plot(kind='bar', color='salmon')
plt.title("Number of Patients by Toxicity Type")
plt.ylabel("Patient Count")
plt.xticks(rotation=15)
plt.tight_layout()
plt.savefig("toxicity_barplot.png")
plt.show()

# === Visualization 2: Any Toxicity Pie Chart ===
labels = ['No Toxicity', 'Any Toxicity']
sizes = [len(tox_df) - tox_df['Any Toxicity'].sum(), tox_df['Any Toxicity'].sum()]
colors = ['#66b3ff', '#ff6666']
plt.figure(figsize=(6,6))
plt.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=140, colors=colors)
plt.title("Overall Toxicity Distribution")
plt.savefig("toxicity_piechart.png")
plt.show()

# === Visualization 3: AUC Distribution Plots ===
plt.figure(figsize=(8,5))
sns.histplot(df["AUC_DOX (ng·min/ml)"], kde=True, color='purple', label='DOX')
plt.axvline(0.04, color='red', linestyle='--', label='DOX Toxicity Threshold')
plt.title("AUC Distribution - Doxorubicin")
plt.xlabel("AUC (ng·min/ml)")
plt.legend()
plt.tight_layout()
plt.savefig("auc_dox_distribution.png")
plt.show()

plt.figure(figsize=(8,5))
sns.histplot(df["AUC_DOXOL (ng·min/ml)"], kde=True, color='green', label='DOXOL')
plt.axvline(0.015, color='red', linestyle='--', label='DOXOL Toxicity Threshold')
plt.title("AUC Distribution - Doxorubicinol")
plt.xlabel("AUC (ng·min/ml)")
plt.legend()
plt.tight_layout()
plt.savefig("auc_doxol_distribution.png")
plt.show()

# === Visualization 4: Correlation Heatmap ===
plt.figure(figsize=(12, 8))
sns.heatmap(final_df.drop(columns=["Patient number"]).corr(), annot=True, fmt=".2f", cmap="coolwarm")
plt.title("Correlation Heatmap: PK & Toxicity Variables")
plt.tight_layout()
plt.savefig("toxicity_correlation_heatmap.png")
plt.show()
