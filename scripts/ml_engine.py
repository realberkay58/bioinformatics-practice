from sklearn.ensemble import RandomForestClassifier
import numpy as np

# TEST FUNCTION: Simulates training data and predicts a new patient's status
@st.cache_data
def run_diagnostic_test(df_length):
    # Simulate training data: 6 samples (3 control, 3 tumor), features = number of genes
    X = np.random.rand(6, df_length) 
    y = np.array([0, 0, 0, 1, 1, 1]) # 0 = Healthy/Control, 1 = Tumor
    
    # Initialize and train the Random Forest Classifier
    clf = RandomForestClassifier(n_estimators=100, random_state=42)
    clf.fit(X, y)
    
    # Simulate a new unseen patient's gene expression data
    new_patient = np.random.rand(1, df_length)
    
    # Make prediction and calculate probability
    prediction = clf.predict(new_patient)
    probability = clf.predict_proba(new_patient)
    
    return prediction[0], np.max(probability)