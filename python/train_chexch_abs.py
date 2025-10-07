import pandas as pd
import xgboost as xgb
import numpy as np
import matplotlib.pyplot as plt

from sklearn.model_selection import train_test_split
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay
from sklearn.metrics import roc_curve, auc
from sklearn.preprocessing import label_binarize
from sklearn.utils.class_weight import compute_class_weight

import convert

def plot_classification_error(history, metric):
    epochs = len(history["validation_0"][metric])
    x_axis = range(epochs)

    plt.figure(figsize=(8, 6))
    plt.plot(x_axis, history["validation_0"][metric], label="Train")
    plt.plot(x_axis, history["validation_1"][metric], label="Test")
    plt.legend()
    plt.ylabel(f"Classification Error ({metric})")
    plt.xlabel("Boosting Round")
    plt.title("XGBoost Classification Error")
    plt.savefig("figs/chexch_abs/classification_error.png", dpi=300, bbox_inches="tight")
    plt.close()

def plot_confusion_matrix(model, X_test, y_test, labels):
    y_pred = model.predict(X_test)
    cm = confusion_matrix(y_test, y_pred)
    disp = ConfusionMatrixDisplay(confusion_matrix=cm, display_labels=labels)
    disp.plot(cmap="Blues")
    plt.title("Confusion Matrix")
    plt.savefig("figs/chexch_abs/confusion_matrix.png", dpi=300, bbox_inches="tight")
    plt.close()

def plot_feature_importance(model):
    xgb.plot_importance(model, importance_type="gain")
    plt.title("Feature Importance by Gain")
    plt.savefig("figs/chexch_abs/feature_importance.png", dpi=300, bbox_inches="tight")
    plt.close()

def plot_roc_curves(model, X_test, y_test, classes):
    """
    Saves ROC curves to 'roc_curves.png'.
    Works for binary (len(classes)==2) and multiclass (len(classes)>2).
    """
    proba = model.predict_proba(X_test)
    class_order = list(model.classes_)  # column order of predict_proba

    if len(classes) == 2:
        # Binary case: plot a single ROC for the positive class (take classes[1] as positive by convention)
        pos_class = classes[1]
        pos_idx = class_order.index(pos_class)

        y_true_bin = (y_test == pos_class).astype(int)
        y_score = proba[:, pos_idx]

        fpr, tpr, _ = roc_curve(y_true_bin, y_score)
        roc_auc = auc(fpr, tpr)

        plt.figure(figsize=(8, 6))
        plt.plot(fpr, tpr, label=f"Class {pos_class} vs. rest (AUC = {roc_auc:.2f})")
        plt.plot([0, 1], [0, 1], "k--")
        plt.xlabel("False Positive Rate")
        plt.ylabel("True Positive Rate")
        plt.title("ROC Curve (Binary)")
        plt.legend()
        plt.savefig("figs/chexch_abs/roc_curves.png", dpi=300, bbox_inches="tight")
        plt.close()

    else:
        # Multiclass one-vs-rest
        y_test_bin = label_binarize(y_test, classes=classes)

        plt.figure(figsize=(8, 6))
        for i, class_id in enumerate(classes):
            # Align probability column with the actual column index for this class
            col_idx = class_order.index(class_id)
            fpr, tpr, _ = roc_curve(y_test_bin[:, i], proba[:, col_idx])
            roc_auc = auc(fpr, tpr)
            plt.plot(fpr, tpr, label=f"Class {class_id} (AUC = {roc_auc:.2f})")

        plt.plot([0, 1], [0, 1], "k--")
        plt.xlabel("False Positive Rate")
        plt.ylabel("True Positive Rate")
        plt.title("ROC Curves by Class")
        plt.legend()
        plt.savefig("figs/chexch_abs/roc_curves.png", dpi=300, bbox_inches="tight")
        plt.close()

if (__name__ == "__main__"):
    # Load DataFrame from pickle file
    df = pd.read_pickle("files/train_chexch_abs_data.pkl")

    # Keep only pion abs 0p, charge exchange, electron showers
    df = df[df["backgroundType"].isin([0, 3, 7])]

    # Map backgroundType to new target: 0 for pion abs 0p, 1 for electron shower, 2 for charge exchange
    df["target"] = df["backgroundType"].apply(lambda x: 2 if x == 7 else (1 if x == 3 else 0))
    
    # Select features and target
    X = df.drop(columns=["backgroundType", "target", "event"])
    y = df["target"]

    # Split into train and test sets
    print("Splitting data into train and test sets")
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.15, random_state=42)

    # Compute class weights
    classes = np.array([0,1,2])
    cw = compute_class_weight("balanced", classes=classes, y=y_train)
    class_to_w = {c:w for c, w in zip(classes, cw)}

    w_train = y_train.map(class_to_w).astype(float).values
    w_test  = y_test.map(class_to_w).astype(float).values
    
    # Three classes
    model = xgb.XGBClassifier(
        num_class=3,
        n_estimators=200,
        objective='multi:softprob',
        eval_metric=['mlogloss', 'merror'],
    )

    print("Starting training")
    model.fit(
        X_train,
        y_train,
        eval_set=[(X_train, y_train), (X_test, y_test)],
        sample_weight_eval_set=[w_train, w_test],
        verbose=20
    )
    evals_result = model.evals_result() 

    # Convert to xml
    variables = []
    dtype_map = {
        'int64': 'I',
        'float64': 'F',
        'float32': 'F',
        'int32': 'I'
    }
    variables = [(col, dtype_map.get(str(dtype), 'F')) for col, dtype in X.dtypes.items()]

    print("Converting model to XML format")
    
    dump = model.get_booster().get_dump()
    per_class_trees = [[] for _ in range(3)]
    for idx, tree in enumerate(dump):
        per_class_trees[idx % 3].append(tree)

    for ci, trees in enumerate(per_class_trees):
        convert.convert_model(
            trees,
            input_variables=variables,
            output_xml=f"model/chexch_abs_model_class_{ci}.xml"
        )

    # Evaluate the model
    accuracy = model.score(X_test, y_test)
    print(f"Test accuracy: {accuracy:.4f}")

    # Test events
    events = df["event"].iloc[:10].values
    test_probs = model.predict_proba(X.iloc[:10])
    print("Test event probabilities:")
    for i, probs in enumerate(test_probs):
        probs_text = " ".join([f"{p:.5f}" for p in probs])
        print(f"Event {events[i]}: Class probabilities: {probs_text}")
        # for j, (col, _) in enumerate(variables):
        #     value = X.iloc[i][col]
        #     print(f"    {col}: {value}")

    plot_classification_error(evals_result, "mlogloss")
    plot_confusion_matrix(model, X_test, y_test, ["abs 0p", "ch. exch.", "electron"])
    plot_feature_importance(model)
    plot_roc_curves(model, X_test, y_test, classes=[0, 1, 2])
