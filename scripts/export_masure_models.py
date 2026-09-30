"""Export the masure-regression scikit-learn models to the .npz files AirfoilAero reads.

    python export_masure_models.py ET_re1e6.pkl ET_re5e6.pkl ET_re2e7.pkl --out DIR
    python export_masure_models.py --fixture DIR

The first form converts the pickles from https://doi.org/10.5281/zenodo.16925758 into
`DIR/ET_re<Re>.npz`; pass `DIR` as `ml_models_dir` to `resolve_aero_geometry`. The
second trains a small model of the same layout on synthetic data and writes it with
the inputs and `predict` outputs the Julia test compares against. Needs numpy and
scikit-learn.

An exported file holds `input_mean` and `input_scale` (the StandardScaler) and, for
output k = 0, 1, 2 (CD, CL, CM), the nodes of all its trees concatenated:
`output{k}_roots`, `_left`, `_right`, `_feature`, `_threshold`, `_value`. Node and
feature indices are 0-based, and `_left` is -1 on a leaf.
"""

import argparse
import pickle
from pathlib import Path

import numpy as np
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.multioutput import MultiOutputRegressor
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


def forest_arrays(forest):
    """Concatenate the nodes of every tree in `forest`, children indexed globally."""
    roots, left, right, feature, threshold, value = [], [], [], [], [], []
    offset = 0
    for estimator in forest.estimators_:
        tree = estimator.tree_
        if tree.n_outputs != 1:
            raise ValueError(f"expected single-output trees, got {tree.n_outputs}")
        is_leaf = tree.children_left == -1
        roots.append(offset)
        left.append(np.where(is_leaf, -1, tree.children_left + offset))
        right.append(np.where(is_leaf, -1, tree.children_right + offset))
        feature.append(tree.feature)
        threshold.append(tree.threshold)
        value.append(tree.value[:, 0, 0])
        offset += tree.node_count
    return {
        "roots": np.asarray(roots, dtype=np.int32),
        "left": np.concatenate(left).astype(np.int32),
        "right": np.concatenate(right).astype(np.int32),
        "feature": np.concatenate(feature).astype(np.int32),
        "threshold": np.concatenate(threshold).astype(np.float64),
        "value": np.concatenate(value).astype(np.float64),
    }


def export_model(model, npz_path):
    """Write a StandardScaler -> MultiOutputRegressor(ExtraTreesRegressor) pipeline."""
    steps = [step for _, step in model.steps] if isinstance(model, Pipeline) else []
    if not (len(steps) == 2 and isinstance(steps[0], StandardScaler)
            and isinstance(steps[1], MultiOutputRegressor)):
        raise ValueError(f"unexpected model layout: {model!r}")
    scaler, multi_output = steps
    if len(multi_output.estimators_) != 3:
        raise ValueError(f"expected 3 outputs, got {len(multi_output.estimators_)}")
    arrays = {"input_mean": scaler.mean_, "input_scale": scaler.scale_}
    for k, forest in enumerate(multi_output.estimators_):
        if not isinstance(forest, ExtraTreesRegressor):
            raise ValueError(f"output {k} is a {type(forest).__name__}")
        for name, array in forest_arrays(forest).items():
            arrays[f"output{k}_{name}"] = array
    np.savez(npz_path, **arrays)


def threshold_rows(model, row):
    """Copies of `row` moved onto each output's first root split, where float32 matters."""
    scaler, multi_output = model.named_steps["scale"], model.named_steps["model"]
    rows = []
    for forest in multi_output.estimators_:
        tree = forest.estimators_[0].tree_
        feature, threshold = tree.feature[0], tree.threshold[0]
        for step in (-1e-9, 0.0, 1e-9):
            moved = row.copy()
            moved[feature] = threshold * (1 + step) * scaler.scale_[feature] + \
                scaler.mean_[feature]
            rows.append(moved)
    return np.array(rows)


def write_fixture(out_dir):
    """Train a small pipeline of the published layout and write it with its predictions."""
    rng = np.random.default_rng(42)
    low = np.array([0.05, 0.1, 0.0, -10.0, 0.1, 0.1, -10.0])
    high = np.array([0.12, 0.6, 0.15, 5.0, 0.4, 0.9, 30.0])
    X = rng.uniform(low, high, size=(400, 7))
    alpha = np.deg2rad(X[:, 6])
    y = np.column_stack([0.02 + 0.3 * alpha**2 + X[:, 0],
                         2 * np.pi * alpha + 5 * X[:, 2],
                         -0.1 - X[:, 2] + 0.01 * X[:, 3]])
    model = Pipeline([
        ("scale", StandardScaler()),
        ("model", MultiOutputRegressor(ExtraTreesRegressor(
            n_estimators=5, max_depth=6, max_features="log2", random_state=42))),
    ]).fit(X, y)
    out_dir.mkdir(parents=True, exist_ok=True)
    export_model(model, out_dir / "ET_re1e6.npz")
    X_test = np.vstack([rng.uniform(low, high, size=(20, 7)), threshold_rows(model, X[0])])
    np.savez(out_dir / "reference.npz", X=X_test, Y=model.predict(X_test))


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("pickles", nargs="*", type=Path)
    parser.add_argument("--out", type=Path, default=Path("."))
    parser.add_argument("--fixture", type=Path)
    args = parser.parse_args()
    if args.fixture is not None:
        write_fixture(args.fixture)
    for pickle_path in args.pickles:
        with open(pickle_path, "rb") as file:
            model = pickle.load(file)
        args.out.mkdir(parents=True, exist_ok=True)
        export_model(model, args.out / pickle_path.with_suffix(".npz").name)


if __name__ == "__main__":
    main()
