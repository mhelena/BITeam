#!/usr/bin/env python3
# Random forest on 1KGP3 PCs assigns a super-population to each study sample (max prob >= threshold),
# with a kNN fallback for samples below it. Usage: python3 05_ancestry_classifier.py --pcs X --labels Y --out-prefix Z
import argparse
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score
from sklearn.neighbors import NearestNeighbors
import joblib
import json
import sys


SUPERPOPS = ['AFR', 'AMR', 'EAS', 'EUR', 'SAS']


def load_pcs(pcs_file: str, n_pcs: int = 8) -> pd.DataFrame:
    df = pd.read_csv(pcs_file, sep='\t')
    pc_cols = [f'PC{i}' for i in range(1, n_pcs + 1)]
    missing = [c for c in pc_cols if c not in df.columns]
    if missing:
        raise ValueError(f"Missing PC columns: {missing}. Available: {list(df.columns)}")
    return df, pc_cols


def load_labels(labels_file: str) -> pd.DataFrame:
    df = pd.read_csv(labels_file, sep='\t')
    col_map = {}
    for c in df.columns:
        cl = c.lower().replace(' ', '_').replace('-', '_')
        if cl in ('iid', 'sample', 'sample_id', 'id'):
            col_map[c] = 'IID'
        elif 'super' in cl and 'pop' in cl:
            col_map[c] = 'SuperPop'
    df.rename(columns=col_map, inplace=True)
    if 'IID' not in df.columns or 'SuperPop' not in df.columns:
        raise ValueError(
            f"Labels file must have IID and SuperPop columns. Found: {list(df.columns)}"
        )
    df = df[df['SuperPop'].isin(SUPERPOPS)].copy()
    return df


def main():
    parser = argparse.ArgumentParser(description='Random forest ancestry classifier')
    parser.add_argument('--pcs',        required=True, help='All PC scores TSV (reference + study)')
    parser.add_argument('--labels',     required=True, help='1KGP3 super-population labels TSV')
    parser.add_argument('--out-prefix', required=True, help='Output prefix')
    parser.add_argument('--n-pcs',      type=int, default=8,   help='Number of PCs to use (default: 8)')
    parser.add_argument('--n-trees',    type=int, default=400,  help='Number of RF trees (default: 400)')
    parser.add_argument('--threshold',  type=float, default=0.8, help='Min probability for assignment (default: 0.8)')
    parser.add_argument('--random-state', type=int, default=42)
    parser.add_argument('--save-model', action='store_true', help='Save trained RF model')
    parser.add_argument('--knn-k',        type=int,   default=50,  help='k for kNN fallback (default: 50)')
    parser.add_argument('--knn-threshold', type=float, default=0.7, help='Min vote fraction for kNN assignment (default: 0.7)')
    args = parser.parse_args()

    print(f"PCs: {args.pcs}")
    print(f"Labels: {args.labels}")
    print(f"Output: {args.out_prefix}")
    print(f"n_pcs={args.n_pcs} n_trees={args.n_trees} threshold={args.threshold} "
          f"knn_k={args.knn_k} knn_threshold={args.knn_threshold}")

    all_pcs_df, pc_cols = load_pcs(args.pcs, args.n_pcs)
    labels_df = load_labels(args.labels)

    ref_pcs  = all_pcs_df[all_pcs_df['source'] == 'reference_1KGP3'].copy()
    study_pcs = all_pcs_df[all_pcs_df['source'] == 'study'].copy()
    print(f"Reference samples: {len(ref_pcs)}, study samples: {len(study_pcs)}")

    ref_labeled = pd.merge(ref_pcs, labels_df[['IID','SuperPop']], on='IID', how='inner')
    print(f"Reference samples with labels: {len(ref_labeled)}")
    print(ref_labeled['SuperPop'].value_counts().to_string())

    if len(ref_labeled) == 0:
        print("ERROR: No reference samples matched with population labels.")
        sys.exit(1)

    X_train = ref_labeled[pc_cols].values
    y_train = ref_labeled['SuperPop'].values

    print(f"Training random forest ({args.n_trees} trees, {args.n_pcs} PCs)")
    rf = RandomForestClassifier(
        n_estimators=args.n_trees,
        random_state=args.random_state,
        n_jobs=-1,
        oob_score=True,
    )
    rf.fit(X_train, y_train)

    oob_error = 1.0 - rf.oob_score_
    print(f"OOB score: {rf.oob_score_:.6f}")
    print(f"OOB error: {oob_error:.6f} (GEL reported 0.0024)")

    # n_jobs=1: running the CV folds in parallel ran out of memory, since each
    # fold already fits a forest with n_jobs=-1.
    cv_scores = cross_val_score(rf, X_train, y_train, cv=5, n_jobs=1)
    print(f"5-fold CV accuracy: {cv_scores.mean():.4f} +/- {cv_scores.std():.4f}")

    ref_proba  = rf.predict_proba(X_train)
    ref_proba_df = pd.DataFrame(ref_proba, columns=rf.classes_)
    ref_proba_df.insert(0, 'IID', ref_labeled['IID'].values)
    ref_proba_df.insert(1, 'TrueLabel', y_train)
    ref_proba_df['MaxProb'] = ref_proba_df[list(rf.classes_)].max(axis=1)
    ref_proba_df['PredLabel'] = ref_proba_df[list(rf.classes_)].idxmax(axis=1)
    ref_proba_df['Assigned'] = ref_proba_df['MaxProb'] >= args.threshold
    ref_accuracy = (ref_proba_df['PredLabel'] == ref_proba_df['TrueLabel']).mean()
    print(f"Reference prediction accuracy: {ref_accuracy:.4f}")

    if len(study_pcs) == 0:
        print("WARNING: No study samples in PC file, check the 'source' column.")
    else:
        X_study = study_pcs[pc_cols].values
        study_proba = rf.predict_proba(X_study)
        study_proba_df = pd.DataFrame(study_proba, columns=rf.classes_)
        study_proba_df.insert(0, 'FID', study_pcs['FID'].values)
        study_proba_df.insert(1, 'IID', study_pcs['IID'].values)

        for col in pc_cols:
            study_proba_df[col] = study_pcs[col].values

        study_proba_df['MaxProb'] = study_proba_df[list(rf.classes_)].max(axis=1)
        study_proba_df['PredLabel'] = study_proba_df[list(rf.classes_)].idxmax(axis=1)
        study_proba_df['AncestryGEL'] = study_proba_df.apply(
            lambda r: r['PredLabel'] if r['MaxProb'] >= args.threshold else 'Unassigned',
            axis=1
        )

        # Study samples pushed off the reference clusters by batch effects fall
        # outside the RF decision regions; a vote among the k nearest reference
        # samples is more robust for those.
        unassigned_mask = study_proba_df['AncestryGEL'] == 'Unassigned'
        if unassigned_mask.any() and args.knn_k > 0:
            knn = NearestNeighbors(n_neighbors=args.knn_k, algorithm='ball_tree')
            knn.fit(X_train)
            X_unassigned = study_pcs.loc[study_pcs.index[unassigned_mask], pc_cols].values
            _, nn_idx = knn.kneighbors(X_unassigned)
            knn_labels = y_train[nn_idx]

            knn_assignments = []
            knn_votes = []
            for row_labels in knn_labels:
                counts = pd.Series(row_labels).value_counts()
                top_pop = counts.index[0]
                top_frac = counts.iloc[0] / args.knn_k
                knn_assignments.append(top_pop if top_frac >= args.knn_threshold else 'Unassigned')
                knn_votes.append(top_frac)

            study_proba_df.loc[unassigned_mask, 'AncestryGEL'] = knn_assignments
            study_proba_df['ClassificationMethod'] = 'RF'
            study_proba_df.loc[unassigned_mask, 'ClassificationMethod'] = [
                'kNN' if a != 'Unassigned' else 'Unassigned'
                for a in knn_assignments
            ]
            n_rescued = sum(a != 'Unassigned' for a in knn_assignments)
            print(f"kNN fallback (k={args.knn_k}, threshold {args.knn_threshold}): "
                  f"{n_rescued} of {unassigned_mask.sum()} RF-unassigned samples assigned")
            for i, (row_idx, _) in enumerate(study_proba_df[unassigned_mask].iterrows()):
                iid = study_proba_df.loc[row_idx, 'IID']
                print(f"  {iid}: {knn_assignments[i]} ({knn_votes[i]*100:.1f}% of {args.knn_k} neighbours)")
        else:
            study_proba_df['ClassificationMethod'] = 'RF'

        print(f"Study samples: {len(study_proba_df)}")
        print(study_proba_df['AncestryGEL'].value_counts().to_string())
        n_unassigned = (study_proba_df['AncestryGEL'] == 'Unassigned').sum()
        pct_unassigned = n_unassigned / len(study_proba_df) * 100
        print(f"Unassigned: {n_unassigned} ({pct_unassigned:.1f}%; GEL reported 7.5%)")

        out_file = f"{args.out_prefix}_ancestry_assignments.tsv"
        study_proba_df.to_csv(out_file, sep='\t', index=False)
        print(f"Written: {out_file}")

        summary = {
            'n_study_samples': int(len(study_proba_df)),
            'n_trees': args.n_trees,
            'n_pcs_used': args.n_pcs,
            'probability_threshold': args.threshold,
            'knn_k': args.knn_k,
            'knn_threshold': args.knn_threshold,
            'oob_score': float(rf.oob_score_),
            'oob_error': float(oob_error),
            'cv_accuracy_mean': float(cv_scores.mean()),
            'cv_accuracy_std': float(cv_scores.std()),
            'ancestry_counts': study_proba_df['AncestryGEL'].value_counts().to_dict(),
            'pct_unassigned': float(pct_unassigned),
            'classification_method_counts': study_proba_df['ClassificationMethod'].value_counts().to_dict(),
        }
        with open(f"{args.out_prefix}_summary.json", 'w') as fh:
            json.dump(summary, fh, indent=2)

    if args.save_model:
        model_file = f"{args.out_prefix}_rf_model.joblib"
        joblib.dump(rf, model_file)
        print(f"Model written: {model_file}")

    ref_proba_df.to_csv(f"{args.out_prefix}_reference_predictions.tsv", sep='\t', index=False)


if __name__ == '__main__':
    main()
