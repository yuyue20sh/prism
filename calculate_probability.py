from pathlib import Path
from itertools import combinations

import pandas as pd


def calculate_probability(enrichment_scores: pd.DataFrame, cluster_names: list, predicted_cluster: str, omics_names: list, omics_weights: list) -> pd.DataFrame:
    """Calculate cluster prediction probability across omics combinations.

    Args:
        enrichment_scores (pd.DataFrame):
        cluster_names (list): Cluster names.
        predicted_cluster (str): Predicted cluster.
        omics_names (list): Omics names.
        omics_weights (list): Weight for omics.

    Returns:
        pd.DataFrame: Predicted probability.

    """
    other_clusters = [i for i in cluster_names if i != predicted_cluster]
    weights = {k: v for k, v in zip(omics_names, omics_weights)}

    # probability of each omics
    cluster_probability = pd.DataFrame(index=enrichment_scores.index.tolist(), columns=omics_names)
    for omics, weight in zip(omics_names, omics_weights):
        denominator = enrichment_scores[['%s_%s' % (omics, cluster) for cluster in other_clusters]].max(axis=1)
        numerator = enrichment_scores['%s_%s' % (omics, predicted_cluster)]
        probability = numerator / denominator
        cluster_probability[omics] = probability.loc[cluster_probability.index]

    omics_groups = []
    for i in range(1, len(omics_names) + 1):
        omics_groups.extend(combinations(omics_names, i))
    for group in omics_groups:
        group_name = 'combination_%s' % ('_'.join(group))
        cluster_probability[group_name] = 0
        for name in group:
            cluster_probability[group_name] = cluster_probability[group_name] + weights[name] * cluster_probability[name]

    return cluster_probability


if __name__ == '__main__':

    enrichment_scores_dir = '/lustre/home/acct-medsjia/medsjia-wuzr/yy/prcc_multiomics/results/enrichment_scores/'
    save_dir = '/lustre/home/acct-medsjia/medsjia-wuzr/yy/prcc_multiomics/results/probability/nmf1/'

    predicted_cluster = 'nmf1'
    omics_weights = [1, 1, 1]

    # load data
    enrichment_scores = None
    cluster_names = []
    omics_names = []
    for f in Path(enrichment_scores_dir).glob('*/*.csv'):
        cluster = f.stem
        omics = f.parent.stem
        cluster_names.append(cluster)
        omics_names.append(omics)
        df = pd.read_csv(str(f), sep=',', index_col=0, header=0)
        if enrichment_scores is None:
            enrichment_scores = df['NES'].to_frame()
            enrichment_scores.columns = ['%s_%s' % (omics, cluster)]
        else:
            sample_names = df.index.intersection(enrichment_scores.index).tolist()
            enrichment_scores = enrichment_scores.loc[sample_names]
            enrichment_scores['%s_%s' % (omics, cluster)] = df.loc[sample_names, 'NES']
    cluster_names = list(set(cluster_names))
    omics_names = list(set(omics_names))

    # calculate probability
    probability = calculate_probability(
        enrichment_scores,
        cluster_names=cluster_names, predicted_cluster=predicted_cluster,
        omics_names=omics_names, omics_weights=omics_weights
    )

    Path(save_dir).mkdir(parents=True, exist_ok=True)
    probability.to_csv('%s/probability.csv' % save_dir, sep=',', index=True, header=True)

    print('Done!')
