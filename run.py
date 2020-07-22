#!/usr/bin/env python
import sys
import os
import json

import pandas as pd
import numpy as np

import networkx as nx

from scipy.stats import spearmanr, pearsonr
from sklearn.preprocessing import scale
from sklearn.linear_model import *


def build_database_level_structure(G, interaction_databases_path, df_expr):
    '''
    Add edges to the network G from the
    user-specified interaction databases.
    '''
    # Iterate though all available interaction databases
    for filename in os.listdir(interaction_databases_path):
        if not filename.endswith(".tsv"):
            continue

        df = pd.read_csv(
            "{}/{}".format(interaction_databases_path, filename),
            sep="\t"
        )
        # Read node types from the header
        node_type1, node_type2 = df.columns
        for i, (n1, n2) in df.iterrows():
            if n1 not in G.nodes or n2 not in G.nodes :
                continue

            x1 = df_expr.loc[n1].to_numpy()
            x2 = df_expr.loc[n2].to_numpy()

            G.add_edge(n1, n2, **{
                "Type": "{}_{}".format(node_type1, node_type2),
                "Source": filename,
                "Spearman correlation": spearmanr(x1, x2)[0]
            })

            # Add attributes to nodes
            for n, x, node_type in [(n1, x1, node_type1), (n2, x2, node_type2)]:
                G.nodes[n]["Median expression"] = np.median(x)
                # This field will be useful for some programs (e.g. yED)
                G.nodes[n]["Label"] = n
                # If gene was already mentioned as TF do not mark it as a Gene
                if G.nodes[n].get("Type") == "TF" and node_type == "Gene":
                    continue

                G.nodes[n]["Type"] = node_type

    G.remove_edges_from(nx.selfloop_edges(G))
    G.remove_nodes_from(list(nx.isolates(G)))


def filter_non_correlated_edges(G, cutoff_percentile, sign_constraints):
    '''
    For each interaction type and each interaction direction
    remove edges with low absolute value of correlation.
    The threshold is defined in terms of distribution percentile.
    Also remove edges with correlation signs not matching
    constraints provided.
    '''
    # First, calculate thresholds for each interaction type and direction
    edge_types = {G.edges[e]["Type"] for e in G.edges}
    corr_thresholds = {}
    for i, edge_type in enumerate(sorted(edge_types)):
        corrs = np.array([
            G.edges[e]["Spearman correlation"] for e in G.edges
            if G.edges[e]["Type"] == edge_type
        ])

        corr_thresholds[(edge_type, "+")] = np.percentile(corrs[corrs > 0], cutoff_percentile)
        corr_thresholds[(edge_type, "-")] = np.percentile(np.abs(corrs[corrs < 0]), cutoff_percentile)

    edges_to_remove = []
    for e in G.edges:
        type_ = G.edges[e]["Type"]
        source = G.edges[e]["Source"]
        corr = G.edges[e]["Spearman correlation"]
        sign = "+" if corr > 0 else "-"

        # Edge is removed if it's correlation is below the threshold or it
        # violates a sign constraint from the config
        if (
            abs(corr) < corr_thresholds[(type_, sign)] or
            (source in sign_constraints and sign_constraints[source] != sign)
        ):
            edges_to_remove.append(e)

    G.remove_edges_from(edges_to_remove)
    G.remove_nodes_from(list(nx.isolates(G)))


def calculate_scores(G, df_expr):
    '''
    Calculate interaction score for each node and
    incoming score for each edge
    '''
    for n in G.nodes:
        regulators = {n1: data for n1, n2, data in G.in_edges(nbunch=[n], data=True)}
        if not regulators:
            continue

        # Regulator nodes are predictors
        X = scale(np.array([df_expr.loc[r] for r in sorted(regulators)]).transpose())
        # Current node is a dependent variable
        y = scale(np.array(df_expr.loc[n]))

        model = RidgeCV(normalize=False)
        model.fit(X, y)

        R2 = model.score(X, y)
        G.nodes[n]["Incoming score"] = R2

        for r, c in zip(sorted(regulators), model.coef_):
            G.edges[(r, n)]["Interaction score"] = abs(c)


def extract_core(G, incoming_score_threshold, interaction_score_cutoff_percentile):
    '''
    Extract the core from the interaction network, i.e.
    subgraph composed of highly regulated nodes (incoming score)
    and highly regulating nodes (interaction score)
    '''
    # Identify nodes with high incoming scores
    core_nodes = [n for n in G.nodes if G.nodes[n].get("Incoming score", 0) >= incoming_score_threshold]
    # Add nodes from adjacent incoming edges
    core_nodes += [n1 for n in core_nodes for n1, n2 in G.in_edges(nbunch=[n])]
    core = G.subgraph(core_nodes).copy()

    # Remove edges with low interaction scores
    interaction_score_threshold = np.percentile(
        [core.edges[e]["Interaction score"] for e in core.edges],
        interaction_score_cutoff_percentile
    )
    edges_to_remove = [e for e in core.edges if core.edges[e]["Interaction score"] < interaction_score_threshold]
    core.remove_edges_from(edges_to_remove)
    core.remove_nodes_from(list(nx.isolates(core)))

    return core


def save_graph(G, out_path, filename):
    if not os.path.exists(out_path):
        os.mkdir(out_path)

    nx.write_graphml(G, "{}/{}.graphml".format(out_path, filename))


def node_report(G, out_file):
    print(file=out_file)
    print("Number of nodes: {}".format(len(G.nodes)), file=out_file)
    print("Node types:", file=out_file)
    node_types = [G.nodes[n]["Type"] for n in G.nodes]
    for node_type in np.unique(node_types):
        count = node_types.count(node_type)
        print("\t- {} nodes with type = {}".format(count, node_type), file=out_file)


def edge_report(G, out_file):
    print(file=out_file)
    print("Number of edges: {}".format(len(G.edges)), file=out_file)
    print("Edge types:", file=out_file)
    edge_types = [G.edges[e]["Type"] for e in G.edges]
    for edge_type in np.unique(edge_types):
        count = edge_types.count(edge_type)
        print("\t- {} edges with type = {}".format(count, edge_type.replace("_", " -> ")), file=out_file)


if __name__ == "__main__":
    '''
    Load input data from provided config file
    '''

    if len(sys.argv) < 2:
        print("Usage: python run.py /path/to/config.json", file=sys.stderr)
        sys.exit(1)

    # Load the config
    try:
        config = json.load(open(sys.argv[1], "r"))
    except:
        print("Please provide a valid JSON config file", file=sys.stderr)
        raise

    # Load gene and miRNA expression tables
    try:
        df_gene = pd.read_csv(
            config["gene_expression_table_path"],
            sep="\t",
            index_col=0
        )
        df_miRNA = pd.read_csv(
            config["miRNA_expression_table_path"],
            sep="\t",
            index_col=0
        )
        df_expr = pd.concat([df_gene, df_miRNA], axis=0)
    except:
        print("Please provide valid gene/miRNA expression tables", file=sys.stderr)
        raise

    # Create output folder if not exists
    try:
        if not os.path.exists(config["output_path"]):
            os.mkdir(config["output_path"])
    except:
        print("Please provide valid output directory")
        raise

    '''
    Build the network
    '''

    # Start building the network
    G = nx.DiGraph()
    # Build network on top of nodes from expression table
    G.add_nodes_from(df_expr.index.to_list())
    # Load edges from provided databases
    build_database_level_structure(G, config["interaction_databases_path"], df_expr)
    # Remove edges with low correlations
    filter_non_correlated_edges(G, config["Spearman_correlation_cutoff_percentile"], config["sign_constraints"])
    # Calculate interaction and incoming scores using ridge regression
    calculate_scores(G, df_expr)
    # Extract subgraph with high scores
    G = extract_core(
        G,
        config["incoming_score_threshold"],
        config["interaction_score_cutoff_percentile"]
    )

    save_graph(G, "{}/networks".format(config["output_path"]), "network")

    '''
    Analyze network and create a report
    '''

    report_file = open("{}/report.txt".format(config["output_path"]), "w")

    # Summary on nodes and edges
    print("Network (summary)", file=report_file)
    node_report(G, out_file=report_file)
    edge_report(G, out_file=report_file)
    print("\n" + "*"*17 + "\n", file=report_file)

    for direction, func in [("out", G.out_degree), ("in", G.in_degree)]:
        print("Distribution of {}-degrees\n".format(direction), file=report_file)
        degrees = {n: func(nbunch=n) for n in G.nodes}
        for deg in np.unique(list(degrees.values())):
            count = list(degrees.values()).count(deg)
            print("\t- {} nodes with {}-degree = {}".format(count, direction, deg), file=report_file)

        print("\nTop-10 {}-degrees\n".format(direction), file=report_file)
        for node, deg in sorted(degrees.items(), key=lambda x: x[1], reverse=True)[0:10]:
            print("\t- {}, {}-degree = {}".format(node, direction, deg), file=report_file)

        print("\n" + "*"*17 + "\n", file=report_file)


    # Summary of weakly connected components
    components = sorted(nx.weakly_connected_components(G), key=lambda comp: len(comp), reverse=True)
    for i, comp in enumerate(components):
        comp = G.subgraph(comp)
        save_graph(comp, "{}/networks".format(config["output_path"]), "weakly_conn_component_{}".format(i + 1))

        print("Weakly connected component #{} (summary)".format(i + 1), file=report_file)
        node_report(comp, out_file=report_file)
        edge_report(comp, out_file=report_file)
        print("\n" + "*"*17 + "\n", file=report_file)

    # Summary of non-trivial strongly connected components
    components = sorted(nx.strongly_connected_components(G), key=lambda comp: len(comp), reverse=True)
    for i, comp in enumerate(components):
        if len(comp) == 1:
            continue

        comp = G.subgraph(comp)
        save_graph(comp, "{}/networks".format(config["output_path"]), "strongly_conn_component_{}".format(i + 1))

        print("Strongly connected component #{} (summary)".format(i + 1), file=report_file)
        node_report(comp, out_file=report_file)
        edge_report(comp, out_file=report_file)
        print("\n" + "*"*17 + "\n", file=report_file)

    # Calculate and export node statistics
    columns = ["Type", "Incoming score", "Median expression"]
    df = pd.DataFrame(
        [[G.nodes[n].get(c) for c in columns] for n in G.nodes],
        columns=columns
    )
    df.index = G.nodes
    df.index.name = "Node"
    df = df.sort_index()
    # In- and out- degrees, centrality measures
    df["In-degree"] = [G.in_degree(nbunch=n) for n in df.index]
    df["Out-degree"] = [G.out_degree(nbunch=n) for n in df.index]
    closeness_centrality = nx.algorithms.centrality.closeness_centrality(G)
    betweenness_centrality = nx.algorithms.centrality.betweenness_centrality(G)
    df["Closeness centrality"] = [closeness_centrality[n] for n in df.index]
    df["Betweenness centrality"] = [betweenness_centrality[n] for n in df.index]

    df.to_csv("{}/node_stats.csv".format(config["output_path"]))

    # Calculate and export edge statistics
    columns = ["Type", "Source", "Interaction score", "Spearman correlation"]
    df = pd.DataFrame(
        [[G.edges[e].get(c) for c in columns] for e in G.edges],
        columns=columns
    )
    df["Type"] = [t.replace("_", " -> ") for t in df["Type"]]
    df["From"] = [e[0] for e in G.edges]
    df["To"] = [e[1] for e in G.edges]
    # Reorder columns
    df = df[["From", "To"] + columns]

    betweenness_centrality = nx.algorithms.centrality.edge_betweenness_centrality(G)
    df["Betweenness centrality"] = [betweenness_centrality[e] for e in G.edges]

    df = df.sort_values("From")
    df.to_csv("{}/edge_stats.csv".format(config["output_path"]), index=None)
