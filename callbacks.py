import dash
from dash import callback, dcc, html, Input, Output, State
import plotly.graph_objects as go
import pandas as pd
import numpy as np
from app_context import * 
from plotly.subplots import make_subplots
import plotly.express as px
import plotly.colors
import matplotlib.colors as mcolors
from plotly.colors import sample_colorscale
import matplotlib.pyplot as plt
from numpy import inf
import colorcet as cc
import re
import utils
from utils import row_to_binary, row_to_trinary, binary_to_decimal, trinary_to_decimal, uppercase_genes, update_bar_plots, parse_gene_input, pool_genes_in_df, get_gene_label, get_tri_thresholds, pooled_gene_names
import logging, time

################################################################################################################
##### Callback decorator
################################################################################################################

@callback(
    [Output('sunburst1', 'figure'),
     Output('sunburst2', 'figure'),
     # Output('gene-bar-plots', 'children'),
     Output('store-fig-bin', 'data'),
     Output('store-fig-trin', 'data'),
     Output('store-genes', 'data'),
     # Output('store-expression-div', 'data'),
     # Output('store-expression-class', 'data'),
     Output('store-expression-neuronal-cells', 'data'),
     Output('store-expression-nonneuronal-cells', 'data')],
     # Output('store-expression-subclass', 'data'),
     # Output('store-expression-supertype', 'data')],
    [Input('update-button', 'n_clicks'),
     Input('toggle-trinary', 'value')], 
    [State('gene-input', 'value'),
     State('store-fig-bin', 'data'),
     State('store-fig-trin', 'data'),
     State('store-genes', 'data'),
     # State('store-expression-div', 'data'),
     # State('store-expression-class', 'data'),
     State('store-expression-neuronal-cells', 'data'),
     State('store-expression-nonneuronal-cells', 'data')]
)

################################################################################################################
##### Actual callback
################################################################################################################

                                   
def update_sunburst(n_clicks, toggle_value, gene_input, stored_fig_bin, stored_fig_trin, store_genes, stored_expression_neuronal_cells, stored_expression_nonneuronal_cells):
    if toggle_value is None:
        toggle_value = ['binary']
    
    ctx = dash.callback_context  # Identifies what triggered the callback
    # bar_plots = []

    if not ctx.triggered:  # If nothing has triggered the function yet
        return go.Figure(), go.Figure(), go.Figure().to_dict(), go.Figure().to_dict(), [], pd.DataFrame().to_dict(), pd.DataFrame().to_dict()

    trigger_id = ctx.triggered[0]['prop_id'].split('.')[0]  # Get the trigger ID
    
    # If button is clicked, update sunburst and store new figures
    if trigger_id == "update-button":
        genes_to_test = parse_gene_input(gene_input)
        # genes_to_test = [gene.upper() for gene in genes_to_test]
        genes_to_test = uppercase_genes(genes_to_test)
        num_genes = sum(len(g) if isinstance(g, list) else 1 for g in genes_to_test)
        
        print(f"Updating plots for genes: {genes_to_test}")

        neuron_cells_with_genes_local = pool_genes_in_df(neuron_cells_with_genes, genes_to_test)
        print('Pooled (if needed) neuronal genes.')
        print(neuron_cells_with_genes_local.columns)
        nonneuron_cells_with_genes_local = pool_genes_in_df(nonneuron_cells_with_genes, genes_to_test)
        print('Pooled (if needed) nonneuronal genes.')
        genes_to_test = pooled_gene_names(genes_to_test)
        print(genes_to_test)

        neuron_cells_with_genes_local[genes_to_test] = neuron_cells_with_genes_local[genes_to_test].apply(pd.to_numeric, errors='coerce')
        nonneuron_cells_with_genes_local[genes_to_test] = nonneuron_cells_with_genes_local[genes_to_test].apply(pd.to_numeric, errors='coerce')
        
        ################################################################################################################
        ##### Sunburst 1 - expression (% of expressed supertypes in subclasses for given genes)
        ################################################################################################################

        # NEURONAL
        # Step 0: Create a boolean column indicating if any gene > 2
        any_gene_above_2 = (
            neuron_cells_with_genes_local[genes_to_test] >= 2
        ).any(axis=1)
        
        # Attach this to the main DataFrame
        neuron_cells_with_genes_local = neuron_cells_with_genes_local.copy()
        neuron_cells_with_genes_local["any_gene_above_2"] = any_gene_above_2
        
        # Step 1 & 2: Group and compute proportion
        cluster_stats_neuron = (
            neuron_cells_with_genes_local
            .groupby(["supercluster", "cluster"])['any_gene_above_2']
            .mean()  # proportion of True
            .mul(100)  # convert to percent
            .reset_index(name="expr_proportion")
        )
        
        # # Step 3: Add counts
        # cluster_counts = (
        #     neuron_cells_with_genes_local
        #     .groupby(["supercluster", "cluster"])
        #     .size()
        #     .reset_index(name="cell_count")
        # )

        # Step 3: Sum cell_count from subclusters to clusters
        cluster_counts = (
            neuron_cells_with_genes_local
            .groupby(["supercluster", "cluster"])['cell_count']
            .sum()
            .reset_index()
        )
        
        cluster_stats_neuron = cluster_stats_neuron.merge(cluster_counts, on=["supercluster", "cluster"])

        # NON-NEURONAL
        
        # Step 0: Create a boolean column indicating if any gene > 2
        any_gene_above_2 = (
            nonneuron_cells_with_genes_local[genes_to_test] >= 2
        ).any(axis=1)
        
        # Attach this to the main DataFrame
        nonneuron_cells_with_genes_local = nonneuron_cells_with_genes_local.copy()
        nonneuron_cells_with_genes_local["any_gene_above_2"] = any_gene_above_2
        
        # Step 1 & 2: Group and compute proportion > 2
        cluster_stats_nonneuron = (
            nonneuron_cells_with_genes_local
            .groupby(["supercluster", "cluster"])['any_gene_above_2']
            .mean()  # proportion of True
            .mul(100)  # convert to percent
            .reset_index(name="expr_proportion")
        )
        
        # Step 3: Add counts (for size)
        # cluster_counts = (
        #     nonneuron_cells_with_genes_local
        #     .groupby(["supercluster", "cluster"])
        #     .size()
        #     .reset_index(name="cell_count")
        # )

        # Step 3: Sum cell_count from subclusters to clusters
        cluster_counts = (
            nonneuron_cells_with_genes_local
            .groupby(["supercluster", "cluster"])['cell_count']
            .sum()
            .reset_index()
        )
        
        cluster_stats_nonneuron = cluster_stats_nonneuron.merge(cluster_counts, on=["supercluster", "cluster"])
        
        # Create subplot layout with 1 row and 2 columns
        fig = make_subplots(rows=1, cols=2, 
                    specs=[[{"type": "domain"}, {"type": "domain"}]],  # "domain" for sunburst
                    subplot_titles=[f"Neuronal Cells {len(neuron_cells_with_genes_local):,}", f"Non-Neuronal Cells {len(nonneuron_cells_with_genes_local):,}"])  # Add titles

        # Create first sunburst (Neuronal)
        fig1 = px.sunburst(
            cluster_stats_neuron,
            #path=["parent", "id"],
            path=["supercluster", "cluster"],
            values="cell_count",
            color="expr_proportion",
            #color_continuous_scale="Blues",
            #range_color=[0, 100],
            #title=f"Expression Proportion for {gene_col}"
        )
        
        # Create second sunburst (Non-Neuronal)
        fig2 = px.sunburst(
            cluster_stats_nonneuron,
            path=["supercluster", "cluster"],
            values="cell_count",
            color="expr_proportion",
            #color_continuous_scale="Blues",
            #range_color=[0, 100],
            #title=f"Expression Proportion for {gene_col}"
        )

        # Add both sunburst plots to the figure
        fig.add_trace(fig1.data[0], row=1, col=1)  # Add first chart
        fig.add_trace(fig2.data[0], row=1, col=2)  # Add second chart

        # Adjust layout
        fig.update_layout(
            width=1200,  # Increase width (default is ~700)
            height=800,  # Increase height (default is ~450)
            title_text=f"{genes_to_test}",
            title_x=0.5,
            coloraxis=dict(
                    colorscale='blues',  # Set colorscale
                    cmin=0,  # Ensure minimum value is 0
                    cmax=100  # Ensure maximum value is 1
            ),
            showlegend=True,
        )

        fig.update_traces(marker=dict(line=dict(width=0.1, color='black')))

        print('Expression plots done.')


        ################################################################################################################
        ##### Sunburst 2 - binary
        ################################################################################################################
        
        avg_expr_neuron = neuron_cells_with_genes_local.groupby(["supercluster", "cluster"])[genes_to_test].mean()
        avg_expr_nonneuron = nonneuron_cells_with_genes_local.groupby(["supercluster", "cluster"])[genes_to_test].mean()

        # Neuronal
        binary_expr_neuron = avg_expr_neuron.copy()
        for gene in genes_to_test:
            threshold = neuronal_gene_thresholds.get(gene, [2])[0]  # fallback to 2
            binary_expr_neuron[gene] = (avg_expr_neuron[gene] >= threshold).astype(int)

        # Non-neuronal
        binary_expr_nonneuron = avg_expr_nonneuron.copy()
        for gene in genes_to_test:
            threshold = nonneuronal_gene_thresholds.get(gene, [2])[0]  # fallback to 2
            binary_expr_nonneuron[gene] = (avg_expr_nonneuron[gene] >= threshold).astype(int)

        binary_expr_neuron['binary_signature'] = binary_expr_neuron.apply(row_to_binary, axis=1)
        binary_expr_nonneuron['binary_signature'] = binary_expr_nonneuron.apply(row_to_binary, axis=1)
        
        binary_expr_neuron['Cluster'] = binary_expr_neuron['binary_signature'].apply(binary_to_decimal)
        binary_expr_nonneuron['Cluster'] = binary_expr_nonneuron['binary_signature'].apply(binary_to_decimal)

        binary_expr_neuron = binary_expr_neuron.reset_index()
        binary_expr_nonneuron = binary_expr_nonneuron.reset_index()
        
        cluster_stats_neuron = cluster_stats_neuron.merge(
            binary_expr_neuron[['supercluster', 'cluster', 'binary_signature', 'Cluster']],
            on=['supercluster', 'cluster'],
            how='left'
        )
        
        cluster_stats_nonneuron = cluster_stats_nonneuron.merge(
            binary_expr_nonneuron[['supercluster', 'cluster', 'binary_signature', 'Cluster']],
            on=['supercluster', 'cluster'],
            how='left'
        )
        

        #Create a colormap
        
        # Define a fixed set of distinct colors 
        #distinct_colors = plotly.colors.qualitative.Safe + plotly.colors.qualitative.Dark24
        #distinct_colors = cc.glasbey[50:306]

        distinct_colors = [mcolors.to_hex(rgb) for rgb in cc.glasbey_hv[:256]]

        # print(type(cc.glasbey[0]))         # likely: <class 'str'>, e.g. '#e60049'
        # print(type(cc.glasbey_hv[0]))      # likely: <class 'tuple'>, e.g. (0.8, 0.4, 0.2)
        
        #distinct_colors = cc.glasbey_hv[:256]
        #print(f'Colors available {len(distinct_colors)}')
        
        combined_clusters = pd.concat([cluster_stats_neuron, cluster_stats_nonneuron], ignore_index=True)
        unique_signatures = combined_clusters['binary_signature'].dropna().unique()
        num_unique_signatures = len(unique_signatures)
        
        if num_unique_signatures > len(distinct_colors):
            raise ValueError(f"Not enough distinct colors for {num_unique_signatures} unique binary signatures.")
        
        possible_binary_signatures = 2**len(genes_to_test)
        
        # Create a mapping dictionary for binary_signature → color
        binary_signature_to_color = {
            sig: distinct_colors[i]
            for i, sig in enumerate(sorted(unique_signatures))  # optional: sort for reproducibility
        }
        
        color_map = {
            **binary_signature_to_color,
            **{d: "white" for d in cluster_stats_neuron["supercluster"].dropna().unique()},
            **{c: "white" for c in cluster_stats_neuron["cluster"].dropna().unique()},
            # 'CB Glut': "white",
            '(?)': 'white'
        }

        #Plot

        # Create subplot layout with 1 row and 2 columns
        fig_bin = make_subplots(
            rows=1, cols=2, 
            specs=[[{"type": "domain"}, {"type": "domain"}]],  # "domain" for sunburst
            subplot_titles=[f"Neuronal Cells {len(neuron_cells_with_genes):,}", f"Non-Neuronal Cells {len(nonneuron_cells_with_genes):,}"]
        )
        
        # Create first sunburst (Neuronal)
        fig1 = px.sunburst(
            cluster_stats_neuron,
            path=["supercluster", "cluster"],  # Keep hierarchy
            values="cell_count",
            color="binary_signature",  # Color by binary signature
            color_discrete_map=color_map  # Use predefined colors
        )
        
        # Create second sunburst (Non-Neuronal)
        fig2 = px.sunburst(
            cluster_stats_nonneuron,
            path=["supercluster", "cluster"],  # Keep hierarchy
            values="cell_count",
            color="binary_signature",  # Color by binary signature
            color_discrete_map=color_map  # Use predefined colors
        )
        
        # Add both sunburst plots to the figure
        fig_bin.add_trace(fig1.data[0], row=1, col=1)  # Add first chart
        fig_bin.add_trace(fig2.data[0], row=1, col=2)  # Add second chart
        
        # **Manually Add a Custom Legend for Binary Signatures**
        legend_traces = []
        for binary_signature, color in binary_signature_to_color.items():
            legend_traces.append(go.Scatter(
                x=[None], y=[None],  # No actual data points, just legend markers
                mode='markers',
                marker=dict(size=12, color=color),
                name=str(binary_signature),  # Label the legend with binary signature
                showlegend=True
            ))
        
        # Adjust layout
        fig_bin.update_layout(
            paper_bgcolor='white',  # Removes the grid-like background
            plot_bgcolor='white',   # Ensures no grid in the plot area
            width=1200,  # Increase width
            height=800,  # Increase height
            xaxis=dict(showticklabels=False, showgrid=False, zeroline=False, title=""),
            yaxis=dict(showticklabels=False, showgrid=False, zeroline=False, title=""),
            title_text=f"{genes_to_test}",
            title_x=0.5,
            showlegend=True,  # Ensure legend is displayed
            #margin=dict(t=50, b=10, l=10, r=10),
            legend=dict(title=f"Binary Clusters<br>{num_unique_signatures}/{possible_binary_signatures}", font=dict(size=12))
        )
        
        # Add legend traces to figure
        for trace in legend_traces:
            fig_bin.add_trace(trace)
        
        fig_bin.update_traces(marker=dict(line=dict(width=0.1, color='black')))

        print('Binary expression plots done.')

        ################################################################################################################
        ##### Sunburst 2 - trinary
        ################################################################################################################

        # Neuronal
        trinary_expr_neuron = avg_expr_neuron.copy()
        for gene in genes_to_test:
            lower, upper = neuronal_gene_thresholds.get(gene, [2, 7])
            trinary_expr_neuron[gene] = avg_expr_neuron[gene].apply(
                lambda x: 0 if x < lower else (1 if x < upper else 2)
            )
        
        # Non-neuronal
        trinary_expr_nonneuron = avg_expr_nonneuron.copy()
        for gene in genes_to_test:
            lower, upper = nonneuronal_gene_thresholds.get(gene, [2, 7])
            trinary_expr_nonneuron[gene] = avg_expr_nonneuron[gene].apply(
                lambda x: 0 if x < lower else (1 if x < upper else 2)
            )
        
        trinary_expr_neuron['trinary_signature'] = trinary_expr_neuron.apply(row_to_trinary, axis=1)
        trinary_expr_nonneuron['trinary_signature'] = trinary_expr_nonneuron.apply(row_to_trinary, axis=1)
        
        trinary_expr_neuron['Cluster'] = trinary_expr_neuron['trinary_signature'].apply(trinary_to_decimal)
        trinary_expr_nonneuron['Cluster'] = trinary_expr_nonneuron['trinary_signature'].apply(trinary_to_decimal)
        
        trinary_expr_neuron = trinary_expr_neuron.reset_index()
        trinary_expr_nonneuron = trinary_expr_nonneuron.reset_index()
        
        cluster_stats_neuron = cluster_stats_neuron.merge(
            trinary_expr_neuron[['supercluster', 'cluster', 'trinary_signature', 'Cluster']],
            on=['supercluster', 'cluster'],
            how='left'
        )
        
        cluster_stats_nonneuron = cluster_stats_nonneuron.merge(
            trinary_expr_nonneuron[['supercluster', 'cluster', 'trinary_signature', 'Cluster']],
            on=['supercluster', 'cluster'],
            how='left'
        )
        
        #distinct_colors = cc.glasbey[:256]
        combined_clusters = pd.concat([cluster_stats_neuron, cluster_stats_nonneuron], ignore_index=True)
        unique_signatures = combined_clusters['trinary_signature'].dropna().unique()
        num_unique_signatures = len(unique_signatures)
        possible_trinary_signatures = 3 ** len(genes_to_test)
        
        if num_unique_signatures > len(distinct_colors):
            raise ValueError(f"Not enough distinct colors for {num_unique_signatures} trinary signatures.")
        
        trinary_signature_to_color = {
            sig: distinct_colors[i]
            for i, sig in enumerate(sorted(unique_signatures))
        }
        
        color_map = {
            **trinary_signature_to_color,
            **{d: "white" for d in cluster_stats_neuron["supercluster"].dropna().unique()},
            **{c: "white" for c in cluster_stats_neuron["cluster"].dropna().unique()},
            '(?)': 'white'
        }


        #Plot
        
        # Create subplot layout with 1 row and 2 columns
        fig_trin = make_subplots(
            rows=1, cols=2, 
            specs=[[{"type": "domain"}, {"type": "domain"}]],  # "domain" for sunburst
            subplot_titles=[f"Neuronal Cells {len(neuron_cells_with_genes_local):,}", f"Non-Neuronal Cells {len(nonneuron_cells_with_genes_local):,}"]
        )
        
        # Create first sunburst (Neuronal)
        fig1 = px.sunburst(
            cluster_stats_neuron,
            path=["supercluster", "cluster"],  # Keep hierarchy
            values="cell_count",
            color="trinary_signature",  # Color by binary signature
            color_discrete_map=color_map  # Use predefined colors
        )
        
        # Create second sunburst (Non-Neuronal)
        fig2 = px.sunburst(
            cluster_stats_nonneuron,
            path=["supercluster", "cluster"],  # Keep hierarchy
            values="cell_count",
            color="trinary_signature",  # Color by binary signature
            color_discrete_map=color_map  # Use predefined colors
        )
        
        # Add both sunburst plots to the figure
        fig_trin.add_trace(fig1.data[0], row=1, col=1)  # Add first chart
        fig_trin.add_trace(fig2.data[0], row=1, col=2)  # Add second chart
        
        # **Manually Add a Custom Legend for Binary Signatures**
        legend_traces = []
        for trinary_signature, color in trinary_signature_to_color.items():
            legend_traces.append(go.Scatter(
                x=[None], y=[None],  # No actual data points, just legend markers
                mode='markers',
                marker=dict(size=12, color=color),
                name=str(trinary_signature),  # Label the legend with binary signature
                showlegend=True
            ))
        
        # Adjust layout
        fig_trin.update_layout(
            paper_bgcolor='white',  # Removes the grid-like background
            plot_bgcolor='white',   # Ensures no grid in the plot area
            width=1200,  # Increase width
            height=800,  # Increase height
            xaxis=dict(showticklabels=False, showgrid=False, zeroline=False, title=""),
            yaxis=dict(showticklabels=False, showgrid=False, zeroline=False, title=""),
            title_text=f"{genes_to_test}",
            title_x=0.5,
            showlegend=True,  # Ensure legend is displayed
            #margin=dict(t=50, b=10, l=10, r=10),
            legend=dict(title=f"Trinary Clusters<br>{num_unique_signatures}/{possible_trinary_signatures}", font=dict(size=12))
        )
        
        # Add legend traces to figure
        for trace in legend_traces:
            fig_trin.add_trace(trace)
        
        fig_trin.update_traces(marker=dict(line=dict(width=0.1, color='black')))

        print('Trinary expression plots done.')

        ################################################################################################################
        ##### Bar plots - expression level in classes
        ################################################################################################################

        # unique_divisions = avg_expression_genesAll_div_df.index
        # unique_classes = avg_expression_genesAll_class_df.index

        # # Create a color map for divisions
        # division_colors = plt.cm.jet(np.linspace(0, 1, len(unique_divisions)))
        # division_colors_hex = [mcolors.to_hex(color) for color in division_colors]
        # division_colors_dict = {division: division_colors_hex[i] for i, division in enumerate(unique_divisions)}
        t = time.perf_counter()
        fig_bin_dict = fig_bin.to_dict()
        logging.info("fig_bin → dict %.3fs", time.perf_counter() - t)
        
        return (
            fig,
            fig_bin,
            # update_bar_plots(
            #     toggle_value,
            #     genes_to_test,
            #     avg_expression_div_df,
            #     avg_expression_class_df,
            #     expression_cells_pooled_df,
            #     gene_thresholds,
            #     binary_threshold=2,
            #     trinary_threshold_low=2,
            #     trinary_threshold_high=7
            # ), 
            fig_bin.to_dict(),
            fig_trin.to_dict(),
            genes_to_test,
            # avg_expression_div_df.to_dict(),     
            # avg_expression_class_df.to_dict(),
            # expression_cells_pooled_df.to_dict(), 
            neuron_cells_with_genes_local.to_dict(),
            nonneuron_cells_with_genes_local.to_dict()
        )
    
    ########################################################################

    # If toggle is switched, retrieve stored figures and update `sunburst2`
    elif trigger_id == "toggle-trinary":
        if stored_fig_bin is None or stored_fig_trin is None:
            return dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update

        fig_bin = go.Figure(stored_fig_bin)
        
        fig_trin = go.Figure(stored_fig_trin)
        fig_to_use = fig_trin if 'trinary' in toggle_value else fig_bin
        genes_to_test = store_genes

        neuron_cells_with_genes_local =  pd.DataFrame(stored_expression_neuronal_cells)
        nonneuron_cells_with_genes_local = pd.DataFrame(stored_expression_nonneuronal_cells)


        return (dash.no_update, 
                fig_to_use, 
                # update_bar_plots(
                #     toggle_value, 
                #     genes_to_test, 
                #     avg_expression_div_df, 
                #     avg_expression_class_df,
                #     expression_cells_pooled_df,
                #     gene_thresholds, 
                #     binary_threshold=2, 
                #     trinary_threshold_low=2, 
                #     trinary_threshold_high=7
                # ), 
                dash.no_update,
                dash.no_update,
                dash.no_update,
                dash.no_update,
                # dash.no_update,
                dash.no_update
        )

    ################################################################################################################
    ################################################################################################################
