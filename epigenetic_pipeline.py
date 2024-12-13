# coding: utf-8 # <- This is an encoding declaration

import requests
import pandas as pd
import json
import plotly.express as px
import panel as pn

'''
This script implements a KEGG-Enrichr Gene Regulatory Network Explorer, 
allowing users to input a KEGG pathway ID, retrieve associated genes, 
identify interacting transcription factors, and visualize the results 
in an interactive dashboard. Results are saved as a CSV file for further analysis.
'''

# Enable Panel extensions
pn.extension("plotly")


def get_genes_for_pathway(pathway_id, progress_log):
    """
    Retrieves a list of genes for a given pathway from the KEGG database.

    Args:
        pathway_id (str): The KEGG pathway ID (e.g., "hsa00010" for Glycolysis in humans).
        progress_log (pn.widgets.TextAreaInput): The progress log widget for updating progress.

    Returns:
        list: A list of gene identifiers associated with the pathway.
    """
    base_url = "https://rest.kegg.jp/link/hsa"
    url = f"{base_url}/{pathway_id}"

    try:
        response = requests.get(url)
        response.raise_for_status()  # Raise exception for HTTP errors

        # Extract gene IDs from the response removing the "hsa:" prefix
        lines = response.text.strip().split("\n")
        gene_ids = [line.split("\t")[1].replace("hsa:", "")
                    for line in lines if line]
        progress_log.value += f"Number of genes found for the pathway {pathway_id}: {len(gene_ids)}\n"
        return gene_ids

    except requests.exceptions.RequestException as e:
        progress_log.value += f"Error retrieving genes: {e}\n"
        return []


def get_gene_names(gene_ids, progress_log):
    """
    Retrieves gene names for the given gene IDs from KEGG.

    Args:
        gene_ids (list): A list of gene IDs.
        progress_log (pn.widgets.TextAreaInput): The progress log widget for updating progress.

    Returns:
        dict: A dictionary mapping gene IDs to gene names.
    """
    base_url = "https://rest.kegg.jp/get"
    gene_names = {}

    for gene_id in gene_ids:
        try:
            url = f"{base_url}/hsa:{gene_id}"
            response = requests.get(url)
            response.raise_for_status()

            gene_name = "Unknown"
            for line in response.text.split("\n"):
                if line.startswith("NAME"):
                    gene_name = line.split(
                        "NAME")[1].strip().replace("(RefSeq) ", "")
                    break

            gene_names[gene_id] = gene_name

        except requests.exceptions.RequestException as e:
            progress_log.value += f"An error occurred while retrieving gene name for {gene_id}: {e}\n"
            gene_names[gene_id] = "Unknown"

    return gene_names


def get_gene_symbols(gene_ids, progress_log):
    """
    Retrieves gene symbols for the given gene IDs from KEGG.

    Args:
        gene_ids (list): A list of gene IDs.
        progress_log (pn.widgets.TextAreaInput): The progress log widget for updating progress.

    Returns:
        dict: A dictionary mapping gene IDs to gene symbols.
    """
    base_url = "https://rest.kegg.jp/get"
    gene_symbols = {}

    for gene_id in gene_ids:
        try:
            url = f"{base_url}/hsa:{gene_id}"
            response = requests.get(url)
            response.raise_for_status()

            gene_symbol = "Unknown"
            for line in response.text.split("\n"):
                if line.startswith("SYMBOL"):
                    raw_symbol = line.split("      ")[1].split(",")[0].strip()
                    # Remove spaces and special characters
                    gene_symbol = ''.join(filter(str.isalnum, raw_symbol))
                    break

            gene_symbols[gene_id] = gene_symbol

        except requests.exceptions.RequestException as e:
            progress_log.value += f"An error occurred while retrieving gene symbol for {gene_id}: {e}\n"
            gene_symbols[gene_id] = "Unknown"

    return gene_symbols


def get_gene_ids(genes, progress_log, batch_size=15):
    """
    Retrieves Enrichr user list IDs for given gene symbols by processing them in batches.

    Args:
        genes (list): A list of gene symbols.
        progress_log (pn.widgets.TextAreaInput): The progress log widget for updating progress.
        batch_size (int): The number of genes to process in each batch (default is 15).

    Returns:
        list: A list of Enrichr user list IDs corresponding to the gene symbols.

    Raises:
        Exception: If the gene list cannot be analyzed after multiple attempts.
    """

    ENRICHR_URL = 'https://maayanlab.cloud/Enrichr/addList'
    gene_ids = []

    progress_log.value += "Batch processing of the retrieved genes initialized:\n"

    # Process genes in batches
    for batch_start in range(0, len(genes), batch_size):
        batch = genes[batch_start:batch_start + batch_size]
        progress_log.value += f"Processing batch: {batch}\n"

        # Now get ID for each gene in the batch
        for gene in batch:
            genes_str = '\n'.join(list(gene))
            description = 'Query gene list'
            payload = {
                'list': (None, genes_str),
                'description': (None, description)
            }

            # Get ID with retries
            response = requests.post(ENRICHR_URL, files=payload)
            if not response.ok:
                # Try a couple more times since NCBI might not return
                # info on the first call
                for i in range(10):
                    response = requests.post(ENRICHR_URL, files=payload)
                    if response.ok:
                        break
                    elif i >= 9:  # Log error on the final attempt
                        progress_log.value += f"Error analyzing gene: {gene} after 10 attempts\n"
                        raise Exception(
                            'Error analyzing gene list for gene: ', gene)

            data = json.loads(response.text)
            gene_ids.append(data['userListId'])

    return gene_ids


def get_tf_interactions(pathway_id, progress_log):
    """
    Retrieves interacting transcription factors (TFs) for genes in a given KEGG pathway using the ENCODE_TF_ChIP-seq_2015 library from Enrichr.

    Parameters:
    - pathway_id: The KEGG pathway ID (e.g., "hsa00010" for Glycolysis in humans).
    - progress_log (pn.widgets.TextAreaInput): The progress log widget for updating progress.

    Returns:
    - A pandas DataFrame containing gene names, symbols, interacting TFs, and their counts.
    """
    # Retrieve gene IDs for the pathway
    gene_ids_kegg = get_genes_for_pathway(pathway_id, progress_log)

    # Get gene names and symbols for the retrieved gene IDs
    progress_log.value += "Getting necessary information for all the genes retrieved...\n"
    gene_names_dict = get_gene_names(gene_ids_kegg, progress_log)
    gene_symbols_dict = get_gene_symbols(gene_ids_kegg, progress_log)

    # Prepare the gene names and symbols in lists
    gene_names = list(gene_names_dict.values())
    gene_symbols = list(gene_symbols_dict.values())

    # Get gene IDs from Enrichr
    gene_ids_enrichr = get_gene_ids(gene_symbols, progress_log)

    # Enrichr URL and query string for ENCODE_TF_ChIP-seq_2015 library
    ENRICHR_URL = 'https://maayanlab.cloud/Enrichr/enrich'
    query_string = '?userListId=%s&backgroundType=ENCODE_TF_ChIP-seq_2015'

    # Initialize a dictionary to store results
    results = {
        "Gene Name": [],
        "Gene Symbol": [],
        "Interacting TFs": [],
        "Count of Interacting TFs": []
    }

    # Loop through each gene to find interacting TFs
    for gene_name, gene_symbol, gene_id in zip(gene_names, gene_symbols, gene_ids_enrichr):
        progress_log.value += f"Retrieving interacting TFs for {gene_symbol} with id: {gene_id}\n"
        response = requests.get(ENRICHR_URL + query_string % gene_id)
        if not response.ok:
            raise Exception(
                f"Error fetching enrichment results for gene: {gene_symbol}")
        data = json.loads(response.text)
        tf_list = []
        if "ENCODE_TF_ChIP-seq_2015" in data:
            tf_info = data["ENCODE_TF_ChIP-seq_2015"]
            for entry in tf_info:
                tf_name = entry[1].split()[0]
                tf_list.append(tf_name)
        unique_tfs = list(set(tf_list))
        results["Gene Name"].append(gene_name)
        results["Gene Symbol"].append(gene_symbol)
        results["Interacting TFs"].append(", ".join(unique_tfs))
        results["Count of Interacting TFs"].append(len(unique_tfs))

    # Create a DataFrame from the results dictionary
    results_df = pd.DataFrame(results)
    results_df = results_df.reset_index(drop=True)
    results_df.index += 1
    return results_df

# Panel Dashboard Implementation


def dashboard():
    dashboard_title = pn.pane.HTML("""
    <h1 style="text-align:center;">KEGG-Enrichr Gene Regulatory Network Explorer: Pathway-to-Transcription Factor Analysis</h1>
    """)

    progress_log = pn.widgets.TextAreaInput(
        name='Progress Log', value='', height=80, width=1450)
    pathway_input = pn.widgets.TextInput(
        name='Enter KEGG Pathway ID:', placeholder='e.g., hsa00010')
    run_button = pn.widgets.Button(name='Run Analysis', button_type='primary')

    @pn.depends(run_button.param.clicks, watch=True)
    def run_analysis(_):
        if pathway_input.value:
            try:
                progress_log.value += f"Fetching genes for pathway {pathway_input.value}...\n"
                results_df = get_tf_interactions(
                    pathway_input.value, progress_log)
                progress_log.value += "Analysis complete.\n"

                # Save the results to a CSV file
                file_name = f"{pathway_input.value}_genes_with_interacting_TFs.csv"

                # Save without index column
                results_df.to_csv(file_name, index=False)

                progress_log.value += f"File saved as {file_name}\n"

                # Plot Results
                fig = px.bar(
                    results_df,
                    x='Gene Symbol',
                    y='Count of Interacting TFs',
                    orientation='v',
                    title=f'Interacting Transcription Factors (TFs) for every gene associated with pathway {pathway_input.value}',
                    labels={'Count of Interacting TFs': 'Count',
                            'Gene Symbol': 'Gene Symbols'}
                )
                plot_pane.object = fig

                # Display Results
                results_pane.object = results_df

            except Exception as e:
                progress_log.value += f"Error: {e}\n"
        else:
            progress_log.value += "Please enter a valid KEGG Pathway ID.\n"

    plot_pane = pn.pane.Plotly(height=370, width=1450)
    results_pane = pn.pane.DataFrame(width=1450, height=235)

    dashboard_layout = pn.Column(
        dashboard_title,
        pn.Row(pathway_input, run_button),
        progress_log,
        plot_pane,
        results_pane
    )

    return dashboard_layout


dashboard().servable()


def main():
    print("Starting KEGG-Enrichr Gene Regulatory Network Explorer")
    dashboard_instance = dashboard()
    print("Dashboard initialized. Results will be saved in a csv file when analysis is run.")
    dashboard_instance.show()


if __name__ == "__main__":
    main()
