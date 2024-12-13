# KEGG-Enrichr Gene Regulatory Network Explorer

## Overview
This Python application provides an innovative dashboard for exploring gene regulatory networks by seamlessly integrating KEGG pathway data with transcription factor interaction insights. By leveraging the ENCODE_TF_ChIP-seq_2015 library from Enrichr, the tool enables researchers to map complex gene-transcription factor relationships within specific biological pathways. 

The dashboard retrieves gene information from KEGG and cross-references transcription factor interactions using the ENCODE_TF_ChIP-seq_2015 dataset, offering a comprehensive view of potential regulatory mechanisms in biological pathways.

### Features

- Retrieves genes for a specified human KEGG pathway
- Fetches gene names and symbols from KEGG
- Finds interacting transcription factors using ENCODE_TF_ChIP-seq_2015 library from Enrichr
- Visualizes results with interactive Plotly bar chart and generates exportable result (`.csv` file) for further analysis
- Provides a user-friendly interface using Panel


## Prerequisites

### Python Version
- Python 3.12.4

### Required Libraries
- requests (v2.32.3)
- pandas (v2.2.3)
- plotly (v5.24.1)
- panel (v1.5.4)

## Installation

1. Download the code file: ```epigentic_pipeline.py```
2. Create and activate a virtual environment:
    ```bash
    python -m venv venv
    source venv/bin/activate  # On Windows: venv\Scripts\activate
    ```
3. Install dependencies:
    ```bash
    pip install requests==2.32.3 pandas==2.2.3 plotly==5.24.1 panel==1.5.4
    ```

## Usage

1. Run the dashboard:
    ```bash
    panel serve epigenetic_pipeline.py
    ```
2. Open a web browser and navigate to web browser in the output. Example url - http://localhost:5006/epigenetic_pipeline
3. Enter a valid human KEGG Pathway ID (e.g., "hsa00010" for Glycolysis) in the input field.
4. Click the ```Run Analysis``` button to start the analysis.

### Input

- **Valid Human KEGG Pathway ID** (must start with "hsa" followed by a five-digit number, e.g., hsa00010 for Glycolysis)

### Output

- **Progress log** displaying analysis steps and any errors encountered.
- **Interactive bar plot** showing interacting transcription factors for each gene in the pathway
- **Table** displaying detailed results, including gene names, symbols, interacting transcription factors, and their counts.
- **CSV file** saved in the current directory: ```{pathway_id}_genes_with_interacting_TFs.csv```

## Troubleshooting

- Ensure stable internet connection
- Check KEGG pathway ID format
- Verify library versions match the specified requirements

## Limitations

- May take some time depending on the number of genes in the pathway
- Relies on external API services (KEGG, Enrichr)
- Rate-limited by external API call restrictions
- Limited to human pathways ('hsa' prefix)

## Acknowledgments

- KEGG Database
- Enrichr
- ENCODE Project
- Panel and Plotly for visualization tools