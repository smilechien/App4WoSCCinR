# App4WoSCCinR

A web-based Shiny application for streamlined bibliometric analysis using **Web of Science Core Collection (WoSCC)** metadata.

## Overview

**App4WoSCCinR** is designed for users who want to perform bibliometric analysis from standard WoSCC exports without writing code. The app converts WoS metadata into structured intermediate tables, generates a summary report across major metadata domains, and provides interactive visualizations for exploratory analysis.

The platform supports a workflow from **data upload or URL loading** to **processed tables, summary reports, network-style visualizations, and downloadable outputs**.

## Features

- Upload **WoSCC Excel files** (`.xls` / `.xlsx`)
- Load data from a **URL** (`.xls`, `.xlsx`, `.csv`, `.txt`)
- Optional **Demo** mode for quick testing
- Built-in preprocessing pipeline:
  - `wide32`
  - `long32`
  - `long16`
  - `metatable`
- Summary report across 10 metadata domains
- Dominance analysis using **Absolute Advantage Coefficient (AAC)**
- Interactive and report-style outputs, including:
  - Summary report
  - AAC panel
  - Network visualization
  - Chord diagram
  - Kano plot
  - Lotka’s law analysis
  - Most-cited article table
  - HTML report preview/download
- Download outputs as:
  - **PNG**
  - **ZIP**

## Supported Input

### Local upload
Upload a WoSCC export file in:

- `.xls`
- `.xlsx`

### URL input
Load a supported file directly from a URL:

- `.xls`
- `.xlsx`
- `.csv`
- `.txt`

### Demo
Run the built-in demo dataset for testing the interface and outputs.

## Access / CMC

The app includes a **CMC / trial code** field.

- A valid **10-digit CMC** may be required for full access
- The value `test` can be used for trial access
- URL-based runs are also subject to CMC validation

> Replace this section if your public GitHub version does not require access control.

## Workflow

1. Enter a valid **CMC / trial code** if needed
2. Upload a local WoSCC file or provide a URL
3. Click **Run WoS (Upload)** or **Run URL**
4. Review the generated outputs in the app tabs
5. Download the summary figure or full ZIP results

## Output Files

The app generates and exposes processed files such as:

- `woswide32.csv`
- `woslong32.csv`
- `woslong16.csv`
- `wosmetatable.csv`
- `summary_report.csv`
- `summary_report.png`

## Main Analysis Panels

Depending on the run, the app provides tabs such as:

- **ReadMe**
- **Summary**
- **AAC**
- **woswide32**
- **woslong32**
- **woslong16**
- **metatable**
- **Network**
- **Chord**
- **Author**
- **Kano**
- **Lotka**
- **Most cited**
- **Report**
- **IP**

These panels support both quick reporting and deeper bibliometric exploration.

## Typical Use Cases

App4WoSCCinR is useful for:

- Researchers preparing a quick bibliometric overview
- Educators demonstrating metadata-based bibliometric methods
- Analysts comparing publication patterns across journals, authors, institutions, countries, or keywords
- Users who need a more transparent and reproducible workflow than manual spreadsheet processing

## Installation

Clone the repository and run the app in **R**:

```r
shiny::runApp()
```

Or specify the app directory:

```r
shiny::runApp("path/to/App4WoSCCinR")
```

## Requirements

Recommended R packages include:

- `shiny`
- `DT`
- `dplyr`
- `stringr`
- `tidyr`
- `readxl`
- `ggplot2`
- `htmltools`
- `base64enc`
- `igraph`
- `visNetwork`

Additional optional modules may be used if present in the app directory.

## Example Citation

Chien TW, Yang S, Chou W, Wang WC. *App4WoSCCinR: A Web-Based Bibliometric Analysis Platform for the Web of Science Core Collection*. 2026.

## Availability

- **Web app**: `https://smilechien.shinyapps.io/woscc/`
- **GitHub**: `https://github.com/<your-account>/<your-repo>`

## Notes

- This app is intended for **WoSCC-based bibliometric workflows**
- For best compatibility, use standard WoS export formats
- If access control is not needed in your GitHub release, you can simplify or remove the CMC section
