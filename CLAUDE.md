# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

This is an R/Quarto research project focused on causal effect estimation of change from baseline measurements. The project demonstrates treatment effect analysis using simulated data and various statistical approaches for both experimental and observational study designs.

## Development Commands

### Building and Rendering
- **Render main document**: `quarto render Baseline_Change.qmd`
- **Preview with live reload**: `quarto preview`
- **Render to GitHub Pages**: `quarto publish gh-pages` (output goes to docs/ directory)

### R Project
- Open in RStudio: `open baselinechange.Rproj`
- Install required packages by running the setup chunk in any .qmd file

## Project Structure

### Core Analysis Files
- `Baseline_Change.qmd` - Main analysis document demonstrating causal effect estimation methods
- `Claude-Change-Baseline.html` - Alternative analysis approach 
- `Claude-Censoring-Demo.qmd` - Censoring analysis demonstration

### R Scripts
- `causaleffect.R` - Causal effect identification using the causaleffect package and DAGs
- `dosearch breskin.R` - DoSearch implementation for causal identification  
- `Breskin el al sim censoring.R` - Simulation code for censoring scenarios

### Configuration
- `quarto.yml` - Quarto project configuration (website type, cosmo theme)
- `baselinechange.Rproj` - RStudio project settings
- `references.bib` - Bibliography for citations

## Key Dependencies

The project primarily uses:
- **Quarto** for literate programming and document generation
- **R packages**: dagitty, ggdag, tidyverse, DT, gt, gtsummary, knitr, kableExtra, causaleffect, igraph
- **Causal inference tools**: dagitty for DAG analysis, causaleffect for identification

## Architecture Notes

The project follows a research notebook structure where:
1. Multiple QMD files explore different aspects of baseline change analysis
2. R scripts contain supporting functions and specific analyses
3. DAGs (Directed Acyclic Graphs) are used extensively for causal reasoning
4. Simulated datasets demonstrate various confounding and censoring scenarios
5. Output is published as HTML to GitHub Pages

## Publication

The rendered analysis is published at: https://tjohnson250.github.io/baselinechange/Baseline_Change.html