#!/bin/bash

#INPUT="/data/baca/users/sc1238/datasets/scripts/for_alexis/scripts_redo/resultsFINAL/results_H1.0_20perTSS/h1_tss_profile_profiles.csv"
#INPUT="/data/baca/users/sc1238/datasets/scripts/for_alexis/scripts_redo/resultsFINAL/results_H1.0_50perTSS/h1_tss_profile_profiles.csv"
#INPUT="/data/baca/users/sc1238/datasets/scripts/for_alexis/scripts_redo/scriptsFINAL/resultsFINAL/results_H1.3/h1_tss_profile_profiles.csv"
#INPUT="/data/baca/users/sc1238/datasets/scripts/for_alexis/scripts_redo/scriptsFINAL/resultsFINAL/results_H1.3_V2/h1_tss_profile_profiles.csv"
#INPUT="/data/baca/users/sc1238/datasets/scripts/for_alexis/scripts_redo/scriptsFINAL/resultsFINAL/results_H1.0/h1_tss_profile_profiles.csv"
#INPUT="/data/baca/users/sc1238/datasets/scripts/for_alexis/scripts_redo/scriptsFINAL/resultsFINAL/results_H1.4/h1_tss_profile_profiles.csv"

INPUT="./resultsFINAL/results_H1.2.tcgaGtexGenes/h1_tss_profile_profiles.csv"
INPUT="./resultsFINAL/results_H1.3.tcgaGtexGenes/h1_tss_profile_profiles.csv"
INPUT="./resultsFINAL/results_H1.0.tcgaGtexGenes/h1_tss_profile_profiles.csv"
INPUT="./resultsFINAL/results_H1.4.tcgaGtexGenes/h1_tss_profile_profiles.csv"
INPUT="./resultsFINAL/results_H1.5.tcgaGtexGenes/h1_tss_profile_profiles.csv"
INPUT="./resultsFINAL/results_H1.X.tcgaGtexGenes/h1_tss_profile_profiles.csv"

INPUT="./resultsFINAL/results_H1.3.eRNA/h1_tss_profile_profiles.csv"
INPUT="./resultsFINAL/results_H1.0.eRNA/h1_tss_profile_profiles.csv"
INPUT="./resultsFINAL/results_H1.2.eRNA/h1_tss_profile_profiles.csv"
INPUT="./resultsFINAL/results_H1.4.eRNA/h1_tss_profile_profiles.csv"
INPUT="./resultsFINAL/results_H1.5.eRNA/h1_tss_profile_profiles.csv"
INPUT="./resultsFINAL/results_H1.X.eRNA/h1_tss_profile_profiles.csv"

#BASE_DIR="/data/baca/users/sc1238/datasets/scripts/for_alexis/scripts_redo/resultsFINAL/smoothing_comparison/H1.0.20perTss"
#BASE_DIR="/data/baca/users/sc1238/datasets/scripts/for_alexis/scripts_redo/resultsFINAL/smoothing_comparison/H1.3.50perTss"
#BASE_DIR="/data/baca/users/sc1238/datasets/scripts/for_alexis/scripts_redo/scriptsFINAL/resultsFINAL/results_H1.3/smoothing_comparison/H1.3"
#BASE_DIR="/data/baca/users/sc1238/datasets/scripts/for_alexis/scripts_redo/scriptsFINAL/resultsFINAL/results_H1.3_V2/smoothing_comparison/H1.3"
#BASE_DIR="/data/baca/users/sc1238/datasets/scripts/for_alexis/scripts_redo/scriptsFINAL/resultsFINAL/results_H1.0/smoothing_comparison/H1.0"
#BASE_DIR="/data/baca/users/sc1238/datasets/scripts/for_alexis/scripts_redo/scriptsFINAL/resultsFINAL/results_H1.4/smoothing_comparison/H1.0"

BASE_DIR="./resultsFINAL/results_H1.2.tcgaGtexGenes/smoothing_comparison"
BASE_DIR="./resultsFINAL/results_H1.3.tcgaGtexGenes/smoothing_comparison"
BASE_DIR="./resultsFINAL/results_H1.0.tcgaGtexGenes/smoothing_comparison"
BASE_DIR="./resultsFINAL/results_H1.4.tcgaGtexGenes/smoothing_comparison"
BASE_DIR="./resultsFINAL/results_H1.5.tcgaGtexGenes/smoothing_comparison"
BASE_DIR="./resultsFINAL/results_H1.X.tcgaGtexGenes/smoothing_comparison"

BASE_DIR="./resultsFINAL/results_H1.3.eRNA/smoothing_comparison"
BASE_DIR="./resultsFINAL/results_H1.0.eRNA/smoothing_comparison"
BASE_DIR="./resultsFINAL/results_H1.2.eRNA/smoothing_comparison"
BASE_DIR="./resultsFINAL/results_H1.4.eRNA/smoothing_comparison"
BASE_DIR="./resultsFINAL/results_H1.5.eRNA/smoothing_comparison"
BASE_DIR="./resultsFINAL/results_H1.X.eRNA/smoothing_comparison"

# Default parameter values
CONDA_ENV="genometools_env"

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate $CONDA_ENV

# Create base directory
mkdir -p "$BASE_DIR"

echo "====================================================="
echo "BATCH PROCESSING TSS PROFILE SMOOTHING"
echo "Input: $INPUT"
echo "Base output directory: $BASE_DIR"
echo "====================================================="

# Create base directory
mkdir -p "$BASE_DIR"
echo "Created base directory: $BASE_DIR"

# Try different methods and parameters
for method in savgol gaussian moving_avg; do
  echo ""
  echo "====================================================="
  echo "PROCESSING METHOD: $method"
  echo "====================================================="
  
  # Create method directory
  METHOD_DIR="$BASE_DIR/$method"
  mkdir -p "$METHOD_DIR"
  echo "Created method directory: $METHOD_DIR"
  
  if [ "$method" == "gaussian" ]; then
    # Test different sigma values for Gaussian
    for sigma in 1 2 3 5 8 15; do
      echo ""
      echo "-------------------------------------------------"
      echo "  Processing $method with sigma = $sigma"
      echo "-------------------------------------------------"
      
      OUTPUT_DIR="$METHOD_DIR/sigma_$sigma"
      mkdir -p "$OUTPUT_DIR"
      echo "  Created output directory: $OUTPUT_DIR"
      
      python smooth_tss_profiles.py \
        --input "$INPUT" \
        --output-dir "$OUTPUT_DIR" \
        --prefix "h1_${method}_${sigma}" \
        --method "$method" \
        --smooth-sigma "$sigma"
    done
  else
    # Test different window sizes for Savgol and Moving Average
    for window in 11 21 31 51 101 251 501; do
      echo ""
      echo "-------------------------------------------------"
      echo "  Processing $method with window = $window"
      echo "-------------------------------------------------"
      
      OUTPUT_DIR="$METHOD_DIR/window_$window"
      mkdir -p "$OUTPUT_DIR"
      echo "  Created output directory: $OUTPUT_DIR"
      
      if [ "$method" == "savgol" ]; then
        # For Savgol, also vary polynomial order
        for poly in 2 3 4; do
          echo ""
          echo "  > Using polynomial order = $poly"
          
          POLY_DIR="$OUTPUT_DIR/poly_$poly"
          mkdir -p "$POLY_DIR"
          echo "    Created polynomial directory: $POLY_DIR"
          
          python smooth_tss_profiles.py \
            --input "$INPUT" \
            --output-dir "$POLY_DIR" \
            --prefix "h1_${method}_w${window}_p${poly}" \
            --method "$method" \
            --smooth-window "$window" \
            --smooth-poly "$poly"
        done
      else
        # For moving average, just use the window
        python smooth_tss_profiles.py \
          --input "$INPUT" \
          --output-dir "$OUTPUT_DIR" \
          --prefix "h1_${method}_w${window}" \
          --method "$method" \
          --smooth-window "$window"
      fi
    done
  fi
done

echo ""
echo "====================================================="
echo "ALL SMOOTHING VARIATIONS COMPLETE"
echo "Results in: $BASE_DIR"
echo "====================================================="
