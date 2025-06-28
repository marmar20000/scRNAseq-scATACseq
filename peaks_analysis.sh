#!/bin/bash

# Define clusters to process
clusters=("Stem.I" "Stem.II" "Stem.III" "Progenitor.I")

for cluster in "${clusters[@]}"; do
  echo "Processing cluster: $cluster"

  # Define BED input files
  KO_BED="${cluster}_KO-reproduciblePeaks.gr.bed"
  WT_BED="${cluster}_WT-reproduciblePeaks.gr.bed"

  # Check BED files exist
  if [[ ! -f "$KO_BED" || ! -f "$WT_BED" ]]; then
    echo "⚠️ BED files missing for $cluster. Skipping."
    continue
  fi

  # Find corresponding bigWig files dynamically
  KO_BW=$(ls ${cluster}_KO-*.bw 2>/dev/null | head -n 1)
  WT_BW=$(ls ${cluster}_WT-*.bw 2>/dev/null | head -n 1)

  if [[ -z "$KO_BW" || -z "$WT_BW" ]]; then
    echo "⚠️ bigWig files missing for $cluster. Skipping."
    continue
  fi

  echo "Using KO bigWig: $KO_BW"
  echo "Using WT bigWig: $WT_BW"

  # Sort and deduplicate BEDs
  KO_SORTED="${cluster}_KO_sorted.bed"
  WT_SORTED="${cluster}_WT_sorted.bed"

  sort -k1,1 -k2,2n "$KO_BED" | uniq > "$KO_SORTED"
  sort -k1,1 -k2,2n "$WT_BED" | uniq > "$WT_SORTED"

  # Find common and unique peaks
  bedtools intersect -a "$KO_SORTED" -b "$WT_SORTED" -u > "${cluster}_common.bed"
  bedtools subtract -a "$KO_SORTED" -b "$WT_SORTED" -A > "${cluster}_unique_KO.bed"
  bedtools subtract -a "$WT_SORTED" -b "$KO_SORTED" -A > "${cluster}_unique_WT.bed"


  # Generate computeMatrix
  computeMatrix reference-point \
    --referencePoint center \
    -b 2000 -a 2000 \
    -R "${cluster}_common.bed" "${cluster}_unique_KO.bed" "${cluster}_unique_WT.bed" \
    -S "$KO_BW" "$WT_BW" \
    --skipZeros \
    --missingDataAsZero \
    --sortRegions descend \
    --outFileName "${cluster}_2000_matrix.gz" \
    --outFileSortedRegions "${cluster}_sorted_regions.bed"\
    -p max/2

  # Plot heatmap
  plotHeatmap \
    -m "${cluster}_2000_matrix.gz" \
    -out "${cluster}_heatmap.png" \
    --colorMap RdBu \
    --regionsLabel "Common" "KO-unique" "WT-unique" \
    --samplesLabel "KO" "WT" \
    --refPointLabel "Peak Center" \
    --zMin 0\
    --heatmapWidth 9

#    -p max/2

  echo "✅ Done with $cluster"
done
