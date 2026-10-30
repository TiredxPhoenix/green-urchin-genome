#!/bin/bash

echo "Scanning user scripts and updating data paths..."

# Find python files, explicitly excluding third-party tool folders like gmes_linux_64_4
PYTHON_FILES=$(find . -maxdepth 2 -name "*.py" ! -name "fix_paths.sh" ! -path "*/gmes_linux_64_4/*")

for py in $PYTHON_FILES; do
    # Loop through actual files in the workspace
    for f in *.tsv *.gff *.gff3 *.fna *.fa *.fai *.bam *.bai *.mmi; do
        [ -e "$f" ] || continue
        
        case "$f" in
            *.tsv|*.gff|*.gff3) folder="annotations/" ;;
            *.fna|*.fa|*.fai)   folder="genomes/" ;;
            *.bam|*.bai)        folder="alignments/" ;;
            *.mmi)              folder="indices/" ;;
            *)                  continue ;;
        esac
        
        # Safe sed replacement supporting both GNU and BSD sed
        if [[ "$OSTYPE" == "darwin"* ]]; then
            sed -i '' "s|\"\(.*/\)?$f\"|\"${folder}$f\"|g" "$py"
            sed -i '' "s|'\(.*/\)?$f'|'${folder}$f'|g" "$py"
        else
            sed -i "s|\"\(.*/\)?$f\"|\"${folder}$f\"|g" "$py"
            sed -i "s|'\(.*/\)?$f'|'${folder}$f'|g" "$py"
        fi
    done
    echo "Checked: $py"
done

echo "Path updating complete!"
