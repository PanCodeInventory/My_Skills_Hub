#!/bin/bash
# Convert BED/narrowPeak format to RIdeogram label format

set -e

if [ $# -lt 2 ]; then
    echo "Usage: $0 <input.bed> <output.txt> [color_mode]"
    echo "  color_mode: signal (default), binary, or fixed"
    exit 1
fi

INPUT=$1
OUTPUT=$2
COLOR_MODE=${3:-signal}

echo "Converting $INPUT to RIdeogram label format..."

COLS=$(head -1 "$INPUT" | awk '{print NF}')

if [ "$COLS" -ge 10 ]; then
    echo "Detected narrowPeak format (10+ columns)"
    
    case $COLOR_MODE in
        binary)
            awk -v OFS='\t' 'BEGIN {
                print "Type", "Shape", "Chr", "Start", "End", "color"
            }
            {
                chrom = $1
                gsub("chr", "", chrom)
                start = $2
                end = $3
                qval = $9
                
                if (qval > 10) {
                    color = "d73027"
                } else {
                    color = "fee08b"
                }
                
                print "marker", "box", chrom, start, end, color
            }' "$INPUT" > "$OUTPUT"
            ;;
        
        fixed)
            awk -v OFS='\t' 'BEGIN {
                print "Type", "Shape", "Chr", "Start", "End", "color"
            }
            {
                chrom = $1
                gsub("chr", "", chrom)
                start = $2
                end = $3
                
                print "marker", "box", chrom, start, end, "ff7f00"
            }' "$INPUT" > "$OUTPUT"
            ;;
        
        *)
            awk -v OFS='\t' 'BEGIN {
                print "Type", "Shape", "Chr", "Start", "End", "color"
            }
            {
                chrom = $1
                gsub("chr", "", chrom)
                start = $2
                end = $3
                signal = $7
                
                if (signal > 100) {
                    color = "d73027"
                } else if (signal > 50) {
                    color = "fd8d3c"
                } else {
                    color = "4575b4"
                }
                
                print "marker", "box", chrom, start, end, color
            }' "$INPUT" > "$OUTPUT"
            ;;
    esac
    
elif [ "$COLS" -ge 6 ]; then
    echo "Detected BED format (6+ columns)"
    
    awk -v OFS='\t' 'BEGIN {
        print "Type", "Shape", "Chr", "Start", "End", "color"
    }
    {
        chrom = $1
        gsub("chr", "", chrom)
        start = $2
        end = $3
        score = $5
        
        if (score > 100) {
            color = "d73027"
        } else {
            color = "ff7f00"
        }
        
        print "marker", "box", chrom, start, end, color
    }' "$INPUT" > "$OUTPUT"
    
else
    echo "Error: Unrecognized format (expected 6+ columns, got $COLS)"
    exit 1
fi

LINES=$(wc -l < "$OUTPUT")
echo "Converted $((LINES - 1)) peaks to $OUTPUT"
echo "Output format: Type, Shape, Chr, Start, End, color"
