# DNA Melting Temperature Analyzer

Calculate and visualize the melting temperature (Tm) of DNA sequences.

## Background
Tm is the temperature at which half of a DNA duplex dissociates. It depends on:
- Sequence length  
- GC content  
- Ion concentration  

**Formulas used:**
- Simple: `Tm = 4(G + C) + 2(A + T)` °C  
- Salt-corrected: `Tm = 81.5 + 16.6*log10([Na+]) + 0.41*(%GC) – 600/length`

## Features
1. **Tm Calculation**
   - Input: DNA sequence  
   - Output: Melting temperature  
2. **Sliding Window Analysis**
   - Window: 9 bases, step: 1  
   - Generates a vector `p` with Tm values  
   - Plots Tm along the sequence  
3. **Threshold Visualization**
   - Set cutoff values  
   - Horizontal bars displayed on the plot  
4. **GUI Interface**
   - Enter or load DNA sequences  
   - View plots and Tm values interactively  

## Technologies
- Python 3.x  
- Tkinter  
- Matplotlib  


