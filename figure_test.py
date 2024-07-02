import os
from Bio import SeqIO
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# Create a figure and axes
fig, ax = plt.subplots(figsize=(10, len(os.listdir('outputfiles/sliced_gbk')) * 2))
ax.set_xlim(0, 2000)  # Adjust as needed based on your sequences
ax.set_ylim(0, len(os.listdir('outputfiles/sliced_gbk')))

y_pos = 0  # Initial y-position for the first sequence

# Iterate over each GenBank file in the directory
for filename in os.listdir('outputfiles/sliced_gbk'):
    if filename.endswith('.gbk'):
        filepath = os.path.join('outputfiles/sliced_gbk', filename)
        record = SeqIO.read(filepath, 'genbank')
        
        # Plot a line for the sequence
        ax.plot([0, len(record)], [y_pos, y_pos], 'k-', lw=2)
        
        # Display the ACCESSION of the parent sequence
        ax.text(-100, y_pos, record.annotations['accessions'][0], verticalalignment='center')
        
        # Add CDS features as filled blue boxes
        for feature in record.features:
            if feature.type == "CDS":
                start = feature.location.start
                end = feature.location.end
                width = end - start
                rect = patches.Rectangle((start, y_pos - 0.1), width, 0.2, linewidth=1, edgecolor='blue', facecolor='blue')
                ax.add_patch(rect)
                
                # Display protein ID in white text inside the box
                protein_id = feature.qualifiers.get('protein_id', [''])[0]
                ax.text(start + width / 2, y_pos, protein_id, horizontalalignment='center', verticalalignment='center', color='white', fontsize=6)
        
        y_pos += 1  # Increment y-position for the next sequence

# Remove y-axis and ticks
ax.set_yticklabels([])
ax.set_yticks([])
ax.set_xlabel('Sequence Position')
ax.set_title('GenBank Sequences and Features')

# Save the figure
plt.tight_layout()
plt.savefig('sequences_features.png', dpi=300)