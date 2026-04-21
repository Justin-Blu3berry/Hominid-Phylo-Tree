Last edited: 04/21/2026 12:00a. EST
Current status: three sections left to write in their entireties
# What Went Right *NOT DONE*
- Specific technical and process-related successes (e.g., design decisions that paid off, effective testing strategies, fruitful data choices).
- What were the personal learning successes that you experienced?
# What Went Wrong or Was Hard *NOT DONE*
- Concrete challenges or missteps (e.g., underestimated complexity, unstable training, bad initial data).
- What would you do differently if you started over?
- What were your personal learning obstacles during this project? How did you attempt to address them, and how would you fix them if you couldn't?
# Algorithmic Lessons *NOT DONE*
- Insights about your chosen algorithm/class in real use.
- How well does it fit the problem?
- Tradeoffs between accuracy, complexity, and implementation complexity.
- Any surprises compared to how it was presented in the lecture or text?
# Future Directions
- Adding features that were cut for time
    - `main.py` checking if a gene doesn't have a fasta in `/data/` 
    - `main.py` checking if a species is missing from one or more of the files
    - having `main.py` attempt to download the missing data identified in the two steps above
- Testing another NCBI database 
    - nucleotides was convenient because of the format of the data it gave back, but its searches were weirdly inflexible
    - Gene might be more appropriate than nucleotide
    - Due to lack of experience with Entrez, I want to spend some time playing around with it to compare the handles it gets with `efetch` from different databases
- Comparisons with other tools as an evaluation method
- Making the construction of the "mean" matrix more robust against any species missing species by figuring out what row/column each species correlates to on each distance matrix. This would allow a species that is missing from one gene's distance matrix to ignore that one gene's distance matrix when calculating the mean for each of its interspecies distances. 
# Generative AI Disclosure (If Used)
As before, include tool/version, prompts, and a transparent account of how you used generative AI, if applicable.