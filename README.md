[![Review Assignment Due Date](https://classroom.github.com/assets/deadline-readme-button-22041afd0340ce965d47ae6ef1cefeee28c7c493a6346c4f15d667ab976d596c.svg)](https://classroom.github.com/a/pWmxMLzQ)
# pro3.global|local-Aln
Khoo Kim Jun 110304027

## Description

* Write a Python script to perform a global or local alignment.
* Creating your own script, i.e. hw3.py.
* In this program, library Biostrings is only used to parse input fasta file.
* Packages you can use: numpy, pandas, Bio
* You should write a program with a function named alignment, ie.
```

def alignment(input_path, score_path, output_path, aln, gap):
    import pandas as pd
    import numpy as np
    from Bio import SeqIO

    # Read sequences from input
    sequences = list(SeqIO.parse(input_path, "fasta"))
    if len(sequences) != 2:
        raise ValueError("The input file must contain exactly two sequences.")

    seq1_id, seq1 = sequences[0].id, str(sequences[0].seq)
    seq2_id, seq2 = sequences[1].id, str(sequences[1].seq)

    # Load scoring matrix
    score_matrix_df = pd.read_csv(score_path, sep='\\s+', skiprows=9, index_col=0)

    # Initialize DP matrix
    n = len(seq1) + 1
    m = len(seq2) + 1
    dp = np.zeros((n, m), dtype=int)

    # Initialize DP matrix for global alignment
    if aln == "global":
        for i in range(n):
            dp[i][0] = i * gap
        for j in range(m):
            dp[0][j] = j * gap

    if aln == "local":
        dp[0, :] = 0
        dp[:, 0] = 0
    
    max_score = 0
    max_positions = []

    # Fill the DP matrix
    for i in range(1, n):
        for j in range(1, m):
            match = dp[i - 1][j - 1] + score_matrix_df.loc[seq1[i - 1], seq2[j - 1]]
            delete = dp[i - 1][j] + gap
            insert = dp[i][j - 1] + gap
            if aln == "global":
                dp[i][j] = max(match, delete, insert)
            elif aln == "local":
                dp[i][j] = max(0, match, delete, insert)
                if dp[i][j] > max_score:
                    max_score = dp[i][j]
                    max_positions = [(i, j)]
                elif dp[i][j] == max_score:
                    max_positions.append((i, j))

    # Traceback
    def trace_back(i, j):
        aligned_seq1 = []
        aligned_seq2 = []

        while i > 0 and j > 0:
            current_score = dp[i][j]
            if aln == "local" and current_score == 0:

                # Include the next character from both sequences, if valid
                if i > 0 and j > 0:
                    aligned_seq1.append(seq1[i - 1])
                    aligned_seq2.append(seq2[j - 1])
                break
                
            if  current_score == dp[i - 1][j - 1] + score_matrix_df.loc[seq1[i - 1], seq2[j - 1]]:
                aligned_seq1.append(seq1[i - 1])
                aligned_seq2.append(seq2[j - 1])
                i -= 1
                j -= 1
                
            elif current_score == dp[i - 1][j] + gap:
                aligned_seq1.append(seq1[i - 1])
                aligned_seq2.append("-")
                i -= 1
                
            elif current_score == dp[i][j - 1] + gap:
                aligned_seq1.append("-")
                aligned_seq2.append(seq2[j - 1])
                j -= 1

         # Include the last valid character if the traceback stopped prematurely
        if aln == "local" and dp[i][j] > 0:
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append(seq2[j - 1])

        
        aligned_seq1 = ''.join(reversed(aligned_seq1))
        aligned_seq2 = ''.join(reversed(aligned_seq2))
        return aligned_seq1, aligned_seq2

    alignments = []
    if aln == "global":
        alignments.append(trace_back(n - 1, m - 1))
    elif aln == "local":
        for i, j in max_positions:
            alignments.append(trace_back(i, j))

    # Write output to file
    with open(output_path, "w") as f:
        if aln == "global":
            aligned_seq1, aligned_seq2 = alignments[0]
            f.write(f">{seq1_id}\n{aligned_seq1}\n")
            f.write(f">{seq2_id}\n{aligned_seq2}\n")
        elif aln == "local":
            for aligned_seq1, aligned_seq2 in alignments:
                f.write(f">{seq1_id}\n{aligned_seq1}\n")
                f.write(f">{seq2_id}\n{aligned_seq2}\n")

```
* If there is more than one local alignment with the same highest score, you should output local alignments with the maximum length. 
* If there is more than one local alignment with the same highest score, you should output those local alignments in string sequential order according to protein1 and then protein2, i.e., 
  ```
  >protein1
  local alignment1
  >protein2
  local alignment1
  >protein1
  local alignment2
  >protein2
  local alignment2
  ```
## Parameters

* input: .fasta file (ex. test_global.fasta)
* score: score file (ex. pam250.txt)
* aln: global|local
* gap: gap score
* output: .fasta file

## Files

* hw3_ref.py: You can start from this reference code, and try to write your own comment in English.
* pam100.txt
* pam250.txt
* test_global.fasta
* result_global.fasta: You should output your alignment in FASTA format.
* test_local.fasta
* result_local.fasta
## Command

Executing your code with the following command.


```Python
alignment("examples/test_global.fasta", "examples/pam250.txt", "examples/result_global.fasta", "global", -10)
alignment("examples/test_local.fasta", "examples/pam100.txt", "examples/result_local.fasta", "local", -10)
```

## Evaluation

10 testing data(5 public, 5 private)

The correct answer gets 10 points for each testing data.



### Penalty

* High code similarity to others: YOUR SCORE = 0

## References
* ChatGPT, respond to my prompt, on November 24, 2024.
* Conversation link: https://chatgpt.com/share/6742d5b3-08b8-8013-81e8-bcb902b897bd
* Below are some websites that I have copied the code and feed the GPT:
* https://stackoverflow.com/questions/36596389/sequence-alignment-algorithm-with-a-group-of-characters-instead-of-one-character
* https://github.com/biopython/biopython/issues/955






