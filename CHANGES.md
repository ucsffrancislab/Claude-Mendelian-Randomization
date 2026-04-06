# Audit-Driven Changes to MR Pipeline
# Applied to: ucsffrancislab/Claude-Mendelian-Randomization

## Summary of Changes

Six targeted fixes were applied to run_bidirectional_mr.R.
No changes were needed to run_hla_sensitivity.R, format_icvf_for_mr.R,
or mr_postprocessing.R (though mr_postprocessing.R should be updated
separately to consume the new output files).

All changes are marked with "# AUDIT FIX" comments in the code.


## Change 1: Fix PGS001458 SNP ID format detection (local_clump function)

ISSUE:  PGS001458 (cingulum cingulate gyrus R) was silently dropped because
        local_clump() checked only the FIRST SNP to decide format. If that
        SNP happened to be chr:pos:ref:alt instead of rsID, the function
        treated ALL SNPs as chr:pos — and the chr:pos:ref:alt format didn't
        match the .bim file's chr:pos lookup, causing 0/7170 to map.

FIX:    Check the MAJORITY of SNPs, not just the first one.

OLD:    is_rsid <- grepl("^rs", original_snps[1])
NEW:    n_rsid <- sum(grepl("^rs", original_snps))
        is_rsid <- n_rsid > length(original_snps) * 0.5

IMPACT: PGS001458 will now be included in forward MR (was completely missing).
        This is the bilateral counterpart of our strongest result (PGS001457).


## Change 2: MR-PRESSO NbDistribution increased from 1,000 to 10,000

ISSUE:  MR-PRESSO failed for 16/17 significant tests with "Not enough elements
        to compute empirical P-values, increase NbDistribution". The default
        1,000 permutations was insufficient for the number of instruments.

FIX:    NbDistribution = 10000

IMPACT: MR-PRESSO will now produce results for all tests. Runtime increases
        ~10x for the MR-PRESSO step (but this runs sequentially anyway).


## Change 3: F-statistic calculation added to run_single_mr()

ISSUE:  No instrument strength assessment was performed. Weak instruments
        (F < 10) bias MR estimates.

FIX:    After harmonization, calculate F = (beta.exposure/se.exposure)^2
        for each instrument. Report mean, median, min, and count of weak
        instruments (F < 10) per tract.

NEW OUTPUT: forward_f_statistics.tsv, reverse_f_statistics.tsv


## Change 4: Reverse MR Steiger results now saved to file

ISSUE:  run_single_mr() already computed Steiger for reverse MR, and
        collect_results() already collected it, but the save block for
        reverse results did NOT write reverse_steiger.tsv.

FIX:    Added fwrite for rev$steiger to the reverse save block.

NEW OUTPUT: reverse_steiger.tsv (was missing entirely)


## Change 5: F-statistics collected in collect_results()

ISSUE:  The new f_stats field from run_single_mr() (Change 3) needed to be
        collected alongside mr, het, pleio, and steiger results.

FIX:    Added fstat_list to collect_results() with the same pattern as
        the other result types.


## Change 6: F-statistics saved for both directions

FIX:    Added fwrite calls for forward_f_statistics.tsv and
        reverse_f_statistics.tsv after the other save calls.


## Files NOT Modified

format_icvf_for_mr.R:   No changes needed. The BIG40 stats33k files have
                          rsIDs in the SNP column. The root cause of the
                          PGS001458 issue was in local_clump(), not here.

run_hla_sensitivity.R:   No changes needed. HLA exclusion works correctly.
                          NOTE: If you want HLA sensitivity to also benefit
                          from the local_clump fix, copy the updated
                          local_clump() function from run_bidirectional_mr.R.

mr_postprocessing.R:     No changes applied here, but SHOULD be updated to:
                          - Read and display reverse_steiger.tsv
                          - Read and display forward/reverse_f_statistics.tsv
                          - Update confidence ratings to use F-stat info


## New Output Files (after running the updated script)

  forward_f_statistics.tsv    - F-stats per forward MR tract x subtype
  reverse_f_statistics.tsv    - F-stats per reverse MR tract x subtype
  reverse_steiger.tsv         - Steiger directionality for reverse MR
  forward_mrpresso.tsv        - MR-PRESSO (now with NbDist=10000)
