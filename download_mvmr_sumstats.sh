#!/bin/bash
# =============================================================================
# Download FA, MD, ISOVF, OD GWAS Summary Statistics from Oxford BIG40
# For multivariable MR alongside ICVF
# =============================================================================
# Source: Smith et al. (2021) Nature Neuroscience
# N = 33,224, GRCh37, 17.1M SNPs per file
# Each file is ~344 MB compressed
# Total: 72 files × ~344 MB = ~24.8 GB
# =============================================================================

OUTPUT_DIR="mvmr_gwas_sumstats"
mkdir -p "${OUTPUT_DIR}"

BASE_URL="https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k"

echo "=============================================="
echo "Downloading FA, MD, ISOVF, OD summary stats"
echo "for 18 white matter tracts (72 files)"
echo "=============================================="

# --- FA (Fractional Anisotropy) ---
echo ""
echo "=== FA (Fractional Anisotropy) - 18 tracts ==="

echo "  Downloading FA Middle_cerebellar_peduncle (TBSS, IDP 1614)..."
curl -s -o "${OUTPUT_DIR}/FA_TBSS_Middle_cerebellar_peduncle_IDP1614.txt.gz" "https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/1614.txt.gz"

echo "  Downloading FA Genu_of_corpus_callosum (TBSS, IDP 1616)..."
curl -s -o "${OUTPUT_DIR}/FA_TBSS_Genu_of_corpus_callosum_IDP1616.txt.gz" "https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/1616.txt.gz"

echo "  Downloading FA Body_of_corpus_callosum (TBSS, IDP 1617)..."
curl -s -o "${OUTPUT_DIR}/FA_TBSS_Body_of_corpus_callosum_IDP1617.txt.gz" "https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/1617.txt.gz"

echo "  Downloading FA Cerebral_peduncle_R (TBSS, IDP 1628)..."
curl -s -o "${OUTPUT_DIR}/FA_TBSS_Cerebral_peduncle_R_IDP1628.txt.gz" "https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/1628.txt.gz"

echo "  Downloading FA Posterior_limb_of_internal_capsule_L (TBSS, IDP 1633)..."
curl -s -o "${OUTPUT_DIR}/FA_TBSS_Posterior_limb_of_internal_capsule_L_IDP1633.txt.gz" "https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/1633.txt.gz"

echo "  Downloading FA Retrolenticular_part_of_internal_capsule_R (TBSS, IDP 1634)..."
curl -s -o "${OUTPUT_DIR}/FA_TBSS_Retrolenticular_part_of_internal_capsule_R_IDP1634.txt.gz" "https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/1634.txt.gz"

echo "  Downloading FA Retrolenticular_part_of_internal_capsule_L (TBSS, IDP 1635)..."
curl -s -o "${OUTPUT_DIR}/FA_TBSS_Retrolenticular_part_of_internal_capsule_L_IDP1635.txt.gz" "https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/1635.txt.gz"

echo "  Downloading FA Superior_corona_radiata_R (TBSS, IDP 1638)..."
curl -s -o "${OUTPUT_DIR}/FA_TBSS_Superior_corona_radiata_R_IDP1638.txt.gz" "https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/1638.txt.gz"

echo "  Downloading FA Superior_corona_radiata_L (TBSS, IDP 1639)..."
curl -s -o "${OUTPUT_DIR}/FA_TBSS_Superior_corona_radiata_L_IDP1639.txt.gz" "https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/1639.txt.gz"

echo "  Downloading FA Sagittal_stratum_R (TBSS, IDP 1644)..."
curl -s -o "${OUTPUT_DIR}/FA_TBSS_Sagittal_stratum_R_IDP1644.txt.gz" "https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/1644.txt.gz"

echo "  Downloading FA Sagittal_stratum_L (TBSS, IDP 1645)..."
curl -s -o "${OUTPUT_DIR}/FA_TBSS_Sagittal_stratum_L_IDP1645.txt.gz" "https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/1645.txt.gz"

echo "  Downloading FA Cingulum_cingulate_gyrus_R (TBSS, IDP 1648)..."
curl -s -o "${OUTPUT_DIR}/FA_TBSS_Cingulum_cingulate_gyrus_R_IDP1648.txt.gz" "https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/1648.txt.gz"

echo "  Downloading FA Cingulum_cingulate_gyrus_L (TBSS, IDP 1649)..."
curl -s -o "${OUTPUT_DIR}/FA_TBSS_Cingulum_cingulate_gyrus_L_IDP1649.txt.gz" "https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/1649.txt.gz"

echo "  Downloading FA Cingulum_hippocampus_R (TBSS, IDP 1650)..."
curl -s -o "${OUTPUT_DIR}/FA_TBSS_Cingulum_hippocampus_R_IDP1650.txt.gz" "https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/1650.txt.gz"

echo "  Downloading FA Cingulum_hippocampus_L (TBSS, IDP 1651)..."
curl -s -o "${OUTPUT_DIR}/FA_TBSS_Cingulum_hippocampus_L_IDP1651.txt.gz" "https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/1651.txt.gz"

echo "  Downloading FA ar_r (ProbtrackX, IDP 1789)..."
curl -s -o "${OUTPUT_DIR}/FA_ProbtrackX_ar_r_IDP1789.txt.gz" "https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/1789.txt.gz"

echo "  Downloading FA cgh_r (ProbtrackX, IDP 1795)..."
curl -s -o "${OUTPUT_DIR}/FA_ProbtrackX_cgh_r_IDP1795.txt.gz" "https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/1795.txt.gz"

echo "  Downloading FA fma (ProbtrackX, IDP 1798)..."
curl -s -o "${OUTPUT_DIR}/FA_ProbtrackX_fma_IDP1798.txt.gz" "https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/1798.txt.gz"

echo ""
echo "=== MD (Mean Diffusivity) - 18 tracts ==="

echo "  Downloading MD Middle_cerebellar_peduncle (TBSS, IDP 1662)..."
curl -s -o "${OUTPUT_DIR}/MD_TBSS_Middle_cerebellar_peduncle_IDP1662.txt.gz" "https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/1662.txt.gz"

echo "  Downloading MD Genu_of_corpus_callosum (TBSS, IDP 1664)..."
curl -s -o "${OUTPUT_DIR}/MD_TBSS_Genu_of_corpus_callosum_IDP1664.txt.gz" "https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/1664.txt.gz"

echo "  Downloading MD Body_of_corpus_callosum (TBSS, IDP 1665)..."
curl -s -o "${OUTPUT_DIR}/MD_TBSS_Body_of_corpus_callosum_IDP1665.txt.gz" "https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/1665.txt.gz"

echo "  Downloading MD Cerebral_peduncle_R (TBSS, IDP 1676)..."
curl -s -o "${OUTPUT_DIR}/MD_TBSS_Cerebral_peduncle_R_IDP1676.txt.gz" "https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/1676.txt.gz"

echo "  Downloading MD Posterior_limb_of_internal_capsule_L (TBSS, IDP 1681)..."
curl -s -o "${OUTPUT_DIR}/MD_TBSS_Posterior_limb_of_internal_capsule_L_IDP1681.txt.gz" "https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/1681.txt.gz"

echo "  Downloading MD Retrolenticular_part_of_internal_capsule_R (TBSS, IDP 1682)..."
curl -s -o "${OUTPUT_DIR}/MD_TBSS_Retrolenticular_part_of_internal_capsule_R_IDP1682.txt.gz" "https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/1682.txt.gz"

echo "  Downloading MD Retrolenticular_part_of_internal_capsule_L (TBSS, IDP 1683)..."
curl -s -o "${OUTPUT_DIR}/MD_TBSS_Retrolenticular_part_of_internal_capsule_L_IDP1683.txt.gz" "https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/1683.txt.gz"

echo "  Downloading MD Superior_corona_radiata_R (TBSS, IDP 1686)..."
curl -s -o "${OUTPUT_DIR}/MD_TBSS_Superior_corona_radiata_R_IDP1686.txt.gz" "https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/1686.txt.gz"

echo "  Downloading MD Superior_corona_radiata_L (TBSS, IDP 1687)..."
curl -s -o "${OUTPUT_DIR}/MD_TBSS_Superior_corona_radiata_L_IDP1687.txt.gz" "https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/1687.txt.gz"

echo "  Downloading MD Sagittal_stratum_R (TBSS, IDP 1692)..."
curl -s -o "${OUTPUT_DIR}/MD_TBSS_Sagittal_stratum_R_IDP1692.txt.gz" "https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/1692.txt.gz"

echo "  Downloading MD Sagittal_stratum_L (TBSS, IDP 1693)..."
curl -s -o "${OUTPUT_DIR}/MD_TBSS_Sagittal_stratum_L_IDP1693.txt.gz" "https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/1693.txt.gz"

echo "  Downloading MD Cingulum_cingulate_gyrus_R (TBSS, IDP 1696)..."
curl -s -o "${OUTPUT_DIR}/MD_TBSS_Cingulum_cingulate_gyrus_R_IDP1696.txt.gz" "https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/1696.txt.gz"

echo "  Downloading MD Cingulum_cingulate_gyrus_L (TBSS, IDP 1697)..."
curl -s -o "${OUTPUT_DIR}/MD_TBSS_Cingulum_cingulate_gyrus_L_IDP1697.txt.gz" "https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/1697.txt.gz"

echo "  Downloading MD Cingulum_hippocampus_R (TBSS, IDP 1698)..."
curl -s -o "${OUTPUT_DIR}/MD_TBSS_Cingulum_hippocampus_R_IDP1698.txt.gz" "https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/1698.txt.gz"

echo "  Downloading MD Cingulum_hippocampus_L (TBSS, IDP 1699)..."
curl -s -o "${OUTPUT_DIR}/MD_TBSS_Cingulum_hippocampus_L_IDP1699.txt.gz" "https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/1699.txt.gz"

echo "  Downloading MD ar_r (ProbtrackX, IDP 1816)..."
curl -s -o "${OUTPUT_DIR}/MD_ProbtrackX_ar_r_IDP1816.txt.gz" "https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/1816.txt.gz"

echo "  Downloading MD cgh_r (ProbtrackX, IDP 1822)..."
curl -s -o "${OUTPUT_DIR}/MD_ProbtrackX_cgh_r_IDP1822.txt.gz" "https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/1822.txt.gz"

echo "  Downloading MD fma (ProbtrackX, IDP 1825)..."
curl -s -o "${OUTPUT_DIR}/MD_ProbtrackX_fma_IDP1825.txt.gz" "https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/1825.txt.gz"

echo ""
echo "=== ISOVF (Isotropic Volume Fraction) - 18 tracts ==="

echo "  Downloading ISOVF Middle_cerebellar_peduncle (TBSS, IDP 1950)..."
curl -s -o "${OUTPUT_DIR}/ISOVF_TBSS_Middle_cerebellar_peduncle_IDP1950.txt.gz" "https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/1950.txt.gz"

echo "  Downloading ISOVF Genu_of_corpus_callosum (TBSS, IDP 1952)..."
curl -s -o "${OUTPUT_DIR}/ISOVF_TBSS_Genu_of_corpus_callosum_IDP1952.txt.gz" "https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/1952.txt.gz"

echo "  Downloading ISOVF Body_of_corpus_callosum (TBSS, IDP 1953)..."
curl -s -o "${OUTPUT_DIR}/ISOVF_TBSS_Body_of_corpus_callosum_IDP1953.txt.gz" "https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/1953.txt.gz"

echo "  Downloading ISOVF Cerebral_peduncle_R (TBSS, IDP 1964)..."
curl -s -o "${OUTPUT_DIR}/ISOVF_TBSS_Cerebral_peduncle_R_IDP1964.txt.gz" "https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/1964.txt.gz"

echo "  Downloading ISOVF Posterior_limb_of_internal_capsule_L (TBSS, IDP 1969)..."
curl -s -o "${OUTPUT_DIR}/ISOVF_TBSS_Posterior_limb_of_internal_capsule_L_IDP1969.txt.gz" "https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/1969.txt.gz"

echo "  Downloading ISOVF Retrolenticular_part_of_internal_capsule_R (TBSS, IDP 1970)..."
curl -s -o "${OUTPUT_DIR}/ISOVF_TBSS_Retrolenticular_part_of_internal_capsule_R_IDP1970.txt.gz" "https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/1970.txt.gz"

echo "  Downloading ISOVF Retrolenticular_part_of_internal_capsule_L (TBSS, IDP 1971)..."
curl -s -o "${OUTPUT_DIR}/ISOVF_TBSS_Retrolenticular_part_of_internal_capsule_L_IDP1971.txt.gz" "https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/1971.txt.gz"

echo "  Downloading ISOVF Superior_corona_radiata_R (TBSS, IDP 1974)..."
curl -s -o "${OUTPUT_DIR}/ISOVF_TBSS_Superior_corona_radiata_R_IDP1974.txt.gz" "https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/1974.txt.gz"

echo "  Downloading ISOVF Superior_corona_radiata_L (TBSS, IDP 1975)..."
curl -s -o "${OUTPUT_DIR}/ISOVF_TBSS_Superior_corona_radiata_L_IDP1975.txt.gz" "https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/1975.txt.gz"

echo "  Downloading ISOVF ar_r (ProbtrackX, IDP 1978)..."
curl -s -o "${OUTPUT_DIR}/ISOVF_ProbtrackX_ar_r_IDP1978.txt.gz" "https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/1978.txt.gz"

echo "  Downloading ISOVF Sagittal_stratum_R (TBSS, IDP 1980)..."
curl -s -o "${OUTPUT_DIR}/ISOVF_TBSS_Sagittal_stratum_R_IDP1980.txt.gz" "https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/1980.txt.gz"

echo "  Downloading ISOVF Sagittal_stratum_L (TBSS, IDP 1981)..."
curl -s -o "${OUTPUT_DIR}/ISOVF_TBSS_Sagittal_stratum_L_IDP1981.txt.gz" "https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/1981.txt.gz"

echo "  Downloading ISOVF Cingulum_cingulate_gyrus_R (TBSS, IDP 1984)..."
curl -s -o "${OUTPUT_DIR}/ISOVF_TBSS_Cingulum_cingulate_gyrus_R_IDP1984.txt.gz" "https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/1984.txt.gz"

echo "  Downloading ISOVF cgh_r (ProbtrackX, IDP 1984)..."
curl -s -o "${OUTPUT_DIR}/ISOVF_ProbtrackX_cgh_r_IDP1984.txt.gz" "https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/1984.txt.gz"

echo "  Downloading ISOVF Cingulum_cingulate_gyrus_L (TBSS, IDP 1985)..."
curl -s -o "${OUTPUT_DIR}/ISOVF_TBSS_Cingulum_cingulate_gyrus_L_IDP1985.txt.gz" "https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/1985.txt.gz"

echo "  Downloading ISOVF Cingulum_hippocampus_R (TBSS, IDP 1986)..."
curl -s -o "${OUTPUT_DIR}/ISOVF_TBSS_Cingulum_hippocampus_R_IDP1986.txt.gz" "https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/1986.txt.gz"

echo "  Downloading ISOVF Cingulum_hippocampus_L (TBSS, IDP 1987)..."
curl -s -o "${OUTPUT_DIR}/ISOVF_TBSS_Cingulum_hippocampus_L_IDP1987.txt.gz" "https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/1987.txt.gz"

echo "  Downloading ISOVF fma (ProbtrackX, IDP 1987)..."
curl -s -o "${OUTPUT_DIR}/ISOVF_ProbtrackX_fma_IDP1987.txt.gz" "https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/1987.txt.gz"

echo ""
echo "=== OD (Orientation Dispersion) - 18 tracts ==="

echo "  Downloading OD Middle_cerebellar_peduncle (TBSS, IDP 1998)..."
curl -s -o "${OUTPUT_DIR}/OD_TBSS_Middle_cerebellar_peduncle_IDP1998.txt.gz" "https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/1998.txt.gz"

echo "  Downloading OD Genu_of_corpus_callosum (TBSS, IDP 2000)..."
curl -s -o "${OUTPUT_DIR}/OD_TBSS_Genu_of_corpus_callosum_IDP2000.txt.gz" "https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/2000.txt.gz"

echo "  Downloading OD Body_of_corpus_callosum (TBSS, IDP 2001)..."
curl -s -o "${OUTPUT_DIR}/OD_TBSS_Body_of_corpus_callosum_IDP2001.txt.gz" "https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/2001.txt.gz"

echo "  Downloading OD ar_r (ProbtrackX, IDP 2005)..."
curl -s -o "${OUTPUT_DIR}/OD_ProbtrackX_ar_r_IDP2005.txt.gz" "https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/2005.txt.gz"

echo "  Downloading OD cgh_r (ProbtrackX, IDP 2011)..."
curl -s -o "${OUTPUT_DIR}/OD_ProbtrackX_cgh_r_IDP2011.txt.gz" "https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/2011.txt.gz"

echo "  Downloading OD Cerebral_peduncle_R (TBSS, IDP 2012)..."
curl -s -o "${OUTPUT_DIR}/OD_TBSS_Cerebral_peduncle_R_IDP2012.txt.gz" "https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/2012.txt.gz"

echo "  Downloading OD fma (ProbtrackX, IDP 2014)..."
curl -s -o "${OUTPUT_DIR}/OD_ProbtrackX_fma_IDP2014.txt.gz" "https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/2014.txt.gz"

echo "  Downloading OD Posterior_limb_of_internal_capsule_L (TBSS, IDP 2017)..."
curl -s -o "${OUTPUT_DIR}/OD_TBSS_Posterior_limb_of_internal_capsule_L_IDP2017.txt.gz" "https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/2017.txt.gz"

echo "  Downloading OD Retrolenticular_part_of_internal_capsule_R (TBSS, IDP 2018)..."
curl -s -o "${OUTPUT_DIR}/OD_TBSS_Retrolenticular_part_of_internal_capsule_R_IDP2018.txt.gz" "https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/2018.txt.gz"

echo "  Downloading OD Retrolenticular_part_of_internal_capsule_L (TBSS, IDP 2019)..."
curl -s -o "${OUTPUT_DIR}/OD_TBSS_Retrolenticular_part_of_internal_capsule_L_IDP2019.txt.gz" "https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/2019.txt.gz"

echo "  Downloading OD Superior_corona_radiata_R (TBSS, IDP 2022)..."
curl -s -o "${OUTPUT_DIR}/OD_TBSS_Superior_corona_radiata_R_IDP2022.txt.gz" "https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/2022.txt.gz"

echo "  Downloading OD Superior_corona_radiata_L (TBSS, IDP 2023)..."
curl -s -o "${OUTPUT_DIR}/OD_TBSS_Superior_corona_radiata_L_IDP2023.txt.gz" "https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/2023.txt.gz"

echo "  Downloading OD Sagittal_stratum_R (TBSS, IDP 2028)..."
curl -s -o "${OUTPUT_DIR}/OD_TBSS_Sagittal_stratum_R_IDP2028.txt.gz" "https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/2028.txt.gz"

echo "  Downloading OD Sagittal_stratum_L (TBSS, IDP 2029)..."
curl -s -o "${OUTPUT_DIR}/OD_TBSS_Sagittal_stratum_L_IDP2029.txt.gz" "https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/2029.txt.gz"

echo "  Downloading OD Cingulum_cingulate_gyrus_R (TBSS, IDP 2032)..."
curl -s -o "${OUTPUT_DIR}/OD_TBSS_Cingulum_cingulate_gyrus_R_IDP2032.txt.gz" "https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/2032.txt.gz"

echo "  Downloading OD Cingulum_cingulate_gyrus_L (TBSS, IDP 2033)..."
curl -s -o "${OUTPUT_DIR}/OD_TBSS_Cingulum_cingulate_gyrus_L_IDP2033.txt.gz" "https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/2033.txt.gz"

echo "  Downloading OD Cingulum_hippocampus_R (TBSS, IDP 2034)..."
curl -s -o "${OUTPUT_DIR}/OD_TBSS_Cingulum_hippocampus_R_IDP2034.txt.gz" "https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/2034.txt.gz"

echo "  Downloading OD Cingulum_hippocampus_L (TBSS, IDP 2035)..."
curl -s -o "${OUTPUT_DIR}/OD_TBSS_Cingulum_hippocampus_L_IDP2035.txt.gz" "https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/2035.txt.gz"

echo ""
echo "=============================================="
echo "Download complete!"
echo "Files saved to: ${OUTPUT_DIR}/"
echo "Total files: $(ls -1 ${OUTPUT_DIR}/*.txt.gz 2>/dev/null | wc -l)"
echo "Total size: $(du -sh ${OUTPUT_DIR}/ | cut -f1)"
echo "=============================================="
