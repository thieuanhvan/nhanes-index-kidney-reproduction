"""M1 builder, 2ecb3dff61 (trunk fat-to-muscle ratio and DKD; source withheld).
Spec: trunk-FMR per 0.1 unit, DKD (UACR >= 30 or eGFR < 60 in diabetes) vs all other adults, fully adjusted model,
NHANES 1999-2006 + 2011-2018 with DXA.

NOT TESTABLE (rule 9), decided 15/11/2026 before any build was executed:
- The exposure needs DXA trunk fat and trunk lean mass. Only the 2011-2018 DXA files (DXX_G..DXX_J) are available here;
  the 1999-2006 DXA multiple-imputation files could not be obtained from the CDC public path tried, and downloads are
  not allowed at this stage, so 4 of the paper's 8 cycles have no exposure.
- The missing cycles are not a random half: 2011-2018 whole-body DXA covers ages 8-59 only (checked: age range of rows
  with DXXTRFAT is 8 to 59 in 2011-2012 and 2017-2018), so every participant aged 60+ in the paper's 17,859 comes from
  1999-2006. DKD is concentrated in older adults (the paper's DKD group has mean age 56), so a 2011-2018-only build
  would analyse a different population (adults 20-59) with a different weight pooling, not the paper's specification.
- The 1999-2006 DXA release is 5 multiple imputations; the paper does not say how they were combined, a second
  unresolved point even if the files were present.
A builder restricted to 2011-2018 would therefore not reconstruct the first-combination specification.
"""

REC_ID = "2ecb3dff61"


def build(lib):
    return {"not_testable": "exposure (DXA trunk fat / trunk lean) absent for 1999-2006 (DXA multiple-imputation files not "
                            "obtainable at the CDC public path tried); 2011-2018 DXA covers ages 8-59 only, so all "
                            "participants aged 60+ of the paper's sample come from the missing cycles and the analytic "
                            "population cannot be reconstructed"}
