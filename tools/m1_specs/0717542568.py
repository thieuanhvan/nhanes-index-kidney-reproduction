"""M1 builder, 0717542568 (source withheld). Serum copper/zinc ratio, CKD, Model 2.
Reported estimate: copper/zinc < 1.2 category, OR 0.57 (Results text; Table 2 referenced but absent from the full text).
Spec ambiguity 'reference category': reading_1 = Methods/Results 'Q1 set as referent' (lowest quartile).
Under reading_1 the '< 1.2' category is the lowest ratio category (the next is '1.20-1.46'), that is Q1, the reference
itself, so the reported OR 0.57 has no defined contrast; the other cut points and the reference category bounds are only
in Table 2, which is not available. Reading_2 (higher-ratio reference) would be buildable as a sensitivity combination:
cycles 2011-2016, age >= 20, CUSEZN subsample weight WTSA2YR / 3, Cockcroft-Gault style eGFR as written
((140 - age) x weight x 1.23 / Scr men, 1.03 women; Scr in umol/L implied by the constants), CKD = eGFR < 60 or UACR >= 30,
Model 2 covariates age, gender, ethnicity (4 groups), education (3), household income (3 classes), BMI, diabetes, smoking.
"""

REC_ID = "0717542568"


def build(lib):
    return {"not_testable": "specification undefined under reading_1: the reported '< 1.2' copper/zinc category is the lowest "
                            "ratio category, i.e. Q1, which the paper sets as the reference (Methods and Results), so its OR 0.57 "
                            "has no defined contrast; cut points and reference bounds are only in Table 2, absent from the full text"}
