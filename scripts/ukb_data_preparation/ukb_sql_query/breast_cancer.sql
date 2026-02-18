.header on
.mode csv
.output BreastCancer.csv
CREATE TEMP TABLE Age_Diagnosis AS
SELECT birth.sample_id AS sample_id,
    age.array AS array,
    strftime('%Y', age.pheno) - CAST(birth.pheno AS INT) AS age,
    "ICD10" AS ICD
FROM f34 birth
    JOIN f41280 age ON age.sample_id = birth.sample_id
UNION ALL
SELECT birth.sample_id AS sample_id,
    age.array AS array,
    strftime('%Y', age.pheno) - CAST(birth.pheno AS INT) AS age,
    "ICD9" AS ICD
FROM f34 birth
    JOIN f41281 age ON age.sample_id = birth.sample_id
UNION ALL
SELECT sample_id,
    CAST(pheno AS INT) AS age,
    instance AS array,
    "Self" AS ICD
FROM f21003
UNION ALL
SELECT death.sample_id AS sample_id,
    CAST(death.pheno AS INT) AS age,
    0 AS array,
    "Death" AS ICD
FROM f40001 death
GROUP BY death.sample_id
HAVING COUNT(DISTINCT death.pheno) = 1; -- don't allow multiple age of death


CREATE TEMP TABLE Cases AS
SELECT ICD.sample_id AS sample_id,
    1 AS Pheno,
    MIN(ICD.age) AS Age
FROM (
        SELECT icd10.sample_id AS sample_id,
            MIN(age.age) AS age
        FROM f41204 icd10 -- Diagnoses - secondary ICD10
            JOIN Age_Diagnosis age ON icd10.sample_id = age.sample_id
            AND icd10.array = age.array
            AND age.ICD = "ICD10"
        WHERE icd10.pheno LIKE '"C50_"'
        GROUP BY icd10.sample_id 
        -------------------------------------------------------
        UNION
        SELECT icd10.sample_id AS sample_id,
            MIN(age.age) AS age
        FROM f41202 icd10 -- Diagnoses - main ICD10
            JOIN Age_Diagnosis age ON icd10.sample_id = age.sample_id
            AND icd10.array = age.array
            AND age.ICD = "ICD10"
        WHERE icd10.pheno LIKE '"C50_"'
        GROUP BY icd10.sample_id 
        -------------------------------------------------------
        UNION
        SELECT icd10.sample_id AS sample_id,
            MIN(age.age) AS age
        FROM f41270 icd10 -- Diagnoses - ICD10
            JOIN Age_Diagnosis age ON icd10.sample_id = age.sample_id
            AND icd10.array = age.array
            AND age.ICD = "ICD10"
        WHERE icd10.pheno LIKE '"C50_"'
        GROUP BY icd10.sample_id 
        -------------------------------------------------------
        UNION
        SELECT icd10.sample_id AS sample_id,
            MIN(age.age) AS age
        FROM f40002 icd10 -- Contributory (secondary) causes of death: ICD10
            JOIN Age_Diagnosis age ON icd10.sample_id = age.sample_id
            AND age.ICD = "Death"
        WHERE icd10.pheno LIKE '"C50_"'
        GROUP BY icd10.sample_id 
        -------------------------------------------------------
        UNION
        SELECT icd10.sample_id AS sample_id,
            MIN(age.age) AS age
        FROM f40001 icd10 -- Underlying (primary) cause of death: ICD10
            JOIN Age_Diagnosis age ON icd10.sample_id = age.sample_id
            AND age.ICD = "Death"
        WHERE icd10.pheno LIKE '"C50_"'
        GROUP BY icd10.sample_id 
        -------------------------------------------------------
        UNION
        SELECT icd9.sample_id AS sample_id,
            MIN(age.age) AS age
        FROM f41271 icd9 -- Diagnoses - ICD9
            JOIN Age_Diagnosis age ON icd9.sample_id = age.sample_id
            AND icd9.array = age.array
            AND age.ICD = "ICD9"
        WHERE icd9.pheno LIKE '"174%"'
        GROUP BY icd9.sample_id 
        -------------------------------------------------------
        UNION
        SELECT icd9.sample_id AS sample_id,
            MIN(age.age) AS age
        FROM f41205 icd9 -- Diagnoses - secondary ICD9
            JOIN Age_Diagnosis age ON icd9.sample_id = age.sample_id
            AND icd9.array = age.array
            AND age.ICD = "ICD9"
        WHERE icd9.pheno LIKE '"174%"'
        GROUP BY icd9.sample_id 
        -------------------------------------------------------
        UNION
        SELECT icd9.sample_id AS sample_id,
            MIN(age.age) AS age
        FROM f41203 icd9 -- Diagnoses - main ICD9
            JOIN Age_Diagnosis age ON icd9.sample_id = age.sample_id
            AND icd9.array = age.array
            AND age.ICD = "ICD9"
        WHERE icd9.pheno LIKE '"174%"'
        GROUP BY icd9.sample_id
        -------------------------------------------------------
        UNION
        SELECT self.sample_id AS sample_id,
            MIN(age.age) AS age
        FROM f20001 self -- Cancer code, self-reported
            JOIN Age_Diagnosis age ON self.sample_id = age.sample_id
            AND self.instance = age.array
            AND age.ICD = "Self"
        WHERE self.pheno = 1002 
        GROUP BY self.sample_id 
    ) AS ICD
GROUP BY ICD.sample_id;

CREATE TEMP TABLE Controls AS
SELECT s.sample_id AS sample_id,
    0 AS Pheno,
    MAX(age.Pheno) AS Age
FROM Participant s
    JOIN f21003 age ON age.sample_id = s.sample_id
WHERE s.withdrawn = 0
    AND s.sample_id NOT IN (
        SELECT sample_id
        FROM f20001 -- Cancer code, self-reported
        -------------------------------------------------------
        UNION ALL
        SELECT sample_id
        FROM Cases -- Breast Cancer Cases
        -------------------------------------------------------
        UNION ALL
        SELECT sample_id
        FROM f41204 -- Diagnoses - secondary ICD10
        WHERE pheno LIKE '"C%"'
            OR pheno LIKE '"D0%"'
            OR pheno LIKE '"D37%"'
            OR pheno LIKE '"D38%"'
            OR pheno LIKE '"D39%"'
            OR pheno LIKE '"D4%"' 
        -------------------------------------------------------
        UNION ALL
        SELECT sample_id
        FROM f41202 -- Diagnoses - main ICD10
        WHERE pheno LIKE '"C%"'
            OR pheno LIKE '"D0%"'
            OR pheno LIKE '"D37%"'
            OR pheno LIKE '"D38%"'
            OR pheno LIKE '"D39%"'
            OR pheno LIKE '"D4%"' 
        -------------------------------------------------------
        UNION ALL
        SELECT sample_id
        FROM f41270 -- Diagnoses - ICD10
        WHERE pheno LIKE '"C%"'
            OR pheno LIKE '"D0%"'
            OR pheno LIKE '"D37%"'
            OR pheno LIKE '"D38%"'
            OR pheno LIKE '"D39%"'
            OR pheno LIKE '"D4%"' 
        -------------------------------------------------------
        UNION ALL
        SELECT sample_id
        FROM f40002 -- Contributory (secondary) causes of death: ICD10
        WHERE pheno LIKE '"C%"'
            OR pheno LIKE '"D0%"'
            OR pheno LIKE '"D37%"'
            OR pheno LIKE '"D38%"'
            OR pheno LIKE '"D39%"'
            OR pheno LIKE '"D4%"' 
        -------------------------------------------------------
        UNION ALL
        SELECT sample_id
        FROM f40001 -- Underlying (primary) cause of death: ICD10
        WHERE pheno LIKE '"C%"'
            OR pheno LIKE '"D0%"'
            OR pheno LIKE '"D37%"'
            OR pheno LIKE '"D38%"'
            OR pheno LIKE '"D39%"'
            OR pheno LIKE '"D4%"' 
        -------------------------------------------------------
        UNION ALL
        SELECT sample_id
        FROM f41271 -- Diagnoses - ICD9
        WHERE pheno LIKE '"1%"'
            OR pheno LIKE '"20%"'
            OR pheno LIKE '"23%"' 
        -------------------------------------------------------
        UNION ALL
        SELECT sample_id
        FROM f41205 -- Diagnoses - secondary ICD9
        WHERE pheno LIKE '"1%"'
            OR pheno LIKE '"20%"'
            OR pheno LIKE '"23%"' 
        -------------------------------------------------------
        UNION ALL
        SELECT sample_id
        FROM f41203 -- Diagnoses - main ICD9
        WHERE pheno LIKE '"1%"'
            OR pheno LIKE '"20%"'
            OR pheno LIKE '"23%"'
    )
GROUP BY s.sample_id;


INSERT INTO Cases(sample_id, pheno, age)
SELECT *
FROM Controls;


/* Now obtain the required data structure*/
SELECT s.sample_id as FID,
    s.sample_id as IID,
    bc.age as Age,
    sex.Pheno as Sex,
    centre.Pheno as Centre,
    bc.pheno as BreastCancer
FROM Participant s
    JOIN f31 sex ON s.sample_id = sex.sample_id
    AND sex.instance = 0
    AND sex.pheno = 0 -- only use female
    JOIN Cases bc ON s.sample_id = bc.sample_id
    JOIN f54 centre ON s.sample_id = centre.sample_id
    AND centre.instance = 0
WHERE s.withdrawn = 0
GROUP BY s.sample_id;
.quit