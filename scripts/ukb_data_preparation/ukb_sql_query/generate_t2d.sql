.header on
.mode csv
.output T2D.csv
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
SELECT death.sample_id AS sample_id,
    CAST(death.pheno AS INT) AS age,
    0 AS array,
    "Death" AS ICD
FROM f40001 death
GROUP BY death.sample_id
HAVING COUNT(DISTINCT death.pheno) = 1; -- don't allow multiple age of death

CREATE TEMP Table T2D AS
SELECT info.sample_id AS sample_id,
    1 AS Pheno,
    MIN(info.age) AS Age
FROM(
        SELECT icd10.sample_id AS sample_id,
            MIN(age.age) AS age
        FROM f41204 icd10 -- Diagnoses - secondary ICD10
            JOIN Age_Diagnosis age ON icd10.sample_id = age.sample_id
            AND icd10.array = age.array
            AND age.ICD = "ICD10"
        WHERE icd10.pheno LIKE '"E11%"'
        GROUP BY icd10.sample_id 
        -------------------------------------------------------
        UNION
        SELECT icd10.sample_id AS sample_id,
            MIN(age.age) AS age
        FROM f41202 icd10 -- Diagnoses - main ICD10
            JOIN Age_Diagnosis age ON icd10.sample_id = age.sample_id
            AND icd10.array = age.array
            AND age.ICD = "ICD10"
        WHERE icd10.pheno LIKE '"E11%"'
        GROUP BY icd10.sample_id 
        -------------------------------------------------------
        UNION
        SELECT icd10.sample_id AS sample_id,
            MIN(age.age) AS age
        FROM f41270 icd10 -- Diagnoses - ICD10
            JOIN Age_Diagnosis age ON icd10.sample_id = age.sample_id
            AND icd10.array = age.array
            AND age.ICD = "ICD10"
        WHERE icd10.pheno LIKE '"E11%"'
        GROUP BY icd10.sample_id 
        -------------------------------------------------------
        UNION
        SELECT icd10.sample_id AS sample_id,
            MIN(age.age) AS age
        FROM f40002 icd10 -- Contributory (secondary) causes of death: ICD10
            JOIN Age_Diagnosis age ON icd10.sample_id = age.sample_id
            AND age.ICD = "Death"
        WHERE icd10.pheno LIKE '"E11%"'
        GROUP BY icd10.sample_id 
        -------------------------------------------------------
        UNION
        SELECT icd10.sample_id AS sample_id,
            MIN(age.age) AS age
        FROM f40001 icd10 -- Underlying (primary) cause of death: ICD10
            JOIN Age_Diagnosis age ON icd10.sample_id = age.sample_id
            AND age.ICD = "Death"
        WHERE icd10.pheno LIKE '"E11%"'
        GROUP BY icd10.sample_id 
        -------------------------------------------------------
        UNION
        SELECT icd9.sample_id AS sample_id,
            MIN(age.age) AS age
        FROM f41271 icd9 -- Diagnoses - ICD9
            JOIN Age_Diagnosis age ON icd9.sample_id = age.sample_id
            AND icd9.array = age.array
            AND age.ICD = "ICD9"
        WHERE icd9.pheno LIKE '"250_"'
            OR icd9.pheno LIKE '"250_0"'
            OR icd9.pheno LIKE '"250_2"'
        GROUP BY icd9.sample_id 
        -------------------------------------------------------
        UNION
        SELECT icd9.sample_id AS sample_id,
            MIN(age.age) AS age
        FROM f41205 icd9 -- Diagnoses - secondary ICD9
            JOIN Age_Diagnosis age ON icd9.sample_id = age.sample_id
            AND icd9.array = age.array
            AND age.ICD = "ICD9"
        WHERE icd9.pheno LIKE '"250_"'
            OR icd9.pheno LIKE '"250_0"'
            OR icd9.pheno LIKE '"250_2"'
        GROUP BY icd9.sample_id 
        -------------------------------------------------------
        UNION
        SELECT icd9.sample_id AS sample_id,
            MIN(age.age) AS age
        FROM f41203 icd9 -- Diagnoses - main ICD9
            JOIN Age_Diagnosis age ON icd9.sample_id = age.sample_id
            AND icd9.array = age.array
            AND age.ICD = "ICD9"
        WHERE icd9.pheno LIKE '"250_"'
            OR icd9.pheno LIKE '"250_0"'
            OR icd9.pheno LIKE '"250_2"'
        GROUP BY icd9.sample_id 
        
    ) info
    GROUP BY info.sample_id;


CREATE TEMP TABLE exclude AS
SELECT info.sample_id AS sample_id,
    1 AS Pheno,
    MIN(info.age) AS Age
FROM(
        SELECT icd10.sample_id AS sample_id,
            MIN(age.age) AS age
        FROM f41204 icd10 -- Diagnoses - secondary ICD10
            JOIN Age_Diagnosis age ON icd10.sample_id = age.sample_id
            AND icd10.array = age.array
            AND age.ICD = "ICD10"
        WHERE icd10.pheno LIKE '"E10%"'  
            OR icd10.pheno LIKE '"O24%"'
        -------------------------------------------------------
        UNION ALL
        SELECT icd10.sample_id AS sample_id,
            MIN(age.age) AS age
        FROM f41202 icd10 -- Diagnoses - main ICD10
            JOIN Age_Diagnosis age ON icd10.sample_id = age.sample_id
            AND icd10.array = age.array
            AND age.ICD = "ICD10"
        WHERE icd10.pheno LIKE '"E10%"'  
            OR icd10.pheno LIKE '"O24%"'
        -------------------------------------------------------
        UNION ALL
        SELECT icd10.sample_id AS sample_id,
            MIN(age.age) AS age
        FROM f41270 icd10 -- Diagnoses - ICD10
            JOIN Age_Diagnosis age ON icd10.sample_id = age.sample_id
            AND icd10.array = age.array
            AND age.ICD = "ICD10"
        WHERE icd10.pheno LIKE '"E10%"'  
            OR icd10.pheno LIKE '"O24%"'
        -------------------------------------------------------
        UNION ALL
        SELECT icd10.sample_id AS sample_id,
            MIN(age.age) AS age
        FROM f40002 icd10 -- Contributory (secondary) causes of death: ICD10
            JOIN Age_Diagnosis age ON icd10.sample_id = age.sample_id
            AND age.ICD = "Death"
        WHERE icd10.pheno LIKE '"E10%"'  
            OR icd10.pheno LIKE '"O24%"'
        -------------------------------------------------------
        UNION ALL
        SELECT icd10.sample_id AS sample_id,
            MIN(age.age) AS age
        FROM f40001 icd10 -- Underlying (primary) cause of death: ICD10
            JOIN Age_Diagnosis age ON icd10.sample_id = age.sample_id
            AND age.ICD = "Death"
        WHERE icd10.pheno LIKE '"E10%"'  
            OR icd10.pheno LIKE '"O24%"'
        -------------------------------------------------------
        UNION ALL
        SELECT icd9.sample_id AS sample_id,
            MIN(age.age) AS age
        FROM f41271 icd9 -- Diagnoses - ICD9
            JOIN Age_Diagnosis age ON icd9.sample_id = age.sample_id
            AND icd9.array = age.array
            AND age.ICD = "ICD9"
        WHERE icd9.pheno LIKE '"250_1%"'
            OR icd9.pheno LIKE '"250_3"' 
            OR icd9.pheno LIKE '"648_"'
            -------------------------------------------------------
        UNION ALL
        SELECT icd9.sample_id AS sample_id,
            MIN(age.age) AS age
        FROM f41205 icd9 -- Diagnoses - secondary ICD9
            JOIN Age_Diagnosis age ON icd9.sample_id = age.sample_id
            AND icd9.array = age.array
            AND age.ICD = "ICD9"
        WHERE icd9.pheno LIKE '"250_1%"'
            OR icd9.pheno LIKE '"250_3"' 
            OR icd9.pheno LIKE '"648_"'
            -------------------------------------------------------
        UNION ALL
        SELECT icd9.sample_id AS sample_id,
            MIN(age.age) AS age
        FROM f41203 icd9 -- Diagnoses - main ICD9
            JOIN Age_Diagnosis age ON icd9.sample_id = age.sample_id
            AND icd9.array = age.array
            AND age.ICD = "ICD9"
        WHERE icd9.pheno LIKE '"250_1%"'
            OR icd9.pheno LIKE '"250_3"'
            OR icd9.pheno LIKE '"648_"'
    ) AS info
GROUP BY info.sample_id;

-- Remove T1D and T2D from controls
CREATE TEMP TABLE Controls AS
SELECT s.sample_id AS sample_id,
    0 AS Pheno,
    MAX(age.Pheno) AS Age
FROM Participant s
    JOIN f21003 age ON age.sample_id = s.sample_id
WHERE s.withdrawn = 0
    AND s.sample_id NOT IN (
        -------------------------------------------------------
        SELECT sample_id
        FROM T2D -- T2D cases
        -------------------------------------------------------
        UNION ALL
        SELECT sample_id
        FROM exclude -- T1D and Gestational diabetes
    )
GROUP BY s.sample_id;


INSERT INTO T2D(sample_id, pheno, age)
SELECT *
FROM Controls;


-- There seems to be only 6 T1D without T2D?

/* Now obtain the required data structure*/
SELECT s.sample_id as FID,
    s.sample_id as IID,
    t2d.age as Age,
    sex.Pheno as Sex,
    centre.Pheno as Centre,
    t2d.pheno as T2D
FROM Participant s
    JOIN f31 sex ON s.sample_id = sex.sample_id
    AND sex.instance = 0
    JOIN f21001 bmi ON s.sample_id = bmi.sample_id
    AND bmi.instance = 0
    JOIN T2D t2d ON s.sample_id = t2d.sample_id 
    AND t2d.age > 35 -- apparently T2D is adult onset
    JOIN f54 centre ON s.sample_id = centre.sample_id
    AND centre.instance = 0
WHERE s.withdrawn = 0
    AND s.sample_id NOT IN (
        SELECT sample_id
        FROM exclude -- T1D and Gestational diabetes
    )
GROUP BY s.sample_id;
.quit