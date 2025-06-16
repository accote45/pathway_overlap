.mode csv
.header on 
.output CAD.csv

/* Extract all CAD Cases and the instance when they got CAD diagnosis */
CREATE TEMP TABLE CAD_Patient AS
SELECT DISTINCT sample_id,
    MIN(age) as age,
    1 AS pheno
FROM 
(
    SELECT DISTINCT icd10.sample_id AS sample_id,
        /* Define age as the earliest when sample has CAD (therefore MIN)*/
            MIN(strftime('%Y', age.pheno) - CAST(birth.pheno AS INT)) AS age
    FROM    f41270 icd10 
            LEFT JOIN f41280 age ON
                icd10.sample_id = age.sample_id 
                AND icd10.array = age.array /* Need to match by array instead of instance */
            LEFT JOIN f34 birth ON
                icd10.sample_id = birth.sample_id
    WHERE   icd10.pheno LIKE '"I21_"' 
            OR icd10.pheno LIKE '"I22_"' 
            OR icd10.pheno LIKE '"I23_"' 
            OR icd10.pheno LIKE '"I24_"' 
            OR icd10.pheno LIKE '"I252"'
    GROUP BY icd10.sample_id
    UNION
    SELECT DISTINCT icd10.sample_id AS sample_id,
            MIN(strftime('%Y', age.pheno) - CAST(birth.pheno AS INT)) AS age
    FROM    f41202 icd10 
            LEFT JOIN f41280 age ON
                icd10.sample_id = age.sample_id 
                AND icd10.array = age.array
            LEFT JOIN f34 birth ON
                icd10.sample_id = birth.sample_id
    WHERE   icd10.pheno LIKE '"I21_"' 
            OR icd10.pheno LIKE '"I22_"' 
            OR icd10.pheno LIKE '"I23_"' 
            OR icd10.pheno LIKE '"I24_"' 
            OR icd10.pheno LIKE '"I252"'
    GROUP BY icd10.sample_id
    UNION
    SELECT DISTINCT icd10.sample_id AS sample_id,
            MIN(strftime('%Y', age.pheno) - CAST(birth.pheno AS INT)) AS age
    FROM    f41204 icd10 
            LEFT JOIN f41280 age ON
                icd10.sample_id = age.sample_id 
                AND icd10.array = age.array
            LEFT JOIN f34 birth ON
                icd10.sample_id = birth.sample_id
    WHERE   icd10.pheno LIKE '"I21_"' 
            OR icd10.pheno LIKE '"I22_"' 
            OR icd10.pheno LIKE '"I23_"' 
            OR icd10.pheno LIKE '"I24_"' 
            OR icd10.pheno LIKE '"I252"'
    GROUP BY icd10.sample_id
    UNION 
    SELECT DISTINCT icd9.sample_id AS sample_id,
            MIN(strftime('%Y', age.pheno)-CAST(birth.pheno AS INT)) AS age
    FROM    f41203 icd9
            LEFT JOIN f41281 age ON 
                icd9.sample_id = age.sample_id
                AND icd9.array = age.array
            LEFT JOIN f34 birth ON 
                icd9.sample_id = birth.sample_id
    WHERE   icd9.pheno LIKE '"410_"' 
            OR icd9.pheno LIKE '"411_"' 
            OR icd9.pheno LIKE '"412_"' 
    GROUP BY icd9.sample_id
    UNION 
    SELECT DISTINCT icd9.sample_id AS sample_id,
            MIN(strftime('%Y', age.pheno)-CAST(birth.pheno AS INT)) AS age
    FROM    f41205 icd9
            LEFT JOIN f41281 age ON 
                icd9.sample_id = age.sample_id
                AND icd9.array = age.array
            LEFT JOIN f34 birth ON 
                icd9.sample_id = birth.sample_id
    WHERE   icd9.pheno LIKE '"410_"' 
            OR icd9.pheno LIKE '"411_"' 
            OR icd9.pheno LIKE '"412_"' 
    GROUP BY icd9.sample_id
    UNION 
    SELECT DISTINCT icd9.sample_id AS sample_id,
            MIN(strftime('%Y', age.pheno)-CAST(birth.pheno AS INT)) AS age
    FROM    f41271 icd9
            LEFT JOIN f41281 age ON 
                icd9.sample_id = age.sample_id
                AND icd9.array = age.array
            LEFT JOIN f34 birth ON 
                icd9.sample_id = birth.sample_id
    WHERE   icd9.pheno LIKE '"410_"' 
            OR icd9.pheno LIKE '"411_"' 
            OR icd9.pheno LIKE '"412_"' 
    GROUP BY icd9.sample_id
    UNION 
    SELECT DISTINCT death.sample_id AS sample_id,
            MIN(CAST(age.pheno AS INT)) AS age
    FROM    f40001 death
            LEFT JOIN f40007 age ON 
                death.sample_id = age.sample_id
                AND death.instance = age.instance
    WHERE   death.pheno LIKE '"I21_"' 
            OR death.pheno LIKE '"I22_"' 
            OR death.pheno LIKE '"I23_"' 
            OR death.pheno LIKE '"I24_"' 
            OR death.pheno LIKE '"I252"'
    GROUP BY death.sample_id
    UNION 
    SELECT DISTINCT death.sample_id AS sample_id,
            MIN(CAST(age.pheno AS INT)) AS age
    FROM    f40002 death
            LEFT JOIN f40007 age ON 
                death.sample_id = age.sample_id
                AND death.instance = age.instance
    WHERE   death.pheno LIKE '"I21_"' 
            OR death.pheno LIKE '"I22_"' 
            OR death.pheno LIKE '"I23_"' 
            OR death.pheno LIKE '"I24_"' 
            OR death.pheno LIKE '"I252"'
    GROUP BY death.sample_id
    HAVING COUNT(DISTINCT age.pheno) = 1
    UNION 
    SELECT DISTINCT heart.sample_id AS sample_id,
            MIN(CAST(age.pheno AS INT)) AS age
    FROM    f6150 heart
            LEFT JOIN f21003 age ON
                heart.sample_id = age.sample_id 
                AND heart.instance = age.instance
    WHERE   heart.pheno =1
    GROUP BY heart.sample_id
    UNION 
    SELECT DISTINCT illness.sample_id AS sample_id,
            MIN(CAST(age.pheno AS INT)) AS age
    FROM    f20002 illness
            LEFT JOIN f21003 age ON
                illness.sample_id = age.sample_id
                AND illness.instance = age.instance
    WHERE   illness.pheno = 1075
    GROUP BY illness.sample_id
    UNION
    SELECT DISTINCT operation.sample_id AS sample_id,
            MIN(CAST(age.pheno AS INT)) AS age
    FROM    f20004 operation
            LEFT JOIN f21003 age ON
                operation.sample_id = age.sample_id
                AND operation.instance = age.instance
    WHERE   operation.pheno in (1070, 1095, 1523) 
    GROUP BY operation.sample_id
    UNION
    SELECT DISTINCT operation.sample_id AS sample_id,
            MIN(strftime('%Y', age.pheno)-CAST(birth.pheno AS INT)) AS age
    FROM    f41200 operation
            LEFT JOIN f41260 age ON 
                operation.sample_id = age.sample_id
                AND operation.array = age.array
            LEFT JOIN f34 birth ON 
                operation.sample_id = birth.sample_id
    WHERE   operation.pheno LIKE '"K40_"' 
            OR operation.pheno LIKE '"K41_"' 
            OR operation.pheno LIKE '"K42_"' 
            OR operation.pheno LIKE '"K43_"' 
            OR operation.pheno LIKE '"K44_"' 
            OR operation.pheno LIKE '"K45_"' 
            OR operation.pheno LIKE '"K46_"' 
            OR operation.pheno LIKE '"K49_"' 
            OR operation.pheno LIKE '"K501"' 
            OR operation.pheno LIKE '"K75_"'
    GROUP BY operation.sample_id
    UNION
    SELECT DISTINCT operation.sample_id AS sample_id,
            MIN(strftime('%Y', age.pheno)-CAST(birth.pheno AS INT)) AS age
    FROM    f41272 operation
            LEFT JOIN f41282 age ON 
                operation.sample_id = age.sample_id
                AND operation.array = age.array
            LEFT JOIN f34 birth ON 
                operation.sample_id = birth.sample_id
    WHERE   operation.pheno LIKE '"K40_"' 
            OR operation.pheno LIKE '"K41_"' 
            OR operation.pheno LIKE '"K42_"' 
            OR operation.pheno LIKE '"K43_"' 
            OR operation.pheno LIKE '"K44_"' 
            OR operation.pheno LIKE '"K45_"' 
            OR operation.pheno LIKE '"K46_"' 
            OR operation.pheno LIKE '"K49_"' 
            OR operation.pheno LIKE '"K501"' 
            OR operation.pheno LIKE '"K75_"'
    GROUP BY operation.sample_id
) AS subquery
GROUP BY sample_id;


/* Get all Control, defined as sample who never got CAD */
CREATE TEMP TABLE Healthy AS 
    SELECT DISTINCT Participant.sample_id
FROM Participant
EXCEPT 
    SELECT DISTINCT CAD_Patient.sample_id
    FROM CAD_Patient;

/* Append to the CAD data structure */
INSERT INTO CAD_Patient(sample_id, pheno, age) 
    SELECT  Healthy.sample_id, 
            0 AS pheno,
            /* Use the latest age of the healthy patients */
            MAX(CAST(age.pheno AS INT)) AS age
    FROM    Healthy
            LEFT JOIN f21003 age ON
            Healthy.sample_id=age.sample_id
    GROUP BY Healthy.sample_id;

/* Now obtain the required data structure*/
SELECT  s.sample_id as FID, 
        s.sample_id as IID,
        cad.age as Age,
        sex.Pheno as Sex,
        centre.Pheno as Centre,
        cad.pheno as CAD
FROM    Participant s 
        LEFT JOIN f31 sex ON
            s.sample_id=sex.sample_id 
            AND sex.instance = 0
        LEFT JOIN CAD_Patient cad ON 
            s.sample_id=cad.sample_id
        LEFT JOIN f54 centre ON 
            s.sample_id=centre.sample_id 
            AND centre.instance = 0
group by s.sample_id;    
.quit



