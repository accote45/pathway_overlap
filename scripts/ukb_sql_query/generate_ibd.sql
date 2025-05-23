.mode csv
.output IBD.csv
.header on 
-- Extract IBD samples based on self report data
CREATE TEMP TABLE IBD_SELF_REPORTED AS
SELECT subquery.sample_id AS sample_id,
    MAX(
        CASE
            WHEN subquery.pheno = 1461 THEN subquery.age
            ELSE NULL
        END
    ) AS IBD,
    MAX(
        CASE
            WHEN subquery.pheno = 1462 THEN subquery.age
            ELSE NULL
        END
    ) AS Crohns,
    MAX(
        CASE
            WHEN subquery.pheno = 1463 THEN subquery.age
            ELSE NULL
        END
    ) AS Ulcerative
FROM (
        SELECT illness.sample_id AS sample_id,
            illness.pheno AS pheno,
            MIN(age.pheno) AS age
        FROM f20002 AS illness
            LEFT JOIN f20009 AS age ON illness.sample_id = age.sample_id
            AND illness.instance = age.instance
        WHERE illness.pheno IN (1461, 1462, 1463)
        GROUP BY illness.sample_id,
            illness.pheno
    ) AS subquery
GROUP BY subquery.sample_id;
-- Extract IBD samples using ICD codes
-- currently there's some problem where we are not picking up comorbid cases
-- also, there seems to be more IBD cases than we'd otherwise expected
CREATE TEMP TABLE IBD_ICD AS
SELECT subquery.sample_id AS sample_id,
    MAX(
        CASE
            WHEN subquery.pheno = "IBD" THEN subquery.age
            ELSE NULL
        END
    ) AS IBD,
    MAX(
        CASE
            WHEN subquery.pheno = "Crohns" THEN subquery.age
            ELSE NULL
        END
    ) AS Crohns,
    (
        CASE
            WHEN subquery.pheno = "Ulcerative" THEN subquery.age
        END
    ) AS Ulcerative
FROM (
        SELECT DISTINCT sub.sample_id AS sample_id,
            sub.pheno AS pheno,
            sub.age AS age
        FROM(
                SELECT DISTINCT icd10.sample_id AS sample_id,
                    /* Define age as the earliest when sample has CAD (therefore MIN)*/
                    MIN(
                        strftime('%Y', age.pheno) - CAST(birth.pheno AS INT)
                    ) AS age,
                    (
                        CASE
                            WHEN icd10.pheno LIKE '"K51_"' THEN "Ulcerative"
                            ELSE CASE
                                WHEN icd10.pheno LIKE '"K50_"' THEN "Crohns"
                                ELSE NULL
                            END
                        END
                    ) AS pheno
                FROM f41270 icd10
                    LEFT JOIN f41280 age ON icd10.sample_id = age.sample_id
                    AND icd10.array = age.array
                    /* Need to match by array instead of instance */
                    LEFT JOIN f34 birth ON icd10.sample_id = birth.sample_id
                WHERE icd10.pheno LIKE '"K51_"'
                    OR icd10.pheno LIKE '"K50_"'
                GROUP BY icd10.sample_id,
                    icd10.pheno
                UNION
                SELECT DISTINCT icd10.sample_id AS sample_id,
                    /* Define age as the earliest when sample has CAD (therefore MIN)*/
                    MIN(
                        strftime('%Y', age.pheno) - CAST(birth.pheno AS INT)
                    ) AS age,
                    (
                        CASE
                            WHEN icd10.pheno LIKE '"K51_"' THEN "Ulcerative"
                            ELSE CASE
                                WHEN icd10.pheno LIKE '"K50_"' THEN "Crohns"
                                ELSE NULL
                            END
                        END
                    ) AS pheno
                FROM f41202 icd10
                    LEFT JOIN f41280 age ON icd10.sample_id = age.sample_id
                    AND icd10.array = age.array
                    /* Need to match by array instead of instance */
                    LEFT JOIN f34 birth ON icd10.sample_id = birth.sample_id
                WHERE icd10.pheno LIKE '"K51_"'
                    OR icd10.pheno LIKE '"K50_"'
                GROUP BY icd10.sample_id,
                    icd10.pheno
                UNION
                SELECT DISTINCT icd10.sample_id AS sample_id,
                    /* Define age as the earliest when sample has CAD (therefore MIN)*/
                    MIN(
                        strftime('%Y', age.pheno) - CAST(birth.pheno AS INT)
                    ) AS age,
                    (
                        CASE
                            WHEN icd10.pheno LIKE '"K51_"' THEN "Ulcerative"
                            ELSE CASE
                                WHEN icd10.pheno LIKE '"K50_"' THEN "Crohns"
                                ELSE NULL
                            END
                        END
                    ) AS pheno
                FROM f41204 icd10
                    LEFT JOIN f41280 age ON icd10.sample_id = age.sample_id
                    AND icd10.array = age.array
                    /* Need to match by array instead of instance */
                    LEFT JOIN f34 birth ON icd10.sample_id = birth.sample_id
                WHERE icd10.pheno LIKE '"K51_"'
                    OR icd10.pheno LIKE '"K50_"'
                GROUP BY icd10.sample_id,
                    icd10.pheno
                UNION
                SELECT DISTINCT icd9.sample_id AS sample_id,
                    /* Define age as the earliest when sample has CAD (therefore MIN)*/
                    MIN(
                        strftime('%Y', age.pheno) - CAST(birth.pheno AS INT)
                    ) AS age,
                    (
                        CASE
                            WHEN icd9.pheno LIKE '"556_"' THEN "Ulcerative"
                            ELSE CASE
                                WHEN icd9.pheno LIKE '"555_"' THEN "Crohns"
                                ELSE NULL
                            END
                        END
                    ) AS pheno
                FROM f41203 icd9
                    LEFT JOIN f41281 age ON icd9.sample_id = age.sample_id
                    AND icd9.array = age.array
                    /* Need to match by array instead of instance */
                    LEFT JOIN f34 birth ON icd9.sample_id = birth.sample_id
                WHERE icd9.pheno LIKE '"555_"'
                    OR icd9.pheno LIKE '"556_"'
                GROUP BY icd9.sample_id,
                    icd9.pheno
                UNION
                SELECT DISTINCT icd9.sample_id AS sample_id,
                    /* Define age as the earliest when sample has CAD (therefore MIN)*/
                    MIN(
                        strftime('%Y', age.pheno) - CAST(birth.pheno AS INT)
                    ) AS age,
                    (
                        CASE
                            WHEN icd9.pheno LIKE '"556_"' THEN "Ulcerative"
                            ELSE CASE
                                WHEN icd9.pheno LIKE '"555_"' THEN "Crohns"
                                ELSE NULL
                            END
                        END
                    ) AS pheno
                FROM f41205 icd9
                    LEFT JOIN f41281 age ON icd9.sample_id = age.sample_id
                    AND icd9.array = age.array
                    /* Need to match by array instead of instance */
                    LEFT JOIN f34 birth ON icd9.sample_id = birth.sample_id
                WHERE icd9.pheno LIKE '"555_"'
                    OR icd9.pheno LIKE '"556_"'
                GROUP BY icd9.sample_id,
                    icd9.pheno
                UNION
                SELECT DISTINCT icd9.sample_id AS sample_id,
                    /* Define age as the earliest when sample has CAD (therefore MIN)*/
                    MIN(
                        strftime('%Y', age.pheno) - CAST(birth.pheno AS INT)
                    ) AS age,
                    (
                        CASE
                            WHEN icd9.pheno LIKE '"556_"' THEN "Ulcerative"
                            ELSE CASE
                                WHEN icd9.pheno LIKE '"555_"' THEN "Crohns"
                                ELSE NULL
                            END
                        END
                    ) AS pheno
                FROM f41271 icd9
                    LEFT JOIN f41281 age ON icd9.sample_id = age.sample_id
                    AND icd9.array = age.array
                    /* Need to match by array instead of instance */
                    LEFT JOIN f34 birth ON icd9.sample_id = birth.sample_id
                WHERE icd9.pheno LIKE '"555_"'
                    OR icd9.pheno LIKE '"556_"'
                GROUP BY icd9.sample_id,
                    icd9.pheno
                UNION
                SELECT DISTINCT death.sample_id AS sample_id,
                    MIN(CAST(age.pheno AS INT)) AS age,
                    (
                        CASE
                            WHEN death.pheno LIKE '"K51_"' THEN "Ulcerative"
                            ELSE CASE
                                WHEN death.pheno LIKE '"K50_"' THEN "Crohns"
                                ELSE NULL
                            END
                        END
                    ) AS pheno
                FROM f40001 death
                    LEFT JOIN f40007 age ON death.sample_id = age.sample_id
                    AND death.instance = age.instance
                WHERE death.pheno LIKE '"K51_"'
                    OR death.pheno LIKE '"K52_"'
                GROUP BY death.sample_id,
                    death.pheno
                UNION
                SELECT DISTINCT death.sample_id AS sample_id,
                    MIN(CAST(age.pheno AS INT)) AS age,
                    (
                        CASE
                            WHEN death.pheno LIKE '"K51_"' THEN "Ulcerative"
                            ELSE CASE
                                WHEN death.pheno LIKE '"K50_"' THEN "Crohns"
                                ELSE NULL
                            END
                        END
                    ) AS pheno
                FROM f40002 death
                    LEFT JOIN f40007 age ON death.sample_id = age.sample_id
                    AND death.instance = age.instance
                WHERE death.pheno LIKE '"K51_"'
                    OR death.pheno LIKE '"K50_"'
                GROUP BY death.sample_id,
                    death.pheno
            ) AS sub
        GROUP BY sub.sample_id,
            sub.pheno
    ) AS subquery
GROUP BY subquery.sample_id;
-- Merge ICD and self-report results
CREATE TEMP TABLE IBD_PATIENT AS
SELECT sub.sample_id AS sample_id,
    MIN(sub.IBD) AS IBD,
    MIN(sub.Crohns) AS Crohns,
    MIN(sub.Ulcerative) AS Ulcerative
FROM(
        SELECT IBD_ICD.sample_id AS sample_id,
            IBD_ICD.IBD AS IBD,
            IBD_ICD.Crohns AS Crohns,
            IBD_ICD.Ulcerative AS Ulcerative
        FROM IBD_ICD
        UNION
        SELECT IBD_SELF_REPORTED.sample_id AS sample_id,
            IBD_SELF_REPORTED.IBD AS IBD,
            IBD_SELF_REPORTED.Crohns AS Crohns,
            IBD_SELF_REPORTED.Ulcerative AS Ulcerative
        FROM IBD_SELF_REPORTED
    ) AS sub
GROUP BY sub.sample_id;
-- Generate final case control data
SELECT s.sample_id as FID,
    s.sample_id as IID,
    sex.Pheno as Sex,
    centre.Pheno as Centre,
    age.Pheno as Age,
    ibd.IBD as IBD,
    ibd.Crohns as Crohns,
    ibd.Ulcerative as Ulcerative
FROM Participant s
    LEFT JOIN f31 sex ON s.sample_id = sex.sample_id
    AND sex.instance = 0
    LEFT JOIN IBD_PATIENT ibd ON s.sample_id = ibd.sample_id
    LEFT JOIN f54 centre ON s.sample_id = centre.sample_id
    AND centre.instance = 0
    LEFT JOIN f21003 age ON s.sample_id = age.sample_id
    AND age.instance = 0;
.quit